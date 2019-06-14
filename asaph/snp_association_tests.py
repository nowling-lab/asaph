"""
Copyright 2017 Ronald J. Nowling

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import argparse
from collections import Counter
from collections import defaultdict
import itertools
import random
import os
import sys

import numpy as np
import numpy.linalg as LA
from sklearn.linear_model import SGDClassifier
from sklearn.preprocessing import OneHotEncoder

from asaph.ml import estimate_lr_iter
from asaph.ml import likelihood_ratio_test
from asaph.ml import upsample_features
from asaph.newioutils import deserialize
from asaph.newioutils import PROJECT_SUMMARY_FLNAME
from asaph.newioutils import read_features
from asaph.newioutils import read_populations
from asaph.newioutils import serialize

def make_labels(sample_labels, sample_populations):
    populations = []
    for name in sample_labels:
        if name in sample_populations:
            populations.append(sample_populations[name])
        else:
            raise Exception("Population unknown for sample '%s'" % name)

    return populations

def run_lrtest_gt_dep(data_model, project_summary, args, stats_dir, class_labels):
    genotypes = data_model.feature_matrix
    
    n_iter = estimate_lr_iter(len(class_labels))

    lr = SGDClassifier(penalty="l2",
                       loss="log",
                       n_iter = n_iter,
                       fit_intercept = False)

    encoder = OneHotEncoder(sparse=False)
    
    flname = "snp_lrtests_gt.tsv"    
    n_pops = len(set(class_labels))
    
    with open(os.path.join(stats_dir, flname), "w") as fl:
        next_output = 1

        headers = ["chrom", "pos", "snp_p_value"]
        fl.write("\t".join(headers))
        fl.write("\n")
        
        for i, pair in enumerate(data_model.snp_feature_map.iteritems()):
            pos_label, feature_idx = pair
            chrom, pos = pos_label

            N_COPIES = 3
            pops, snp_genotypes = upsample_features(class_labels,
                                                    genotypes[:, feature_idx])
            
            # as we're using the genotypes as the labels,
            # they need to be one dimensional
            snp_genotypes = snp_genotypes.argmax(axis=1)

            p_value = 1.0
            if len(set(snp_genotypes)) > 1:
                # likewise, the pops need to 2D and one-hot encoded
                pops = pops.reshape(-1, 1)
                pops = encoder.fit_transform(pops)

                # since we make multiple copies of the original samples,
                # we need to scale the log loss so that it is correct for
                # the original sample size
                p_value = likelihood_ratio_test(pops,
                                                snp_genotypes,
                                                lr,
                                                g_scaling_factor = 1.0 / N_COPIES)
            
            if i == next_output:
                print i, "Position", pos_label, "has p-value", p_value
                next_output *= 2

            fl.write("\t".join([chrom, pos, "%.2E" % p_value]))
            fl.write("\n")


def run_lrtest_pop_dep(features, project_summary, args, stats_dir, class_labels):
    n_iter = estimate_lr_iter(len(class_labels))
    labels = np.array(class_labels)

    fit_intercept = False
    if args.intercept == "free-parameter":
        fit_intercept = True
    
    lr = SGDClassifier(penalty="l2",
                       loss="log",
                       n_iter = n_iter,
                       fit_intercept = fit_intercept)
    
    flname = "snp_lrtests_pop.tsv"
    with open(os.path.join(stats_dir, flname), "w") as fl:
        next_output = 1
        for i, pair in enumerate(features.snp_feature_map.iteritems()):
            snp_label, feature_idx = pair
            chrom, pos = snp_label

            snp_features = features.feature_matrix[:, feature_idx]

            if args.adjustment != "none":
                labels, snp_features = upsample_features(labels,
                                                         snp_features)

            # remove columns that are all zeros since these
            # aren't true degrees of freedom.  prevents
            # under-estimating significance
            if args.remove_empty_columns:
                mask = np.all(snp_features == 0.,
                              axis=0)
                snp_features = snp_features[:, ~mask]

            set_intercept_to_class_prob = False
            if args.intercept == "class-probabilities":
                set_intercept_to_class_prob = True

            scaling_factor = 1.0
            if args.adjustment == "training-set":
                snp_features = (snp_features,
                                features.feature_matrix[:, feature_idx])
                labels = (labels,
                          np.array(features.class_labels))
            elif args.adjustment == "scaling-factor":
                scaling_factor = 1.0 / 3.0

            p_value = likelihood_ratio_test(snp_features,
                                            labels,
                                            lr,
                                            set_intercept=set_intercept_to_class_prob,
                                            g_scaling_factor=scaling_factor)

            if i == next_output:
                print i, "SNP", snp_label, "has p-value", p_value
                next_output *= 2

            fl.write("\t".join([chrom, pos, "%.2E" % p_value]))
            fl.write("\n")
            
def parseargs():
    parser = argparse.ArgumentParser(description="Asaph - Likelihood Ratio Test")

    parser.add_argument("--workdir", type=str, help="Work directory", required=True)

    parser.add_argument("--populations", type=str, help="Populations file", required=True)

    parser.add_argument("--intercept",
                        type=str,
                        default="class-probabilities",
                        choices=["class-probabilities",
                                 "none",
                                 "free-parameter"])

    parser.add_argument("--adjustment",
                        type=str,
                        default="training-set",
                        choices=["training-set",
                                 "scaling-factor",
                                 "none"])

    parser.add_argument("--remove-empty-columns",
                        action="store_true")

    parser.add_argument("--dependent-variable",
                        type=str,
                        default="population",
                        choices=["genotype",
                                 "population"])

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    populations, sample_populations = read_populations(args.populations)

    n_pops = len(populations)
    if n_pops < 2:
        raise Exception("Need at least 2 populations!")
    
    if not os.path.exists(args.workdir):
        print "Work directory '%s' does not exist." % args.workdir
        sys.exit(1)

    stats_dir = os.path.join(args.workdir, "statistics")
    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    project_summary = deserialize(os.path.join(args.workdir, PROJECT_SUMMARY_FLNAME))
    
    features = read_features(args.workdir)

    class_labels = make_labels(features.sample_labels,
                               sample_populations)
    
    if args.dependent_variable == "genotype":
        run_lrtest_gt_dep(features,
                          project_summary,
                          args,
                          stats_dir,
                          class_labels)
    else:
        run_lrtest_pop_dep(features,
                           project_summary,
                           args,
                           stats_dir,
                           class_labels)
