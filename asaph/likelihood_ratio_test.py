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
from collections import defaultdict
import itertools
import random
import os
import sys

import numpy as np
import numpy.linalg as LA
from sklearn.linear_model import SGDClassifier

from asaph.ml import estimate_lr_iter
from asaph.ml import likelihood_ratio_test
from asaph.models import COUNTS_FEATURE_TYPE
from asaph.models import CATEGORIES_FEATURE_TYPE
from asaph.newioutils import read_features
from asaph.newioutils import deserialize
from asaph.newioutils import PROJECT_SUMMARY_FLNAME
from asaph.newioutils import serialize

def generate_training_set(feature_type, labels, features, sample_unknown_genotype):
    n_samples, n_features = features.shape

    if feature_type == COUNTS_FEATURE_TYPE:
        if n_features != 1:
            raise ValueError, "Feature matrix for counts feature type must have 1 column"
    elif feature_type == CATEGORIES_FEATURE_TYPE:
        if n_features != 3:
            raise ValueError, "Feature matrix for categories feature type must have 3 columns"
    else:
        raise Exception, "LR Test does not support feature type %s" % feature_type

    # we make 3 copies so we can impute each unknown genotype
    # with each of the 3 known genotypes
    N_COPIES = 3
    training_labels = np.zeros(N_COPIES * n_samples)
    training_features = np.zeros((N_COPIES * n_samples, n_features))

    for i in xrange(n_samples):
        for j in xrange(N_COPIES):
            idx = N_COPIES * i + j
            training_labels[idx] = labels[i]

            if not sample_unknown_genotype[i]:
                training_features[idx, :] = features[i, :]
            elif feature_type == CATEGORIES_FEATURE_TYPE:
                training_features[idx, j] = 1.
            elif feature_type == COUNTS_FEATURE_TYPE:
                training_features[idx, 0] = j

    return training_labels, training_features

        
def run_likelihood_ratio_tests(features, project_summary, args, stats_dir):
    if len(set(features.class_labels)) != 2:
        raise ValueError, "LR Test currently only supports 2 populations."

    n_iter = estimate_lr_iter(len(features.class_labels))

    fit_intercept = False
    if args.intercept == "free-parameter":
        fit_intercept = True
    
    lr = SGDClassifier(penalty="l2",
                       loss="log",
                       n_iter = n_iter,
                       fit_intercept = fit_intercept)
    
    flname = "snp_likelihood_ratio_tests.tsv"
    with open(os.path.join(stats_dir, flname), "w") as fl:
        next_output = 1
        for i, pair in enumerate(features.snp_feature_map.iteritems()):
            snp_label, feature_idx = pair
            chrom, pos = snp_label

            labels = np.array(features.class_labels)
            snp_features = features.feature_matrix[:, feature_idx]

            training_labels, training_features = generate_training_set(project_summary.feature_encoding,
                                                                       labels,
                                                                       snp_features,
                                                                       features.unknown_genotypes[snp_label])

            set_intercept_to_class_prob = False
            if args.intercept == "class-probabilities":
                set_intercept_to_class_prob = True
            
            p_value = likelihood_ratio_test((training_features, snp_features),
                                            (training_labels, labels),
                                            lr,
                                            set_intercept=set_intercept_to_class_prob)
            
            if i == next_output:
                print i, "SNP", snp_label, "has p-value", p_value
                next_output *= 2

            fl.write("\t".join([chrom, pos, str(p_value)]))
            fl.write("\n")
            
def parseargs():
    parser = argparse.ArgumentParser(description="Asaph - Likelihood Ratio Test")

    parser.add_argument("--workdir", type=str, help="Work directory", required=True)

    parser.add_argument("--intercept",
                        type=str,
                        default="class-probabilities",
                        choices=["class-probabilities",
                                 "none",
                                 "free-parameter"])

    parser.add_argument("--training-set",
                        type=str,
                        default="adjusted-for-unknown",
                        choices=["adjusted-for-unknown",
                                 "unadjusted"])

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    if not os.path.exists(args.workdir):
        print "Work directory '%s' does not exist." % args.workdir
        sys.exit(1)

    stats_dir = os.path.join(args.workdir, "statistics")
    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    project_summary = deserialize(os.path.join(args.workdir, PROJECT_SUMMARY_FLNAME))
    
    features = read_features(args.workdir)

    run_likelihood_ratio_tests(features,
                               project_summary,
                               args,
                               stats_dir)
