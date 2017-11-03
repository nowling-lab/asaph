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
from asaph.newioutils import read_features
from asaph.newioutils import deserialize
from asaph.newioutils import PROJECT_SUMMARY_FLNAME
from asaph.newioutils import serialize

def generate_training_set(labels, features):
    # we make 3 copies so we can impute each unknown genotype
    # with each of the 3 known genotypes
    N_COPIES = 3
    training_labels = np.zeros(N_COPIES * labels.shape[0])
    training_features = np.zeros((N_COPIES * features.shape[0],
                                  features.shape[1]))

    for i in xrange(features.shape[0]):
        for j in xrange(N_COPIES):
            idx = N_COPIES * i + j
            training_labels[idx] = labels[i]
            
            if features[i, :].sum() == 1.:
                training_features[idx, :] = features[i, :]
            else:
                training_features[idx, j] = 1.

    return training_labels, training_features

        
def run_likelihood_ratio_tests(features, stats_dir):
    n_iter = estimate_lr_iter(len(features.class_labels))
    lr = SGDClassifier(penalty="l2",
                       loss="log",
                       n_iter = n_iter,
                       fit_intercept = False)
    
    flname = "snp_likelihood_ratio_tests.tsv"
    with open(os.path.join(stats_dir, flname), "w") as fl:
        next_output = 1
        for i, pair in enumerate(features.snp_feature_map.iteritems()):
            snp_label, feature_idx = pair
            chrom, pos = snp_label

            labels = np.array(features.class_labels)
            snp_features = features.feature_matrix[:, feature_idx]

            training_labels, training_features = generate_training_set(labels,
                                                                       snp_features)
            
            p_value = likelihood_ratio_test((training_features, snp_features),
                                            (training_labels, labels),
                                            lr,
                                            set_intercept=True)
            
            if i == next_output:
                print i, "SNP", snp_label, "has p-value", p_value
                next_output *= 2

            fl.write("\t".join([chrom, pos, str(p_value)]))
            fl.write("\n")
            
def parseargs():
    parser = argparse.ArgumentParser(description="Asaph - Likelihood Ratio Test")

    parser.add_argument("--workdir", type=str, help="Work directory", required=True)
    
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
                               stats_dir)
