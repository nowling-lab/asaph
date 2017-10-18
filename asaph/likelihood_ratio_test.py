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

gts = [(2., 0.), (0., 2.), (1., 1.), (1., 0.), (0., 1.), (0., 0.)]

subsets = { (0., 0.) : [(0., 0.), (0., 1.), (1., 0.), (1., 1.), (2., 0.), (0., 2.)],
            (1., 0.) : [(1., 0.), (2., 0.), (1., 1.)],
            (0., 1.) : [(0., 1.), (0., 2.), (1., 1.)],
            (1., 1.) : [(1., 1.)],
            (2., 0.) : [(2., 0.)],
            (0., 2.) : [(0., 2.)]
}

class ProbabilitySolver(object):
    def __init__(self):
        self.expected_prob = None
        
    def fit(self, X, y):
        assert X.shape[1] == 2, "X must have only two columns"

        gt_counts = defaultdict(lambda: {"total": 0, "class_1": 0})
        for r in xrange(X.shape[0]):
            gt_counts[tuple(X[r, :])]["total"] += 1
            if y[r] == 1.:
                gt_counts[tuple(X[r, :])]["class_1"] += 1.

        self.expected_prob = dict()
        for i, gt in enumerate(gts):
            total = 0.
            class_1 = 0.
            prob = 0.
            for subset in subsets[gt]:
                total += gt_counts[subset]["total"]
                class_1 += gt_counts[subset]["class_1"]
            if total != 0.:
                prob = class_1 / total
            self.expected_prob[gt] = prob

    def predict_proba(self, X):
        pred_prob = np.zeros((X.shape[0], 2))
        for r in xrange(X.shape[0]):
            prob = self.expected_prob[tuple(X[r, :])] 
            pred_prob[r, 0] = 1. - prob
            pred_prob[r, 1] = prob

        return pred_prob

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
    
    flname = "snp_likelihood_ratio_tests.txt"
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
