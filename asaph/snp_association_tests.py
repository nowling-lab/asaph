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
from sklearn.preprocessing import OneHotEncoder

from asaph.ml import estimate_lr_iter
from asaph.ml import likelihood_ratio_test
from asaph.newioutils import read_features
from asaph.newioutils import deserialize
from asaph.newioutils import PROJECT_SUMMARY_FLNAME
from asaph.newioutils import serialize

OUTPUT_DIR = "statistics"
OUTPUT_FLNAME = "snp_association_tests.tsv"

def prepare_model_variables(n_copies, testing_variables, null_variables):
    N_GENOTYPES = 3
    n_samples = testing_variables.shape[0]
    testing_features = np.zeros((n_copies * n_samples,
                                 testing_variables.shape[1]))
    if null_variables is not None:
        null_features = np.zeros((n_copies * n_samples,
                                  null_variables.shape[1]))
    for i in xrange(n_samples):
        for j in xrange(n_copies):
            idx = n_copies * i + j
            testing_features[idx, :] = testing_variables[i, :]
            if null_variables is not None:
                null_features[idx, :] = null_variables[i, :]

    if null_variables is not None:
        testing_features = np.hstack([testing_features,
                                      null_features])
    else:
        null_features = None

    return testing_features, null_features

def prepare_class_labels(n_copies, genotypes, class_labels=None):
    """
    Converts genotype categories to class labels and imputes unknown genotypes.
    If class_labels is None, 
    """

    n_samples = genotypes.shape[0]
    n_genotypes = genotypes.shape[1]
    N_GENOTYPES = 3
    if n_genotypes != N_GENOTYPES:
        raise Exception("Must be using genotype categories!")

    if n_samples != genotypes.shape[0]:
        raise Exception("Genotype and feature matrices do not have same number of samples!")

    if class_labels is None:
        class_labels = np.zeros(n_copies * n_samples)

    for i in xrange(n_samples):
        gt = None
        for j in xrange(N_GENOTYPES):
            if genotypes[i, j] == 1.0:
                gt = j
        
        for j in xrange(n_copies):
            idx = n_copies * i + j

            if gt is None:
                class_labels[idx] = j
            else:
                class_labels[idx] = gt

    return class_labels

def select_samples(data_model, selected_sample_ids):
    sample_indices = dict([(sample_name, idx) for idx, sample_name in
                               enumerate(data_model.sample_labels)])
        
    selected_indices = []
    for sample_id in selected_sample_ids:
        if sample_id not in sample_indices:
            raise Exception("Unknown sample id '%s', known ids are: %s" %
                                (sample_id, ",".join(data_model.sample_labels)))
        selected_indices.append(sample_indices[sample_id])

    return selected_indices

def run_likelihood_ratio_tests(args):
    if not os.path.exists(args.workdir):
        print "Work directory '%s' does not exist." % args.workdir
        sys.exit(1)

    stats_dir = os.path.join(args.workdir, OUTPUT_DIR)
    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    project_summary = deserialize(os.path.join(args.workdir, PROJECT_SUMMARY_FLNAME))
    
    data_model = read_features(args.workdir)
    genotypes = data_model.feature_matrix
    
    n_iter = estimate_lr_iter(len(data_model.class_labels))

    lr = SGDClassifier(penalty="l2",
                       loss="log",
                       n_iter = n_iter,
                       fit_intercept = False)

    testing_variables = np.array(data_model.class_labels).reshape(-1, 1)
    null_variables = None

    N_COPIES = 3
    class_labels = None
    testing_features, null_features = prepare_model_variables(N_COPIES,
                                                              testing_variables,
                                                              null_variables)
    
    with open(os.path.join(stats_dir, OUTPUT_FLNAME), "w") as fl:
        next_output = 1
        for i, pair in enumerate(data_model.snp_feature_map.iteritems()):
            pos_label, feature_idx = pair
            chrom, pos = pos_label

            pos_genotypes = genotypes[:, feature_idx]

            class_labels = prepare_class_labels(N_COPIES,
                                                pos_genotypes,
                                                class_labels)

            # since we make multiple copies of the original samples,
            # we need to scale the log loss so that it is correct for
            # the original sample size

            p_value = likelihood_ratio_test(testing_features,
                                            class_labels,
                                            lr,
                                            g_scaling_factor = 1.0 / N_COPIES)
            
            if i == next_output:
                print i, "Position", pos_label, "has p-value", p_value
                next_output *= 2

            fl.write("\t".join([chrom, pos, str(p_value)]))
            fl.write("\n")
            
def parseargs():
    parser = argparse.ArgumentParser(description="Asaph - Single SNP Association Tests")

    parser.add_argument("--workdir", type=str, help="Work directory", required=True)

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    run_likelihood_ratio_tests(args)
