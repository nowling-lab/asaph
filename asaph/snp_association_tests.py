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
from asaph.models import COUNTS_FEATURE_TYPE
from asaph.models import CATEGORIES_FEATURE_TYPE
from asaph.newioutils import read_features
from asaph.newioutils import deserialize
from asaph.newioutils import PROJECT_SUMMARY_FLNAME
from asaph.newioutils import serialize

OUTPUT_DIR = "statistics"
OUTPUT_FLNAME = "snp_association_tests.tsv"

NUMERICAL_TYPE = "n"
CATEGORICAL_TYPE = "c"
SAMPLE_IDENTIFIER = "id"
ALLOWED_TYPES = set([NUMERICAL_TYPE,
                     CATEGORICAL_TYPE,
                     SAMPLE_IDENTIFIER])

def parse_variables_file(flname):
    with open(flname) as fl:
        field_types = []
        for field_type in next(fl).split():
            if field_type not in ALLOWED_TYPES:
                raise Exception("Unknown field type '%s'" % field_type)
            field_types.append(field_type)

        # ignore sample identifier
        data = [[] for _ in field_types[1:]]
        sample_identifiers = []
        for ln in fl:
            fields = ln.split()
            if len(fields) != len(field_types):
                raise Exception("Found %s fields but have I have headers for %s fields" % (len(fields), len(field_types)))

            # first column is the sample identifier
            sample_identifiers.append(fields[0])
            
            row = [float(v) if t == NUMERICAL_TYPE else v
                   for t, v in zip(field_types[1:], fields[1:])]

            for i, v in enumerate(row):
                data[i].append(v)

    matrices = []
    for field_type, column in zip(field_types, data):
        column = np.array(column).reshape(-1, 1)
        if field_type == NUMERICAL_TYPE:
            matrices.append(column)
        else:
            encoder = OneHotEncoder(sparse=False)
            matrix = encoder.fit_transform(column)
            matrices.append(matrix)

    features = np.hstack(matrices)

    return sample_identifiers, features
            

def prepare_model_variables(genotypes, testing_variables, null_variables):
    """
    Converts genotype categories to class labels and imputed unknown genotypes.
    Imputation involves duplicating samples, so we create copies of the rows of
    the feature matrix to align.
    """

    n_samples = genotypes.shape[0]
    n_genotypes = genotypes.shape[1]
    N_GENOTYPES = 3
    if n_genotypes != N_GENOTYPES:
        raise Exception("Must be using genotype categories!")

    if n_samples != genotypes.shape[0]:
        raise Exception("Genotype and feature matrices do not have same number of samples!")

    N_COPIES = 3
    class_labels = np.zeros(N_COPIES * n_samples)
    testing_features = np.zeros((N_COPIES * n_samples,
                                 testing_variables.shape[1]))
    if null_variables is not None:
        null_features = np.zeros((N_COPIES * n_samples,
                                  null_variables.shape[1]))
    for i in xrange(n_samples):
        gt = None
        for j in xrange(N_GENOTYPES):
            if genotypes[i, j] == 1.0:
                gt = j
        
        for j in xrange(N_COPIES):
            idx = N_COPIES * i + j
            testing_features[idx, :] = testing_variables[i, :]
            if null_variables is not None:
                null_features[idx, :] = null_variables[i, :]

            if gt is None:
                class_labels[idx] = j
            else:
                class_labels[idx] = gt

    if null_variables is not None:
        testing_features = np.hstack([testing_features,
                                      null_features])
    else:
        null_features = None

    return N_COPIES, class_labels, testing_features, null_features

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
    if args.variables_fl:
        selected_samples, null_variables = parse_variables_file(args.variables_fl)
        sample_indices = dict([(sample_name, idx) for idx, sample_name in
                               enumerate(data_model.sample_labels)])
        
        selected_indices = []
        for sample_id in selected_samples:
            if sample_id not in sample_indices:
                raise Exception("Unknown sample id '%s', known ids are: %s" %
                                (sample_id, ",".join(data_model.sample_labels)))
            selected_indices.append(sample_indices[sample_id])

        # select subset and re-order
        genotypes = genotypes[selected_indices, :]
    
    with open(os.path.join(stats_dir, OUTPUT_FLNAME), "w") as fl:
        next_output = 1
        for i, pair in enumerate(data_model.snp_feature_map.iteritems()):
            pos_label, feature_idx = pair
            chrom, pos = pos_label

            pos_genotypes = genotypes[:, feature_idx]

            n_copies, class_labels, testing_features, null_features = prepare_model_variables(
                pos_genotypes,
                testing_variables,
                null_variables)

            # since we make multiple copies of the original samples,
            # we need to scale the log loss so that it is correct for
            # the original sample size

            p_value = likelihood_ratio_test(testing_features,
                                            class_labels,
                                            lr,
                                            features_null = null_features,
                                            g_scaling_factor = 1.0 / n_copies)
            
            if i == next_output:
                print i, "Position", pos_label, "has p-value", p_value
                next_output *= 2

            fl.write("\t".join([chrom, pos, str(p_value)]))
            fl.write("\n")
            
def parseargs():
    parser = argparse.ArgumentParser(description="Asaph - Single SNP Association Tests")

    parser.add_argument("--workdir", type=str, help="Work directory", required=True)

    parser.add_argument("--variables-fl",
                        type=str,
                        help="Specify labels and variables for the null model")

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    run_likelihood_ratio_tests(args)
