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

from asaph.newioutils import read_features
from asaph.cramers_v import cramers_v

from asaph.newioutils import deserialize
from asaph.newioutils import PROJECT_SUMMARY_FLNAME
from asaph.newioutils import serialize

onehot_to_ints = np.array([1, 2, 3])

def pair_association(pair1, pair2, features, model_type):
    snp_label1, feature_idx1 = pair1
    snp_label2, feature_idx2 = pair2

    if model_type == "regular" or model_type == "permutation":
        # convert one-hot encoding back into integer encoding
        snp_features1 = np.dot(features.feature_matrix[:, feature_idx1],
                               onehot_to_ints)
        snp_features2 = np.dot(features.feature_matrix[:, feature_idx2],
                               onehot_to_ints)
        if model_type == "permutation":
            np.random.shuffle(snp_features1)
            np.random.shuffle(snp_features2)
    else:
        n_samples = features.feature_matrix.shape[0]
        snp_features1 = np.random.randint(1, 4, n_samples)
        snp_features2 = np.random.randint(1, 4, n_samples)
    
    v = cramers_v(snp_features1,
                  snp_features2)

    return snp_label1, snp_label2, v

def sample_pairs(n_samples, items):
    sampled = set()
    while len(sampled) < n_samples:
        key1, value1 = random.choice(items)
        key2, value2 = random.choice(items)

        if key1 == key2:
            continue
        
        pair = frozenset([key1, key2])
        if pair in sampled:
            continue

        sampled.add(pair)
        yield (key1, value1), (key2, value2)

def snp_pair_associations(workdir, n_samples, model_type):
    if not os.path.exists(workdir):
        raise Exception, "workdir '%s' does not exist." % workdir

    stats_dir = os.path.join(workdir, "statistics")
    if not os.path.exists(stats_dir):
        os.makedirs(stats_dir)

    project_summary = deserialize(os.path.join(workdir, PROJECT_SUMMARY_FLNAME))
    if project_summary.feature_encoding != "categories":
        print "Cramer's V only works with the 'categories' feature encoding."
        sys.exit(1)
    
    features = read_features(workdir)
    print features.feature_matrix.shape

    next_output = 1
    if n_samples is not None:
        pairs = sample_pairs(n_samples, features.snp_feature_map.items())
    else:
        pairs = itertools.combinations(features.snp_feature_map.iteritems(), 2)

    flname = "snp_pairwise_associations.txt"
    if model_type == "permutation":
        flname = "snp_pairwise_associations_permutation.txt"
    elif model_type == "uniform-random":
        flname = "snp_pairwise_associations_uniform-random.txt"
        
    with open(os.path.join(stats_dir, flname), "w") as fl:
        for i, (pair1, pair2) in enumerate(pairs):
            snp_label1, snp_label2, v = pair_association(pair1,
                                                         pair2,
                                                         features,
                                                         model_type)
            
            if i == next_output:
                print "Pair", i, "has an association of", v
                next_output *= 2

            chrom1, pos1 = snp_label1
            chrom2, pos2 = snp_label2
            fl.write("\t".join([chrom1, pos1, chrom2, pos2, str(v)]))
            fl.write("\n")
        
def parseargs():
    parser = argparse.ArgumentParser(description="Asaph - Cramer's V")

    parser.add_argument("--workdir", type=str, help="Work directory", required=True)

    parser.add_argument("--samples",
                        type=int,
                        help="Number of samples")

    parser.add_argument("--model-type",
                        choices = ["regular",
                                   "uniform-random",
                                   "permutation"],
                        default = "regular")

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    snp_pair_associations(args.workdir,
                          args.samples,
                          args.model_type)
