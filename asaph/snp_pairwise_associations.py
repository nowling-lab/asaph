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
import os
import sys

import numpy as np

from asaph.newioutils import read_features
from asaph.cramers_v import cramers_v

from asaph.newioutils import deserialize
from asaph.newioutils import PROJECT_SUMMARY_FLNAME
from asaph.newioutils import serialize

onehot_to_ints = np.array([1, 2, 3])

def cramers_v_pair(pair1, pair2, features):
    snp_label1, feature_idx1 = pair1
    snp_label2, feature_idx2 = pair2
    # convert one-hot encoding back into integer encoding
    snp_features1 = np.dot(features.feature_matrix[:, feature_idx1],
                           onehot_to_ints)
    snp_features2 = np.dot(features.feature_matrix[:, feature_idx2],
                           onehot_to_ints)
    
    v = cramers_v(snp_features1,
                  snp_features2)

    return snp_label1, snp_label2, v

def snp_pair_associations(workdir):
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
    snp_vs = []
    pairs = itertools.combinations(features.snp_feature_map.iteritems(), 2)
    with open(os.path.join(stats_dir, "snp_pairwise_associations.txt"), "w") as fl:
        for i, (pair1, pair2) in enumerate(pairs):
            snp_label1, snp_label2, v = cramers_v_pair(pair1, pair2, features)
            
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

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    snp_pair_associations(args.workdir)
