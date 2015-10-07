"""
Copyright 2015 Ronald J. Nowling

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

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rankdata
import sys
import os
from collections import defaultdict

from ioutils import *

def rank_features(cutoff, feature_importances1, feature_importances2):
    sorted_indices1 = np.argsort(-1.0 * feature_importances1)
    sorted_indices2 = np.argsort(-1.0 * feature_importances2)

    sorted_importances1 = feature_importances1[sorted_indices1]
    sorted_importances2 = feature_importances2[sorted_indices2]

    nonzero_indices1 = sorted_indices1[sorted_importances1 != 0.0]
    nonzero_indices2 = sorted_indices2[sorted_importances2 != 0.0]

    nonzero_importances1 = sorted_importances1[sorted_importances1 != 0.0]
    nonzero_importances2 = sorted_importances2[sorted_importances2 != 0.0]

    return nonzero_indices1[:cutoff], nonzero_importances1[:cutoff], \
        nonzero_indices2[:cutoff], nonzero_importances2[:cutoff]

def plot_errors(basename, tree_sizes, common_features):
    if not os.path.exists(basename + ".figures"):
        os.makedirs(basename + ".figures")
    
    plt.clf()
    plt.hold(True)
    plt.grid(True)
    colors = ["r.-", "g.-", "b.-", "k.-"]
    for i, threshold in enumerate(sorted(common_features.keys())):
        normalized_common = common_features[threshold]
        color = colors[i]
        plt.semilogx(tree_sizes, normalized_common, color, \
            label="Top %s SNPs" % threshold)
    
    plt.xlabel("Number of Trees", fontsize=16)
    plt.ylabel("Common SNPs (%)", fontsize=16)
    plt.title("Ranking Convergence Analysis", fontsize=16)
    plt.legend(loc="lower right")
    plt.ylim([0, 100])
    plt.xlim([min(tree_sizes), max(tree_sizes)])

    plt.savefig(basename + ".figures/ranking_convergence_analysis.pdf", DPI=200)

if __name__ == "__main__":
    basename = sys.argv[1]
    top_n_features = [int(i) for i in sys.argv[2].split(",")]

    model_tree_counts = list_model_trees(basename)

    common_features = defaultdict(list)
    for n_trees in model_tree_counts:
        feature_importances1, feature_importances2 = \
                        deserialize_feature_importances(basename, n_trees)

        for threshold in top_n_features:
            indices1, importances1, indices2, importances2 = \
                    rank_features(threshold, feature_importances1, \
                                  feature_importances2)

            common_indices = set(indices1).intersection(indices2)
            normalized_common = 100.0 * len(common_indices) / threshold
            common_features[threshold].append(normalized_common)

    plot_errors(basename, model_tree_counts, common_features)
    
    
