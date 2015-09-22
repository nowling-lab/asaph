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

import sys

from ioutils import *
import numpy as np
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt

def plot_explained_var_ratios(outbase, pca):
    plt.plot(pca.explained_variance_ratio_, "b.-")
    plt.xlabel("Component", fontsize=16)
    plt.ylabel("Explained Variance Ratio", fontsize=16)
    plt.ylim([0.0, 1.0])
    plt.savefig(outbase + "_explained_var_ratio.png", DPI=200)

def plot_projection(outbase, projected, pc_idx1, pc_idx2, group_rows, leftover_rows):
    colors = ["r", "g", "b", "c", "m", "y"]
    plt.clf()
    
    for c, (group_name, rows) in enumerate(group_rows):
        plt.scatter(projected[rows, pc_idx1], projected[rows, pc_idx2], color=colors[c], \
                    label=group_name)

    if len(leftover_rows) > 0:
        plt.scatter(projected[leftover_rows, pc_idx1], projected[leftover_rows, pc_idx2], \
                    color="k", label="Other")

    plt.xlabel("Principal Component %s" % (pc_idx1 + 1), fontsize=16)
    plt.ylabel("Principal Component %s" % (pc_idx2 + 1), fontsize=16)
    plt.legend()
    plt.savefig(outbase + "_pca_proj_%s_%s.png" % \
                (pc_idx1 + 1, pc_idx2 + 1), DPI=200)

def round_even(n):
    return n - n % 2
    
if __name__ == "__main__":
    inbase = sys.argv[1]
    outbase = sys.argv[2]

    group_rows, leftover_rows = get_group_indices(inbase)
    
    feature_matrix = open_feature_matrix(inbase)
        
    pca = PCA(whiten=True)
    projected = pca.fit_transform(feature_matrix)

    plot_explained_var_ratios(outbase, pca)

    n_components = round_even(feature_matrix.shape[0])
    for pc_idx in xrange(0, n_components, 2):
        plot_projection(outbase, projected, pc_idx, pc_idx + 1, group_rows, \
                        leftover_rows)
    

