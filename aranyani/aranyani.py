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

import argparse
from collections import defaultdict
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata

from sklearn.ensemble import RandomForestClassifier

from utils import *

import itertools

def read_features(basename):
    class_labels, selected_indices, selected_individuals = select_groups(basename, True)
    feature_labels = read_features(basename, True)
    feature_matrix = open_feature_matrix(basename, True)

    features = Features(feature_matrix, feature_labels, class_labels)

    return features

def plot_errors(basename, tree_sizes, common_features):
    if not os.path.exists(basename + os.sep + "figures"):
        os.makedirs(basename + os.sep + "figures")
    
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

    plt.savefig(basename + os.sep + "figures" + os.sep + \
                "ranking_convergence_analysis.pdf", DPI=200)

    plt.savefig(basename + os.sep + "figures" + os.sep + \
                "ranking_convergence_analysis.eps", DPI=200)
    
def tree_sweep(args):
    if args["tree_list"] is None:
        print "Need to specify list of tree sizes for tree-sweep mode."
        sys.exit(1)

    if args["snp_list"] is None:
        print "Need to specify list of SNP counts for tree-sweep mode."
        sys.exit(1)

    features = read_features(args["data"])
    
    print features.feature_matrix.shape

    common_features = defaultdict(list)
    for n_trees in args["tree_list"]:
        print "Training forest of", n_trees, "trees"
        rf1 = features.train_rf(n_trees)
        rf2 = features.train_rf(n_trees)

        for threshold in args["snp_list"]:
            snps1 = features.snp_importances(rf1).rank()
            snps2 = features.snp_importances(rf2).rank()

            common_indices = set(indices1).intersection(indices2)
            normalized_common = 100.0 * len(common_indices) / threshold
            common_features[threshold].append(normalized_common)

    plot_errors(args["data"], args["tree_list"], common_features)

    
def select_snps(args):
    if args["trees"] is None:
        print "Need to specify number of trees for selecting SNPs."
        sys.exit(1)

    if args["snps"] is None:
        print "Need to specify number of SNPs to select."
        sys.exit(1)

    labels, selected_indices, selected_individuals = select_groups(args["data"], True)
    feature_labels = read_features(args["data"], True)
    feature_matrix = open_feature_matrix(args["data"], True)

    rf = RandomForestClassifier(n_estimators=args["trees"])
    rf.fit(feature_matrix, labels)

    snp_importances = calc_snp_importances(feature_labels, rf.feature_importances_)
    ranked_snp_importances = rank_features(args["snps"], snp_importances)

    print ranked_snp_importances

    #for i, label in enumerate(feature_labels):
    #    if position in snp_indices:
    #        print 
    

def parseargs():
    parser = argparse.ArgumentParser(description="Pipeline for SNP analysis")

    parser.add_argument("--mode", required=True,
                        choices=["tree-sweep", "select-snps"],
                        help="Run mode")

    parser.add_argument("--data", required=True,
                        help="Data Directory")

    tree_group = parser.add_mutually_exclusive_group()
    tree_group.add_argument("--trees", type=int,
                        help="Size of Random Forest to use in selecting or validating SNPs.")
    tree_group.add_argument("--tree-list", type=int, nargs="+",
                        help="Sizes of Random Forests for tree parameter sweep.")

    snp_group = parser.add_mutually_exclusive_group()
    snp_group.add_argument("--snps", type=int)
    snp_group.add_argument("--snp-list", type=int, nargs="+")
    
    return vars(parser.parse_args())

if __name__ == "__main__":
    args = parseargs()

    if args["mode"] == "tree-sweep":
        tree_sweep(args)
    elif args["mode"] == "select-snps":
        select_snps(args)
    else:
        print "Invalid mode", args["mode"]
        sys.exit(1)
