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
import os
import sys

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

import numpy as np

from asaph.ml import ConstrainedBaggingRandomForest
from asaph.newioutils import read_features
from asaph.newioutils import read_rf_snps
from asaph.newioutils import write_rf_snps
    
def train_model(args):
    workdir = args["workdir"]

    n_trees = args["trees"]
    if n_trees is None:
        print "Number of trees must be specified for training"
        sys.exit(1)

    n_resamples = args["resamples"]
    if n_resamples is None:
        print "Number of additional samples must be specified for training"
        sys.exit(1)

    features = read_features(workdir)

    rf = ConstrainedBaggingRandomForest(n_trees,
                                        n_resamples)
    feature_importances = rf.feature_importances(features.feature_matrix,
                                                 features.class_labels)
    snp_importances = features.rank_snps(feature_importances)
    write_rf_snps(workdir, snp_importances, n_trees, "model1")

    rf = ConstrainedBaggingRandomForest(n_trees,
                                        n_resamples)
    feature_importances = rf.feature_importances(features.feature_matrix,
                                                 features.class_labels)
    snp_importances = features.rank_snps(feature_importances)
    write_rf_snps(workdir, snp_importances, n_trees, "model2")

def analyze_rankings(args):
    workdir = args["workdir"]

    figures_dir = os.path.join(workdir, "figures")

    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    all_snps = read_rf_snps(workdir)
    ordered_trees = sorted(all_snps.keys())

    thresholds = [0.05, 0.1, 0.25, 0.5] 
    
    common_feature_counts = []
    snp1_feature_counts = []
    snp2_feature_counts = []
    common_feature_threshold_percentages = defaultdict(list)
    for n_trees in ordered_trees:
        snps1, snps2 = all_snps[n_trees]

        common_feature_counts.append(snps1.count_intersection(snps2))
        snp1_feature_counts.append(len(snps1))
        snp2_feature_counts.append(len(snps2))
            
        for threshold in thresholds:
            n = max(1, int(threshold * min(len(snps1), len(snps2))))
            percentage = 100.0 * float(snps1.take(n).count_intersection(snps2.take(n))) \
                         / float(n)
            common_feature_threshold_percentages[threshold].append(percentage)

    plt.clf()
    plt.hold(True)
    plt.grid(True)
    plt.semilogx(ordered_trees, common_feature_counts, "k.-", label="Common")
    plt.semilogx(ordered_trees, snp1_feature_counts, "c.-", label="Model 1")
    plt.semilogx(ordered_trees, snp2_feature_counts, "m.-", label="Model 2")
    plt.xlabel("Number of Trees", fontsize=16)
    plt.ylabel("SNPs (Count)", fontsize=16)
    plt.legend(loc="lower right")
    plt.ylim([0, max(max(common_feature_counts), max(snp1_feature_counts), max(snp2_feature_counts)) + 10])
    plt.xlim([min(ordered_trees), max(ordered_trees)])

    plt.savefig(os.path.join(figures_dir, "snp_counts.png"), DPI=200) 
    plt.savefig(os.path.join(figures_dir, "snp_counts.pdf"), DPI=200)

    plt.clf()
    plt.hold(True)
    plt.grid(True)
    colors = ["r.-", "g.-", "b.-", "m.-", "c.-"]
    for i, threshold in enumerate(thresholds):
        c = colors[i]
        label = str(int(100.0 * threshold))
        plt.semilogx(ordered_trees, common_feature_threshold_percentages[threshold],
                     c, label="Top %s%%" % label)
    plt.xlabel("Number of Trees", fontsize=16)
    plt.ylabel("Common SNPs (%)", fontsize=16)
    plt.legend(loc="lower right")
    plt.ylim([0, 100])
    plt.xlim([min(ordered_trees), max(ordered_trees)])

    plt.savefig(os.path.join(figures_dir, "common_snps.png"), DPI=200) 
    plt.savefig(os.path.join(figures_dir, "common_snps.pdf"), DPI=200)

def output_rankings(args):
    workdir = args["workdir"]

    n_trees = args["trees"]
    if n_trees is None:
        print "Number of trees must be specified for outputting ranks"
        sys.exit(1)

    ranks_flname = args["ranks_file"]
    if ranks_flname is None:
        print "Output filename must be specified"
        sys.exit(1)

    all_models = read_rf_snps(workdir)
    if n_trees not in all_models:
        print "No model with %s trees. The available models have %s trees" \
            % (n_trees, sorted(all_models.keys()))
        sys.exit(1)

    snps1, snps2 = all_models[n_trees]

    fl = open(ranks_flname, "w")
    for i in xrange(len(snps1)):
        chrom, pos = snps1.labels[i]
        importance = snps1.importances[i]
        fl.write("%s\t%s\t%s\n" % (chrom, pos, importance))

    fl.close()

def parseargs():
    parser = argparse.ArgumentParser(description="Asaph - Random Forests")

    parser.add_argument("--workdir", type=str, help="Work directory", required=True)

    subparsers = parser.add_subparsers(dest="mode")
    train_parser = subparsers.add_parser("train",
                                         help="Train RF model")
    train_parser.add_argument("--trees",
                              type=int,
                              help="Number of trees in Random Forest")
    train_parser.add_argument("--resamples",
                              type=int,
                              help="Number of additional samples")

    analyze_parser = subparsers.add_parser("analyze-rankings",
                                           help="Analyze rankings")

    output_parser = subparsers.add_parser("output-rankings",
                                          help="Output rankings")
    output_parser.add_argument("--trees",
                               type=int,
                               help="Number of trees in Random Forest")

    output_parser.add_argument("--ranks-file",
                               type=str,
                               help="Output file for SNP ranks")

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    if args.mode == "train":
        train_model(vars(args))
    elif args.mode == "analyze-rankings":
        analyze_rankings(vars(args))
    elif args.mode == "output-rankings":
        output_rankings(vars(args))
    else:
        print "Unknown mode '%s'" % args["mode"]
        sys.exit(1)
