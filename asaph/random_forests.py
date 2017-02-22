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

from asaph.analysis import write_similarity_curves
from asaph.analysis import plot_sampled_snps_curves
from asaph.analysis import plot_similarity_curves
from asaph.analysis import sampled_snps_curves
from asaph.analysis import similarity_curves
from asaph.analysis import histogram_sparse_to_dense
from asaph.analysis import plot_feature_histogram
from asaph.ml import ConstrainedBaggingRandomForest
from asaph.newioutils import read_features
from asaph.newioutils import read_rf_snps
from asaph.newioutils import write_rf_snps

from asaph.newioutils import deserialize
from asaph.newioutils import PROJECT_SUMMARY_FLNAME
from asaph.newioutils import serialize

def write_interactions(basedir, n_trees, used_feature_sets):
    model_dir = os.path.join(basedir, "models", "rf", str(n_trees))
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    flname = os.path.join(model_dir, "interactions")
    serialize(flname, used_feature_sets)

def train_model(args):
    workdir = args["workdir"]

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)
    
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
                                        n_resamples,
                                        args["batch_size"])
    feature_importances, used_feature_counts, used_feature_sets = \
                                               rf.feature_importances(
                                                   features.feature_matrix,
                                                   features.class_labels,
                                                   statistics=args["statistics"],
                                                   interactions=args["interactions"])

    snp_importances = features.rank_snps(feature_importances)
    write_rf_snps(workdir, snp_importances, n_trees, "model1")
    
    if args["statistics"]:
        dense = histogram_sparse_to_dense(used_feature_counts)
        flname = os.path.join(figures_dir, "features_used_histogram_rf_%s_trees.png" \
                              % args["trees"])
        plot_feature_histogram(flname, dense)

    if args["interactions"]:
        write_interactions(workdir, n_trees, used_feature_sets)

    rf = ConstrainedBaggingRandomForest(n_trees,
                                        n_resamples,
                                        args["batch_size"])
    feature_importances, _, _ = rf.feature_importances(features.feature_matrix,
                                                    features.class_labels,
                                                    statistics=False)
    snp_importances = features.rank_snps(feature_importances)
    write_rf_snps(workdir, snp_importances, n_trees, "model2")

def analyze_rankings(args):
    workdir = args["workdir"]

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    analysis_dir = os.path.join(workdir, "analysis")
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    project_summary = deserialize(os.path.join(workdir, PROJECT_SUMMARY_FLNAME))
    all_snps = read_rf_snps(workdir)

    n_models, common_feature_percentages = similarity_curves(args["thresholds"],
                                                             all_snps,
                                                             project_summary.filtered_positions)

    analysis_flname = os.path.join(analysis_dir,
                                   "snp_ranking_overlaps_rf.tsv")
    write_similarity_curves(analysis_flname,
                            args["thresholds"],
                            n_models,
                            common_feature_percentages)

    flname_base = os.path.join(figures_dir, "snp_ranking_overlaps_rf")
    plot_similarity_curves(flname_base,
                           args["thresholds"],
                           n_models,
                           common_feature_percentages)

    n_models, common_feature_counts, snp1_feature_counts, \
        snp2_feature_counts = sampled_snps_curves(all_snps)
    flname_base = os.path.join(figures_dir, "snp_counts")
    plot_sampled_snps_curves(flname_base,
                             n_models,
                             common_feature_counts,
                             snp1_feature_counts,
                             snp2_feature_counts)

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
                              default=-1,
                              type=int,
                              help="Number of additional samples")

    train_parser.add_argument("--batch-size",
                              default=100,
                              type=int,
                              help="Number of trees to train in each batch. Trade off between speed and memory usage.")

    train_parser.add_argument("--statistics",
                              action="store_true",
                              help="Record statistics on forest")

    train_parser.add_argument("--interactions",
                              action="store_true",
                              help="Record feature interactions")
    
    analyze_parser = subparsers.add_parser("analyze-rankings",
                                           help="Analyze rankings")

    analyze_parser.add_argument("--thresholds",
                                type=float,
                                nargs="+",
                                default=[0.0001,0.001,0.01,0.1],
                                help="Thresholds for similarity curves")

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
