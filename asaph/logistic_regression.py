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
import glob
import os
import sys

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

import numpy as np

from asaph.analysis import write_similarity_curves
from asaph.analysis import plot_similarity_curves
from asaph.analysis import similarity_curves
from asaph.models import ProjectSummary
from asaph.ml import LogisticRegressionEnsemble
from asaph.newioutils import read_features
from asaph.newioutils import serialize
from asaph.newioutils import deserialize
from asaph.newioutils import PROJECT_SUMMARY_FLNAME

def write_snps(basedir, snps, method, n_models, model_id):
    model_dir = os.path.join(basedir, "models", "lr-" + method, str(n_models))
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    flname = os.path.join(model_dir, model_id)
    serialize(flname, snps)

def read_snps(basedir, method):
    model_dir = os.path.join(basedir, "models", "lr-" + method)
    if not os.path.exists(model_dir):
        return dict()

    model_dirs = glob.glob(os.path.join(model_dir, "*"))

    models = defaultdict(list)
    for model_dir in model_dirs:
        model_flnames = glob.glob(os.path.join(model_dir, "*"))
        for flname in model_flnames:
            snps = deserialize(flname)
            n_models = int(os.path.basename(os.path.dirname(flname)))
            models[n_models].append(snps)

    return models

def analyze_rankings(args):
    workdir = args.workdir

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    analysis_dir = os.path.join(workdir, "analysis")
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    project_summary = deserialize(os.path.join(workdir, PROJECT_SUMMARY_FLNAME))
    all_snps = read_snps(workdir, args.method)

    n_models, common_feature_percentages = similarity_curves(args.thresholds,
                                                             all_snps,
                                                             project_summary.filtered_positions)

    analysis_flname = os.path.join(analysis_dir,
                                   "snp_ranking_overlaps_" + args.method + ".tsv")
    write_similarity_curves(analysis_flname,
                            args.thresholds,
                            n_models,
                            common_feature_percentages)

    flname_base = os.path.join(figures_dir, "snp_ranking_overlaps_" + args.method)
    plot_similarity_curves(flname_base,
                           args.thresholds,
                           n_models,
                           common_feature_percentages)

def train(args):
    workdir = args.workdir

    features = read_features(workdir)

    print "Training Ensemble 1"
    lr1 = LogisticRegressionEnsemble(args.n_models,
                                     args.method,
                                     args.batch_size,
                                     bagging=args.bagging)
    feature_importances = lr1.feature_importances(features.feature_matrix,
                                                  features.class_labels)
    snp_importances = features.rank_snps(feature_importances)
    write_snps(workdir, snp_importances, args.method, args.n_models, "1")

    print "Training ensemble 2"
    lr2 = LogisticRegressionEnsemble(args.n_models,
                                     args.method,
                                     args.batch_size,
                                     bagging=args.bagging)
    feature_importances = lr2.feature_importances(features.feature_matrix,
                                                  features.class_labels)
    snp_importances = features.rank_snps(feature_importances)
    write_snps(workdir, snp_importances, args.method, args.n_models, "2")

def output_rankings(args):
    workdir = args.workdir

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    rankings_dir = os.path.join(workdir, "rankings")
    if not os.path.exists(rankings_dir):
        os.makedirs(rankings_dir)

    snp_models = read_snps(workdir, args.method)

    if len(snp_models) == 0:
        raise Exception, "Model type '%s' hasn't been trained yet." % args.method

    if args.n_models not in snp_models:
        raise Exception, "No ensemble with '%s' models has been trained yet." \
            % args.n_models

    model = snp_models[args.n_models][0]

    ranks_flname = os.path.join(rankings_dir,
                                "rankings_lr_%s_%s.tsv" \
                                % (args.method, args.n_models))

    with open(ranks_flname, "w") as fl:
        for i in xrange(len(model)):
            chrom, pos = model.labels[i]
            importance = model.importances[i]
            fl.write("%s\t%s\t%s\n" % (chrom, pos, importance))

    fig_flname = os.path.join(figures_dir, "lr_weights_%s_%s.png" % \
                              (args.method, args.n_models))
    plt.clf()
    plt.grid(True)
    plt.plot(model.importances, "m.-")
    plt.xlabel("SNPs (ordered)", fontsize=16)
    plt.ylabel("Weight", fontsize=16)
    plt.savefig(fig_flname, DPI=300)


def parseargs():
    method_choices = ["sgd-l2", "sgd-en", "sag-l2", "asgd-l2"]
    default_method = "sgd-l2"

    parser = argparse.ArgumentParser(description="Asaph - Logistic Regression")

    parser.add_argument("--workdir",
                        type=str,
                        required=True,
                        help="Work directory")

    subparsers = parser.add_subparsers(dest="mode")
    train_parser = subparsers.add_parser("train",
                                         help="Train Logistic Regression model")
    train_parser.add_argument("--method",
                              choices=method_choices,
                              default=default_method,
                              help="LR algorithm to use")

    train_parser.add_argument("--bagging",
                              action="store_true",
                              help="Enable bagging")

    train_parser.add_argument("--n-models",
                              type=int,
                              required=True,
                              help="Number of models to train in LR ensemble")

    train_parser.add_argument("--batch-size",
                              type=int,
                              default=100,
                              help="Number of models to train in each batch.")

    analyze_parser = subparsers.add_parser("analyze-rankings",
                                           help="Analyze rankings")

    analyze_parser.add_argument("--method",
                              choices=method_choices,
                              default=default_method,
                              help="LR algorithm to use")

    analyze_parser.add_argument("--thresholds",
                                type=float,
                                nargs="+",
                                default=[0.0001,0.001,0.01,0.1],
                                help="Thresholds for similarity curves")

    output_parser = subparsers.add_parser("output-rankings",
                                          help="Output rankings and plots")
    output_parser.add_argument("--method",
                               choices=method_choices,
                               default=default_method,
                               help="LR algorithm to use")

    output_parser.add_argument("--n-models",
                              type=int,
                              required=True,
                              help="Pick LR ensemble with this many models.")

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    if args.mode == "train":
        train(args)
    elif args.mode == "analyze-rankings":
        analyze_rankings(args)
    elif args.mode == "output-rankings":
        output_rankings(args)
    else:
        print "Unknown mode '%s'" % args.mode
        sys.exit(1)
