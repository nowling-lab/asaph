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
from sklearn.linear_model import SGDClassifier

from asaph.newioutils import read_features
from asaph.newioutils import serialize
from asaph.newioutils import deserialize

def write_snps(basedir, snps, model_id):
    model_dir = os.path.join(basedir, "models", "lr")
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    flname = os.path.join(model_dir, model_id)
    serialize(flname, snps)
    
def read_snps(basedir):
    model_dir = os.path.join(basedir, "models", "lr")
    if not os.path.exists(model_dir):
        return dict()

    models = defaultdict(list)
    model_flnames = glob.glob(os.path.join(model_dir, "*"))
    for flname in model_flnames:
        snps = deserialize(flname)
        model_id = os.path.basename(flname)
        models[model_id] = snps

    return models


def train(args):
    workdir = args.workdir

    if args.method == "sgd-l2":
        lr = SGDClassifier(loss="log",
                           penalty="l2")
    elif args.method  == "sgd-en":
        lr = SGDClassifier(loss="log",
                           penalty="elasticnet")
    else:
        raise Exception, "Unknown method '%s'" % args.method

    features = read_features(workdir)
    lr.fit(features.feature_matrix,
           features.class_labels)
    snp_importances = features.rank_snps(lr.coef_[0, :])
    write_snps(workdir, snp_importances, args.method)

def rankings(args):
    workdir = args.workdir

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)
    
    snp_models = read_snps(workdir)

    if args.method not in snp_models:
        raise Exception, "Model type '%s' hasn't been trained yet." % args.method

    model = snp_models[args.method]

    ranks_flname = os.path.join(workdir, "rankings_lr_%s.tsv" % args.method)
    with open(ranks_flname, "w") as fl:
        for i in xrange(len(model)):
            chrom, pos = model.labels[i]
            importance = model.importances[i]
            fl.write("%s\t%s\t%s\n" % (chrom, pos, importance))

    fig_flname = os.path.join(figures_dir, "lr_weights_%s.png" % args.method)
    plt.clf()
    plt.grid(True)
    plt.plot(model.importances, "m.-")
    plt.xlabel("SNPs (ordered)", fontsize=16)
    plt.ylabel("Weight", fontsize=16)
    plt.savefig(fig_flname, DPI=300)
    

def parseargs():
    parser = argparse.ArgumentParser(description="Asaph - Logistic Regression")

    parser.add_argument("--workdir",
                        type=str,
                        required=True,
                        help="Work directory")

    subparsers = parser.add_subparsers(dest="mode")
    train_parser = subparsers.add_parser("train",
                                         help="Train Logistic Regression model")
    train_parser.add_argument("--method",
                              choices=["sgd-l2", "sgd-en"],
                              required=True,
                              help="LR algorithm to use")
    
    output_parser = subparsers.add_parser("rankings",
                                          help="Output rankings and plots")
    output_parser.add_argument("--method",
                               choices=["sgd-l2", "sgd-en"],
                               required=True,
                               help="LR algorithm to use")

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    if args.mode == "train":
        train(args)
    elif args.mode == "rankings":
        rankings(args)
    else:
        print "Unknown mode '%s'" % args.mode
        sys.exit(1)
