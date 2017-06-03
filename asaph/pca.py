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
import os
import sys

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

from asaph.ml import mcfadden_r2
from asaph.ml import PCA
from asaph.newioutils import read_features

def explained_variance_analysis(args):
    workdir = args.workdir

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    features = read_features(workdir)

    pca = PCA(args.n_components)
    explained_variance_ratios = pca.explained_variance(features.feature_matrix)

    fig_flname = os.path.join(figures_dir,
                              "pca_explained_variance_ratios.png")

    plt.clf()
    plt.grid(True)
    plt.plot(explained_variance_ratios, "m.-")
    plt.xlabel("Principal Component", fontsize=16)
    plt.ylabel("Explained Variance Ratio", fontsize=16)
    plt.ylim([0., 1.])
    plt.savefig(fig_flname,
                DPI=300)


def coefficient_of_determination(args):
    workdir = args.workdir

    features = read_features(workdir)

    pca = PCA(args.n_components)
    projected = pca.transform(features.feature_matrix)
    r2 = mcfadden_r2(projected,
                     features.class_labels)

    print "McFadden r**2", r2

def parseargs():
    parser = argparse.ArgumentParser(description="Asaph - PCA")

    parser.add_argument("--workdir",
                        type=str,
                        required=True,
                        help="Work directory")

    subparsers = parser.add_subparsers(dest="mode")
    eva_parser = subparsers.add_parser("explained-variance-analysis",
                                       help="Compute explained variances of PCs")
    eva_parser.add_argument("--n-components",
                            type=int,
                            required=True,
                            help="Number of PCs to compute")

    r2_parser = subparsers.add_parser("mcfadden-r2",
                                      help="Compute McFadden's R2 of PC-projected features under a LR model")
    
    r2_parser.add_argument("--n-components",
                            type=int,
                            required=True,
                            help="Number of PCs to compute")

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    if args.mode == "explained-variance-analysis":
        explained_variance_analysis(args)
    elif args.mode == "mcfadden-r2":
        coefficient_of_determination(args)
    else:
        print "Unknown mode '%s'" % args.mode
        sys.exit(1)
