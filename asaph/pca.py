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
from itertools import tee, izip
import os
import sys

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import numpy as np

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
    selected = projected[:, args.model_components]
    r2 = mcfadden_r2(selected,
                     features.class_labels)

    print "McFadden r**2", r2

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def plot_projections(args):
    if len(args.pairs) != 2:
        print "Error: PCs must be provided in pairs of 2"
        sys.exit(1)

    workdir = args.workdir

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)
    
    features = read_features(workdir)

    pca = PCA(args.n_components)
    projected = pca.transform(features.feature_matrix)

    labels = np.array(features.class_labels)
    pop1 = labels == 0.
    pop2 = labels == 1.
    pop1_name = "Pop 1"
    pop2_name = "Pop 2"

    for p1, p2 in pairwise(args.pairs):
        fig_flname = os.path.join(figures_dir,
                                  "pca_projection_%s_%s.png" % (str(p1), str(p2)))
        plt.clf()
        plt.grid(True)
        plt.scatter(projected[pop1, p1],
                    projected[pop1, p2],
                    color="m",
                    edgecolor="k",
                    alpha=0.7,
                    label=pop1_name)
        plt.scatter(projected[pop2, p1],
                    projected[pop2, p2],
                    color="c",
                    edgecolor="k",
                    alpha=0.7,
                    label=pop2_name)
        plt.xlabel("Principal Component %s" % p1, fontsize=16)
        plt.ylabel("Principal Component %s" % p2, fontsize=16)
        plt.legend()
        plt.savefig(fig_flname,
                    DPI=300)

def output_coordinates(args):
    workdir = args.workdir

    features = read_features(workdir)

    pca = PCA(args.n_components)
    projected = pca.transform(features.feature_matrix)
    selected = projected[:, args.selected_components]
    
    labels = np.array(features.class_labels)
    pop1 = labels == 0.
    pop2 = labels == 1.
    pop1_name = "Pop 1"
    pop2_name = "Pop 2"

    with open(args.output_fl, "w") as fl:
        headers = ["sample", "population"]
        headers.extend(map(str, args.selected_components))
        fl.write("\t".join(headers))
        fl.write("\n")

        for i in xrange(len(features.sample_labels)):
            sample = features.sample_labels[i]
            pop = str(features.class_labels[i])
            line = [sample, pop]
            line.extend(map(str, selected[i, :]))
            fl.write("\t".join(line))
            fl.write("\n")
    
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

    r2_parser.add_argument("--model-components",
                           type=int,
                           nargs="+",
                           help="Components to use in model")

    output_parser = subparsers.add_parser("output-coordinates",
                                      help="Output PC projected coordinates")
    
    output_parser.add_argument("--n-components",
                               type=int,
                               required=True,
                               help="Number of PCs to compute")
    
    output_parser.add_argument("--selected-components",
                               type=int,
                               nargs="+",
                               help="Components to output")

    output_parser.add_argument("--output-fl",
                               type=str,
                               required=True,
                               help="Output file")
    
    plot_parser = subparsers.add_parser("plot-projections",
                                        help="Plot samples on principal coordinates")

    plot_parser.add_argument("--n-components",
                             type=int,
                             required=True,
                             help="Number of PCs to compute")

    plot_parser.add_argument("--pairs",
                             type=int,
                             nargs="+",
                             required=True,
                             help="Pairs of PCs to plot")

    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    if args.mode == "explained-variance-analysis":
        explained_variance_analysis(args)
    elif args.mode == "mcfadden-r2":
        coefficient_of_determination(args)
    elif args.mode == "plot-projections":
        plot_projections(args)
    elif args.mode == "output-coordinates":
        output_coordinates(args)
    else:
        print "Unknown mode '%s'" % args.mode
        sys.exit(1)
