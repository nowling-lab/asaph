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

from aranyani.newioutils import read_features
from aranyani.newioutils import read_snps
from aranyani.newioutils import write_snps
from aranyani.vcf import convert

def import_vcf(args):
    vcf_flname = args["vcf"]
    if vcf_flname is None:
        print "VCF file must be specified for import"
        sys.exit(1)

    groups_flname = args["groups"]
    if groups_flname is None:
        print "Groups file must be specified for import"
        sys.exit(1)

    workdir = args["workdir"]
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    convert(groups_flname, vcf_flname, workdir)

def train_model(args):
    workdir = args["workdir"]

    n_trees = args["trees"]
    if n_trees is None:
        print "Number of trees must be specified for training"
        sys.exit(1)

    features = read_features(workdir)

    rf1 = features.train_rf(n_trees)
    rf2 = features.train_rf(n_trees)

    snp_importances1 = features.snp_importances(rf1).rank()
    snp_importances2 = features.snp_importances(rf2).rank()

    write_snps(workdir, snp_importances1, "model1")
    write_snps(workdir, snp_importances2, "model2")

def analyze_rankings(args):
    workdir = args["workdir"]

    figures_dir = os.path.join(workdir, "figures")

    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    all_snps = read_snps(workdir)
    ordered_trees = sorted(all_snps.keys())

    thresholds = [0.01, 0.05, 0.1, 0.25, 0.5] 
    
    common_feature_counts = []
    snp1_feature_counts = []
    snp2_feature_counts = []
    common_feature_threshold_percentages = defaultdict(list)
    for n_trees in ordered_trees:
        snps1, snps2 = all_snps[n_trees]

        common_feature_counts.append(snps1.count_intersection(snps2))
        snp1_feature_counts.append(len(snps1))
        snp2_feature_counts.append(len(snps2))
            
        #print n_trees, len(snps1), len(snps2), snps1.count_intersection(snps2)

        for threshold in thresholds:
            n = max(1, int(threshold * min(len(snps1), len(snps2))))
            percentage = 100.0 * float(snps1.take(n).count_intersection(snps2.take(n))) \
                         / float(n)
            common_feature_threshold_percentages[threshold].append(percentage)

            print n_trees, n, threshold, percentage

    plt.clf()
    plt.hold(True)
    plt.grid(True)
    plt.semilogx(ordered_trees, common_feature_counts, "k.-", label="Common")
    plt.semilogx(ordered_trees, snp1_feature_counts, "c.-", label="Model 1")
    plt.semilogx(ordered_trees, snp2_feature_counts, "m.-", label="Model 2")
    plt.xlabel("Number of Trees", fontsize=16)
    plt.ylabel("SNPs (Count)", fontsize=16)
    plt.legend(loc="upper left")
    plt.ylim([0, max(common_feature_counts) + 10])
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
    plt.ylabel("Common SNPs (Percentage)", fontsize=16)
    plt.legend(loc="upper left")
    plt.ylim([0, 100])
    plt.xlim([min(ordered_trees), max(ordered_trees)])

    plt.savefig(os.path.join(figures_dir, "common_snps.png"), DPI=200) 
    plt.savefig(os.path.join(figures_dir, "common_snps.pdf"), DPI=200)
    

def parseargs():
    parser = argparse.ArgumentParser(description="Aranyani")

    parser.add_argument("--mode", required=True,
                        choices=["import",
                                 "train",
                                 "analyze-rankings"],
                        help="Operating mode")

    parser.add_argument("--vcf", type=str, help="VCF file to import")
    parser.add_argument("--groups", type=str, help="Groups file to import")
    parser.add_argument("--workdir", type=str, help="Work directory", required=True)

    parser.add_argument("--trees", type=int, help="Number of trees in Random Forest")

    parser.add_argument("--thresholds", type=int, nargs="+",
                        help="Number of SNPs to use in ranking convergence analysis")

    return vars(parser.parse_args())

if __name__ == "__main__":
    args = parseargs()

    if args["mode"] == "import":
        import_vcf(args)
    elif args["mode"] == "train":
        train_model(args)
    elif args["mode"] == "analyze-rankings":
        analyze_rankings(args)
    else:
        print "Unknown mode '%s'" % args["mode"]
        sys.exit(1)
