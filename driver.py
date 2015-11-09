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
import os
import sys

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
    
    snp_thresholds = args["thresholds"]
    if snp_thresholds is None:
        print "Need to specify at least one threshold for ranking convergence analysis."
        sys.exit(1)

    plot_fl = args["plot_file"]
    if plot_fl is None:
        print "Need to specify plot filename."
        sys.exit(1)
    
    all_snps = read_snps(workdir)
    ordered_trees = sorted(all_snps.keys())

    common_features_percentages = defaultdict(list)
    for threshold in snp_thresholds:
        for n_trees in ordered_trees:
            snps = all_snps[n_trees]
            snps1 = snps[0].take(threshold)
            snps2 = snps[1].take(threshold)

            percentage = 100.0 * float(snps1.count_intersection(snps2)) / float(threshold)
            common_feature_percentages[threshold].append(percentage)
            
            print threshold, n_trees, percentage

    
    plt.clf()
    plt.hold(True)
    plt.grid(True)
    colors = ["r.-", "g.-", "b.-", "k.-", "c.-", "m.-"]
    for i, threshold in enumerate(sorted(common_features_percentages.keys())):
        normalized_common = common_features_percentages[threshold]
        color = colors[i]
        plt.semilogx(ordered_trees, normalized_common, color, \
            label="Top %s SNPs" % threshold)
    
    plt.xlabel("Number of Trees", fontsize=16)
    plt.ylabel("Common SNPs (%)", fontsize=16)
    plt.title("SNP Ranking Convergence Analysis", fontsize=16)
    plt.legend(loc="lower right")
    plt.ylim([0, 100])
    plt.xlim([min(ordered_trees), max(ordered_trees)])

    plt.savefig(plot_fl, DPI=200) 

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

    parser.add_argument("--plot-file", type=str,
                        help="File to write plot to. Extensions supported include .PNG, .PDF, and .EPS")

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
