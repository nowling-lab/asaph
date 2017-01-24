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

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

import numpy as np

from asaph.newioutils import read_features
from asaph.newioutils import read_snps
from asaph.newioutils import write_snps
from asaph.vcf import convert as convert_vcf
from asaph.fasta import convert as convert_fasta

def import_fasta(args):
    fasta_flname = args["fasta"]
    if fasta_flname is None:
        print "FASTA file must be specified for import"
        sys.exit(1)

    groups_flname = args["groups"]
    if groups_flname is None:
        print "Groups file must be specified for import"
        sys.exit(1)

    seq_type = args["seq"]
    if seq_type is None:
        print "Sequence type must be specified for import"
        sys.exit(1)

    workdir = args["workdir"]
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    convert_fasta(groups_flname, fasta_flname, seq_type, workdir)

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

    convert_vcf(groups_flname, vcf_flname, workdir, args["compress"])

def pca(args):
    workdir = args["workdir"]

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    n_pcs = args["n_pcs"]
    if n_pcs is None:
        print "Number of pcs must be specified for PCA"
        sys.exit(1)
    
    features = read_features(workdir)
    proj, explained_variance_ratios = features.svd(n_pcs)

    plt.clf()
    plt.plot(explained_variance_ratios, "b.-")
    plt.xlabel("PC", fontsize=16)
    plt.ylabel("Explained Variance Ratio", fontsize=16)
    plt.savefig(os.path.join(figures_dir, "pca_explained_variance_ratios.png"), DPI=200)
    
    for i in xrange(n_pcs - 1):
        for j in xrange(i + 1, n_pcs):
            plt.clf()
            plt.scatter(proj[:, i], proj[:, j])
            plt.xlabel("PC " + str(i), fontsize=16)
            plt.ylabel("PC " + str(j), fontsize=16)
            plt.savefig(os.path.join(figures_dir, "pca_%s_%s.png" % (i, j)), DPI=200)
    
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

    snp_importances1 = features.snp_importances(n_trees, n_resamples).rank()
    snp_importances2 = features.snp_importances(n_trees, n_resamples).rank()

    write_snps(workdir, snp_importances1, "model1")
    write_snps(workdir, snp_importances2, "model2")

def analyze_rankings(args):
    workdir = args["workdir"]

    figures_dir = os.path.join(workdir, "figures")

    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    all_snps = read_snps(workdir)
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

    all_models = read_snps(workdir)
    if n_trees not in all_models:
        print "No model with %s trees. The available models have %s trees" \
            % (n_trees, sorted(all_models.keys()))
        sys.exit(1)

    snps1, snps2 = all_models[n_trees]

    fl = open(ranks_flname, "w")
    for i in xrange(len(snps1)):
        chrom, pos = snps1.labels[i]
        importance = snps1.importances[i]
        fd = snps1.fixed_differences[i]
        missing = snps1.missing_data[i]
        fl.write("%s\t%s\t%s\t%s\t%s\n" % (chrom, pos, importance, fd, missing))

    fl.close()

def parseargs():
    parser = argparse.ArgumentParser(description="Aranyani")

    parser.add_argument("--mode", required=True,
                        choices=["import",
                                 "train",
                                 "analyze-rankings",
                                 "output-rankings",
                                 "pca"],
                        help="Operating mode")

    parser.add_argument("--seq", choices=["DNA", "AA"],
                        help="Sequence type for FASTA input")

    parser.add_argument("--compress", action="store_true")
    parser.add_argument("--vcf", type=str, help="VCF file to import")
    parser.add_argument("--fasta", type=str, help="FASTA file to import")
    parser.add_argument("--groups", type=str, help="Groups file to import")
    parser.add_argument("--workdir", type=str, help="Work directory", required=True)

    parser.add_argument("--trees", type=int, help="Number of trees in Random Forest")
    parser.add_argument("--resamples", type=int, help="Number of additional samples")
    parser.add_argument("--n-pcs", type=int, help="Number of PCs to use in PCA")

    parser.add_argument("--ranks-file", type=str,
                        help="Output file for SNP ranks")

    return vars(parser.parse_args())

if __name__ == "__main__":
    args = parseargs()

    if args["mode"] == "import":
        if args["vcf"] is not None:
            import_vcf(args)
        elif args["fasta"] is not None:
            import_fasta(args)
        else:
            print "Need to specify VCF or FASTA file for import"
            sys.exit(1)
    elif args["mode"] == "train":
        train_model(args)
    elif args["mode"] == "analyze-rankings":
        analyze_rankings(args)
    elif args["mode"] == "output-rankings":
        output_rankings(args)
    elif args["mode"] == "pca":
        pca(args)
    else:
        print "Unknown mode '%s'" % args["mode"]
        sys.exit(1)
