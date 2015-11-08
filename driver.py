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

def parseargs():
    parser = argparse.ArgumentParser(description="Aranyani")

    parser.add_argument("--mode", required=True,
                        choices=["import",
                                 "train"],
                        help="Operating mode")

    parser.add_argument("--vcf", type=str, help="VCF file to import")
    parser.add_argument("--groups", type=str, help="Groups file to import")
    parser.add_argument("--workdir", type=str, help="Work directory", required=True)

    parser.add_argument("--trees", type=int, help="Number of trees in Random Forest")

    return vars(parser.parse_args())

if __name__ == "__main__":
    args = parseargs()

    if args["mode"] == "import":
        import_vcf(args)
    elif args["mode"] == "train":
        train_model(args)
    else:
        print "Unknown mode '%s'" % args["mode"]
        sys.exit(1)
