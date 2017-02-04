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
from asaph.vcf import convert as convert_vcf
from asaph.fasta import convert as convert_fasta

def import_fasta(args):
    fasta_flname = args["fasta"]
    if fasta_flname is None:
        print "FASTA file must be specified for import"
        sys.exit(1)

    pops_flname = args["populations"]
    if pops_flname is None:
        print "Populations file must be specified for import"
        sys.exit(1)

    seq_type = args["seq"]
    if seq_type is None:
        print "Sequence type must be specified for import"
        sys.exit(1)

    workdir = args["workdir"]
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    convert_fasta(pops_flname, fasta_flname, seq_type, workdir)

def import_vcf(args):
    vcf_flname = args["vcf"]
    if vcf_flname is None:
        print "VCF file must be specified for import"
        sys.exit(1)

    pops_flname = args["populations"]
    if pops_flname is None:
        print "Populations file must be specified for import"
        sys.exit(1)

    workdir = args["workdir"]
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    convert_vcf(pops_flname, vcf_flname, workdir, args["compress"], args["feature_type"])


def parseargs():
    parser = argparse.ArgumentParser(description="Asaph")

    parser.add_argument("--seq",
                        choices=["DNA", "AA"],
                        default="AA",
                        help="Sequence type for FASTA input")

    parser.add_argument("--compress", action="store_true")

    parser.add_argument("--feature-type",
                        choices=["categories", "counts"],
                        default="counts",
                        help="Feature representation to use")
    
    format_group = parser.add_mutually_exclusive_group(required=True)
    format_group.add_argument("--vcf", type=str, help="VCF file to import")
    format_group.add_argument("--fasta", type=str, help="FASTA file to import")
    
    parser.add_argument("--populations",
                        type=str,
                        help="Populations file to import")
    
    parser.add_argument("--workdir",
                        type=str,
                        help="Work directory",
                        required=True)

    return vars(parser.parse_args())

if __name__ == "__main__":
    args = parseargs()

    if args["vcf"] is not None:
        import_vcf(args)
    elif args["fasta"] is not None:
        import_fasta(args)
    else:
        print "Need to specify VCF or FASTA file for import"
        sys.exit(1)
