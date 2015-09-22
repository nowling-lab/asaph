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

import csv
import sys

import numpy as np
import pandas as pd

DEFAULT_COLUMNS = {'CHROM' : 0, 'POS' : 1, 'ID' : 2, 'REF' : 3, 'ALT' : 4, 'QUAL' : 5, 'FILTER' : 6, 'INFO' : 7, 'FORMAT' : 8}

def vcf_line_to_seq(ln):
    cols = ln.split()
    ref_seq = cols[DEFAULT_COLUMNS["REF"]]
    alt_seq = cols[DEFAULT_COLUMNS["ALT"]]

    snps1 = []
    snps2 = []
    for i, col in enumerate(cols[len(DEFAULT_COLUMNS):]):
        tag1, tag2 = col.split(":")[0].split("/")

        nucl1 = "X"
        if tag1 == "0":
            nucl1 = ref_seq
        elif tag1 == "1":
            nucl1 = alt_seq

        snps1.append(nucl1)

        nucl2 = "X"
        if tag2 == "0":
            nucl2 = ref_seq
        elif tag2 == "1":
            nucl2 = alt_seq

        snps2.append(nucl2)

    return cols[DEFAULT_COLUMNS["CHROM"]], cols[DEFAULT_COLUMNS["POS"]], tuple(snps1), tuple(snps2)

def extract_features(seq):
    unique_nucl = set(seq)
    unique_nucl.discard("X")
    return { nucl : offset for offset, nucl in enumerate(unique_nucl) }

def stream_vcf_fl(flname):
    with open(flname) as fl:
        for ln in fl:
            if ln.startswith("#CHROM"):
                column_names = ln[1:].strip().split()
                break

        individual_ids = column_names[len(DEFAULT_COLUMNS):]

        yield individual_ids

        for ln in fl:
            if not ln.startswith("#"):
                chrom, pos, snps1, snps2 = vcf_line_to_seq(ln)
                yield chrom, pos, snps1, snps2

def read_dimensions(inflname):
    # Data columns 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
    # individuals are columns after data columns
    # variants are rows

    # read file once to get dimensions
    gen = stream_vcf_fl(inflname)
    individual_ids = next(gen)        
    n_individuals = len(individual_ids)

    n_features = 0
    for _, _, snps1, snps2 in gen:
        n_features += len(extract_features(snps1))
        n_features += len(extract_features(snps2))

    return individual_ids, n_individuals, n_features

def write_individuals(outbase, individual_ids):
    individuals_fl = open(outbase + ".individuals", "w")
    for indiv in individual_ids:
        individuals_fl.write(indiv)
        individuals_fl.write("\n")
    individuals_fl.close()

def append_feature_names(features_fl, feat_map, chrom, pos, diploid):
    for nucl, idx in sorted(feat_map.iteritems(), key=lambda t: t[1]):
        features_fl.write(chrom)
        features_fl.write(",")
        features_fl.write(pos)
        features_fl.write(",")
        features_fl.write(diploid)
        features_fl.write(",")
        features_fl.write(nucl)
        features_fl.write("\n")

def convert(inflname, outbase):
    individual_ids, n_individuals, n_features = read_dimensions(inflname)
    print n_individuals, n_features

    write_individuals(outbase, individual_ids)
    
    feature_matrix = np.memmap(filename=outbase + ".npy", dtype="float32", \
                               mode="w+", shape=(n_individuals, n_features))

    features_fl = open(outbase + ".features", "w")

    feature_idx = 0
    gen = stream_vcf_fl(inflname)
    # discard ids
    next(gen)
        
    for chrom, pos, snps1, snps2 in gen:
        for diploid, snps in enumerate([snps1, snps2]):
            feat_map = extract_features(snps)
            for indiv, nucl in enumerate(snps):
                if nucl != "X":
                    feature_matrix[indiv, feature_idx + feat_map[nucl]] = 1.0

            feature_idx += len(feat_map)
            append_feature_names(features_fl, feat_map, chrom, pos, str(diploid))

    # close
    del feature_matrix
    features_fl.close()                

if __name__ == "__main__":
    outbase = sys.argv[1]
    vcf_flnames = sys.argv[2]

    convert(vcf_flnames, outbase)
