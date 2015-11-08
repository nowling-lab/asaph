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

from exceptions import NotImplementedError
import os
import sys

import numpy as np

from newioutils import *

DEFAULT_COLUMNS = {'CHROM' : 0, 'POS' : 1, 'ID' : 2, 'REF' : 3, 'ALT' : 4, 'QUAL' : 5, 'FILTER' : 6, 'INFO' : 7, 'FORMAT' : 8}

def read_groups(flname):
    fl = open(flname)
    groups = dict()
    for group_id, ln in enumerate(fl):
        cols = ln.strip().split(",")
        for ident in cols[1:]:
            groups[ident] = group_id

    return groups

def vcf_line_to_seq(ln):
    cols = ln.split()
    ref_seq = cols[DEFAULT_COLUMNS["REF"]]
    alt_seq = cols[DEFAULT_COLUMNS["ALT"]]

    snps1 = []
    snps2 = []
    for i, col in enumerate(cols[len(DEFAULT_COLUMNS):]):
        tag1, tag2 = col.split(":")[0].split("/")

        if tag1 == ".":
            nucl1 = "X"
        elif tag1 == "0":
            nucl1 = ref_seq
        elif tag1 == "1":
            nucl1 = alt_seq
        else:
            raise NotImplementedError("Support for non-biallelic SNPs not implemented. Found '" + tag1 + "'")

        snps1.append(nucl1)

        if tag2 == ".":
            nucl2 = "X"
        elif tag2 == "0":
            nucl2 = ref_seq
        elif tag2 == "1":
            nucl2 = alt_seq
        else:
            raise NotImplementedError("Support for non-biallelic SNPs not implemented. Found '" + tag1 + "'")

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

def read_dimensions(inflname, groups):
    # Data columns 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
    # individuals are columns after data columns
    # variants are rows

    # read file once to get dimensions
    gen = stream_vcf_fl(inflname)
    all_ids = list(next(gen))
    individual_ids = [ident for ident in all_ids
                      if ident in groups.keys()]
    individual_idx = set([i for i, ident in enumerate(all_ids)
                          if ident in groups.keys()])
    
    n_individuals = len(individual_idx)

    n_features = 0
    for _, _, haploid1, haploid2 in gen:
        filtered_haploid1 = [value for i, value in enumerate(haploid1)
                             if i in individual_idx]
        filtered_haploid2 = [value for i, value in enumerate(haploid2)
                             if i in individual_idx]

        # ignore positions where we only have 1 value
        if len(extract_features(filtered_haploid1)) >= 2:
            n_features += len(extract_features(filtered_haploid1))
        if len(extract_features(filtered_haploid2)) >= 2:
            n_features += len(extract_features(filtered_haploid2))

    return individual_ids, n_individuals, n_features

def convert(groups_flname, vcf_flname, outbase):
    groups = read_groups(groups_flname)
    
    individual_ids, n_individuals, n_features = read_dimensions(vcf_flname, groups)

    flname = os.path.join(outbase, FEATURE_MATRIX_FLNAME)
    feature_matrix = create_feature_matrix(flname, n_individuals, n_features)

    feature_idx = 0

    gen = stream_vcf_fl(vcf_flname)
    all_ids = list(next(gen))
    individual_idx = set([i for i, ident in enumerate(all_ids)
                          if ident in groups.keys()])

    feature_labels = [None] * n_features
    for chrom, pos, haploid1, haploid2 in gen:
        for haploid_idx, snps in enumerate([haploid1, haploid2]):
            filtered_snps = [value for i, value in enumerate(snps)
                             if i in individual_idx]
            feature_offsets = extract_features(filtered_snps)

            # ignore positions where we only have 1 value
            if len(feature_offsets) < 2:
                continue
            
            for indiv, nucl in enumerate(filtered_snps):
                if nucl != "X":
                    feature_matrix[indiv, feature_idx + feature_offsets[nucl]] = 1.0

            for value, offset in feature_offsets.items():
                feature_labels[feature_idx + offset] = (chrom, pos, haploid_idx, value)
                
            feature_idx += len(feature_offsets)

    del feature_matrix

    to_json(os.path.join(outbase, FEATURE_LABELS_FLNAME), feature_labels)
    to_json(os.path.join(outbase, SAMPLE_LABELS_FLNAME), individual_ids)

    class_labels = [None] * n_individuals
    for idx, ident in enumerate(individual_ids):
        class_labels[idx] = groups[ident]

    to_json(os.path.join(outbase, CLASS_LABELS_FLNAME), class_labels)

    
    
