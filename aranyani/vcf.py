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

def extract_features(chrom, pos, snps):
    n_individuals = len(snps)

    all_nucl = set()
    
    for nucl_dict in snps:
        all_nucl.update(nucl_dict.keys())

    features = dict()
    for nucl in all_nucl:
        features[(chrom, pos, nucl)] = np.zeros((n_individuals))

    for idx, nucl_dict in enumerate(snps):
        for nucl, count in nucl_dict.items():
            features[(chrom, pos, nucl)][idx] = count

    return features

def vcf_line_to_seq(ln, kept_individual_idx):
    cols = ln.strip().split()
    ref_seq = cols[DEFAULT_COLUMNS["REF"]]
    alt_seq = cols[DEFAULT_COLUMNS["ALT"]]

    snps = []
    for i, col in enumerate(cols[len(DEFAULT_COLUMNS):]):
        if i not in kept_individual_idx:
            continue
        
        tags = col.split(":")[0].split("/")

        nucleotides = defaultdict(int)

        for tag in tags:
            if tag == ".":
                # ignore unknown nucleotides
                pass
            elif tag == "0":
                nucleotides[ref_seq] += 1
            elif tag == "1":
                nucleotides[alt_seq] += 1
            else:
                raise NotImplementedError("Support for non-biallelic SNPs not implemented. Found '" + tag + "'")

        snps.append(dict(nucleotides))

    return extract_features(cols[DEFAULT_COLUMNS["CHROM"]], cols[DEFAULT_COLUMNS["POS"]], \
                            snps)

def stream_vcf_fl(flname, kept_individuals):
    with open(flname) as fl:
        for ln in fl:
            if ln.startswith("#CHROM"):
                column_names = ln[1:].strip().split()
                break

        all_ids = column_names[len(DEFAULT_COLUMNS):]

        individual_ids = [ident for ident in all_ids
                          if ident in kept_individuals]
        
        individual_idx = list([i for i, ident in enumerate(all_ids)
                              if ident in kept_individuals])

        yield individual_ids, individual_idx

        for ln in fl:
            if not ln.startswith("#"):
                yield vcf_line_to_seq(ln, individual_idx)

def read_dimensions(inflname, groups):
    # Data columns 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
    # individuals are columns after data columns
    # variants are rows

    # read file once to get dimensions
    gen = stream_vcf_fl(inflname, groups.keys())
    kept_ids, kept_idx = list(next(gen))

    n_individuals = len(kept_ids)

    n_features = 0
    for features in gen:
        # all nucleotides are missing or only one genotype
        if len(features) <= 1:
            continue

        n_features += len(features)

    return kept_ids, n_individuals, n_features

def is_fixed_difference(snp_features, class_labels):
    n_individuals = len(class_labels)

    individuals_without_missing = []
    for idx in xrange(n_individuals):
        nucl_count = 0
        for key, features in snp_features.items():
            nucl_count += features[idx]
        if nucl_count == 2:
            individuals_without_missing.append(idx)

    class_entries = defaultdict(set)
    for idx in individuals_without_missing:
        nucl_counts = []
        for key, features in snp_features.items():
            nucl_counts.append(features[idx])
        class_entries[class_labels[idx]].add(tuple(nucl_counts))

    all_class_entries = reduce(lambda x, y: x.intersection(y), class_entries.values())

    # fixed differences, any missing
    return len(all_class_entries) == 0, \
        len(individuals_without_missing) != len(class_labels)
    

def convert(groups_flname, vcf_flname, outbase):
    groups = read_groups(groups_flname)
    
    individual_ids, n_individuals, n_features = read_dimensions(vcf_flname, groups)

    print n_individuals, "individuals", n_features, "features"
    # 4 bytes / float, 1 MB = 1024 * 1024 bytes
    print "Estimated file size (MB):", 4 * n_individuals * n_features/ float(1024 * 1024)

    flname = os.path.join(outbase, FEATURE_MATRIX_FLNAME)
    feature_matrix = create_feature_matrix(flname, n_individuals, n_features)

    feature_idx = 0

    gen = stream_vcf_fl(vcf_flname, groups.keys())
    individual_ids, individual_idx = list(next(gen))

    class_labels = [None] * n_individuals
    for idx, ident in enumerate(individual_ids):
        class_labels[idx] = groups[ident]

    
    feature_labels = [None] * n_features
    fixed_differences = []
    missing = []
    column_idx = 0
    for snp_features in gen:
        # all nucleotides are missing or only one genotype
        if len(snp_features) <= 1:
            continue

        fd, missing_data = is_fixed_difference(snp_features, class_labels)
        
        for label, column in snp_features.iteritems():
            feature_labels[column_idx] = label
            feature_matrix[:, column_idx] = column
            column_idx += 1
            fixed_differences.append(fd)
            missing.append(missing_data)

    # close and save
    del feature_matrix

    to_json(os.path.join(outbase, FEATURE_LABELS_FLNAME), feature_labels)
    to_json(os.path.join(outbase, SAMPLE_LABELS_FLNAME), individual_ids)
    to_json(os.path.join(outbase, CLASS_LABELS_FLNAME), class_labels)
    to_json(os.path.join(outbase, FIXED_DIFFERENCES_FLNAME), fixed_differences)
    to_json(os.path.join(outbase, MISSING_DATA_FLNAME), missing)

    
    
