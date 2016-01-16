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

UNKNOWN_GENOTYPE = ("X", "X")

def read_groups(flname):
    fl = open(flname)
    groups = dict()
    for group_id, ln in enumerate(fl):
        cols = ln.strip().split(",")
        for ident in cols[1:]:
            groups[ident] = group_id

    return groups

def extract_features(triplet):
    label, snps, genotypes = triplet
    chrom, pos = label
    n_individuals = len(snps)

    features = dict()
    for genotype in genotypes:
        features[(chrom, pos, genotype)] = np.zeros(n_individuals)

    for idx, genotype in enumerate(snps):
        if genotype != UNKNOWN_GENOTYPE:
            features[(chrom, pos, genotype)][idx] = 1

    return ((chrom, pos), [(key, tuple(value)) for key, value in features.iteritems()])

def parse_vcf_line(ln):
    cols = ln.strip().split()
    ref_seq = cols[DEFAULT_COLUMNS["REF"]]
    alt_seq = cols[DEFAULT_COLUMNS["ALT"]]

    snps = []
    for i, col in enumerate(cols[len(DEFAULT_COLUMNS):]):
        left, right = col.split(":")[0].split("/")

        homo1 = (ref_seq, ref_seq)
        homo2 = (alt_seq, alt_seq)
        hetero = tuple(sorted((ref_seq, alt_seq)))

        if left == "0" and right == "0":
            snps.append(homo1)
        elif left == "1" and right == "1":
            snps.append(homo2)
        elif (left == "0" and right == "1") or \
             (left == "1" and right == "0"):
            snps.append(hetero)
        # even if only one of the two alleles is unknown,
        # consider both unknown
        else:
            snps.append(UNKNOWN_GENOTYPE)

    return ((cols[DEFAULT_COLUMNS["CHROM"]], cols[DEFAULT_COLUMNS["POS"]]),
            snps,
            (homo1, homo2, hetero))

class SelectIndividuals(object):
    def __init__(self, kept_individual_idx):
        self.kept_individual_idx = kept_individual_idx

    def __call__(self, triplets):
        label, snps, genotypes = triplets
        selected = [snps[idx] for idx in self.kept_individual_idx]
        return (label, selected, genotypes)

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
                yield ln

class ImputeUnknown(object):
    def __init__(self, class_labels, threshold):
        self.class_labels = class_labels
        self.threshold = threshold

    def __call__(self, triplet):
        label, snps, genotypes = triplet
        n_individuals = len(self.class_labels)
        
        class_entries = defaultdict(lambda: defaultdict(int))
        for idx in xrange(n_individuals):
            genotype = snps[idx]
            class_entries[self.class_labels[idx]][genotype] += 1

        class_modes = dict()
        for class_label, genotypes in class_entries.iteritems():
            class_size = float(sum(genotypes.values()))
            _, mode = max([(count, genotype) for genotype, count in genotypes.items()
                           if genotype != UNKNOWN_GENOTYPE])
            mode_ratio = float(genotypes[mode]) / class_size
            class_modes[class_label] = (mode, mode_ratio)

        imputed_snps = []
        for idx, genotype in enumerate(snps):
            class_label = self.class_labels[idx]
            mode, ratio = class_modes[class_label]
            if genotype == UNKNOWN_GENOTYPE and ratio > self.threshold:
                imputed_snps.append(mode)
            else:
                imputed_snps.append(genotype)

        return (label, imputed_snps, genotypes)

class FilterUnknown(object):
    def __init__(self, class_labels):
        self.class_labels = class_labels

    def __call__(self, triplet):
        label, snps, genotypes = triplet
        n_individuals = len(self.class_labels)
        
        class_entries = defaultdict(set)
        for idx in xrange(n_individuals):
            genotype = snps[idx]
            class_entries[self.class_labels[idx]].add(genotype)

        entire_class_unknown = False
        for class_label, genotypes in class_entries.items():
            if genotypes == set([UNKNOWN_GENOTYPE]):
                entire_class_unknown = True

        all_known_genotypes = set(snps)
        all_known_genotypes.discard(UNKNOWN_GENOTYPE)
        single_known_genotype = False
        if len(all_known_genotypes) <= 1:
            single_known_genotype = True

        return not(single_known_genotype or entire_class_unknown)

class AnnotateTrivial(object):
    def __init__(self, class_labels):
        self.class_labels = class_labels
        self.trivial_snps = dict()
        self.unknown_genotypes = dict()

    def __call__(self, triplet):
        snp_label, snps, genotypes = triplet
        n_individuals = len(self.class_labels)

        individuals_without_missing = []
        for idx in xrange(n_individuals):
            if snps[idx] != UNKNOWN_GENOTYPE:
                individuals_without_missing.append(idx)

        class_entries = defaultdict(set)
        for idx in individuals_without_missing:
            genotype = snps[idx]
            class_entries[self.class_labels[idx]].add(genotype)

        all_class_entries = reduce(lambda x, y: x.intersection(y), class_entries.values())

        # fixed differences, any missing
        self.trivial_snps[snp_label] = len(all_class_entries) == 0
        self.unknown_genotypes[snp_label] = len(individuals_without_missing) != len(self.class_labels)

        return triplet
    

def convert(groups_flname, vcf_flname, outbase, compress, impute_threshold):
    groups = read_groups(groups_flname)
    
    gen = stream_vcf_fl(vcf_flname, groups.keys())
    individual_ids, individual_idx = list(next(gen))

    class_labels = [groups[ident] for ident in individual_ids]
    annotation = AnnotateTrivial(class_labels)
    filter_unknown = FilterUnknown(class_labels)

    parsed_lines = map(parse_vcf_line, gen)
    selected_individuals = map(SelectIndividuals(individual_idx), parsed_lines)
    # remove SNPs with < 2 known genotypes
    known_features = filter(lambda triplet: filter_unknown(triplet), selected_individuals)
    if impute_threshold != None:
        impute_unknown = ImputeUnknown(class_labels, impute_threshold)
        known_features = map(impute_unknown, known_features)
    annotated_features = map(annotation, known_features)
    extracted_features = map(extract_features, annotated_features)
    
    feature_labels = []
    column_idx = 0
    feature_column = dict()
    feature_columns = []
    for pairs in extracted_features:
        snp_label, labeled_columns = pairs

        if compress:
            for label, column in labeled_columns:
                if column not in feature_column:
                    feature_column[column] = column_idx
                    feature_labels.append([label])
                    feature_columns.append(column)
                    column_idx += 1
                else:
                    feature_column_idx = feature_column[column]
                    feature_labels[feature_column_idx].append(label)
        else:
            for label, column in labeled_columns:
                feature_labels.append([label])
                feature_columns.append(column)
                column_idx += 1

    # need to transpose, otherwise we get (n_features, n_individuals) instead
    feature_matrix = np.array(feature_columns).T

    print feature_matrix.shape[0], "individuals", feature_matrix.shape[1], "features"

    np.save(os.path.join(outbase, FEATURE_MATRIX_FLNAME), feature_matrix)
    to_json(os.path.join(outbase, FEATURE_LABELS_FLNAME), feature_labels)
    to_json(os.path.join(outbase, SAMPLE_LABELS_FLNAME), individual_ids)
    to_json(os.path.join(outbase, CLASS_LABELS_FLNAME), class_labels)
    to_json(os.path.join(outbase, FIXED_DIFFERENCES_FLNAME), annotation.trivial_snps)
    to_json(os.path.join(outbase, MISSING_DATA_FLNAME), annotation.unknown_genotypes)

    
    
