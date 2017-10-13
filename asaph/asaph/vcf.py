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

import gzip
import os
import sys

import numpy as np

from newioutils import *

from models import ProjectSummary

DEFAULT_COLUMNS = {'CHROM' : 0, 'POS' : 1, 'ID' : 2, 'REF' : 3, 'ALT' : 4, 'QUAL' : 5, 'FILTER' : 6, 'INFO' : 7, 'FORMAT' : 8}

UNKNOWN_GENOTYPE = (0, 0)

def read_groups(flname):
    fl = open(flname)
    groups = dict()
    group_names = dict()
    for group_id, ln in enumerate(fl):
        cols = ln.strip().split(",")
        for ident in cols[1:]:
            groups[ident] = group_id
            group_names[group_id] = cols[0]

    return groups, group_names

## Parsing
def parse_vcf_line(ln, individual_names):
    """
    Takes a string containing a line of a VCF file
    and a list of individual names.

    Returns a pair of (variant_label, alleles, individual_genotypes).

    variant_label is a pair of (chromosome, position)
    alelles is a pair of (ref_seq, alt_seq)
    individual_genotypes is a dictionary of individual ids to pairs
    of (ref_count, alt_count)
    """

    cols = ln.strip().split()
    # TODO: Allow for more than 1 alternative sequence
    alleles = (cols[DEFAULT_COLUMNS["REF"]],
               cols[DEFAULT_COLUMNS["ALT"]])

    individual_genotypes = {}
    for i, col in enumerate(cols[len(DEFAULT_COLUMNS):]):
        genotype_pair = col.split(":")[0]

        ref_count = 0
        alt_count = 0

        # avoid caring whether / or | is used as separator by indexing
        # ignore unknown genotype (.)
        if genotype_pair[0] == "0":
            ref_count += 1
        elif genotype_pair[0] == "1":
            alt_count += 1

        if genotype_pair[2] == "0":
            ref_count += 1
        elif genotype_pair[2] == "1":
            alt_count += 1

        individual_genotypes[individual_names[i]] = (ref_count, alt_count)

    variant_label = (cols[DEFAULT_COLUMNS["CHROM"]], cols[DEFAULT_COLUMNS["POS"]])
    return (variant_label, alleles, individual_genotypes)

class VCFStreamer(object):
    def __init__(self, flname, compressed):
        self.flname = flname
        self.individual_names = None
        self.compressed = compressed
        self.positions_read = 0

    def __open__(self):
        if self.compressed:
            with gzip.open(self.flname) as fl:
                for ln in fl:
                    yield ln
        else:
            with open(self.flname) as fl:
                for ln in fl:
                    yield ln

    def __iter__(self):
        stream = self.__open__()
        for ln in stream:
            if ln.startswith("#CHROM"):
                column_names = ln[1:].strip().split()
                self.individual_names = column_names[len(DEFAULT_COLUMNS):]
                break
            elif ln.startswith("#"):
                continue

        for ln in stream:
            if not ln.startswith("#"):
                self.positions_read += 1
                yield parse_vcf_line(ln, self.individual_names)

## Filters
def select_individuals(stream, individual_ids):
    """
    Filter stream of pairs, only keeping genotypes
    for individuals with ids in individual_ids
    """
    kept_individuals = set(individual_ids)

    for label, alleles, genotypes in stream:
        selected = { name : genotypes[name] for name in kept_individuals }
        yield (label, alleles, selected)

def filter_invariants(min_percentage, stream):
    """
    Filter out variants where the least-frequently occurring allele occurs less than some threshold.

    0 <= min_percentage < 1
    """
    for label, alleles, genotypes in stream:
        total_ref_count = 0
        total_alt_count = 0
        for sample_ref_count, sample_alt_count in genotypes.itervalues():
            total_ref_count += sample_ref_count
            total_alt_count += sample_alt_count

        # all SNPs have unknown genotypes
        if total_ref_count == 0 and total_alt_count == 0:
            continue
            
        min_count = min(total_ref_count,
                        total_alt_count)

        fraction = min_count / float(total_ref_count + total_alt_count)
        
        if fraction >= min_percentage:
            yield (label, alleles, genotypes)
            
def filter_unknown(class_labels, stream):
    """
    Filter out variants where 1 or more classes contains all unknown genotypes.
    """
    for label, variant_alleles, genotypes in stream:
        class_entries = defaultdict(set)
        for idx, allele_counts in genotypes.iteritems():
            class_entries[class_labels[idx]].add(allele_counts)

        entire_class_unknown = False
        for class_label, genotypes in class_entries.items():
            if genotypes == set([UNKNOWN_GENOTYPE]):
                entire_class_unknown = True

        if not entire_class_unknown:
            yield (label, variant_alleles, genotypes)

## Feature extraction
class CountFeaturesExtractor(object):
    def __init__(self, stream, individual_names):
        self.stream = stream
        self.name_to_row = dict()
        self.rows_to_names = []
        for idx, name in enumerate(individual_names):
            self.name_to_row[name] = idx
            self.rows_to_names.append(name)

    def __iter__(self):
        for variant_label, alleles, genotypes in self.stream:
            chrom, pos = variant_label
            ref_column = [0.] * len(self.name_to_row)
            alt_column = [0.] * len(self.name_to_row)
            for name, allele_counts in genotypes.items():
                row_idx = self.name_to_row[name]
                ref_column[row_idx] = allele_counts[0]
                alt_column[row_idx] = allele_counts[1]
            yield (chrom, pos, alleles[0]), tuple(ref_column)
            yield (chrom, pos, alleles[1]), tuple(alt_column)

class CategoricalFeaturesExtractor(object):
    def __init__(self, stream, individual_names):
        self.stream = stream
        self.name_to_row = dict()
        self.rows_to_names = []
        for idx, name in enumerate(individual_names):
            self.name_to_row[name] = idx
            self.rows_to_names.append(name)

    def __iter__(self):
        for variant_label, alleles, genotypes in self.stream:
            chrom, pos = variant_label
            homo1_column = [0.] * len(self.name_to_row)
            homo2_column = [0.] * len(self.name_to_row)
            het_column = [0.] * len(self.name_to_row)
            for name, allele_counts in genotypes.items():
                row_idx = self.name_to_row[name]
                if allele_counts == (2, 0):
                    homo1_column[row_idx] = 1
                elif allele_counts == (0, 2):
                    homo2_column[row_idx] = 1
                elif allele_counts == (1, 1):
                    het_column[row_idx] = 1
                # ignore partially unknown e.g., A/X
            yield (chrom, pos, (alleles[0], alleles[0])), tuple(homo1_column)
            yield (chrom, pos, (alleles[1], alleles[1])), tuple(homo2_column)
            yield (chrom, pos, (alleles[0], alleles[1])), tuple(het_column)

class StreamCounter(object):
    def __init__(self, stream):
        self.count = 0
        self.stream = stream

    def __iter__(self):
        for item in self.stream:
            self.count += 1
            yield item


def convert(groups_flname, vcf_flname, outbase, compress, feature_type, compressed_vcf, allele_min_freq_threshold):
    # dictionary of individual ids to population ids
    populations, population_names = read_groups(groups_flname)

    # returns triplets of (variant_label, variant_alleles, genotype_counts)
    stream = VCFStreamer(vcf_flname, compressed_vcf)
    selected_individuals = select_individuals(stream,
                                              populations.keys())

    # remove SNPs with least-frequently occurring alleles less than a threshold
    variants = filter_invariants(allele_min_freq_threshold,
                                 selected_individuals)

    filtered_positions_counter = StreamCounter(variants)

    # extract features
    if feature_type == "counts":
        extractor = CountFeaturesExtractor(filtered_positions_counter,
                                           populations.keys())
    elif feature_type == "categories":
        extractor = CategoricalFeaturesExtractor(filtered_positions_counter,
                                                 populations.keys())
    else:
        raise Exception, "Unknown feature type: %s" % feature_type

    snp_features = defaultdict(list)
    column_idx = 0
    col_dict = dict()
    feature_columns = []
    for feature_idx, ((chrom, pos, _), column) in enumerate(extractor):
        if compress:
            if column not in col_dict:
                col_dict[column] = column_idx
                snp_features[(chrom, pos)].append(column_idx)
                feature_columns.append(column)
                column_idx += 1
            else:
                feature_column_idx = col_dict[column]
                snp_features[(chrom, pos)].append(feature_column_idx)
        else:
            snp_features[(chrom, pos)].append(column_idx)
            feature_columns.append(column)
            column_idx += 1

    # need to transpose, otherwise we get (n_features, n_individuals) instead
    feature_matrix = np.array(feature_columns).T

    print feature_matrix.shape[0], "individuals"
    if compress:
        print (feature_idx + 1), "features before compression"
        print feature_matrix.shape[1], "features after compression"
    else:
        print feature_matrix.shape[1], "features"

    class_labels = [populations[ident] for ident in extractor.rows_to_names]

    project_summary = ProjectSummary(original_positions = stream.positions_read,
                                     filtered_positions = filtered_positions_counter.count,
                                     n_features = feature_matrix.shape[1],
                                     feature_encoding = feature_type,
                                     compressed = compress,
                                     n_samples = len(class_labels),
                                     population_names = population_names)

    np.save(os.path.join(outbase, FEATURE_MATRIX_FLNAME), feature_matrix)
    serialize(os.path.join(outbase, SAMPLE_LABELS_FLNAME), extractor.rows_to_names)
    serialize(os.path.join(outbase, CLASS_LABELS_FLNAME), class_labels)
    serialize(os.path.join(outbase, SNP_FEATURE_INDICES_FLNAME), snp_features)
    serialize(os.path.join(outbase, PROJECT_SUMMARY_FLNAME), project_summary)
