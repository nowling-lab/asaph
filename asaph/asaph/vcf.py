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

from collections import OrderedDict
import gzip
import os
import sys

import numpy as np

from newioutils import *

from models import ProjectSummary
from models import COUNTS_FEATURE_TYPE
from models import CATEGORIES_FEATURE_TYPE

DEFAULT_COLUMNS = {'CHROM' : 0, 'POS' : 1, 'ID' : 2, 'REF' : 3, 'ALT' : 4, 'QUAL' : 5, 'FILTER' : 6, 'INFO' : 7, 'FORMAT' : 8}

GENOTYPE_OFFSET = 9

UNKNOWN_GENOTYPE = (0, 0)

## Parsing
def parse_vcf_line(ln, kept_pairs):
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

    individual_genotypes = [None] * len(kept_pairs)
    for i, (idx, name) in enumerate(kept_pairs):
        col = cols[GENOTYPE_OFFSET + idx]

        ref_count = 0
        alt_count = 0

        # avoid caring whether / or | is used as separator by indexing
        # ignore unknown genotype (.)
        if col[0] == "0":
            ref_count += 1
        elif col[0] == "1":
            alt_count += 1

        if col[2] == "0":
            ref_count += 1
        elif col[2] == "1":
            alt_count += 1

        individual_genotypes[i] = (ref_count, alt_count)

    variant_label = (cols[DEFAULT_COLUMNS["CHROM"]], cols[DEFAULT_COLUMNS["POS"]])
    return (variant_label, alleles, tuple(individual_genotypes))

class VCFStreamer(object):
    def __init__(self, flname, compressed, kept_individuals=None):
        self.flname = flname
        if kept_individuals:
            self.kept_individuals = set(kept_individuals)
        else:
            self.kept_individuals = None
        self.compressed = compressed
        self.positions_read = 0

        self.stream = self.__open__()
        for ln in self.stream:
            if ln.startswith("#CHROM"):
                column_names = ln[1:].strip().split()
                self.individual_names = column_names[len(DEFAULT_COLUMNS):]
                break
            elif ln.startswith("#"):
                continue

        self.kept_pairs = [(i, name) for i, name in enumerate(self.individual_names)
                           if self.kept_individuals is None \
                           or name in self.kept_individuals]

        self.rows_to_names = [name for name in self.individual_names
                              if self.kept_individuals is None \
                              or name in self.kept_individuals]

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
        for ln in self.stream:
            if not ln.startswith("#"):
                self.positions_read += 1
                yield parse_vcf_line(ln,
                                     self.kept_pairs)

## Filters
def filter_invariants(min_percentage, stream):
    """
    Filter out variants where the least-frequently occurring allele occurs less than some threshold.

    0 <= min_percentage < 1
    """
    for label, alleles, genotypes in stream:
        total_ref_count = 0
        total_alt_count = 0
        for sample_ref_count, sample_alt_count in genotypes:
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

## Feature extraction
class CountFeaturesExtractor(object):
    def __init__(self, stream):
        self.stream = stream

    def __iter__(self):
        for variant_label, alleles, genotypes in self.stream:
            chrom, pos = variant_label
            ref_column = [0.] * len(genotypes)
            alt_column = [0.] * len(genotypes)

            for row_idx, allele_counts in enumerate(genotypes):
                ref_column[row_idx] = allele_counts[0]
                alt_column[row_idx] = allele_counts[1]

            yield (chrom, pos, alleles[0]), tuple(ref_column)
            yield (chrom, pos, alleles[1]), tuple(alt_column)

class CategoricalFeaturesExtractor(object):
    def __init__(self, stream):
        self.stream = stream

    def __iter__(self):
        for variant_label, alleles, genotypes in self.stream:
            chrom, pos = variant_label
            homo1_column = [0.] * len(genotypes)
            homo2_column = [0.] * len(genotypes)
            het_column = [0.] * len(genotypes)

            for row_idx, allele_counts in enumerate(genotypes):
                if allele_counts == (2, 0):
                    homo1_column[row_idx] = 1
                elif allele_counts == (0, 2):
                    homo2_column[row_idx] = 1
                elif allele_counts == (1, 1):
                    het_column[row_idx] = 1

            yield (chrom, pos, (alleles[0] + "/" + alleles[0])), tuple(homo1_column)
            yield (chrom, pos, (alleles[1] + "/" + alleles[1])), tuple(homo2_column)
            yield (chrom, pos, (alleles[0] + "/" + alleles[1])), tuple(het_column)

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
    if groups_flname is not None:
        populations, population_names = read_populations(groups_flname)
        # returns triplets of (variant_label, variant_alleles, genotype_counts)
        stream = VCFStreamer(vcf_flname, compressed_vcf, populations.keys())
    else:
        population_names = None
        stream = VCFStreamer(vcf_flname, compressed_vcf)

    # remove SNPs with least-frequently occurring alleles less than a threshold
    variants = filter_invariants(allele_min_freq_threshold,
                                 stream)

    filtered_positions_counter = StreamCounter(variants)

    # extract features
    if feature_type == COUNTS_FEATURE_TYPE:
        extractor = CountFeaturesExtractor(filtered_positions_counter)
    elif feature_type == CATEGORIES_FEATURE_TYPE:
        extractor = CategoricalFeaturesExtractor(filtered_positions_counter)
    else:
        raise Exception, "Unknown feature type: %s" % feature_type

    snp_features = defaultdict(list)
    snp_genotypes = defaultdict(dict)
    column_idx = 0
    feature_columns = []
    for feature_idx, ((chrom, pos, gt), column) in enumerate(extractor):
        snp_genotypes[(chrom, pos)][column_idx] = gt
        snp_features[(chrom, pos)].append(column_idx)
        feature_columns.append(column)
        column_idx += 1

    # need to transpose, otherwise we get (n_features, n_individuals) instead
    feature_matrix = np.array(feature_columns).T

    print feature_matrix.shape[0], "individuals"
    print feature_matrix.shape[1], "features"

    project_summary = ProjectSummary(original_positions = stream.positions_read,
                                     filtered_positions = filtered_positions_counter.count,
                                     n_features = feature_matrix.shape[1],
                                     feature_encoding = feature_type,
                                     compressed = compress,
                                     n_samples = len(stream.rows_to_names),
                                     population_names = population_names)

    if compress:
        np.savez_compressed(os.path.join(outbase, FEATURE_MATRIX_FLNAME + ".npz"),
                            feature_matrix = feature_matrix)
    else:
        np.save(os.path.join(outbase, FEATURE_MATRIX_FLNAME), feature_matrix)
    serialize(os.path.join(outbase, SAMPLE_LABELS_FLNAME), stream.rows_to_names)
    serialize(os.path.join(outbase, SNP_FEATURE_INDICES_FLNAME), snp_features)
    serialize(os.path.join(outbase, SNP_FEATURE_GENOTYPES_FLNAME), snp_genotypes)
    serialize(os.path.join(outbase, PROJECT_SUMMARY_FLNAME), project_summary)
