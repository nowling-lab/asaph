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
import random
import sys

import numpy as np
from scipy import sparse

from sklearn.feature_extraction import FeatureHasher
from sklearn.random_projection import johnson_lindenstrauss_min_dim as jl_min_dim
from sklearn.random_projection import SparseRandomProjection

from .newioutils import *
from .models import ProjectSummary
from .models import COUNTS_FEATURE_TYPE
from .models import CATEGORIES_FEATURE_TYPE

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
            with gzip.open(self.flname, mode="rt") as fl:
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
                yield parse_vcf_line(ln, self.kept_pairs)

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

class FeatureStringsExtractor(object):
    def __init__(self, stream):
        self.stream = stream

    def __iter__(self):
        for variant_label, alleles, genotypes in self.stream:
            string_features = [None] * len(genotypes)

            homo_ref = ("%s_%s_het_%s" % (variant_label[0],
                                          variant_label[1],
                                          alleles[0]),
                        1)

            homo_alt = ("%s_%s_het_%s" % (variant_label[0],
                                          variant_label[1],
                                          alleles[1]),
                        1)

            het = ("%s_%s_%s_%s" % (variant_label[0],
                                    variant_label[1],
                                    alleles[0],
                                    alleles[1]),
                   1)

            for i, allele_counts in enumerate(genotypes):
                if allele_counts == (2, 0):
                    string_features[i] = homo_ref
                elif allele_counts == (0, 2):
                    string_features[i] = homo_alt
                elif allele_counts == (1, 1):
                    string_features[i] = het

            yield string_features

class StreamCounter(object):
    def __init__(self, stream):
        self.count = 0
        self.stream = stream

    def __iter__(self):
        for item in self.stream:
            self.count += 1
            yield item

class Chunker(object):
    def __init__(self, chunk_size, n_samples, stream):
        self.chunk_size = chunk_size
        self.n_samples = n_samples
        self.buffer = [list() for i in range(n_samples)]
        self.chunk_count = 0
        self.processed_chunks = 0
        self.stream = stream

    def __iter__(self):
        for feature_pairs in self.stream:
            if self.chunk_count >= self.chunk_size:
                self.processed_chunks += 1
                print("Processed chunk", self.processed_chunks)
                yield self.buffer
                self.buffer = [list() for i in range(self.n_samples)]
                self.chunk_count = 0

            self.chunk_count += 1
            for i, pair in enumerate(feature_pairs):
                if pair != None:
                    self.buffer[i].append(pair)

        if self.chunk_count > 0:
            yield self.buffer

class FeatureHashingAccumulator(object):
    def __init__(self, n_features, n_samples):
        self.n_features = n_features
        self.n_samples = n_samples
        self.transformer = FeatureHasher(n_features = n_features,
                                         input_type = "pair")
        self.features = np.zeros((n_samples, n_features), dtype=np.float32)

    def transform(self, stream):
        for chunk in stream:
            self.features += self.transformer.transform(chunk).toarray()
        return self.features

class FullMatrixAccumulator(object):
    def transform(self, stream):
        feature_columns = []
        for (chrom, pos, gt), column in stream:
            feature_columns.append(column)

        # need to transpose, otherwise we get (n_features, n_individuals) instead
        feature_matrix = np.array(feature_columns).T

        return feature_matrix

class ReservoirMatrixAccumulator(object):
    """
    Online sampling of columns using reservoir sampling.
    """
    def __init__(self, n_features):
        self.n_features = n_features

    def transform(self, stream):
        feature_columns = []
        for feature_idx, ((chrom, pos, gt), column) in enumerate(stream):
            if feature_idx < self.n_features:
                feature_columns.append(column)
            else:
                j = random.randint(0, feature_idx + 1)
                if j < self.n_features:
                    feature_columns[j] = column

        # need to transpose, otherwise we get (n_features, n_individuals) instead
        feature_matrix = np.array(feature_columns).T

        return feature_matrix

def convert(vcf_flname, outbase, matrix_type, feature_type, subsample_features, chunk_size, compressed_vcf, allele_min_freq_threshold):
    # dictionary of individual ids to population ids
    stream = VCFStreamer(vcf_flname, compressed_vcf)
    n_samples = len(stream.individual_names)

    # remove SNPs with least-frequently occurring alleles less than a threshold
    variants = filter_invariants(allele_min_freq_threshold,
                                 stream)
    filtered_positions_counter = StreamCounter(variants)

    min_dim = jl_min_dim(n_samples, eps=0.1)
    
    # extract features
    if matrix_type == BAG_OF_WORDS_MATRIX_TYPE and subsample_features != "reservoir":
        if feature_type == COUNTS_FEATURE_TYPE:
            extractor = CountFeaturesExtractor(filtered_positions_counter)
        elif feature_type == CATEGORIES_FEATURE_TYPE:
            extractor = CategoricalFeaturesExtractor(filtered_positions_counter)
        else:
            raise Exception("Unknown feature type: %s" % feature_type)

        accumulator = FullMatrixAccumulator()
        feature_matrix = accumulator.transform(extractor)

        if subsample_features == "column-sampling":
            print(feature_matrix.shape[1], "features")
            print("You specified subsampling features.  Subsampling columns of feature matrix.")
            selected_idx = random.sample(range(feature_matrix.shape[1]), min_dim)
            feature_matrix = feature_matrix[:, selected_idx]
        elif subsample_features == "random-projection":
            print(feature_matrix.shape[1], "features")
            print("You specified subsampling features.  Using sparse random projection.")
            srp = SparseRandomProjection(n_components = min_dim)
            feature_matrix = srp.fit_transform(feature_matrix)

    elif matrix_type == BAG_OF_WORDS_MATRIX_TYPE and subsample_features == "reservoir":
        if feature_type == COUNTS_FEATURE_TYPE:
            extractor = CountFeaturesExtractor(filtered_positions_counter)
        elif feature_type == CATEGORIES_FEATURE_TYPE:
            extractor = CategoricalFeaturesExtractor(filtered_positions_counter)
        else:
            raise Exception("Unknown feature type: %s" % feature_type)

        accumulator = ReservoirMatrixAccumulator(min_dim)
        feature_matrix = accumulator.transform(extractor)
            
    elif matrix_type == HASHED_MATRIX_TYPE:
        string_features = FeatureStringsExtractor(filtered_positions_counter)
        chunker = Chunker(chunk_size, n_samples, string_features)
        accumulator = FeatureHashingAccumulator(min_dim, n_samples)
        feature_matrix = accumulator.transform(chunker)
    else:
        raise Exception("Unknown matrix type: %s" % matrix_type)

    print(feature_matrix.shape[0], "individuals")
    print(stream.positions_read, "variants read")
    print(filtered_positions_counter.count, "variants kept")
    print(feature_matrix.shape[1], "features")

    project_summary = ProjectSummary(original_positions = stream.positions_read,
                                     filtered_positions = filtered_positions_counter.count,
                                     n_features = feature_matrix.shape[1],
                                     n_samples = n_samples,
                                     matrix_type = matrix_type,
                                     feature_type = feature_type)

    np.savez_compressed(os.path.join(outbase, FEATURE_MATRIX_FLNAME + ".npz"),
                        feature_matrix = feature_matrix)
    serialize(os.path.join(outbase, SAMPLE_LABELS_FLNAME), stream.rows_to_names)
    serialize(os.path.join(outbase, PROJECT_SUMMARY_FLNAME), project_summary)
