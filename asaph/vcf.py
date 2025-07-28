"""
This module provides VCF file parsing and streaming functionality for population genetics
analysis, including classes for reading compressed and uncompressed VCF files, filtering
invariant sites, and converting genotype data into standardized formats for downstream
processing.

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

from .newioutils import *
from .models import ProjectSummary

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

        individual_genotypes[i] = (name, (ref_count, alt_count))

    variant_label = (cols[DEFAULT_COLUMNS["CHROM"]], cols[DEFAULT_COLUMNS["POS"]])
    return (variant_label, alleles, tuple(individual_genotypes))

class VCFStreamer:
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
            if ln.startswith("#"):
                continue

        self.kept_pairs = [(i, name) for i, name in enumerate(self.individual_names)
                           if self.kept_individuals is None \
                           or name in self.kept_individuals]

        self.rows_to_names = [name for name in self.individual_names
                              if self.kept_individuals is None \
                              or name in self.kept_individuals]

    def __open__(self):
        if self.compressed:
            with gzip.open(self.flname, mode="rt", encoding="utf-8") as fl:
                yield from fl
        else:
            with open(self.flname, "rt", encoding="utf-8") as fl:
                yield from fl

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
        for _, (sample_ref_count, sample_alt_count) in genotypes:
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

class StreamCounter:
    def __init__(self, stream):
        self.count = 0
        self.stream = stream

    def __iter__(self):
        for item in self.stream:
            self.count += 1
            yield item

def stream_vcf_variants(vcf_flname, compressed_vcf, allele_min_freq_threshold):
    # dictionary of individual ids to population ids
    stream = VCFStreamer(vcf_flname, compressed_vcf)

    # remove SNPs with least-frequently occurring alleles less than a threshold
    variants = filter_invariants(allele_min_freq_threshold,
                                 stream)

    return variants, stream.rows_to_names
