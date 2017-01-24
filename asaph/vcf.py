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
import shelve
import sys

import numpy as np

from newioutils import *

DEFAULT_COLUMNS = {'CHROM' : 0, 'POS' : 1, 'ID' : 2, 'REF' : 3, 'ALT' : 4, 'QUAL' : 5, 'FILTER' : 6, 'INFO' : 7, 'FORMAT' : 8}

UNKNOWN_GENOTYPE = (0, 0)

def read_groups(flname):
    fl = open(flname)
    groups = dict()
    for group_id, ln in enumerate(fl):
        cols = ln.strip().split(",")
        for ident in cols[1:]:
            groups[ident] = group_id

    return groups

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
    def __init__(self, flname):
        self.flname = flname
        self.individual_names = None

    def __iter__(self):
        with open(self.flname) as fl:
            for ln in fl:
                if ln.startswith("#CHROM"):
                    column_names = ln[1:].strip().split()
                    continue
                elif ln.startswith("#"):
                    continue

                self.individual_names = column_names[len(DEFAULT_COLUMNS):]
                
                for ln in fl:
                    if not ln.startswith("#"):
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

def filter_invariants(stream):
    """
    Filter out variants that have < 2 genotypes present
    """
    for label, alleles, genotypes in stream:
        observed_genotypes = set(genotypes.itervalues())
        if len(observed_genotypes) >= 2:
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
        
            
def convert(groups_flname, vcf_flname, outbase, compress):
    # dictionary of individual ids to population ids
    populations = read_groups(groups_flname)

    # returns triplets of (variant_label, variant_alleles, genotype_counts)
    stream = VCFStreamer(vcf_flname)
    selected_individuals = select_individuals(stream,
                                              populations.keys())
    # remove SNPs with < 2 known genotypes
    variants = filter_invariants(selected_individuals)

    # extract features
    extractor = CountFeaturesExtractor(variants, populations.keys())

    feature_names = shelve.open(os.path.join(outbase, FEATURE_LABELS_FLNAME))
    column_idx = 0
    col_dict = dict()
    feature_columns = []
    for col_name, column in extractor:
        if compress:
            if column not in col_dict:
                col_dict[column] = column_idx
                feature_names[str(column_idx)] = [col_name]
                feature_columns.append(column)
                column_idx += 1
            else:
                feature_column_idx = col_dict[column]
                all_names = feature_names[str(feature_column_idx)]
                all_names.append(col_name)
                feature_names[str(feature_column_idx)]
        else:
            feature_names[str(column_idx)].append([col_name])
            feature_columns.append(column)
            column_idx += 1

    # need to transpose, otherwise we get (n_features, n_individuals) instead
    feature_matrix = np.array(feature_columns).T

    print feature_matrix.shape[0], "individuals", feature_matrix.shape[1], "features"

    class_labels = [populations[ident] for ident in extractor.rows_to_names]

    np.save(os.path.join(outbase, FEATURE_MATRIX_FLNAME), feature_matrix)
    feature_names.close()
    to_json(os.path.join(outbase, SAMPLE_LABELS_FLNAME), extractor.rows_to_names)
    to_json(os.path.join(outbase, CLASS_LABELS_FLNAME), class_labels)
    
