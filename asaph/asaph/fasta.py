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

import numpy as np
from scipy import sparse

from newioutils import *
from vcf import UNKNOWN_GENOTYPE

DNA_MAP = {"A" : 0,
           "T" : 1,
           "C" : 2,
           "G" : 3}

AA_MAP = {"A" : 0,
          "C" : 1,
          "D" : 2,
          "E" : 3,
          "F" : 4,
          "G" : 5,
          "H" : 6,
          "I" : 7,
          "K" : 8,
          "L" : 9,
          "M" : 10,
          "N" : 11,
          "P" : 12,
          "Q" : 13,
          "R" : 14,
          "S" : 15,
          "T" : 16,
          "U" : 17,
          "V" : 18,
          "W" : 19,
          "Y" : 20}

def stream_fasta_fl(flname, kept_individuals=None):
    with open(flname) as fl:
        identifier = None
        sequence = ""
        for ln in fl:
            if ln.startswith(">"):
                if identifier != None and \
                   (kept_individuals is None or identifier in kept_individuals):
                    yield identifier, sequence

                identifier = ln[1:].strip()
                sequence = ""
            else:
                sequence += ln.strip()

        yield identifier, sequence

def extract_feature_vector(sequence, char_map, row):
    n = len(sequence) * len(char_map)
    features = []
    feature_labels = []
    seq_labels = [(idx, char) for char, idx in char_map.items()]
    for offset, char in enumerate(sequence):
        if char in char_map:
            idx = offset * len(char_map) + char_map[char]
            features.append((row, idx, 1.0))
        for idx, char in seq_labels:
            feature_labels.append([(offset + idx, char)])
    return features, feature_labels

def convert(groups_flname, fasta_flname, seq_type, outbase):
    char_map = None
    if seq_type == "DNA":
        char_map = DNA_MAP
    elif seq_type == "AA":
        char_map = AA_MAP
    else:
        raise NotImplementedError, "Unknown sequence type (%s)" % seq_type

    if groups_flname:
        groups = read_populations(groups_flname)
        kept_individuals = set(groups.keys())
    else:
        kept_individuals = None

    sample_labels = []
    rows = []
    cols = []
    values = []
    feature_labels = None
    stream = stream_fasta_fl(fasta_flname, kept_individuals)
    for row, (sample_label, sequence) in enumerate(stream):
        sample_features, feature_labels = extract_feature_vector(sequence, char_map, row)
        for r, c, v in sample_features:
            rows.append(r)
            cols.append(c)
            values.append(v)
        sample_labels.append(sample_label)


    feature_matrix = sparse.coo_matrix((values, (rows, cols))).tocsr()

    print feature_matrix.shape[0], "individuals", feature_matrix.shape[1], "features"

    np.savez(os.path.join(outbase, FEATURE_MATRIX_FLNAME),
             data = feature_matrix.data,
             indices = feature_matrix.indices,
             indptr = feature_matrix.indptr,
             shape = feature_matrix.shape)

    serialize(os.path.join(outbase, FEATURE_LABELS_FLNAME), feature_labels)
    serialize(os.path.join(outbase, SAMPLE_LABELS_FLNAME), sample_labels)
    serialize(os.path.join(outbase, FIXED_DIFFERENCES_FLNAME), dict())
    serialize(os.path.join(outbase, MISSING_DATA_FLNAME), dict())
