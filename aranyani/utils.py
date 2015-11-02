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

from collections import defaultdict
import os
import struct

import numpy as np

FORMAT_STRING = "ii"
HEADER_SIZE = struct.calcsize(FORMAT_STRING) # bytes


def open_feature_matrix(basename, filtered=False):
    """
    Opens a memory-mapped Numpy array.  Reads the number of
    rows (individuals) and columns (features) from the header
    of the file and the matrix from the rest of the file. The
    Numpy array is opened in read-only mode.

    Returns a tuple of the number of individuals, number of 
    features, and the Numpy array."
    """

    if not filtered:
        flname = basename + os.sep + "matrix"
    else:
        flname = basename + os.sep + "filtered_matrix"

    # read header
    fl = open(flname)
    header = fl.read(HEADER_SIZE)
    n_individuals, n_features = struct.unpack(FORMAT_STRING, header)
    fl.close()

    feature_matrix = np.memmap(filename=flname, dtype="float32",
                                mode="r", offset=HEADER_SIZE, \
                                shape=(n_individuals, n_features))

    return feature_matrix


def read_groups(inbase):
    groups = dict()
    fl = open(inbase + os.sep + "groups")
    for ln in fl:
        cols = ln.strip().split(",")
        group = cols[0]
        ids = cols[1:]
        groups[group] = ids
    fl.close()
    return groups

def get_group_indices(inbase, filtered):
    groups = read_groups(inbase)
    individuals = read_individuals(inbase, filtered=filtered)

    leftover_rows = set(xrange(len(individuals)))
    groups_rows = []
    for group_name, ids in groups.iteritems():
        group_rows = tuple([i for i, id in enumerate(individuals) if id in ids])
        groups_rows.append((group_name, group_rows))
        leftover_rows.difference_update(group_rows)

    return groups_rows, tuple(leftover_rows)

def read_features(inbase, filtered=False):
    if not filtered:
        flname = inbase + os.sep + "features"
    else:
        flname = inbase + os.sep + "filtered_features"
    
    fl = open(flname)
    features = [line.strip() for line in fl.readlines() if line.strip() != ""]
    fl.close()
    return features

def read_individuals(inbase, filtered=False):
    if not filtered:
        flname = inbase + os.sep + "individuals"
    else:
        flname = inbase + os.sep + "filtered_individuals"
    
    fl = open(flname)
    individuals = [ln.strip() for ln in fl]
    fl.close()
    return individuals

def select_groups(inbase, filtered):
    group_indices, _ = get_group_indices(inbase, filtered=filtered)
    individuals = read_individuals(inbase, filtered=filtered)

    labels = []
    selected_indices = []
    selected_individuals = []
    for label, (group_name, indices) in enumerate(group_indices):
        for j, idx in enumerate(indices):
            labels.append(label)
            selected_indices.append(idx)
            selected_individuals.append(individuals[idx])

    return labels, selected_indices, selected_individuals
