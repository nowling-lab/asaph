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

import numpy as np
import struct

FORMAT_STRING = "ii"
HEADER_SIZE = struct.calcsize(FORMAT_STRING) # bytes

def create_feature_matrix(basename, n_individuals, n_features, filtered=False):
    """
    Creates a memory-mapped Numpy array.  On-disk file format
    stores the matrix dimensions for easy reading later.

    Returns the Numpy array.
    """

    if not filtered:
        flname = basename + ".matrix"
    else:
        flname = basename + ".filtered.matrix"

    # writer header
    fl = open(flname, "w")
    fl.write(struct.pack(FORMAT_STRING, n_individuals, n_features))
    fl.close()

    feature_matrix = np.memmap(filename=flname, dtype="float32", \
                               mode="r+", offset=HEADER_SIZE, \
                               shape=(n_individuals, n_features))

    return feature_matrix

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
        flname = basename + ".matrix"
    else:
        flname = basename + ".filtered.matrix"

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
    fl = open(inbase + ".groups")
    for ln in fl:
        cols = ln.strip().split(",")
        group = cols[0]
        ids = cols[1:]
        groups[group] = ids
    fl.close()
    return groups

def read_individuals(inbase, filtered=False):
    if not filtered:
        flname = inbase + ".individuals"
    else:
        flname = inbase + ".filtered.individuals"
    
    fl = open(flname)
    individuals = [ln.strip() for ln in fl]
    fl.close()
    return individuals

def read_features(inbase, filtered=False):
    if not filtered:
        flname = inbase + ".features"
    else:
        flname = inbase + ".filtered.features"
    
    fl = open(flname)
    features = fl.readlines()
    fl.close()
    return features

def write_features(basename, features, filtered=False):
    if not filtered:
        flname = basename + ".features"
    else:
        flname = basename + ".filtered.features"

    fl = open(flname, "w")
    for feature in features:
        fl.write(feature)
        fl.write("\n")
    fl.close()
    
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

def write_individuals(basename, labels, filtered=False):
    if not filtered:
        flname = basename + ".individuals"
    else:
        flname = basename + ".filtered.individuals"

    fl = open(flname, "w")
    for label in labels:
        fl.write(label)
        fl.write("\n")
    fl.close()

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
