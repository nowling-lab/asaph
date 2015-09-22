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

def create_feature_matrix(basename, n_individuals, n_features):
    """
    Creates a memory-mapped Numpy array.  On-disk file format
    stores the matrix dimensions for easy reading later.

    Returns the Numpy array.
    """
    
    flname = basename + ".matrix"

    # writer header
    fl = open(flname, "w")
    fl.write(struct.pack(FORMAT_STRING, n_individuals, n_features))
    fl.close()

    feature_matrix = np.memmap(filename=flname, dtype="float32", \
                               mode="r+", offset=HEADER_SIZE, \
                               shape=(n_individuals, n_features))

    return feature_matrix

def open_feature_matrix(basename):
    """
    Opens a memory-mapped Numpy array.  Reads the number of
    rows (individuals) and columns (features) from the header
    of the file and the matrix from the rest of the file. The
    Numpy array is opened in read-only mode.

    Returns a tuple of the number of individuals, number of 
    features, and the Numpy array."
    """
    
    flname = basename + ".matrix"

    # read header
    fl = open(flname)
    header = fl.read(HEADER_SIZE)
    n_individuals, n_features = struct.unpack(FORMAT_STRING, header)
    fl.close()

    feature_matrix = np.memmap(filename=flname, dtype="float32",
                                mode="r", offset=HEADER_SIZE, \
                                shape=(n_individuals, n_features))

    return feature_matrix
    
