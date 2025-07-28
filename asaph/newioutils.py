"""
This module provides input/output utilities for reading and writing population genetics data
including population labels, sample names, model serialization, and file format handling for the
Asaph analysis pipeline.

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

import pickle
from collections import OrderedDict
import os
import struct

import numpy as np
from scipy import sparse

from .models import *


SAMPLE_LABELS_FLNAME = "sample_labels"
MODEL_FLNAME = "model"
MODEL_KEY = "pca"
PROJECT_SUMMARY_FLNAME = "project_summary"
PROJECTION_KEY = "projected-coordinates"
FEATURES_FLNAME = "features"
COORDINATES_FLNAME = "pca_coordinates.tsv"

def read_populations(flname):
    """
    Read populations from a file.

    We assume that the file has one population per line.  The population name
    comes first and is followed by a comma.  The name of each sample follows,
    also separated by comma.

    For example:

    pop1,sample_1,sample_2
    pop2,sample_3,sample_4
    """
    with open(flname) as fl:
        groups = OrderedDict()
        group_names = OrderedDict()
        for group_id, ln in enumerate(fl):
            cols = ln.strip().split(",")
            for ident in cols[1:]:
                groups[ident] = group_id
                group_names[group_id] = cols[0]

    return groups, group_names


def serialize(flname, obj):
    with open(flname, "wb") as fl:
        pickle.dump(obj, fl)

def deserialize(flname):
    with open(flname, "rb") as fl:
        obj = pickle.load(fl)

    return obj

def read_sample_names(workdir):
    sample_labels = deserialize(os.path.join(workdir, SAMPLE_LABELS_FLNAME))
    return sample_labels
