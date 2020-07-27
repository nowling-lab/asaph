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

import pickle
from collections import defaultdict
from collections import OrderedDict
import glob
import os
import struct

import numpy as np
from scipy import sparse

from .models import *

SAMPLE_LABELS_FLNAME = "sample_labels"
FEATURE_MATRIX_FLNAME = "feature_matrix.npy"
PROJECT_SUMMARY_FLNAME = "project_summary"

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
    fl = open(flname, "w")
    pickle.dump(obj, fl)
    fl.close()

def deserialize(flname):
    fl = open(flname)
    obj = pickle.load(fl)
    fl.close()

    return obj

def read_features(basename):
    sample_labels = deserialize(os.path.join(basename, SAMPLE_LABELS_FLNAME))
    if os.path.exists(os.path.join(basename, FEATURE_MATRIX_FLNAME + ".npz")):
        with np.load(os.path.join(basename, FEATURE_MATRIX_FLNAME + ".npz")) as loader:
            if "feature_matrix" in loader:
                feature_matrix = loader["feature_matrix"]
            else:
                feature_matrix = sparse.csr_matrix((loader["data"],
                                                    loader["indices"],
                                                    loader["indptr"]),
                                                   shape = loader["shape"])
    else:
        feature_matrix = np.load(os.path.join(basename, FEATURE_MATRIX_FLNAME))

    return Features(feature_matrix,
                    sample_labels)

def write_rf_snps(basedir, snps, n_trees, model_id):
    model_dir = os.path.join(basedir, "models", "rf", str(n_trees))
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    flname = os.path.join(model_dir, model_id)
    serialize(flname, snps)


def read_rf_snps(basedir):
    model_base_dir = os.path.join(basedir, "models", "rf")
    if not os.path.exists(model_base_dir):
        return dict()

    model_tree_dirs = glob.glob(os.path.join(model_base_dir, "*"))

    models = defaultdict(list)
    for model_dir in model_tree_dirs:
        model_flnames = glob.glob(os.path.join(model_dir, "*"))

        for flname in model_flnames:
            if not os.path.basename(flname).startswith("model"):
                continue
            snps = deserialize(flname)
            n_trees = int(os.path.basename(os.path.dirname(flname)))
            models[n_trees].append(snps)

    return models
