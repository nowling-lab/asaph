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

import cPickle
from collections import defaultdict
import glob
import os
import shelve
import struct

import numpy as np
from scipy import sparse

from .models import *

FEATURE_LABELS_FLNAME = "feature_labels"
CLASS_LABELS_FLNAME = "class_labels"
SAMPLE_LABELS_FLNAME = "sample_labels"
FEATURE_MATRIX_FLNAME = "feature_matrix.npy"

def to_json(flname, obj):
    fl = open(flname, "w")
    cPickle.dump(obj, fl)
    fl.close()

def from_json(flname):
    fl = open(flname)
    obj = cPickle.load(fl)
    fl.close()

    return obj

def read_features(basename):
    feature_labels = shelve.open(os.path.join(basename, FEATURE_LABELS_FLNAME))
    class_labels = from_json(os.path.join(basename, CLASS_LABELS_FLNAME))
    sample_labels = from_json(os.path.join(basename, SAMPLE_LABELS_FLNAME))
    if os.path.exists(os.path.join(basename, FEATURE_MATRIX_FLNAME + ".npz")):
        loader = np.load(os.path.join(basename, FEATURE_MATRIX_FLNAME + ".npz"))
        feature_matrix = sparse.csr_matrix((loader["data"],
                                           loader["indices"],
                                           loader["indptr"]),
                                           shape = loader["shape"])
    else:
        feature_matrix = np.load(os.path.join(basename, FEATURE_MATRIX_FLNAME))

    return Features(feature_matrix,
                    feature_labels,
                    class_labels,
                    sample_labels)

def write_snps(basedir, snps, model_id):
    model_dir = os.path.join(basedir, "models", str(snps.n_trees))
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)

    flname = os.path.join(model_dir, model_id)
    to_json(flname, vars(snps))

    
def read_snps(basedir):
    model_base_dir = os.path.join(basedir, "models")
    if not os.path.exists(model_base_dir):
        return dict()

    model_tree_dirs = glob.glob(os.path.join(model_base_dir, "*"))

    models = defaultdict(list)
    for model_dir in model_tree_dirs:
        model_flnames = glob.glob(os.path.join(model_dir, "*"))

        for flname in model_flnames:
            model_data = from_json(flname)

            snps = SNPs(model_data["n_trees"],
                        model_data["labels"],
                        model_data["importances"],
                        model_data["ranked"])

            models[snps.n_trees].append(snps)

    return models
