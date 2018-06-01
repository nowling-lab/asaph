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
import struct

import numpy as np
from scipy import sparse

from .models import *

FEATURE_LABELS_FLNAME = "feature_labels"
CLASS_LABELS_FLNAME = "class_labels"
SAMPLE_LABELS_FLNAME = "sample_labels"
FEATURE_MATRIX_FLNAME = "feature_matrix.npy"
SNP_FEATURE_INDICES_FLNAME = "snp_feature_indices"
SNP_FEATURE_GENOTYPES_FLNAME = "snp_feature_genotypes"
PROJECT_SUMMARY_FLNAME = "project_summary"

def serialize(flname, obj):
    fl = open(flname, "w")
    cPickle.dump(obj, fl)
    fl.close()

def deserialize(flname):
    fl = open(flname)
    obj = cPickle.load(fl)
    fl.close()

    return obj

def read_features(basename):
    snp_features_map = deserialize(os.path.join(basename, SNP_FEATURE_INDICES_FLNAME))
    class_labels = deserialize(os.path.join(basename, CLASS_LABELS_FLNAME))
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

    # to avoid breaking older formats
    genotypes_flname = os.path.join(basename, SNP_FEATURE_GENOTYPES_FLNAME)
    genotypes = None
    if os.path.exists(genotypes_flname):
        genotypes = deserialize(genotypes_flname)
    unknown_genotypes_flname = os.path.join(basename, UNKNOWN_GENOTYPES_FLNAME)
    unknown_genotypes = None
    if os.path.exists(unknown_genotypes_flname):
        unknown_genotypes = deserialize(unknown_genotypes_flname)

    return Features(feature_matrix,
                    snp_features_map,
                    class_labels,
                    sample_labels,
                    genotypes,
                    unknown_genotypes)

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
