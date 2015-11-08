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
import json
import os
import random
import struct
import tempfile
import unittest

import numpy as np

from ..models import *
from ..newioutils import *

def create_features_dataset(basename):
    if not os.path.exists(basename):
        os.mkdir(basename)

    feature_labels = [(1,1,1,"A"),
                      (1,1,1,"T"),
                      (1,1,2,"G"),
                      (1,1,2,"C"),
                      (1,2,1,"A"),
                      (1,2,1,"T"),
                      (1,2,2,"G"),
                      (1,2,2,"C")]

    sample_labels = ["sample1",
                     "sample2",
                     "sample3",
                     "sample4"]

    class_labels = [0, 0, 1, 1]

    feature_matrix = create_feature_matrix(os.path.join(basename, FEATURE_MATRIX_FLNAME),
                                           len(class_labels), len(feature_labels))

    feature_matrix = np.array([[0, 1, 0, 1],
                               [1, 0, 1, 0],
                               [1, 0, 1, 0],
                               [0, 1, 0, 1]])

    del feature_matrix

    to_json(os.path.join(basename, FEATURE_LABELS_FLNAME), feature_labels)
    to_json(os.path.join(basename, SAMPLE_LABELS_FLNAME), sample_labels)
    to_json(os.path.join(basename, CLASS_LABELS_FLNAME), class_labels)

class TestIOUtils(unittest.TestCase):
    def test_read_features(self):
        dirname = tempfile.mkdtemp()

        create_features_dataset(dirname)

        features = read_features(dirname)

        self.assertEqual(features.feature_matrix.shape[0],
                         len(features.class_labels))
        self.assertEqual(features.feature_matrix.shape[1],
                         len(features.feature_labels))
        self.assertEqual(features.feature_matrix.shape[0],
                         len(features.sample_labels))

    def test_write_snps(self):
        dirname = tempfile.mkdtemp()

        labels = [[l] for l in range(10)]
        importances = [random.random() for i in xrange(10)]
        snps = SNPs(10, labels, importances, False)

        model_dir = os.path.join(dirname, "models", "10")
        self.assertFalse(os.path.exists(model_dir))
        
        write_snps(dirname, snps, "1")

        self.assertTrue(os.path.exists(model_dir))
        self.assertTrue(os.path.exists(os.path.join(model_dir, "1")))

    def test_read_snps(self):
        dirname = tempfile.mkdtemp()

        tree_sizes = [10, 25, 50]
        
        labels = [[l] for l in range(10)]
        for n_trees in tree_sizes:
            importances = [random.random() for i in xrange(10)]
            snps = SNPs(n_trees, labels, importances, False)
            write_snps(dirname, snps, "1")
            importances = [random.random() for i in xrange(10)]
            snps = SNPs(n_trees, labels, importances, False)
            write_snps(dirname, snps, "2")

        all_snps = read_snps(dirname)

        self.assertEqual(len(all_snps), len(tree_sizes))
        self.assertItemsEqual(all_snps.keys(), tree_sizes)

        for n_trees in tree_sizes:
            self.assertEqual(len(all_snps[n_trees]), 2)
            
