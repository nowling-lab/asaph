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
import unittest

import numpy as np
from ..models import *

class TestFeatures(unittest.TestCase):
    feature_labels = [[(1,1,"A")],
                      [(1,1,"T")],
                      [(1,2,"G")],
                      [(1,2,"C")]]
    
    def test_snps_labels(self):
        features = Features(None, self.feature_labels, None, None, None, None)
        snp_labels = features.snp_labels()

        self.assertIn((1,1), snp_labels)
        self.assertIn((1,2), snp_labels)
        self.assertIn(0, snp_labels[(1,1)])
        self.assertIn(1, snp_labels[(1,1)])
        self.assertIn(2, snp_labels[(1,2)])
        self.assertIn(3, snp_labels[(1,2)])

    def test_snp_importances(self):
        n_trees = 10
        batch_size = 3
        features = np.array([[0, 1, 0, 1],
                             [1, 0, 1, 0],
                             [1, 0, 1, 0],
                             [0, 1, 0, 1]])
        class_labels = [0, 0, 1, 1]
        sample_labels = [0, 1, 2, 3]
        fixed_differences = {(1, 1) : False, (1, 2) : False}
        missing_data = {(1, 1) : False, (1, 2) : False}
        
        features = Features(features, self.feature_labels, class_labels, sample_labels,
                            fixed_differences, missing_data)

        snps = features.snp_importances(n_trees)

        self.assertEqual(snps.n_trees, n_trees)

class TestSNPs(unittest.TestCase):
    def test_rank(self):
        labels = [[(1, 1)], [(1, 2)], [(1, 3)]]
        importances = np.array([0.5, 0.75, 0.0])
        fixed_differences = {(1, 1) : False, (1, 2) : False}
        missing_data = {(1, 1) : False, (1, 2) : False}

        snps = SNPs(None, labels, importances, False, fixed_differences, missing_data)
        ranked_snps = snps.rank()

        self.assertEqual(ranked_snps.labels, [((1, 2),), ((1, 1),)])
        self.assertAlmostEqual(ranked_snps.importances[0], 0.75)
        self.assertAlmostEqual(ranked_snps.importances[1], 0.5)
        self.assertTrue(ranked_snps.ranked)
        self.assertFalse(snps.ranked)
