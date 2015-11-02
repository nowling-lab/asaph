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
from ..aranyani import *

class MockRF(object):
    def __init__(self, feature_importances):
        self.feature_importances_ = feature_importances

class TestFeatures(unittest.TestCase):
    feature_labels = [(1,1,1,"A"),
                      (1,1,1,"T"),
                      (1,1,2,"G"),
                      (1,1,2,"C")]
    
    def test_snps_labels(self):
        features = Features(None, self.feature_labels)
        snp_labels = features.snp_labels()

        self.assertIn((1,1,1), snp_labels)
        self.assertIn((1,1,2), snp_labels)
        self.assertIn(0, snp_labels[(1,1,1)])
        self.assertIn(1, snp_labels[(1,1,1)])
        self.assertIn(2, snp_labels[(1,1,2)])
        self.assertIn(3, snp_labels[(1,1,2)])

    def test_snp_importances(self):
        features = Features(None, self.feature_labels)
        importances = np.array([0.0, 1.0, 0.5, 1.0])
        rf = MockRF(importances)

        labels, snp_importances = features.snp_importances(rf)

        self.assertEqual(labels[0], (1, 1, 2))
        self.assertEqual(labels[1], (1, 1, 1))
        self.assertAlmostEqual(snp_importances[0], 0.75)
        self.assertAlmostEqual(snp_importances[1], 0.5)
