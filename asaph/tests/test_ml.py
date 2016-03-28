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

import unittest

import numpy as np

from ..ml import ConstrainedBaggingRandomForest

class TestConstrainedBaggingRandomForest(unittest.TestCase):
    features = np.array([[0, 1, 0, 1],
                             [1, 0, 1, 0],
                             [1, 0, 1, 0],
                             [0, 1, 0, 1]])
    class_labels = np.array([0, 0, 1, 1])
    n_trees = 10
    n_resamples = 3

    def test_resample(self):
        rf = ConstrainedBaggingRandomForest(self.n_trees, self.n_resamples)
        X, y = rf._resample(self.features, self.class_labels)

        self.assertEqual(y.shape[0], self.class_labels.shape[0] + self.n_resamples)
        self.assertEqual(X.shape[0], self.features.shape[0] + self.n_resamples)
        self.assertEqual(X.shape[1], self.features.shape[1])

    def test_feature_importances(self):
        rf = ConstrainedBaggingRandomForest(self.n_trees, self.n_resamples)
        fi = rf.feature_importances(self.features, self.class_labels)

        self.assertEqual(fi.shape[0], self.features.shape[1])
        
