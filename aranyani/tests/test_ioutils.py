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
