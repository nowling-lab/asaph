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

from ..fasta import *
from ..newioutils import *

class TestFASTAFunctions(unittest.TestCase):
    def extract_feature_vector(self):
        sequence = "AWY"
        expected_features = [(0, 0, 1.0),
                             (0, 21 + 19, 1.0),
                             (0, 42 + 20, 1.0)]

        features, labels = extract_feature_vector(sequence, AA_MAP, 0)

        self.assertEqual(expected_features, features)
        self.assertEqual(len(labels), len(AA_MAP) * 3)
        
