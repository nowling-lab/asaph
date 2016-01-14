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
from exceptions import NotImplementedError
import json
import os
import struct
import tempfile
import unittest

import numpy as np

from ..vcf import *
from ..newioutils import *

VCF_TEST_FILE = os.path.join("test_data", "test.vcf")
GROUP_TEST_FILE = os.path.join("test_data", "groups")

class TestVCFFunctions(unittest.TestCase):
    def test_parse_vcf_line(self):
        # Data columns 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
        # individuals are columns after data columns
        idx = range(5)
        test_line = "\t".join(["1", "2", ".", "A", "T", "0", "PASS", ".", ".",
                               "0/0:12,0:12:33:0", "0/1:12,0:12:33:0", "1/0:12,0:12:33:0",
                               "1/1:12,0:12:33:0", "./.:12,0:12:33:0"])

        triplet = parse_vcf_line(test_line)

        chrom, pos, snps = triplet
        self.assertEqual(chrom, "1")
        self.assertEqual(pos, "2")                              
        self.assertRaises(NotImplementedError, parse_vcf_line, test_line.replace("0/0", "2/2"))
    
    def test_read_groups(self):
        groups = read_groups(GROUP_TEST_FILE)

        self.assertEquals(len(set(groups.values())), 2)
        self.assertEquals(len(groups.keys()), 16)
    
    def test_convert(self):
        dirname = tempfile.mkdtemp()

        flname = os.path.join(dirname, FEATURE_MATRIX_FLNAME)
        convert(GROUP_TEST_FILE, VCF_TEST_FILE, dirname, False)

        features = read_features(dirname)
        
        self.assertEquals(features.feature_matrix.shape[0], 16)
        self.assertEquals(features.feature_matrix.shape[1], 4)
        self.assertEquals(len(features.feature_labels), 4)
        self.assertEquals(len(features.sample_labels), 16)
        self.assertEquals(len(features.class_labels), 16)
        self.assertEquals(len(set(features.class_labels)), 2)

    def test_convert_compress(self):
        dirname = tempfile.mkdtemp()

        flname = os.path.join(dirname, FEATURE_MATRIX_FLNAME)
        convert(GROUP_TEST_FILE, VCF_TEST_FILE, dirname, True)

        features = read_features(dirname)
        
        self.assertEquals(features.feature_matrix.shape[0], 16)
        self.assertEquals(features.feature_matrix.shape[1], 3)
        self.assertEquals(len(features.feature_labels), 3)
        self.assertEquals(len(features.sample_labels), 16)
        self.assertEquals(len(features.class_labels), 16)
        self.assertEquals(len(set(features.class_labels)), 2)
        
