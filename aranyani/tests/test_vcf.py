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

from ..vcf import *

VCF_TEST_FILE = os.path.join("test_data", "test.vcf")
GROUP_TEST_FILE = os.path.join("test_data", "groups")

class TestVCFFunctions(unittest.TestCase):
    def test_vcf_line_to_seq(self):
        # Data columns 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
        # individuals are columns after data columns
        test_line = "\t".join(["1", "2", ".", "A", "T", "0", "PASS", ".", ".",
                               "0/0:12,0:12:33:0", "0/1:12,0:12:33:0", "1/0:12,0:12:33:0",
                               "1/1:12,0:12:33:0", "./.:12,0:12:33:0"])

        chrom, pos, haploid1, haploid2 = vcf_line_to_seq(test_line)

        self.assertEqual(chrom, "1")
        self.assertEqual(pos, "2")
        self.assertEqual(haploid1, ("A", "A", "T", "T", "X"))
        self.assertEqual(haploid2, ("A", "T", "A", "T", "X"))
    
    def test_read_groups(self):
        groups = read_groups(GROUP_TEST_FILE)

        self.assertIn("_all", groups)
        self.assertEquals(len(groups), 3)
        self.assertEquals(len(groups["group1"]), 8)
        self.assertEquals(len(groups["group2"]), 8)
        self.assertEquals(len(groups["_all"]), 16)
    
    def test_read_dimensions(self):
        sample_ids, n_samples, n_features = read_dimensions(VCF_TEST_FILE)

        print n_samples, n_features
