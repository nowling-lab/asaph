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

import numpy as np

class CountFeaturesExtractor(object):
    def __init__(self, stream):
        self.stream = stream

    def __iter__(self):
        for variant_label, alleles, genotypes in self.stream:
            chrom, pos = variant_label
            ref_column = [0.] * len(genotypes)
            alt_column = [0.] * len(genotypes)

            for row_idx, (sample_name, allele_counts) in enumerate(genotypes):
                ref_column[row_idx] = allele_counts[0]
                alt_column[row_idx] = allele_counts[1]

            yield (chrom, pos, alleles[0]), tuple(ref_column)
            yield (chrom, pos, alleles[1]), tuple(alt_column)

class CategoricalFeaturesExtractor(object):
    def __init__(self, stream):
        self.stream = stream

    def __iter__(self):
        for variant_label, alleles, genotypes in self.stream:
            chrom, pos = variant_label
            homo1_column = [0.] * len(genotypes)
            homo2_column = [0.] * len(genotypes)
            het_column = [0.] * len(genotypes)

            for row_idx, (sample_name, allele_counts) in enumerate(genotypes):
                if allele_counts == (2, 0):
                    homo1_column[row_idx] = 1
                elif allele_counts == (0, 2):
                    homo2_column[row_idx] = 1
                elif allele_counts == (1, 1):
                    het_column[row_idx] = 1

            yield (chrom, pos, (alleles[0] + "/" + alleles[0])), tuple(homo1_column)
            yield (chrom, pos, (alleles[1] + "/" + alleles[1])), tuple(homo2_column)
            yield (chrom, pos, (alleles[0] + "/" + alleles[1])), tuple(het_column)

class FeatureStringsExtractor(object):
    def __init__(self, stream):
        self.stream = stream

    def __iter__(self):
        for variant_label, alleles, genotypes in self.stream:
            string_features = [None] * len(genotypes)

            homo_ref = ("%s_%s_het_%s" % (variant_label[0],
                                          variant_label[1],
                                          alleles[0]),
                        1)

            homo_alt = ("%s_%s_het_%s" % (variant_label[0],
                                          variant_label[1],
                                          alleles[1]),
                        1)

            het = ("%s_%s_%s_%s" % (variant_label[0],
                                    variant_label[1],
                                    alleles[0],
                                    alleles[1]),
                   1)

            for i, (sample_name, allele_counts) in enumerate(genotypes):
                if allele_counts == (2, 0):
                    string_features[i] = (sample_name, homo_ref)
                elif allele_counts == (0, 2):
                    string_features[i] = (sample_name, homo_alt)
                elif allele_counts == (1, 1):
                    string_features[i] = (sample_name, het)
                else:
                    string_features[i] = (sample_name, None)

            yield variant_label, string_features
