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

import heapq
import random

import mmh3

import numpy as np

from sklearn.feature_extraction import FeatureHasher

from .feature_extraction import *

COUNTS_FEATURE_TYPE = "allele-counts"
CATEGORIES_FEATURE_TYPE = "genotype-categories"

RESERVOIR_SAMPLING = "reservoir"
FEATURE_HASHING = "feature-hashing"
BOTTOMK_SKETCHING = "bottom-k"

class FeatureHashingAccumulator(object):
    def __init__(self, n_features, n_samples):
        self.n_features = n_features
        self.n_samples = n_samples

    def transform(self, stream):
        feature_columns = dict()
        for feature_idx, ((chrom, pos, gt), column) in enumerate(stream, start=1):
            feature_name = "{}_{}_{}".format(chrom, pos, gt)

            # this will cause collisions.  that's okay -- we want that.
            hash_ = abs(mmh3.hash(feature_name)) % self.n_features

            # allows for sparse storage
            if hash_ in feature_columns:
                feature_columns[hash_] += np.array(column)
            else:
                feature_columns[hash_] = np.array(column)

            if feature_idx % 10000 == 0:
                print("Chunk", feature_idx // 10000, len(feature_columns))

        # need to transpose, otherwise we get (n_features, n_individuals) instead
        feature_matrix = np.array(list(feature_columns.values())).T

        print(feature_matrix.shape)

        return feature_matrix

class BottomKAccumulator(object):
    """
    Online sampling of columns using bottom-k sketching
    """
    def __init__(self, n_features):
        self.n_features = n_features

    def transform(self, stream):
        feature_columns = []
        for feature_idx, ((chrom, pos, gt), column) in enumerate(stream):
            feature_name = "{}_{}_{}".format(chrom, pos, gt)

            # Python's built-in heap is a min heap, so it
            # will keep the largest elements.  In practice,
            # this probably doesn't matter but we are going
            # to negate the hashes anyway so it keeps the
            # smallest elements.  Note that for signed ints,
            # 0 replaces one of the positive values so we can
            # convert any positive value to negative without an
            # overflow but not the other way around
            # Also, mmh3.hash returns a 32-bit signed int
            hash_ = abs(mmh3.hash(feature_name))

            # we use the feature_idx to break ties since numpy arrays
            # are not comparable (sortable)
            if len(feature_columns) < self.n_features:
                heapq.heappush(feature_columns,
                               (hash_, feature_idx, column))
            else:
                heapq.heapreplace(feature_columns,
                                  (hash_, feature_idx, column))

            if feature_idx % 10000 == 0:
                print("Chunk", feature_idx // 10000, len(feature_columns))

        # drop the hash and feature idx
        feature_columns = [column for _, _, column in feature_columns]

        # need to transpose, otherwise we get (n_features, n_individuals) instead
        feature_matrix = np.array(feature_columns).T

        return feature_matrix
    
class FullMatrixAccumulator(object):
    def transform(self, stream):
        feature_columns = []
        for feature_idx, ((chrom, pos, gt), column) in enumerate(stream, start=1):
            feature_columns.append(column)

            if feature_idx % 10000 == 0:
                print("Chunk", feature_idx // 10000, len(feature_columns))


        # need to transpose, otherwise we get (n_features, n_individuals) instead
        feature_matrix = np.array(feature_columns).T

        return feature_matrix

class ReservoirMatrixAccumulator(object):
    """
    Online sampling of columns using reservoir sampling.
    """
    def __init__(self, n_features):
        self.n_features = n_features

    def transform(self, stream):
        feature_columns = []
        for feature_idx, ((chrom, pos, gt), column) in enumerate(stream):
            if feature_idx < self.n_features:
                feature_columns.append(column)
            else:
                j = random.randint(0, feature_idx + 1)
                if j < self.n_features:
                    feature_columns[j] = column

            if feature_idx % 10000 == 0:
                print("Chunk", feature_idx // 10000, len(feature_columns))
                    
        # need to transpose, otherwise we get (n_features, n_individuals) instead
        feature_matrix = np.array(feature_columns).T

        return feature_matrix

def construct_feature_matrix(variant_stream, n_samples, feature_type, sampling_method, n_dim):
    print("Using feature type:", feature_type)
    if sampling_method is not None:
        print("Using sampling method:", sampling_method)
        print("Using", n_dim, "dimensions")

    if feature_type == COUNTS_FEATURE_TYPE:
        extractor = CountFeaturesExtractor(variant_stream)
    elif feature_type == CATEGORIES_FEATURE_TYPE:
        extractor = CategoricalFeaturesExtractor(variant_stream)
    else:
        raise Exception("Unknown feature type: %s" % feature_type)

    if sampling_method is None:
        accumulator = FullMatrixAccumulator()
    elif sampling_method == RESERVOIR_SAMPLING:
        accumulator = ReservoirMatrixAccumulator(n_dim)
    elif sampling_method == FEATURE_HASHING:
        accumulator = FeatureHashingAccumulator(n_dim, n_samples)
    elif sampling_method == BOTTOMK_SKETCHING:
        accumulator = BottomKAccumulator(n_dim)
    else:
        raise Exception("Sampling method '%s' not implemented" % \
                            sampling_method)

    feature_matrix = accumulator.transform(extractor)

    return feature_matrix

