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

import random

import numpy as np
from scipy import sparse

from sklearn.feature_extraction import FeatureHasher
from sklearn.random_projection import SparseRandomProjection

from .feature_extraction import *

COUNTS_FEATURE_TYPE = "counts"
CATEGORIES_FEATURE_TYPE = "categories"
HASHED_FEATURE_TYPE = "hashed"

class Chunker(object):
    def __init__(self, chunk_size, n_samples, stream):
        self.chunk_size = chunk_size
        self.n_samples = n_samples
        self.buffer = [list() for i in range(n_samples)]
        self.chunk_count = 0
        self.processed_chunks = 0
        self.stream = stream

    def __iter__(self):
        for variant_label, feature_pairs in self.stream:
            if self.chunk_count >= self.chunk_size:
                self.processed_chunks += 1
                print("Processed chunk", self.processed_chunks)
                yield self.buffer
                self.buffer = [list() for i in range(self.n_samples)]
                self.chunk_count = 0

            self.chunk_count += 1
            for i, (sample_name, feature) in enumerate(feature_pairs):
                if feature != None:
                    self.buffer[i].append(feature)

        if self.chunk_count > 0:
            yield self.buffer

class FeatureHashingAccumulator(object):
    def __init__(self, n_features, n_samples):
        self.n_features = n_features
        self.n_samples = n_samples
        self.transformer = FeatureHasher(n_features = n_features,
                                         input_type = "pair")
        self.features = sparse.dok_matrix((n_samples, n_features), dtype=np.float32)

    def transform(self, stream):
        for chunk in stream:
            self.features += self.transformer.transform(chunk)
        return self.features.todense()

class FeatureHashingRandomProjectionAccumulator(object):
    def __init__(self, n_features, n_hashed_features, n_samples):
        self.n_features = n_features
        self.n_samples = n_samples
        
        self.hasher = FeatureHasher(n_features = n_hashed_features,
                                         input_type = "pair")

        # SRP needs the dimensions of the input matrix
        # to generate the projection matrix
        self.random_proj = SparseRandomProjection(n_features)
        dummy_sparse = sparse.csr_matrix((n_samples, n_hashed_features))
        self.random_proj.fit(dummy_sparse)

        self.features = np.zeros((n_samples, n_features), dtype=np.float32)

    def transform(self, stream):
        for chunk in stream:
            hashed_features = self.hasher.transform(chunk)
            self.features += self.random_proj.transform(hashed_features).toarray()
        return self.features
    
class FullMatrixAccumulator(object):
    def transform(self, stream):
        feature_columns = []
        for (chrom, pos, gt), column in stream:
            feature_columns.append(column)

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

        # need to transpose, otherwise we get (n_features, n_individuals) instead
        feature_matrix = np.array(feature_columns).T

        return feature_matrix

def construct_feature_matrix(variant_stream, n_samples, feature_type, subsampling_method, chunk_size, n_dim, n_inner_dim=None):
    print("Using feature type:", feature_type)
    print("Using subsampling method:", subsampling_method)
    if subsampling_method is not None or feature_type == HASHED_FEATURE_TYPE:
        print("Using", n_dim, "dimensions")
    
    # extract features
    if feature_type != HASHED_FEATURE_TYPE:
        if feature_type == COUNTS_FEATURE_TYPE:
            extractor = CountFeaturesExtractor(variant_stream)
        elif feature_type == CATEGORIES_FEATURE_TYPE:
            extractor = CategoricalFeaturesExtractor(variant_stream)
        else:
            raise Exception("Unknown feature type: %s" % feature_type)

        if subsampling_method is None:
            accumulator = FullMatrixAccumulator()
            feature_matrix = accumulator.transform(extractor)
        elif subsampling_method == "reservoir":
            accumulator = ReservoirMatrixAccumulator(n_dim)
            feature_matrix = accumulator.transform(extractor)
        else:
            raise Exception("Only reservoir sampling is supported with count or categorical features.")
    elif feature_type == HASHED_FEATURE_TYPE:
        string_features = FeatureStringsExtractor(variant_stream)
        chunker = Chunker(chunk_size, n_samples, string_features)

        if subsampling_method is None:
            accumulator = FeatureHashingAccumulator(n_dim, n_samples)
        elif subsampling_method == "random-projection":
            accumulator = FeatureHashingRandomProjectionAccumulator(n_dim, n_inner_dim, n_samples)
        else:
            raise Exception("Subsampling method '%s' not implemented for hashed features" % \
                            subsampling_method)
        
        feature_matrix = accumulator.transform(chunker)        
    else:
        s = "Combination of feature and subsampling types '%s' and '%s' not implemented" % \
            (feature_type, subsampling_method)
        raise Exception(s)

    return feature_matrix

