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
import itertools
import random

import numpy as np
from scipy.stats import rankdata

from .ml import ConstrainedBaggingRandomForest

class SNPs(object):
    def __init__(self, labels, importances, ranked):
        # convert to tuple for hashing
        self.labels = map(tuple, labels)
        self.importances = importances
        self.ranked = ranked
        
    def rank(self):
        paired = zip(self.importances, self.labels)
        nonzero = filter(lambda t: t[0] != 0.0, paired)
        sorted_pairs = sorted(nonzero, reverse=True)
        
        sorted_labels = [label for _, label in sorted_pairs]
        sorted_importances = [importance for importance, _ in sorted_pairs]

        # convert to list for serialization
        return SNPs(list(sorted_labels),
                    list(sorted_importances),
                    True)

    def take(self, p1, p2=None):
        if p2 == None:
            start = 0
            end = p1
        else:
            start = p1
            end = p2
            
        if not self.ranked:
            ranked = self.rank()
            return SNPs(ranked.labels[start:end],
                        ranked.importances[start:end],
                        ranked.ranked)
        
        return SNPs(self.labels[start:end],
                    self.importances[start:end],
                    self.ranked)

    def __len__(self):
        return len(self.labels)
    
    def count_intersection(self, other):
        return len(set(self.labels).intersection(other.labels))

    def __repr__(self):
        return "SNPs(n_snps=%s, ranked=%s)" \
            % (len(self), self.ranked)

class Features(object):
    def __init__(self, feature_matrix, feature_names, class_labels, sample_labels):
        self.feature_matrix = feature_matrix
        self.feature_labels = feature_names
        self.class_labels = class_labels
        self.sample_labels = sample_labels

    def snp_labels(self):
        snp_feature_indices = defaultdict(list)
        for feature_idx, feature_labels in self.feature_labels.iteritems():
            feature_idx = int(feature_idx)
            for label in feature_labels:
                # chrom and pos only
                snp_label = label[:2]
                snp_feature_indices[snp_label].append(feature_idx)

        return snp_feature_indices

    def rank_snps(self, feature_importances):
        feature_importances = np.abs(feature_importances)
        
        snp_labels = self.snp_labels()
        snp_importances = [ (np.mean(feature_importances[feature_idx]), snp_label)
                            for snp_label, feature_idx in snp_labels.iteritems() ]

        snp_importances = list(snp_importances)

        labels = [label for _, label in snp_importances]
        importances = np.array([importance for importance, _ in snp_importances])

        return SNPs(labels, importances, False).rank()
        

