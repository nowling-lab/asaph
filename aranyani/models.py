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

import numpy as np
from scipy.stats import rankdata

from sklearn.ensemble import RandomForestClassifier

class SNPs(object):
    def __init__(self, labels, importances):
        self.labels = labels
        self.importances = importances

    def rank(self):
        sorted_indices = np.argsort(-1.0 * self.importances)
        sorted_importances = self.importances[sorted_indices]
        nonzero_indices = sorted_indices[sorted_importances != 0.0]
        nonzero_importances = sorted_importances[sorted_importances != 0.0]

        ranked_labels = [self.labels[idx] for idx in nonzero_indices]
        
        return SNPs(ranked_labels, nonzero_importances)

class Features(object):
    def __init__(self, feature_matrix, feature_labels, class_labels):
        self.feature_matrix = feature_matrix
        self.feature_labels = feature_labels
        self.class_labels = class_labels

    def snp_labels(self):
        snp_feature_indices = defaultdict(list)
        for feature_idx, feature_label in enumerate(self.feature_labels):
            snp_label = feature_label[:3]
            snp_feature_indices[snp_label].append(feature_idx)

        return snp_feature_indices

    def train_rf(self, n_trees):
        rf = RandomForestClassifier(n_estimators=n_trees)
        rf.fit(self.feature_matrix, self.class_labels)
        return rf

    def snp_importances(self, rf):
        feature_importances = rf.feature_importances_
        
        snp_labels = self.snp_labels()

        snp_importances = ( (np.mean(feature_importances[feature_idx]), snp_label)
                            for snp_label, feature_idx in snp_labels.iteritems() )

        snp_importances = sorted(snp_importances, reverse=True)

        labels = [label for importance, label in snp_importances]
        importances = np.array([importance for importance, label in snp_importances])
    
        return SNPs(labels, importances)
