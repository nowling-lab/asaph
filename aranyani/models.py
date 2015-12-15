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

from sklearn.cross_validation import LeaveOneOut
from sklearn.ensemble import RandomForestClassifier

class SNPs(object):
    def __init__(self, n_trees, labels, importances, ranked, fixed_differences,
                 missing_data):
        self.n_trees = n_trees
        # convert to tuple for hashing
        self.labels = map(tuple, labels)
        self.importances = importances
        self.ranked = ranked
        self.fixed_differences = fixed_differences
        self.missing_data = missing_data
        
    def rank(self):
        paired = zip(self.importances, self.labels, self.fixed_differences,
                     self.missing_data)
        nonzero = filter(lambda t: t[0] != 0.0, paired)
        sorted_pairs = sorted(nonzero, reverse=True)
        
        sorted_labels = [label for _, label, _, _ in sorted_pairs]
        sorted_importances = [importance for importance, _, _, _ in sorted_pairs]
        sorted_fixed_differences = [fd for _, _, fd, _ in sorted_pairs]
        sorted_missing_data = [md for _, _, _, md in sorted_pairs]

        # convert to list for serialization
        return SNPs(self.n_trees, list(sorted_labels), list(sorted_importances),
                    True, sorted_fixed_differences, sorted_missing_data)

    def take(self, p1, p2=None):
        if p2 == None:
            start = 0
            end = p1
        else:
            start = p1
            end = p2
            
        if not self.ranked:
            ranked = self.rank()
            return SNPs(ranked.n_trees, ranked.labels[start:end],
                        ranked.importances[start:end], ranked.ranked,
                        ranked.fixed_differences, ranked.missing_data)
        
        return SNPs(self.n_trees, self.labels[start:end],
                    self.importances[start:end], self.ranked,
                    self.fixed_differences, self.missing_data)

    def __len__(self):
        return len(self.labels)
    
    def count_intersection(self, other):
        return len(set(self.labels).intersection(other.labels))

    def __repr__(self):
        return "SNPs(n_trees=%s, n_snps=%s, ranked=%s)" \
            % (self.n_trees, len(self), self.ranked)


class Features(object):
    def __init__(self, feature_matrix, feature_labels, class_labels, sample_labels,
                 fixed_differences, missing_data):
        self.feature_matrix = feature_matrix
        # convert to tuple for hashing
        self.feature_labels = [map(tuple, labels) for labels in feature_labels]
        self.class_labels = class_labels
        self.sample_labels = sample_labels
        self.fixed_differences = fixed_differences
        self.missing_data = missing_data

    def snp_labels(self):
        snp_feature_indices = defaultdict(list)
        for feature_idx, feature_labels in enumerate(self.feature_labels):
            for label in feature_labels:
                # chrom and pos
                snp_label = label[:2]
                snp_feature_indices[snp_label].append(feature_idx)

        return snp_feature_indices

    def train_rf(self, n_trees):
        rf = RandomForestClassifier(n_estimators=n_trees)
        rf.fit(self.feature_matrix, self.class_labels)
        return rf

    def snp_importances(self, n_trees, max_batch_size):
        remaining_trees = n_trees
        feature_importances = np.zeros(self.feature_matrix.shape[1])
        while remaining_trees > 0:
            batch_size = min(remaining_trees, max_batch_size)
            rf = self.train_rf(batch_size)
            feature_importances += batch_size * rf.feature_importances_
            remaining_trees -= batch_size

        feature_importances = feature_importances / n_trees
        
        snp_labels = self.snp_labels()

        # feature_idx is a list of associated features
        # FD and Missing are the same for all features of the same SNP
        snp_importances = ( (np.mean(feature_importances[feature_idx]), snp_label,
                             self.fixed_differences[snp_label],
                             self.missing_data[snp_label])
                            for snp_label, feature_idx in snp_labels.iteritems() )

        snp_importances = sorted(snp_importances, reverse=True)

        labels = [label for _, label, _, _ in snp_importances]
        importances = np.array([importance for importance, _, _, _ in snp_importances])
        fixed_differences = [fd for _, _, fd, _ in snp_importances]
        missing_data = [md for _, _, _, md in snp_importances]
    
        return SNPs(n_trees, labels, importances, False, fixed_differences, \
                    missing_data)

    def validate_rf(self, n_trees):
        score = 0.0
        for train_index, test_index in LeaveOneOut(len(self.class_labels)):
            X_train = self.feature_matrix[train_index]
            X_test = self.feature_matrix[test_index]
            y_train = np.array(self.class_labels)[train_index]
            y_test = np.array(self.class_labels)[test_index]
            
            rf = RandomForestClassifier(n_estimators=n_trees)
            rf.fit(X_train, y_train)
            score += rf.score(X_test, y_test)
            
        return score / len(self.class_labels)

    def select_from_snps(self, snps):
        snp_labels = set(snps.labels)
        selected_indices = []
        selected_feature_labels = []
        selected_fds = []
        selected_missing = []
        
        for i, (chrom, pos, nucleotide) in enumerate(self.feature_labels):
            if (chrom, pos) in snp_labels:
                selected_indices.append(i)
                selected_feature_labels.append(self.feature_labels[i])
                selected_fds.append(self.fixed_differences[i])
                selected_missing.append(self.missing_data[i])

        return Features(self.feature_matrix[:, selected_indices],
                        selected_feature_labels,
                        self.class_labels,
                        self.sample_labels,
                        selected_fds,
                        selected_missing)

