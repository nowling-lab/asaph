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

from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier

class ConstrainedBaggingRandomForest(object):
    """
    Implementation of a Random Forest using constrained
    bagging for generating the training sets for the
    underlying Decision Trees.
    """

    def __init__(self, n_trees, n_resamples, batch_size):
        self.n_trees = n_trees
        self.n_resamples = n_resamples
        self.batch_size = batch_size

    def _resample(self, X, y):
        new_indices = list(xrange(X.shape[0]))
        for i in xrange(self.n_resamples):
            idx = random.randint(0, X.shape[0] - 1)
            new_indices.append(idx)

        X_new = X[new_indices, :]
        y_new = y[new_indices]

        return X_new, y_new

    def feature_importances(self, X, y):
        y = np.array(y)
        feature_importances = np.zeros(X.shape[1])
        if self.n_resamples == -1:
            completed_trees = 0
            batch_size = 100
            while completed_trees < self.n_trees:
                n_classifiers = min(batch_size, self.n_trees - completed_trees)
                rf = RandomForestClassifier(n_estimators=n_classifiers,
                                            n_jobs=-1)
                rf.fit(X, y)
                feature_importances += rf.feature_importances_ * n_classifiers
                completed_trees += n_classifiers
            else:
                for i in xrange(self.n_trees):
                    dt = DecisionTreeClassifier(max_features="sqrt")
                    X_new, y_new = self._resample(X, y)
                    dt.fit(X_new, y_new)
                    feature_importances += dt.feature_importances_

        feature_importances = feature_importances / self.n_trees
        return feature_importances
