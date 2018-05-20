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
import random

import numpy as np
import numpy.linalg as la

from scipy.stats import chi2

from sklearn.decomposition import TruncatedSVD
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import log_loss
from sklearn.preprocessing import LabelEncoder
from sklearn.tree import DecisionTreeClassifier

def estimate_lr_iter(n_samples):
    return max(20,
               int(np.ceil(25000. / n_samples)))

def calculate_intercept_(labels):
    n_class_1 = float(sum(labels))
    n = float(len(labels))
    prob = n_class_1 / n

    # p = 1 / (1 + e^{-b})
    # (1 + e^{-b})p = 1
    # 1 + e^{-b} = 1/p
    # e^{-b} = 1/p - 1
    # -b = log(1/p - 1)
    # b = -log(1/p - 1)

    # clip prob to avoid numerical errors
    prob = min(prob, 0.9999999)
    prob = max(0.0000001, prob)
    
    return - np.log(1. / prob - 1.)

def calculate_null_model(n_samples, intercepts):
    n_classes = len(intercepts)
    null_prob = np.zeros((n_samples, n_classes))
    for i in xrange(n_samples):
        for j in xrange(n_classes):
            null_prob[i, j] = 1.0 / (1.0 + np.exp(-intercepts[j]))

    return null_prob


def calculate_intercept(labels):
    label_ids = sorted(set(labels))

    intercepts = np.zeros(len(label_ids))
    for i, pos_label in enumerate(label_ids):
        binary_labels = [1 if l == pos_label else 0 for l in labels]
        binary_labels = np.array(binary_labels)
        intercepts[i] = calculate_intercept_(binary_labels)

    return intercepts

def likelihood_ratio_test(features_alternate, labels, lr_model, features_null=None, set_intercept=True, g_scaling_factor=1.0):
    if isinstance(features_alternate, tuple) and len(features_alternate) == 2:
        training_features_alternate, testing_features_alternate = features_alternate
        training_labels, testing_labels = labels
    else:
        training_features_alternate = features_alternate
        testing_features_alternate = features_alternate
        training_labels = labels
        testing_labels = labels

    # re-orders labels from smallest to biggest
    # and makes labels consecutive
    # so that intercepts match
    label_encoder = LabelEncoder()
    training_labels = label_encoder.fit_transform(training_labels)
    testing_labels = label_encoder.transform(testing_labels)

    n_classes = len(set(training_labels))

    if set_intercept:
        if n_classes == 2:
            intercept_init = calculate_intercept(training_labels)[-1]
        else:
            intercept_init = calculate_intercept(training_labels)
    else:
        intercept_init = None

    if features_null is not None:
        if isinstance(features_null, tuple) and len(features_null) == 2:
            training_features_null, testing_features_null = features_null
        else:
            training_features_null = features_null
            testing_features_null = features_null

        lr_model.fit(training_features_null,
                     training_labels,
                     intercept_init = intercept_init)
        
        null_prob = lr_model.predict_proba(testing_features_null)
        df = testing_features_alternate.shape[1] - testing_features_null.shape[1]
    else:
        intercepts = calculate_intercept(testing_labels)
        null_prob = calculate_null_model(testing_labels.shape[0],
                                         intercepts)
        df = testing_features_alternate.shape[1]

    lr_model.fit(training_features_alternate,
                 training_labels,
                 intercept_init = intercept_init)
        
    alt_prob = lr_model.predict_proba(testing_features_alternate)

    alt_log_likelihood = -log_loss(testing_labels,
                                   alt_prob,
                                   normalize=False)
    null_log_likelihood = -log_loss(testing_labels,
                                    null_prob,
                                    normalize=False)

    G = g_scaling_factor * 2.0 * (alt_log_likelihood - null_log_likelihood)
    p_value = chi2.sf(G, df)

    return p_value


def null_predict_proba(intercept):
    prob = 1.0 / (1.0 + np.exp(-intercept))
    return prob

def mcfadden_r2(X, y, n_models=250, n_iter=20):
    """
    Computes coefficient of determination using McFadden's method
    """
    y = np.array(y)

    sgd = SGDClassifier(loss="log",
                        penalty="l2",
                        n_iter=n_iter)
    ensemble = BaggingClassifier(sgd,
                                 n_estimators=n_models,
                                 bootstrap=True)
    ensemble.fit(X, y)

    model_probs = ensemble.predict_proba(X)
    model_log_likelihood = -1. * log_loss(y,
                                          model_probs[:, 1],
                                          normalize=False)

    # "null model" with only intercept
    intercept = np.float(np.sum(y)) / y.shape[0]
    null_probs = np.ones(y.shape) * null_predict_proba(intercept)
    null_log_likelihood = -1. * log_loss(y,
                                         null_probs,
                                         normalize=False)
    
    return 1.0 - model_log_likelihood / null_log_likelihood

class PCA(object):
    def __init__(self, n_components):
        self.svd = TruncatedSVD(n_components=n_components)

    def explained_variance(self, features):
        self.svd.fit(features)

        return self.svd.explained_variance_ratio_

    def transform(self, features):
        coordinates = self.svd.fit_transform(features)
        self.components_ = self.svd.components_
        return coordinates


class LogisticRegressionEnsemble(object):
    """
    Implementation of an ensemble of Logistic
    Regression classifiers that supports bagging.
    """

    def __init__(self, n_models, method, batch_size, bagging=True, n_iter=20):
        self.n_models = n_models
        self.bagging = bagging
        self.batch_size = batch_size
        self.method = method
        self.n_iter = n_iter

    def get_base(self, n_samples):
        if self.method == "sag-l2":
            # copied from http://scikit-learn.org/stable/auto_examples/linear_model/plot_sgd_comparison.html#sphx-glr-auto-examples-linear-model-plot-sgd-comparison-py
            return LogisticRegression(solver="sag", tol=1e-1, C=1.e4 / n_samples)
        elif self.method == "sgd-l2":
            return SGDClassifier(loss="log", penalty="l2", n_iter = self.n_iter)
        elif self.method == "sgd-en":
            return SGDClassifier(loss="log", penalty="elasticnet", n_iter = self.n_iter)
        elif self.method == "asgd-l2":
            return SGDClassifier(loss="log", penalty="l2", average=True, n_iter = self.n_iter)
        else:
            raise Exception, "Unknown logistic regression method '%s'" % self.method

    def feature_importances(self, X, y):
        y = np.array(y)
        feature_importances = np.zeros(X.shape[1])
        trained_models = 0
        while trained_models < self.n_models:
            to_train = min(self.batch_size, self.n_models - trained_models)
            ensemble = BaggingClassifier(self.get_base(X.shape[0]),
                                         n_estimators=to_train,
                                         bootstrap=self.bagging)
            ensemble.fit(X, y)
            for model in ensemble.estimators_:
                feature_importances += model.coef_[0] / la.norm(model.coef_)
            trained_models += to_train
            print "Trained %s of %s models" % (trained_models, self.n_models)

        return np.abs(feature_importances / self.n_models)

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

    def feature_importances(self, X, y, statistics=False, interactions=False):
        y = np.array(y)
        feature_importances = np.zeros(X.shape[1])

        used_features_histogram = None
        if statistics:
            used_features_histogram = defaultdict(int)

        used_feature_sets = None
        if interactions:
            used_feature_sets = defaultdict(int)
            
        if self.n_resamples == -1:
            completed_trees = 0
            while completed_trees < self.n_trees:
                n_classifiers = min(self.batch_size, self.n_trees - completed_trees)
                rf = RandomForestClassifier(n_estimators=n_classifiers,
                                            n_jobs=1)
                rf.fit(X, y)
                feature_importances += rf.feature_importances_ * n_classifiers
                if statistics or interactions:
                    for dt in rf.estimators_:
                        tree = dt.tree_
                        used_features = set(tree.feature)
                        # leaves denoted by -2
                        used_features.remove(-2)
                        if statistics:
                            used_features_histogram[len(used_features)] += 1
                        if interactions:
                            used_feature_sets[frozenset(used_features)] += 1
                completed_trees += n_classifiers
                print "Trained", completed_trees, "of", self.n_trees, "trees"
        else:
            for i in xrange(self.n_trees):
                dt = DecisionTreeClassifier(max_features="sqrt")
                X_new, y_new = self._resample(X, y)
                dt.fit(X_new, y_new)
                if statistics or interactions:
                    used_features = set(dt.tree_.feature)
                    used_features.remove(-2)
                    if statistics:
                        used_features_histogram[len(used_features)] += 1
                    if interactions:
                        used_feature_sets[frozenset(used_features)] += 1
                feature_importances += dt.feature_importances_

        feature_importances = feature_importances / self.n_trees
        if interactions:
            used_feature_sets = dict(used_feature_sets)
            
        return feature_importances, used_features_histogram, used_feature_sets
