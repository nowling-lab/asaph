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
from scipy.stats import shapiro
from scipy.stats import ttest_1samp

from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model import SGDRegressor
from sklearn.metrics import log_loss
from sklearn.preprocessing import LabelEncoder

def estimate_lr_iter(n_samples):
    return max(20,
               int(np.ceil(100000. / n_samples)))

def upsample_features(labels, features):
    n_samples, n_features = features.shape

    # we make 1 copy for each variable so we can impute each unknown genotype
    N_COPIES = features.shape[1]
    training_labels = np.zeros(N_COPIES * n_samples)
    training_features = np.zeros((N_COPIES * n_samples, n_features))

    for i in xrange(n_samples):
        gt = None
        if features[i, :].sum() > 0:
            gt = features[i, :].argmax()

        for j in xrange(N_COPIES):
            idx = N_COPIES * i + j
            training_labels[idx] = labels[i]

            if gt is not None:
                training_features[idx, :] = features[i, :]
            else:
                training_features[idx, j] = 1.

    return training_labels, training_features

def lin_reg_log_likelihood(lr, X, y):
    pred_y = lr.predict(X)
    N, n_params = X.shape

    # estimate variance (sigma2)
    avg_y = np.mean(y)
    diff = y - avg_y
    diff2 = np.dot(diff, diff)
    sigma2 = diff2 / (N - 1)

    # residual sum of squares
    error = y - pred_y
    error2 = np.dot(error, error)
    
    log_likelihood = -N * np.log(2. * np.pi * sigma2) / 2. - error2 / (2.0 * sigma2)
    
    return log_likelihood

def lin_reg_lrtest(X, y, n_iter, g_scaling_factor=1.0):
    null_lr = SGDRegressor(fit_intercept = True, max_iter=n_iter)
    null_X = np.zeros((X.shape[0], 1))
    null_lr.fit(null_X,
                y)

    alt_lr = SGDRegressor(fit_intercept = False, max_iter=n_iter)
    alt_lr.fit(X,
               y)

    null_likelihood = lin_reg_log_likelihood(null_lr,
                                             null_X,
                                             y)

    alt_likelihood = lin_reg_log_likelihood(alt_lr,
                                            X,
                                            y)

    G = g_scaling_factor * 2. * (alt_likelihood - null_likelihood)

    p_value = chi2.sf(G, X.shape[1] - 1)

    p_value = max(1e-300, p_value)

    return p_value, alt_lr

def genotype_ttest(X, y):
    flattened = X.argmax(axis=1)

    # default to 1.0
    p_values = np.ones(3)
    for i in range(3):
        in_group = y[flattened == i]

        # need at least 2 samples to do a t-test
        if in_group.shape[0] >= 2:
            _, p_value = ttest_1samp(in_group, 0.0)
            p_values[i] = p_value

    return p_values

def genotype_normality_test(X, y):
    flattened = X.argmax(axis=1)

    # default to 1.0
    p_values = np.ones(3)
    for i in range(3):
        in_group = y[flattened == i]

        # need at least 3 samples to do a shapiro test
        if in_group.shape[0] >= 3:
            _, p_value = shapiro(in_group)
            p_values[i] = p_value

    return p_values

def snp_linreg_pvalues(X, y):
    N_GENOTYPES = 3
    n_iter = estimate_lr_iter(len(y))

    adj_y, adj_X = upsample_features(y, X)
    g_scaling_factor = 1.0 / N_GENOTYPES

    snp_p_value, model = lin_reg_lrtest(adj_X,
                                        adj_y,
                                        n_iter,
                                        g_scaling_factor=g_scaling_factor)

    gt_pred_ys = model.predict(np.eye(N_GENOTYPES))
    gt_ttest_pvalues = genotype_ttest(X, y)
    gt_normality_pvalues = genotype_normality_test(X, y)

    return snp_p_value, gt_ttest_pvalues, gt_normality_pvalues, gt_pred_ys

def likelihood_ratio_test(features_alternate, labels, lr_model, set_intercept=True, g_scaling_factor=1.0):
    if isinstance(features_alternate, tuple) and len(features_alternate) == 2:
        training_features, testing_features = features_alternate
        training_labels, testing_labels = labels
    else:
        training_features = features_alternate
        testing_features = features_alternate
        training_labels = labels
        testing_labels = labels

    n_training_samples = training_features.shape[0]
    n_testing_samples = testing_features.shape[0]
    n_iter = estimate_lr_iter(n_testing_samples)

    # null model
    null_lr = SGDClassifier(loss = "log",
                            fit_intercept = False,
                            max_iter = n_iter)
    null_training_X = np.ones((n_training_samples, 1))
    null_testing_X = np.ones((n_testing_samples, 1))
    null_lr.fit(null_training_X,
                training_labels)
    null_prob = null_lr.predict_proba(null_testing_X)

    intercept_init = None
    if set_intercept:
        intercept_init = null_lr.coef_[:, 0]

    lr_model.fit(training_features,
                 training_labels,
                 intercept_init = intercept_init)
    alt_prob = lr_model.predict_proba(testing_features)
        
    alt_log_likelihood = -log_loss(testing_labels,
                                   alt_prob,
                                   normalize=False)
    null_log_likelihood = -log_loss(testing_labels,
                                    null_prob,
                                    normalize=False)

    G = g_scaling_factor * 2.0 * (alt_log_likelihood - null_log_likelihood)
    
    # both models have intercepts so the intercepts cancel out
    df = training_features.shape[1]
    p_value = chi2.sf(G, df)

    return p_value


def null_predict_proba(intercept):
    prob = 1.0 / (1.0 + np.exp(-intercept))
    return prob
