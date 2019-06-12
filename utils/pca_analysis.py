"""
Copyright 2019 Ronald J. Nowling

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


import argparse
from collections import defaultdict
import os
import sys
import warnings

import numpy as np

from sklearn.cluster import k_means

from scipy.stats import chi2
from scipy.stats import chisquare

from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import log_loss
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import seaborn as sns

# apply seaborn style config
sns.set()

def read_pca_coordinates(flname):
    sample_coordinates = []
    sample_names = []
    with open(flname) as fl:
        # skip header
        next(fl)
        for ln in fl:
            cols = ln.split("\t")

            sample_name = cols[0]
            coordinates = list(map(float, cols[3:]))
            
            sample_names.append(sample_name)
            sample_coordinates.append(coordinates)

    return sample_names, np.array(sample_coordinates)

def read_labels(flname):
    sample_indices = dict()

    with open(flname) as fl:
        for label_idx, ln in enumerate(fl):
            cols = ln.strip().split(",")

            label = cols[0]

            for sample_name in cols[1:]:
                sample_indices[sample_name] = label_idx

    return sample_indices

def pairwise(iterable):
    iterable = iter(iterable)
    try:
        while True:
            a = next(iterable)
            b = next(iterable)
            yield a, b
    except StopIteration:
        pass
    
def plot_projections(coordinates, pairs, dirname):
    if len(args.pairs) % 2 != 0:
        print("Error: PCs must be provided in pairs of 2")
        sys.exit(1)

    for p1, p2 in pairwise(pairs):
        fig_flname = os.path.join(dirname,
                                  "pca_projection_%s_%s.png" % (str(p1), str(p2)))
        plt.clf()
        plt.scatter(coordinates[:, p1 - 1],
                    coordinates[:, p2 - 1],
                    color="m",
                    marker="o",
                    edgecolor="k",
                    alpha=0.7)
        plt.xlabel("Principal Component %s" % p1, fontsize=16)
        plt.ylabel("Principal Component %s" % p2, fontsize=16)
        plt.savefig(fig_flname,
                    DPI=300)

def sweep_clusters(coordinates, components, n_clusters, plot_flname):
    components = list(map(lambda idx: idx - 1, components))
    
    selected = coordinates[:, components]

    inertia_values = []
    for k in n_clusters:
        _, _, inertia = k_means(selected,
                                k,
                                n_jobs=-2)
        print("Clustering with %s clusters: %s" % (k, inertia))
        inertia_values.append(inertia)

    plt.plot(n_clusters,
             inertia_values,
             "k.-")
    plt.xlabel("Number of Clusters", fontsize=16)
    plt.ylabel("Inertia", fontsize=16)

    plt.savefig(plot_flname,
                DPI=300)

def test_clusters(coordinates, components, n_clusters, sample_names, labels):
    components = list(map(lambda idx: idx - 1, components))
    
    selected = coordinates[:, components]

    _, cluster_idx, inertia = k_means(selected,
                                      n_clusters,
                                      n_jobs=-2)

    # calculate observation table
    n_labels = len(set(labels.values()))
    sample_labels = []
    obs_table = np.zeros((n_clusters, n_labels))
    for idx, sample_name in zip(cluster_idx, sample_names):
        label_idx = labels[sample_name]
        obs_table[idx, label_idx] += 1
        sample_labels.append(label_idx)
    sample_labels = np.array(sample_labels)

    # calculate table of expectations
    # that assumes that variables are independent
    exp_table = np.zeros((n_clusters, n_labels))
    obs_row_sums = obs_table.sum(axis=1)
    obs_col_sums = obs_table.sum(axis=0)
    n_samples = np.sum(obs_table)
    for r in range(obs_table.shape[0]):
        for c in range(obs_table.shape[1]):
            exp_table[r, c] = obs_row_sums[r] * obs_col_sums[c] / n_samples

    # scipy's dof
    k = obs_table.shape[0] * obs_table.shape[1] - 1

    # true dof
    dof = (obs_table.shape[0] - 1) * (obs_table.shape[1] - 1)
    ddof = k - dof
    

    _, pvalue = chisquare(obs_table.ravel(),
                          exp_table.ravel(),
                          ddof=ddof)

    print(obs_table)
    print()
    print(exp_table)
    
    print("Chi2 p-value:", pvalue)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        feature_encoder = OneHotEncoder(sparse=False)
        features = feature_encoder.fit_transform(sample_labels.reshape(-1, 1))

        label_encoder = LabelEncoder()
        sample_labels = label_encoder.fit_transform(cluster_idx)
    
        lr = SGDClassifier(loss="log", penalty="l2", max_iter=1000, tol=1e-4)
        lr.fit(features, sample_labels)
        pred_labels = lr.predict(features)
        acc = accuracy_score(sample_labels, pred_labels)
        bal_acc = balanced_accuracy_score(sample_labels, pred_labels)

        print(pred_labels)
        print()
        print(sample_labels)

        print("Classifier accuracy:", acc)
        print("Balanced accuracy:", bal_acc)

def output_clusters(coordinates, components, n_clusters, sample_names, labels_fl):
    components = list(map(lambda idx: idx - 1, components))
    
    selected = coordinates[:, components]

    _, cluster_idx, inertia = k_means(selected,
                                      n_clusters,
                                      n_jobs=-2)

    # group samples by cluster
    populations = defaultdict(set)
    for i, sample_name in enumerate(sample_names):
        cluster_assignment = cluster_idx[i]
        populations[cluster_assignment].add(sample_name)

    with open(labels_fl, "w") as fl:
        for pop_name, samples in populations.items():
            fl.write(str(pop_name))
            for name in samples:
                fl.write(",")
                fl.write(name)
            fl.write("\n")
        

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
    null_lr = SGDClassifier(penalty="l2",
                            loss = "log",
                            max_iter = n_iter * 10.,
                            fit_intercept = False,
                            tol=1e-8)
    
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

def estimate_lr_iter(n_samples):
    return max(20,
               int(np.ceil(1000000. / n_samples)))

def create_class_labels(sample_names, sample_labels):
    str_labels = [sample_labels[name] for name in sample_names]

    encoder = LabelEncoder()
    labels = encoder.fit_transform(str_labels)

    return labels, encoder.classes_
    

def test_labels(coordinates, sample_names, sample_labels):
    class_labels, class_names = create_class_labels(sample_names, sample_labels)

    n_iter = estimate_lr_iter(len(coordinates))
    # we set the intercept to the class ratios in the lr test function
    lr = SGDClassifier(penalty="l2",
                       loss="log",
                       max_iter = n_iter * 10.,
                       tol=1e-8,
                       fit_intercept=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i in range(coordinates.shape[1]):
            features = coordinates[:, i].reshape(-1, 1)
        
            p_value = likelihood_ratio_test(features,
                                            class_labels,
                                            lr,
                                            set_intercept=False)

            lr.fit(features, class_labels)
            pred_labels = lr.predict(features)
            acc = 100. * accuracy_score(class_labels,
                                        pred_labels)

            bal_acc = 100. * balanced_accuracy_score(class_labels,
                                                     pred_labels)

            cm = confusion_matrix(class_labels,
                                  pred_labels)

            print("Component:", (i+1))
            print("p-value: ", p_value)
            print("Accuracy:", acc)
            print("Balanced accuracy:", bal_acc)
            print(cm)
            print()


def parseargs():
    parser = argparse.ArgumentParser()

    parser.add_argument("--coordinates",
                        type=str,
                        required=True)

    subparsers = parser.add_subparsers(dest="mode")

    plot_parser = subparsers.add_parser("plot-projections",
                                        help="Plot PCA projections")
    
    plot_parser.add_argument("--plot-dir",
                             type=str,
                             required=True)

    plot_parser.add_argument("--pairs",
                             nargs="+",
                             type=int,
                             required=True)

    sweep_clusters_parser = subparsers.add_parser("sweep-clusters",
                                                   help="K-Means for a range of clusters")
    
    sweep_clusters_parser.add_argument("--components",
                                       type=int,
                                       nargs="+",
                                       required=True,
                                       help="Components to use in projection")

    sweep_clusters_parser.add_argument("--n-clusters",
                                       type=int,
                                       nargs="+",
                                       required=True,
                                       help="Cluster counts to try")

    sweep_clusters_parser.add_argument("--plot-fl",
                                       type=str,
                                       required=True,
                                       help="Plot filename")

    test_clusters_parser = subparsers.add_parser("test-clusters",
                                                   help="Association tests against clusters")
    
    test_clusters_parser.add_argument("--components",
                                       type=int,
                                       nargs="+",
                                       required=True,
                                       help="Components to use in projection")

    test_clusters_parser.add_argument("--n-clusters",
                                       type=int,
                                       required=True,
                                       help="Number of clusters")
    
    test_clusters_parser.add_argument("--labels-fl",
                                       type=str,
                                       required=True,
                                       help="Labels file")

    out_clusters_parser = subparsers.add_parser("output-clusters",
                                                help="Output cluster labels")

    out_clusters_parser.add_argument("--components",
                                     type=int,
                                     nargs="+",
                                     required=True,
                                     help="Components to use in projection")

    out_clusters_parser.add_argument("--n-clusters",
                                     type=int,
                                     required=True,
                                     help="Number of clusters")
    
    out_clusters_parser.add_argument("--labels-fl",
                                     type=str,
                                     required=True,
                                     help="Labels file to output")
    
    label_test_parser = subparsers.add_parser("test-labels",
                                              help="Run association tests on PCs vs labels")

    label_test_parser.add_argument("--labels-fl",
                                   type=str,
                                   required=True,
                                   help="Labels file")
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    if not os.path.exists(args.coordinates):
        print("Coordinates file path is invalid")
        sys.exit(1)

    sample_names, coordinates = read_pca_coordinates(args.coordinates)

    print(coordinates.shape)

    if args.mode == "plot-projections":
        plot_projections(coordinates,
                         args.pairs,
                         args.plot_dir)
        
    elif args.mode == "sweep-clusters":
        sweep_clusters(coordinates,
                       args.components,
                       args.n_clusters,
                       args.plot_fl)
        
    elif args.mode == "test-clusters":
        labels = read_labels(args.labels_fl)
        test_clusters(coordinates,
                      args.components,
                      args.n_clusters,
                      sample_names,
                      labels)

    elif args.mode == "output-clusters":
        output_clusters(coordinates,
                        args.components,
                        args.n_clusters,
                        sample_names,
                        args.labels_fl)

    elif args.mode == "test-labels":
        labels = read_labels(args.labels_fl)
        test_labels(coordinates,
                    sample_names,
                    labels)
        
    else:
        print("Unknown mode '%s'" % args.mode)
        sys.exit(1)
