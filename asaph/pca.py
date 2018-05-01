"""
Copyright 2017 Ronald J. Nowling

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
from collections import OrderedDict
from itertools import tee, izip
import os
import sys

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import numpy as np
import seaborn

from sklearn.cluster import k_means
from sklearn.decomposition import PCA
from sklearn.externals import joblib
from sklearn.linear_model import SGDClassifier
from sklearn.preprocessing import binarize

from asaph.ml import estimate_lr_iter
from asaph.ml import likelihood_ratio_test
from asaph.newioutils import read_features
from asaph.newioutils import deserialize
from asaph.newioutils import PROJECT_SUMMARY_FLNAME


MODEL_KEY = "model"
PROJECTION_KEY = "projected-coordinates"

def train(args):
    workdir = args.workdir

    models_dir = os.path.join(workdir,
                              "models")
    if not os.path.exists(models_dir):
        os.makedirs(models_dir)

    features = read_features(workdir)
        
    pca = PCA(n_components = args.n_components,
              whiten = True)
    projections = pca.fit_transform(features.feature_matrix)

    model = { MODEL_KEY : pca,
              PROJECTION_KEY : projections}

    model_fl = os.path.join(models_dir,
                            "pca.pkl")
    joblib.dump(model,
                model_fl)

def explained_variance_analysis(args):
    workdir = args.workdir

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    model_fl = os.path.join(workdir,
                            "models",
                            "pca.pkl")

    model = joblib.load(model_fl)

    explained_variance_ratios = model[MODEL_KEY].explained_variance_ratio_

    fig_flname = os.path.join(figures_dir,
                              "pca_explained_variance_ratios.png")

    plt.clf()
    plt.grid(True)
    plt.plot(xrange(1, len(explained_variance_ratios) + 1),
             explained_variance_ratios, "m.-")
    plt.xlabel("Principal Component", fontsize=16)
    plt.ylabel("Explained Variance Ratio", fontsize=16)
    plt.ylim([0., 1.])
    plt.savefig(fig_flname,
                DPI=300)

def min_components_explained_variance(args):
    workdir = args.workdir

    features = read_features(workdir)

    n_components = args.init_n_components
    while True:
        print "Computing PCA with %s components" % n_components
        pca = PCA(n_components = n_components,
                  whiten = True)
        pca.fit(features.feature_matrix)
        explained_variance_ratios = pca.explained_variance_ratio_
        sorted_ratios = np.sort(explained_variance_ratios)[::-1]
        cum_ratios = np.cumsum(sorted_ratios)
        total_explained_variance = cum_ratios[-1]
        if total_explained_variance >= args.explained_variance_threshold:
            break
        n_components *= 2

    needed_components = 0
    achieved_ev_ratio = 0.0
    for i, ev_ratio in enumerate(cum_ratios):
        if ev_ratio >= args.explained_variance_threshold:
            needed_components = i + 1
            achieved_ev_ratio = ev_ratio
            break

    print "Explained-variance threshold of %s surpassed at %s with %s components" % \
        (args.explained_variance_threshold,
         achieved_ev_ratio,
         needed_components)

def pairwise(iterable):
    iterable = iter(iterable)
    try:
        while True:
            a = next(iterable)
            b = next(iterable)
            yield a, b
    except StopIteration:
        pass
    
def plot_projections(args):
    if len(args.pairs) % 2 != 0:
        print "Error: PCs must be provided in pairs of 2"
        sys.exit(1)

    workdir = args.workdir

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    model_fl = os.path.join(workdir,
                            "models",
                            "pca.pkl")
    model = joblib.load(model_fl)
    projected = model[PROJECTION_KEY]
    
    features = read_features(workdir)
    
    project_summary = deserialize(os.path.join(workdir,
                                               PROJECT_SUMMARY_FLNAME))

    all_labels = set(features.class_labels)
    labels = np.array(features.class_labels, dtype=np.int32)
    populations = []
    for l in all_labels:
        pop = labels == l
        pop_name = project_summary.population_names[l]
        populations.append((pop, pop_name))

    for p1, p2 in pairwise(args.pairs):
        fig_flname = os.path.join(figures_dir,
                                  "pca_projection_%s_%s.png" % (str(p1), str(p2)))
        plt.clf()
        plt.grid(True)
        colors = ["m", "c", "k", "r", "g", "b"]
        markers = ["o"] * len(colors) + \
                  ["s"] * len(colors) + \
                  ["+"] * len(colors)
        for idx, (pop_idx, pop_name) in enumerate(populations):
            plt.scatter(projected[pop_idx, p1 - 1],
                        projected[pop_idx, p2 - 1],
                        color=colors[idx % len(colors)],
                        marker=markers[idx % len(markers)],
                        edgecolor="k",
                        alpha=0.7,
                        label=pop_name)
        plt.xlabel("Principal Component %s" % p1, fontsize=16)
        plt.ylabel("Principal Component %s" % p2, fontsize=16)
        if len(all_labels) > 1:
            plt.legend()
        plt.savefig(fig_flname,
                    DPI=300)

def plot_densities(args):
    workdir = args.workdir

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)

    model_fl = os.path.join(workdir,
                            "models",
                            "pca.pkl")
    model = joblib.load(model_fl)
    projected = model[PROJECTION_KEY]
    
    for i in args.components:
        fig_flname = os.path.join(figures_dir,
                                  "pca_density_%s.png" % i)
        plt.clf()
        seaborn.distplot(projected[:, i - 1])
        plt.xlabel("PC %s Coordinate" % i, fontsize=16)
        plt.ylabel("Density", fontsize=16)
        plt.ylim([0.0, 1.0])
        plt.savefig(fig_flname,
                    DPI=300)
        
def output_coordinates(args):
    workdir = args.workdir

    project_summary = deserialize(os.path.join(workdir,
                                               PROJECT_SUMMARY_FLNAME))

    model_fl = os.path.join(workdir,
                            "models",
                            "pca.pkl")
    model = joblib.load(model_fl)    
    projected = model[PROJECTION_KEY]
    selected = projected[:, args.selected_components]

    features = read_features(workdir)

    with open(args.output_fl, "w") as fl:
        headers = ["sample", "population_index", "population_name"]
        headers.extend(map(str, args.selected_components))
        fl.write("\t".join(headers))
        fl.write("\n")

        for i in xrange(len(features.sample_labels)):
            sample = features.sample_labels[i]
            pop_idx = features.class_labels[i]
            pop_name = project_summary.population_names[pop_idx]
            line = [sample, str(pop_idx), pop_name]
            line.extend(map(str, selected[i, :]))
            fl.write("\t".join(line))
            fl.write("\n")

def output_loading_magnitudes(args):
    workdir = args.workdir

    analysis_dir = os.path.join(workdir, "analysis")
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    
    data_model = read_features(workdir)

    model_fl = os.path.join(workdir,
                            "models",
                            "pca.pkl")
    model = joblib.load(model_fl)    
    pca = model[MODEL_KEY]
    components = pca.components_
    selected = components[map(lambda idx: idx - 1, args.components), :]

    output_fl = os.path.join(analysis_dir, "pca_loading_magnitudes.tsv")
    with open(output_fl, "w") as fl:
        header = ["chromosome", "position"]
        header.extend(map(str, args.components))
        fl.write("\t".join(header))
        fl.write("\n")
        for i, pair in enumerate(data_model.snp_feature_map.iteritems()):
            snp_label, feature_idx = pair
            chrom, pos = snp_label

            features = selected[:, feature_idx]
            norms = np.sqrt(np.sum(features**2, axis=1))
            
            fl.write("%s\t%s\t" % (chrom, pos))
            fl.write("\t".join(map(str, norms)))
            fl.write("\n")

def output_loading_factors(args):
    workdir = args.workdir

    analysis_dir = os.path.join(workdir, "analysis")
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    
    data_model = read_features(workdir)

    model_fl = os.path.join(workdir,
                            "models",
                            "pca.pkl")
    model = joblib.load(model_fl)    
    pca = model[MODEL_KEY]
    components = pca.components_
    selected = components[map(lambda idx: idx - 1, args.components), :]

    output_fl = os.path.join(analysis_dir, "pca_loading_factors.tsv")
    with open(output_fl, "w") as fl:
        header = ["chromosome", "position", "dummy variable"]
        header.extend(map(str, args.components))
        fl.write("\t".join(header))
        fl.write("\n")
        for i, pair in enumerate(data_model.snp_feature_map.iteritems()):
            snp_label, feature_idx = pair
            chrom, pos = snp_label

            for j, idx in enumerate(feature_idx):
                features = selected[:, idx]
            
                fl.write("%s\t%s\t%s\t" % (chrom, pos, j))
                fl.write("\t".join(map(str, features)))
                fl.write("\n")

def analyze_weights(args):
    workdir = args.workdir

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)
    
    project_summary = deserialize(os.path.join(workdir,
                                               PROJECT_SUMMARY_FLNAME))

    model_fl = os.path.join(workdir,
                            "models",
                            "pca.pkl")
    model = joblib.load(model_fl)    
    pca = model[MODEL_KEY]
    
    for i, w in enumerate(args.weights):
        plt.clf()
        seaborn.distplot(w * pca.components_[i, :],
                         kde=False)
        plt.xlabel("Feature Weights", fontsize=16)
        plt.ylabel("Count(Features)", fontsize=16)
        plot_fl = os.path.join(figures_dir,
                               "pca_feature_weights_pc%s.png" % (i + 1))
        plt.savefig(plot_fl,
                    DPI=300)

def pop_association_tests(args):
    workdir = args.workdir

    analysis_dir = os.path.join(workdir, "analysis")
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    project_summary = deserialize(os.path.join(workdir,
                                               PROJECT_SUMMARY_FLNAME))
    
    model_fl = os.path.join(workdir,
                            "models",
                            "pca.pkl")
    model = joblib.load(model_fl)    
    projections = model[PROJECTION_KEY]

    data_model = read_features(workdir)

    n_iter = estimate_lr_iter(len(data_model.sample_labels))
    # we set the intercept to the class ratios in the lr test function
    lr = SGDClassifier(penalty="l2",
                       loss="log",
                       n_iter = n_iter,
                       fit_intercept=False)

    pvalues_fl = os.path.join(analysis_dir, "population_pca_association_tests.tsv")
    class_labels = np.array(data_model.class_labels)
    with open(pvalues_fl, "w") as fl:
        for i in xrange(projections.shape[1]):
            features = projections[:, i].reshape(-1, 1)

            p_value = likelihood_ratio_test(features,
                                            class_labels,
                                            lr)
            fl.write("%s\t%s\n" % (i + 1, p_value))

def generate_training_set(features, projections):
    """
    Converts genotype categories to class labels and imputed unknown genotypes.
    Imputation involves duplicating samples, so we create copies of projection
    coordinates to align.
    """

    n_samples = features.shape[0]
    n_features = features.shape[1]
    N_GENOTYPES = 3
    if n_features != N_GENOTYPES:
        raise Exception("Must be using genotype categories!")

    N_COPIES = 3
    class_labels = np.zeros(N_COPIES * n_samples)
    imputed_projections = np.zeros(N_COPIES * n_samples)
    for i in xrange(n_samples):
        gt = None
        for j in xrange(N_GENOTYPES):
            if features[i, j] == 1.0:
                gt = j
        
        for j in xrange(N_COPIES):
            idx = N_COPIES * i + j
            imputed_projections[idx] = projections[i]

            if gt is None:
                class_labels[idx] = j
            else:
                class_labels[idx] = gt

    return N_COPIES, class_labels, imputed_projections
    

def snp_association_tests(args):
    workdir = args.workdir

    analysis_dir = os.path.join(workdir, "analysis")
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    project_summary = deserialize(os.path.join(workdir,
                                               PROJECT_SUMMARY_FLNAME))
    
    model_fl = os.path.join(workdir,
                            "models",
                            "pca.pkl")
    model = joblib.load(model_fl)    
    projections = model[PROJECTION_KEY]

    data_model = read_features(workdir)

    n_iter = estimate_lr_iter(len(data_model.class_labels))
    # we set the intercept to the class ratios in the lr test function
    lr = SGDClassifier(penalty="l2",
                       loss="log",
                       n_iter = n_iter,
                       fit_intercept=False)
    
    n_pcs = projections.shape[0]
    for pc in args.components:
        flname = os.path.join(analysis_dir, "snp_pc_%s_association_tests.tsv" % pc)
        with open(flname, "w") as fl:
            next_output = 1
            for i, pair in enumerate(data_model.snp_feature_map.iteritems()):
                snp_label, feature_idx = pair
                chrom, pos = snp_label

                snp_features = data_model.feature_matrix[:, feature_idx]
                triplet = generate_training_set(snp_features,
                                                projections[:, pc - 1])
                n_copies, class_labels, imputed_projections = triplet

                imputed_projections = imputed_projections.reshape(-1, 1)

                # since we make multiple copies of the original samples,
                # we need to scale the log loss so that it is correct for
                # the original sample size
                try:
                    p_value = likelihood_ratio_test(imputed_projections,
                                                    class_labels,
                                                    lr,
                                                    g_scaling_factor = 1.0 / n_copies)
                # in case of underflow or overflow in a badly-behaving model
                except ValueError:
                    p_value = 1.0
                
                if i == next_output:
                    print i, "SNP", snp_label, "and PC", pc, "has p-value", p_value
                    next_output *= 2

                fl.write("\t".join([chrom, pos, str(p_value)]))
                fl.write("\n")

def sweep_clusters(args):
    workdir = args.workdir

    figures_dir = os.path.join(workdir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)
    
    project_summary = deserialize(os.path.join(workdir,
                                               PROJECT_SUMMARY_FLNAME))

    model_fl = os.path.join(workdir,
                            "models",
                            "pca.pkl")
    model = joblib.load(model_fl)    
    projected = model[PROJECTION_KEY]
    components = map(lambda idx: idx - 1, args.components)
    selected = projected[:, components]

    features = read_features(workdir)

    inertia_values = []
    for k in args.n_clusters:
        print "Clustering with %s clusters" % k
        _, _, inertia = k_means(selected,
                                k,
                                n_jobs=-2)
        inertia_values.append(inertia)


    plt.plot(args.n_clusters,
             inertia_values,
             "k.-")
    plt.xlabel("Number of Clusters", fontsize=16)
    plt.ylabel("Inertia", fontsize=16)

    fig_flname = os.path.join(figures_dir,
                              "cluster_inertia")
    for dim in args.components:
        fig_flname += "_%s" % dim
    fig_flname += ".png"

    plt.savefig(fig_flname,
                DPI=300)

def cluster_samples(args):
    workdir = args.workdir

    analysis_dir = os.path.join(workdir, "analysis")
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    
    project_summary = deserialize(os.path.join(workdir,
                                               PROJECT_SUMMARY_FLNAME))

    model_fl = os.path.join(workdir,
                            "models",
                            "pca.pkl")
    model = joblib.load(model_fl)    
    projected = model[PROJECTION_KEY]
    components = map(lambda idx: idx - 1, args.components)
    selected = projected[:, components]

    features = read_features(workdir)

    _, labels, inertia = k_means(selected,
                                 args.n_clusters,
                                 n_jobs=-2)

    fig_flname = os.path.join(analysis_dir,
                              "clusters_%s.tsv" % args.n_clusters)

    clusters = defaultdict(list)
    for name, cluster in zip(features.sample_labels, labels):
        clusters[cluster].append(name)
        
    with open(fig_flname, "w") as fl:
        for cluster, samples in clusters.iteritems():
            fl.write(str(cluster))
            fl.write(",")
            fl.write(",".join(samples))
            fl.write("\n")
    
def parseargs():
    parser = argparse.ArgumentParser(description="Asaph - PCA")

    parser.add_argument("--workdir",
                        type=str,
                        required=True,
                        help="Work directory")

    subparsers = parser.add_subparsers(dest="mode")
    
    train_parser = subparsers.add_parser("train",
                                         help="Train PCA model")
    
    train_parser.add_argument("--n-components",
                              type=int,
                              required=True,
                              help="Number of PCs to compute")
    
    eva_parser = subparsers.add_parser("explained-variance-analysis",
                                       help="Compute explained variances of PCs")

    output_parser = subparsers.add_parser("output-coordinates",
                                      help="Output PC projected coordinates")
        
    output_parser.add_argument("--selected-components",
                               type=int,
                               nargs="+",
                               help="Components to output")

    output_parser.add_argument("--output-fl",
                               type=str,
                               required=True,
                               help="Output file")
    
    plot_parser = subparsers.add_parser("plot-projections",
                                        help="Plot samples on principal coordinates")

    plot_parser.add_argument("--pairs",
                             type=int,
                             nargs="+",
                             required=True,
                             help="Pairs of PCs to plot")

    density_parser = subparsers.add_parser("plot-densities",
                                           help="Plot densities along PC coordinates")
    density_parser.add_argument("--components",
                                type=int,
                                nargs="+",
                                required=True,
                                help="Components to perform testing on")


    count_parser = subparsers.add_parser("min-components-explained-variance",
                                         help="Find number of components to explain specified percentage of variance")

    count_parser.add_argument("--init-n-components",
                              type=int,
                              required=True,
                              help="Initial number of components to compute")

    count_parser.add_argument("--explained-variance-threshold",
                              type=float,
                              required=True,
                              help="Minimum explained variance")

    weight_analysis_parser = subparsers.add_parser("analyze-weights",
                                                   help="Plot component weight distributions")

    weight_analysis_parser.add_argument("--weights",
                                        type=float,
                                        nargs="+",
                                        required=True,
                                        help="Component weights")

    snp_association_parser = subparsers.add_parser("snp-association-tests",
                                                   help="Run association tests on PCs vs SNPs")
    snp_association_parser.add_argument("--components",
                                        type=int,
                                        nargs="+",
                                        required=True,
                                        help="Components to perform testing on")

    pop_association_parser = subparsers.add_parser("pop-association-tests",
                                                   help="Run association tests on PCs vs population labels")

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

    cluster_samples_parser = subparsers.add_parser("cluster-samples",
                                                   help="Cluster samples with K-Means")
                     
    cluster_samples_parser.add_argument("--components",
                                       type=int,
                                       nargs="+",
                                       required=True,
                                       help="Components to use in projection")

    cluster_samples_parser.add_argument("--n-clusters",
                                       type=int,
                                       required=True,
                                       help="Number of clusters")

    loading_magnitudes_parser = subparsers.add_parser("output-loading-magnitudes",
                                                      help="Output loading magnitudes")

    loading_magnitudes_parser.add_argument("--components",
                                           type=int,
                                           nargs="+",
                                           required=True,
                                           help="Components to output")
    
    loading_factors_parser = subparsers.add_parser("output-loading-factors",
                                                   help="Output loading factors")

    loading_factors_parser.add_argument("--components",
                                        type=int,
                                        nargs="+",
                                        required=True,
                                        help="Components to output")
                     
    return parser.parse_args()

if __name__ == "__main__":
    args = parseargs()

    if args.mode == "train":
       train(args) 
    elif args.mode == "explained-variance-analysis":
        explained_variance_analysis(args)
    elif args.mode == "plot-projections":
        plot_projections(args)
    elif args.mode == "plot-densities":
        plot_densities(args)
    elif args.mode == "output-coordinates":
        output_coordinates(args)
    elif args.mode == "min-components-explained-variance":
        min_components_explained_variance(args)
    elif args.mode == "analyze-weights":
        analyze_weights(args)
    elif args.mode == "output-loading-magnitudes":
        output_loading_magnitudes(args)
    elif args.mode == "output-loading-factors":
        output_loading_factors(args)
    elif args.mode == "pop-association-tests":
        pop_association_tests(args)
    elif args.mode == "snp-association-tests":
        snp_association_tests(args)
    elif args.mode == "sweep-clusters":
        sweep_clusters(args)
    elif args.mode == "cluster-samples":
        cluster_samples(args)
    else:
        print "Unknown mode '%s'" % args.mode
        sys.exit(1)
