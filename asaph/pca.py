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

from sklearn.cluster import k_means
from sklearn.decomposition import NMF
from sklearn.decomposition import PCA
from sklearn.externals import joblib
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import binarize

from asaph.ml import estimate_lr_iter
from asaph.ml import likelihood_ratio_test
from asaph.ml import snp_linreg_pvalues

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

    if args.model_type == "PCA":
        pca = PCA(n_components = args.n_components,
                  whiten = True)
    elif args.model_type == "NMF":
        pca = NMF(n_components = args.n_components)
    else:
        raise Exception("Unknown model type %s" % args.model_type)
        
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

    print explained_variance_ratios

    plt.clf()
    plt.grid(True)
    plt.plot(xrange(1, len(explained_variance_ratios) + 1),
             explained_variance_ratios, "m.-")
    plt.xlabel("Principal Component", fontsize=16)
    plt.ylabel("Explained Variance Ratio", fontsize=16)
    plt.ylim([0., 1.])
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
    selected = projected[:, map(lambda idx: idx - 1, args.selected_components)]

    features = read_features(workdir)

    with open(args.output_fl, "w") as fl:
        headers = ["sample", "population_name"]
        headers.extend(map(str, args.selected_components))
        fl.write("\t".join(headers))
        fl.write("\n")

        for i in xrange(len(features.sample_labels)):
            sample = features.sample_labels[i]
            line = [sample]
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


def snp_logreg_association_tests(args):
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
    
    n_pcs = projections.shape[0]
    for pc in args.components:
        flname = os.path.join(analysis_dir, "snp_pc_%s_logreg_assoc_tests.tsv" % pc)
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

def snp_linreg_association_tests(args):
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

    n_pcs = projections.shape[0]
    for pc in args.components:
        flname = os.path.join(analysis_dir, "snp_pc_%s_linreg_assoc_tests.tsv" % pc)
        with open(flname, "w") as fl:
            next_output = 1
            for i, pair in enumerate(data_model.snp_feature_map.iteritems()):
                snp_label, feature_idx = pair
                chrom, pos = snp_label

                snp_features = data_model.feature_matrix[:, feature_idx]

                # since we make multiple copies of the original samples,
                # we need to scale the log loss so that it is correct for
                # the original sample size
                triplet = snp_linreg_pvalues(snp_features,
                                             projections[:, pc - 1])
                snp_p_value, gt_ttest_pvalues, gt_normality_pvalues, gt_pred_ys = triplet

                if i == next_output:
                    print i, "SNP", snp_label, "and PC", pc, "has p-value", snp_p_value
                    next_output *= 2

                fl.write("\t".join([chrom, pos, str(snp_p_value)]))
                for j in xrange(3):
                    fl.write("\t")
                    fl.write(str(gt_ttest_pvalues[j]))
                    fl.write("\t")
                    fl.write(str(gt_normality_pvalues[j]))
                    fl.write("\t")
                    fl.write(str(gt_pred_ys[j]))
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

    train_parser.add_argument("--model-type",
                              type=str,
                              choices=["PCA",
                                       "NMF"],
                              default="PCA")

    
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
    
    snp_association_parser = subparsers.add_parser("snp-association-tests",
                                                   help="Run association tests on PCs vs SNPs")
    snp_association_parser.add_argument("--components",
                                        type=int,
                                        nargs="+",
                                        required=True,
                                        help="Components to perform testing on")
    
    snp_association_parser.add_argument("--model-type",
                                        type=str,
                                        choices=["logistic",
                                                 "linear"],
                                        required=True,
                                        help="Type of model")

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
    elif args.mode == "output-coordinates":
        output_coordinates(args)
    elif args.mode == "output-loading-magnitudes":
        output_loading_magnitudes(args)
    elif args.mode == "output-loading-factors":
        output_loading_factors(args)
    elif args.mode == "snp-association-tests":
        if args.model_type == "linear":
            snp_linreg_association_tests(args)
        elif args.model_type == "logistic":
            snp_logreg_association_tests(args)
    else:
        print "Unknown mode '%s'" % args.mode
        sys.exit(1)
