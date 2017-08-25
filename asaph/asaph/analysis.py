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

from collections import defaultdict

import numpy as np

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

def sampled_snps_curves(all_models):
    ordered_trees = sorted(all_models.keys())

    common_feature_counts = []
    snp1_feature_counts = []
    snp2_feature_counts = []
    used_models = []
    for n_trees in ordered_trees:
        models = all_models[n_trees]
        if len(models) != 2:
            continue

        snps1, snps2 = models
        used_models.append(n_trees)

        common_feature_counts.append(snps1.count_intersection(snps2))
        snp1_feature_counts.append(len(snps1))
        snp2_feature_counts.append(len(snps2))

    return used_models, \
        common_feature_counts, \
        snp1_feature_counts, \
        snp2_feature_counts

def similarity_within_model(models):
    similarities = []
    snps1, snps2 = models
    min_snps = min(len(snps1), len(snps2))
    for n in xrange(1, min_snps+1):
        percentage = 100.0 * float(snps1.take(n).count_intersection(snps2.take(n))) \
                     / float(n)
        similarities.append(percentage)
        
    return similarities

def plot_similarity_within_model(flname, similarities):
    plt.clf()
    plt.plot(np.arange(len(similarities)) + 1, similarities, "k.-")
    plt.xlabel("Top-Ranked SNPs", fontsize=16)
    plt.ylabel("Jaccard Similarity (%)", fontsize=16)
    plt.xlim([1, len(similarities) + 1])
    plt.ylim([0, 100])
    plt.savefig(flname, DPI=200)

def similarity_curves(thresholds, all_models, universe_size):
    ordered_models = sorted(all_models.keys())

    common_feature_threshold_percentages = defaultdict(list)
    used_models = []
    for n_models in ordered_models:
        models = all_models[n_models]
        if len(models) != 2:
            continue

        snps1, snps2 = models
        used_models.append(n_models)

        min_snps = min(len(snps1), len(snps2))
        for threshold in thresholds:
            cutoff = min(int(threshold * universe_size), min_snps)
            n = max(1, cutoff)
            percentage = 100.0 * float(snps1.take(n).count_intersection(snps2.take(n))) \
                         / float(n)
            common_feature_threshold_percentages[threshold].append(percentage)

    return used_models, common_feature_threshold_percentages

def write_similarity_curves(flname, thresholds, n_models, common_feature_percentages):
    with open(flname, "w") as fl:
        fl.write("threshold\t")
        fl.write("\t".join(map(str, n_models)))
        fl.write("\n")
        for t in thresholds:
            fl.write(str(t))
            fl.write("\t")
            fl.write("\t".join(map(str, common_feature_percentages[t])))
            fl.write("\n")

def plot_similarity_curves(flname_base, thresholds, n_models, common_feature_percentages):
    plt.clf()

    colors = ["r.-", "g.-", "b.-", "m.-", "c.-"]
    for i, threshold in enumerate(thresholds):
        c = colors[i]
        label = str(100.0 * threshold)
        plt.semilogx(n_models, common_feature_percentages[threshold],
                     c, label="Top %s%%" % label)

    plt.grid(True)
    plt.xlabel("Number of Models", fontsize=16)
    plt.ylabel("Overlapping SNPs (%)", fontsize=16)
    plt.legend(loc="lower right")
    plt.ylim([0, 100])
    plt.xlim([min(n_models), max(n_models)])
    plt.yticks(np.linspace(0., 100., num=11))

    plt.savefig(flname_base + ".png", DPI=200)
    plt.savefig(flname_base + ".pdf", DPI=200)

def plot_sampled_snps_curves(flname_base, used_models, common_feature_counts,
        snp1_feature_counts, snp2_feature_counts):

    plt.clf()
    plt.semilogx(used_models, common_feature_counts, "k.-", label="Common")
    plt.semilogx(used_models, snp1_feature_counts, "c.-", label="Model 1")
    plt.semilogx(used_models, snp2_feature_counts, "m.-", label="Model 2")
    plt.grid(True)
    plt.xlabel("Number of Models", fontsize=16)
    plt.ylabel("SNPs (Count)", fontsize=16)
    plt.legend(loc="lower right")
    plt.ylim([0, max(max(common_feature_counts), max(snp1_feature_counts), max(snp2_feature_counts)) + 10])
    plt.xlim([min(used_models), max(used_models)])

    plt.savefig(flname_base + ".png", DPI=200)
    plt.savefig(flname_base + ".pdf", DPI=200)

def histogram_sparse_to_dense(sparse):
    n_entries = max(sparse.keys()) + 1
    dense = [0] * n_entries
    for key, count in sparse.iteritems():
        dense[key] = count
    return dense

def plot_feature_histogram(flname, histogram):
    xs = [0]
    ys = [0]
    for i, count in enumerate(histogram):
        xs.append(i-0.5)
        ys.append(count)
        xs.append(i+0.5)
        ys.append(count)
    xs.append(i+0.5)
    ys.append(0)
    plt.clf()
    plt.plot(xs, ys, "b-")
    plt.grid(True)
    plt.xlabel("Number of Features", fontsize=16)
    plt.ylabel("Number of Trees", fontsize=16)

    plt.savefig(flname, DPI=200)
