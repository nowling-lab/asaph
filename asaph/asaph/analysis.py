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

def similarity_curves(thresholds, all_models):
    ordered_models = sorted(all_models.keys())

    common_feature_threshold_percentages = defaultdict(list)
    used_models = []
    for n_models in ordered_models:
        models = all_models[n_models]
        if len(models) != 2:
            continue

        snps1, snps2 = models
        used_models.append(n_models)

        for threshold in thresholds:
            n = max(1, int(threshold * min(len(snps1), len(snps2))))
            percentage = 100.0 * float(snps1.take(n).count_intersection(snps2.take(n))) \
                         / float(n)
            common_feature_threshold_percentages[threshold].append(percentage)

    return used_models, common_feature_threshold_percentages

def plot_similarity_curves(flname_base, thresholds, n_models, common_feature_percentages):
    plt.clf()

    colors = ["r.-", "g.-", "b.-", "m.-", "c.-"]
    for i, threshold in enumerate(thresholds):
        c = colors[i]
        label = str(int(100.0 * threshold))
        plt.semilogx(n_models, common_feature_percentages[threshold],
                     c, label="Top %s%%" % label)

    plt.grid(True)
    plt.xlabel("Number of Models", fontsize=16)
    plt.ylabel("Overlapping SNPs (%)", fontsize=16)
    plt.legend(loc="lower right")
    plt.ylim([0, 100])
    plt.xlim([min(n_models), max(n_models)])

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
