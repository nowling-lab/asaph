# Genotyping Inversions

Once you've determined that some of the PCs capture inversions, you can proceed to genotyping the samples.  You can cluster samples using their coordinates along the PCs.  The clusters can reveal the genotypes of the samples.  You can use the Manhattan plots from [localizing inversions](localizing-inversions.md) to figured out which PCs capture which inversions.

# Clustering Samples
For each inversion, we can cluster the samples using the coordinates along the PCs that capture the inversion.  The [DBSCAN](https://scikit-learn.org/stable/modules/clustering.html#dbscan) clustering algorithm by separating areas dense in points from those that are not.  It can also identify outliers; in our case, outliers are samples for which we cannot confidently genotype.

For example, if an inversion is captured by PC 1, we can cluster the samples as so:

```bash
$ asaph_genotype \
	cluster-samples \
	--coordinates <workdir>/pca_coordinates.tsv \
	--components 1 \
	--n-clusters 3 \
	--output-labels-fl cluster_labels.pops
```

If an inversion is captured by PCs 1 and 2, we can cluster the samples using both PCs:

```bash
$ asaph-genotype \
	cluster-samples \
	--coordinates <workdir>/pca_coordinates.tsv \
	--components 1 2 \
	--n-clusters 3 \
	--output-labels-fl cluster_labels.pops
```

Samples determined to be outliers will be assigned a cluster label of -1 in the output file.

If you happen to know the genotypes for your samples, you can test the cluster and other labels for agreement:

```bash
$ asaph_genotype \
	evaluate-clusters \
	--clusters-labels-fl cluster_labels.pops \
	--output-labels-fl known_labels.pops
```

# PCA Projection Plots
We can then generate scatter plots for the PCA:

```bash
$ asaph_genotype \
	plot-projections \
	--coordinates pca_coordinates.tsv \
	--pairs 1 2 3 4 \
	--labels-fl cluster_labels.pops \
	--plot-dir pca_plots
```

The two plot files `pca_projection_1_2.png` and `pca_projection_3_4.png` will be created.  The samples will be colored according to their cluster assignment.

You can also plot the projections without cluster labels:

```bash
$ asaph_genotype \
	plot-projections \
	--coordinates pca_coordinates.tsv \
	--pairs 1 2 3 4 \
	--plot-dir pca_plots
```
