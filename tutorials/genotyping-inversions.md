# Genotyping Inversions

Once you've determined that some of the PCs capture inversions, you can proceed to genotyping the samples.  You can cluster samples using their coordinates along the PCs.  The clusters can reveal the genotypes of the samples.  You can use the Manhattan plots from [localizing inversions](localizing-inversions.md) to figured out which PCs capture which inversions.

# Clustering Samples
For each inversion, we can cluster the samples using the coordinates along the PCs that capture the inversion.  The [DBSCAN](https://scikit-learn.org/stable/modules/clustering.html#dbscan) clustering algorithm by separating areas dense in points from those that are not.  It can also identify outliers; in our case, outliers are samples for which we cannot confidently genotype.

For example, if an inversion is captured by PC 1, we can cluster the samples as so:

```bash
$ asaph_clustering \
	cluster-samples-dbscan \
	--coordinates <workdir>/pca_coordinates.tsv \
	--components 1 \
	--output-labels-fl cluster_labels.pops
```

If an inversion is captured by PCs 1 and 2, we can cluster the samples using both PCs:

```bash
$ asaph_clustering \
	--coordinates <workdir>/pca_coordinates.tsv \
	cluster-samples-dbscan \
	--components 1 2 \
	--output-labels-fl cluster_labels.pops
```

Samples determined to be outliers will not be assigned a cluster label in the output file.

# PCA Projection Plots
We can then generate scatter plots for the PCA:

```bash
$ asaph_pca_analysis --coordinates pca_coordinates.tsv \
	plot-projections \
	--pairs 1 2 3 4 \
	--labels-fl cluster_labels.pops \
	--plot-dir pca_plots
```

The two plot files `pca_projection_1_2.png` and `pca_projection_3_4.png` will be created.  The samples will be colored according to their cluster assignment.

# Comparing Cluster Labels to Known Labels
If you happen to know the genotypes for your samples, you can test the cluster and other labels for agreement:

```bash
$ asaph_clustering \
	test-clusters \
	--clusters-labels-fl cluster_labels.pops \
	--output-labels-fl known_labels.pops
```
