# Genotyping Inversions

Once you've determined that some of the PCs capture inversions, you can proceed to genotyping the samples.  You can cluster samples using their coordinates along the PCs.  The clusters can reveal the genotypes of the samples.  You can use the Manhattan plots from [localizing inversions](localizing-inversions.md) to figured out which PCs capture which inversions.

## Unsupervised Genotyping with Clustering
For each inversion, we can cluster the samples using the coordinates along the PCs that capture the inversion.

For example, if an inversion is captured by PC 1, we can cluster the samples as so:

```bash
$ asaph_genotype \
    --workdir <workdir> \
	unsupervised-genotyping \
	--components 1 \
	--n-clusters 3 \
	--predicted-labels-fl predicted_labels.pops
```

If an inversion is captured by PCs 1 and 2, we can cluster the samples using both PCs:

```bash
$ asaph-genotype \
    --workdir <workdir> \
	unsupervised-genotyping \
	--components 1 2 \
	--n-clusters 3 \
	--predicted-labels-fl predicted_labels.pops
```

## Supervised Genotyping
If we have two sets of samples in which one set has known labels and the other does not, we can use supervised learning to predict the labels of the unknown set:

```bash
$ asaph_genotype \
    --workdir <workdir> \
	supervised-genotyping \
	--known-labels-fl known_labels.pops \
	--predicted-labels-fl predicted_labels.pops
```

The model will predict the label for any sample which is in the coordinates file but not the known labels file.  You should include both sets of samples when
performing PCA.

## Evaluating Predictions
If you happen to know the genotypes for your samples, you can test the cluster and other labels for agreement:

```bash
$ asaph_genotype \
	evaluate-predicted-genotypes \
	--predicted-labels-fl predicted_labels.pops \
	--output-labels-fl known_labels.pops
```

## PCA Projection Plots
We can then generate scatter plots for the PCA:

```bash
$ asaph_genotype \
    --workdir <workdir> \
	plot-projections \
	--pairs 1 2 3 4 \
	--labels-fl predicted_labels.pops \
	--plot-dir pca_plots
```

Two plot files `pca_projection_1_2.png` and `pca_projection_3_4.png` will be created in the `<workdir>/plots` directory.  The samples will be colored according to their cluster assignment.

You can also plot the projections without cluster labels:

```bash
$ asaph_genotype \
    --workdir <workdir> \
	plot-projections \
	--pairs 1 2 3 4 \
```
