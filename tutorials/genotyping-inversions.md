# Genotyping Inversions

Once you've determined that some of the PCs capture inversions, you can proceed to genotyping the samples.  You can cluster samples using their coordinates along the PCs.  The clusters can reveal the genotypes of the samples.  You can use the Manhattan plots from [localizing inversions](localizing-inversions.md) to figured out which PCs capture which inversions.

## Unsupervised Genotyping with Clustering
For each inversion, we can cluster the samples using the coordinates along the PCs that capture the inversion.

For example, if an inversion is captured by PC 1, we can cluster the samples as so:

```bash
$ asaph_genotype \
	cluster \
	--workdir <workdir> \
	--components 1 \
	--n-clusters 3 \
	--predicted-labels-fl predicted_labels.pops
```

If an inversion is captured by PCs 1 and 2, we can cluster the samples using both PCs:

```bash
$ asaph-genotype \
	cluster \
	--workdir <workdir> \
	--components 1 2 \
	--n-clusters 3 \
	--predicted-labels-fl predicted_labels.pops
```

## Optimizing Clustering Parameters
If you have known labels, you can find the parameters that optimize the clustering.  Each clustering is scored by how well the cluster labels predict the ground truth labels and vice versa.

```bash
$ asaph-genotype \
	sweep-parameters \
	--workdir <workdir> \
	--components 1 2 3 4 \
	--n-clusters 2 3 4 5 \
	--labels-fl known_labels.pops
```

The output will look like so:

```
Best score: 0.9722845902897717
Best parameters: {'n_clusters': 3, 'components': [1], 'feature_scaling': False}
Best centroids: [[-0.3983714 ]
 [ 3.01494853]
 [ 1.08898097]]
```


## Evaluating Predictions
If you happen to know the genotypes for your samples, you can test the cluster and other labels for agreement:

```bash
$ asaph_genotype \
	evaluate-predicted-genotypes \
	--predicted-labels-fl predicted_labels.pops \
	--output-labels-fl known_labels.pops
```

