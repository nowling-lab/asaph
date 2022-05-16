# Detecting and Localizing Inversions

In the previous [tutorial](pca.md), you prepared, imported, and performed principal component analysis (PCA) on SNP data.  In this tutorial, we will perform association tests between each SNP and principal component (PC).  This will allow us to check each PC to see if it captured an inversion.  If so, we can localize the inversion using the Manhattan plot.

## Association Testing
We can now run single-SNP association tests.  The genotypes of each SNP are tested against the samples' PC coordinates.  The p-values will be written to tab-separated value (TSV) files.  The file will contain four columns (the component, the chromosome, position, and p-value) with one row for each SNP-component pair.

```bash
$ asaph_localize \
    --workdir <workdir> \
    association-tests \
    --components 1 2 \
	--vcf <path/to/vcf>
```

The association test results (p-values) will be written out to the file `<workdir>/pca_associations.tsv`.

## Manhattan Plots
Secondly, we will use manhattan plots to show the p-values of the SNPs across the chromosome.  SNPs marked in orange are statistically significant using a significance threshold of 0.01 with a Bonferroni correction based on the number of SNPs.  Spatial correlation can indicate structural variations such as inversions.  To generate a plot for the association tests against component 1, run the following:

```bash
$ asaph_localize \
    --workdir <workdir> \
    plot \
	--component 1
```

The plots will be written out to files with the prefix `manhattan` in the directory `<workdir>/plots`.

Inversions will be indicated by a step function-like pattern in the Manhattan plot.  Different karyotypes of the same inversion may be captured by separate PCs, so you may see the inversion present in more than one plot.

## Automated Boundary Detection
Asaph includes an algorithm for detecting the boundaries of an inversion.  The chromosome is divide into windows.  The number of statistically significant SNPs in each window is compared to an expected number based on a uniform distribution across the chromosome.  The window p-values are calculated using a binomial test and tested for significance using a threshold of 0.0001 with a Bonferroni correction based on the number of windows. The boundaries are determined from the left coordinate of the leftmost significant window and the right coordinate of the rightmost significant window. This will print out coordinates.

```bash
$ asaph_localize \
    --workdir <workdir> \
	detect-boundaries \
	--component 1
```

You will see output such as:

```
122 of 10000 were significant
Left boundary: 19032733
Right boundary: 30828378
```

You can evaluate the predicted boundaries against known (expected) boundaries like so:

```bash
$ asaph-localize \
    --workdir <workdir> \
	evaluate-boundaries \
	--component 1 \
	--boundaries 19032733 30828378
```

which will result in the following output:

```
122 of 10000 were significant
Left expected boundary: 26758676
Right expected boundary: 31488544

Left predicted boundary: 19032733
Right predicted boundary: 30828378

Recall: 86.0%
Precision: 34.5%
Jaccard: 32.7%
```

You can add these boundaries to your Manhattan plot like so:

```bash
$ asaph_localize \
    --workdir <workdir> \
    plot \
	--component 1 \
	--boundaries 19032733 30828378
```

## What Next?
If any of the PCs appear to capture inversions, we can move on to [predicting sample genotypes](genotyping-inversions.md).
