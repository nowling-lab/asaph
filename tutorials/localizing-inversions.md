# Detecting and Localizing Inversions

In the previous [tutorial](pca.md), you prepared, imported, and performed principal component analysis (PCA) on SNP data.  In this tutorial, we will perform association tests between each SNP and principal component (PC).  This will allow us to check each PC to see if it captured an inversion.  If so, we can localize the inversion using the Manhattan plot.

## Association Testing
We can now run single-SNP association tests.  The genotypes of each SNP are tested against the samples' PC coordinates.  The p-values will be written to tab-separated value (TSV) files.  The file will contain four columns (the component, the chromosome, position, and p-value) with one row for each SNP-component pair.

```bash
$ asaph_detect_and_localize \
    association-tests \
    --pca-coordinates-tsv coordinates.tsv \
    --components 1 2 \
	--vcf <path/to/vcf> \
	--pca-associations-tsv pca_snp_tests.tsv
```
## Manhattan Plots
Secondly, we will use manhattan plots to show the p-values of the SNPs across the chromosome.  Spatial correlation can indicate structural variations such as inversions.  To generate a plot for the association tests against component 1, run the following:

```bash
$ asaph_detect_and_localize \
    plot \
    --pca-associations-tsv pca_snp_tests.tsv \
	--component 1 \
    --plot-fl manhattan_pc_1.png
```

Inversions will be indicated by a step function-like pattern in the Manhattan plot.  Different karyotypes of the same inversion may be captured by separate PCs, so you may see the inversion present in more than one plot.

## What Next?
If any of the PCs appear to capture inversions, we can move on to [predicting sample genotypes](genotyping-inversions.md).
