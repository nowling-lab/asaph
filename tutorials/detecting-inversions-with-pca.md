# Detecting Inversions

In our 2018 ACM-BCB paper [Detecting Chromosomal Inversions from Dense SNPs by Combining PCA and Association Tests](https://dl.acm.org/citation.cfm?id=3233571), we described an approach for detecting large chromosomal inversions from SNPs.  In this tutorial, we describe how to use Asaph to run PCA, the association tests, and interpret the results.  We assume that you already followed the instructions in the tutorial for [installing Asaph](installing-asaph.md).

## Importing Data
To create an Asaph project, we first need to import the data.  We will use feature hashing to construct a very small feature matrix.  A minimal command for importing biallelic SNPs from a VCF file would look like so:

```bash
$ asaph_import --workdir <path/to/workdir> --feature-type hashed --subsampling-method hashing --vcf <path/to/vcf>
```

If you only wish to use a subset of the samples in the VCF file, you can provide a populations file using the `--selected-samples` flag.  The populations file contains one group per line, with the first column indicating the population name, and the sample names separated by commas like so:

```
Population 1,Sample1,Sample2
Population 2,Sample3,Sample4
```

The sample ids have to match the ids in the VCF file.  If a sample name from the VCF file is not present in the populations file, that sample is ignored.

The work directory will be created and contain the resulting Asaph data structures such as a feature matrix.

## Principal Component Analysis (PCA)
In the next step, we perform PCA on the feature matrix. The `asaph_pca` script facilitates this:

```bash
$ asaph_pca --workdir <workdir> \
    train \
    --n-components 10
```

You can plot the "explained variance ratios":

```bash
$ asaph_pca --workdir <workdir> \
    explained-variance-analysis
```

The script will print the explained variance ratios and write a plot to `<workdir>/figures/explained_variance_analysis.png`.

We will output the PCA coordinates for each sample and create scatter plots.  To output the coordinates:

```bash
$ asaph_pca --workdir <workdir> \
	output-coordinates \
	--components 1 2 3 4 \
	--output-fl pca_coordinates.tsv
```

The file will look like so:

```
sample 	1	2	3	4
line_21	-0.5352247448457047	-0.28348205963519324	-0.06379920847846524	-0.22146000350907072
line_26	2.1680548066430783	-1.8192455949345472	0.7504759437957356	-0.21425896651796422
line_28	-0.48820111599841115	-0.25884157351597736	0.05877259706868454	-0.4690312566354758
line_31	-0.36233433754549305	-0.4606695545925513	-1.0556892238848392	-0.21483260734307102
line_32	2.4899268824159315	-2.177898665769419	-0.15243794901799274	-0.13712827848692657
line_38	-0.5689557880617656	-0.34942765856267216	-0.2130572294672218	-0.5187670246143677
```

We can create scatter plots using the built-in PCA analysis script:

```bash
$ asaph_pca_analysis --coordinates pca_coordinates.tsv \
	plot-projections \
	--pairs 1 2 3 4 \
	--plot-dir pca_plots
```

The two plot files `pca_projection_1_2.png` and `pca_projection_3_4.png` will be created.

## Association Tests
We can now run single-SNP association tests.  The genotypes of each SNP are tested against the samples' PC coordinates.  The p-values will be written to tab-separated value (TSV) files.  The file will contain four columns (the component, the chromosome, position, and p-value) with one row for each SNP-component pair.

```bash
$ asaph_pca_association_tests --workdir <workdir> \
    --components 1 2 \
	--vcf <path/to/vcf> \
	--output-tsv pca_snp_tests.tsv
```
## Manhattan Plots
Lastly, we will use manhattan plots to show the p-values of the SNPs across the chromosome.  Spatial correlation can indicate structural variations such as inversions.  To generate a plot for the association tests against component 1, run the following:

```bash
$ asaph_manhattan_plot \
    --input-tsv pca_snp_tests.tsv \
	--component 1 \
    --plot-fl manhattan_pc_1.png
```

Inversions will be indicated by a step function-like pattern in the Manhattan plot.  Different karyotypes of the same inversion may be captured by separate PCs, so you may see the inversion present in more than one plot.

## What Next?
Now that Asaph is installed, check out some of our other [tutorials](README.md).
