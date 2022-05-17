# Import and Principal Component Analysis (PCA) of SNPs

In this tutorial, we describe how to use Asaph to perform principal component analysis (PCA) of SNP data.  This is a necessary preliminary step before moving on to localizing and genotyping inversions.  We assume that you already followed the instructions in the tutorial for [installing Asaph](installing-asaph.md).

## Data Preparation
Asaph assumes that the input VCF only contains biallelic SNPs.  Further, inversion detection works best when samples are all drawn from a single population and SNPs are all from a single chromosome (or chromosome arm).  This avoids confounding factors.  You can use [VCFTools](https://vcftools.github.io/) to filter the SNPs and samples accordingly.

## Import Data and Perform PCA (Principal Component Analysis)
The Asaph analysis pipeline begins with importing data and performing PCA.  Feature hashing (the default setting) is used to construct a very small feature matrix.  A minimal command for importing biallelic SNPs from a VCF file would look like so:

```bash
$ asaph_pca \
	--workdir <workdir> \
	pca \
	--vcf <path/to/vcf>
```

The work directory will be created and contain the resulting Asaph data structures such as the sample PC coordinates and names.

If your inversion is not detected in the later stages, you can adjust the minimum inversion size (as a fraction of SNPs) to change the number of dimensions.  The default is 10%.

```bash
$ asaph_pca \
	--workdir <workdir> \
	pca \
	--vcf <path/to/vcf> \
	--min-inversion-fraction 0.01
```

## Outputing PCA Coordinates
The PCA coordinates for each sample for use in the detection, localization, and genotyping steps will automatically be output to a file named `<workdir>/pca_coordinates.tsv`.  The file will look like so:

```
sample 	1	2	3	4
line_21	-0.5352247448457047	-0.28348205963519324	-0.06379920847846524	-0.22146000350907072
line_26	2.1680548066430783	-1.8192455949345472	0.7504759437957356	-0.21425896651796422
line_28	-0.48820111599841115	-0.25884157351597736	0.05877259706868454	-0.4690312566354758
line_31	-0.36233433754549305	-0.4606695545925513	-1.0556892238848392	-0.21483260734307102
line_32	2.4899268824159315	-2.177898665769419	-0.15243794901799274	-0.13712827848692657
line_38	-0.5689557880617656	-0.34942765856267216	-0.2130572294672218	-0.5187670246143677
```

## PCA Projection Plots
We can then generate scatter plots for the PCA:

```bash
$ asaph_pca \
    --workdir <workdir> \
	plot-projections \
	--pairs 1 2 3 4
```

Two plot files `pca_projection_1_2.png` and `pca_projection_3_4.png` will be created in the `<workdir>/plots` directory.

You can also supply a labels file to color the points by population, genotype, etc..

```bash
$ asaph_genotype \
    --workdir <workdir> \
	plot-projections \
	--pairs 1 2 3 4 \
	--labels-fl predicted_labels.pops
```

The labels file format is as follows:

```
label_1,sample_1,sample_2,sample_3
label_2,sample_4,sample_5,sample_6
label_3,sample_7,sample_8
```

One group per line.  First entry is the label name.  The remaining entries on the line are the sample ids and must match the VCF file.  The labels file may contain sample ids that are not present in the VCF file but the other way around will result in an error.  Entries are separated by commas.

## More Details
Asaph provides two ways (allele counts and genotype categories) of encoding SNPs as features.  In the first approach, a separate column in the feature matrix is created for each allele.  For example, if using biallelic SNPs and a site has "A" and "T" alleles, then two columns will be created.  The columns will store the number of copies of each allele.  For diploid organisms, this means the two columns will have values of (0, 2), (2, 0), (1, 1), or (0, 0).  For the second approach, a separate column is created for each genotype.  If using biallelic SNPs and a site has "A" and "T" alleles, then three columns correspond to "A/A", "A/T", and "T/T" will be created.  These columns are treated as mutually exclusive so only one column will have a 1 for each sample.  If the genotype is unknown for a sample, then all three columns will have values of 0.

Asaph supports three ways (feature hashing, bottom-k sketching, and reservoir sampling) of subsampling variants.  For feature hashing and bottom-k sketching, a string is generated for each column such as "2L\_5453\_A" (allele counts), "2L\_345345\_T" (allele counts), or "2L\_345345\_homo\_0" (genotype categories).  With feature hashing, the number of dimensions is specified ahead of time and columns are mapped to a particular feature as `feature_idx = abs(hash(s)) % n_dimensions`.  This will cause collisions in which the values of multiple variants will be summed up.  This is a form of lossy compression.  For bottom-k sketching, only the `n_dimensions` columns with the smallest hash values are kept.   For reservoir sampling, `n_dimensions` columns are chosen randomly with uniform probability.

By default, Asaph uses allele counts with bottom-k sketching to encode the feature matrix.  The number of dimensions is calculated from the Johnson–Lindenstrauss lemma based on the number of samples detected and the expected minimum fraction of the chromosome that the inversion occupies.  The default assumption is that an inversion occupies 10% of a chromosome.

## What's Next?
Now that the data has been prepared, imported, and PCA was performed, you can move on to [detecting and localizing inversions](localizing-inversions.md).
