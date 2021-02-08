# Import and Principal Component Analysis (PCA) of SNPs

In this tutorial, we describe how to use Asaph to perform principal component analysis (PCA) of SNP data.  This is a necessary preliminary step before moving on to localizing and genotyping inversions.  We assume that you already followed the instructions in the tutorial for [installing Asaph](installing-asaph.md).

## Data Preparation
Asaph assumes that the input VCF only contains biallelic SNPs.  Further, inversion detection works best when samples are all drawn from a single population and SNPs are all from a single chromosome (or chromosome arm).  This avoids confounding factors.  You can use [VCFTools](https://vcftools.github.io/) to filter the SNPs and samples accordingly.

## Importing Data and Perform PCA (Principal Component Analysis)
The Asaph analysis pipeline begins with importing data and performing PCA.  Feature hashing (the default setting) is used to construct a very small feature matrix.  A minimal command for importing biallelic SNPs from a VCF file would look like so:

```bash
$ asaph_pca \
	--workdir <workdir> \
    train \
	--vcf <path/to/vcf>
```

The work directory will be created and contain the resulting Asaph data structures such as the sample PC coordinates and names.

If your inversion is not detected in the later stages, you can adjust the minimum inversion size (as a fraction of SNPs) to change the number of dimensions.  The default is 10%.

```bash
$ asaph_pca \
	--workdir <workdir> \
    train \
	--vcf <path/to/vcf> \
	--min-inversion-fraction 0.01
```

## Outputing PCA Coordiantes
We will output the PCA coordinates for each sample for use in the detection, localization, and genotyping steps.  To output the coordinates:

```bash
$ asaph_pca \
	--workdir <workdir> \
	output-coordinates \
	--components 1 2 3 4 \
	--output-fl <workdir>/pca_coordinates.tsv
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

In this example, the resulting `pca_coordinates.tsv` file will be located in the working directory that you specified.

## What Next?
Now that the data has been prepared, imported, and PCA was performed, you can move on to [detecting and localizing inversions](localizing-inversions.md).
