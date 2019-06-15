# Detecting Inversions with PCA

In our 2019 ACM-BCB paper [Detecting Chromosomal Inversions from Dense SNPs by Combining PCA and Association Tests](https://dl.acm.org/citation.cfm?id=3233571), we described an approach for detecting large chromosomal inversions from SNPs.  In this tutorial, we describe how to use Asaph to run PCA, run the association tests, and make a plot.

## Getting Asaph
To download Asaph, all you need to do is clone the git repository:

```
$ git clone https://github.com/rnowling/asaph.git
```

Asaph requires that [numpy and scipy](http://www.numpy.org/), [matplotlib](http://matplotlib.org/), [scikit-learn](http://scikit-learn.org/stable/), and [seaborn](https://seaborn.pydata.org/index.html) are installed.

## Importing Data
To create an Asaph project, we first need to import the data.  A minimal command for importing biallelic SNPs from a VCF file would look like so:

```
$ bin/import --vcf <path/to/vcf> \
    --compress \
    --populations <path/to/populations_file> \
    --workdir <path/to/workdir>
```

Asaph currently supports encoding SNPs as features in two ways: categories and counts.  In the categorical scheme, each genotype (e.g., A/T, A/A, T/T) is represented as a binary variable. In the counts scheme, each allele (e.g., A, T) is represented as an integer giving the number of copies (e.g., 0, 1, 2) of each allele the sample has.  The default and recommended scheme is to use categories.

A file listing the sample ids for each population must also be specified.  The populations file contains one group per line, with the first column indicating the population name, and the sample names separated by commas like so:

```
Population 1,Sample1,Sample2
Population 2,Sample3,Sample4
```

The sample ids have to match the ids in the VCF file.  If a sample name from the VCF file is not present in the populations file, that sample is ignored.

The work directory will be created and contain the resulting Asaph data structures such as a feature matrix.

## Principal Component Analysis (PCA)
The `pca` script provides several options for training a PCA model, determining how many components to use, and outputting the coordinates.

A PCA model can be trained like so:

```
$ bin/pca --workdir <workdir> \
    train \
    --n-components 10
```

You can calculated the "explained variance ratios" like so:

```
$ bin/pca --workdir <workdir> \
    explained-variance-analysis
```

The script will print the explained variance ratios and write a plot to `<workdir>/figures/explained_variance_analysis.png`.

We can now run single-SNP association tests.  Each SNP is tested against the projected coordinates for each selected PC.  The p-values will be written to tab-separated value (TSV) files under `<workdir>/analysis`.  The file will contain three columns (the chromosome, position, and p-value) with one row per SNP.  The tests from each PC will be written to its own file.

The association tests take a LONG, LONG time to run.  You are best running them in a `screen` or `tmux` session.

```
$ bin/pca --workdir <workdir> \
    snp-association-tests
    --components 1 2 3 4
```

Asaph currently implements both Logistic and Linear Regression models for the association tests.  We recommend using the Logistic Regression-based association tests as this appears to provide greater power.  This is the default.

It may also be useful to write out the projected coordinates from the PCA for downstream analysis.  They will be written as a tab-separated value (TSV) file, one sample per row.  The first column will be the sample name and the other columns will be the coordinates along the components that you selected.

```
$ bin/pca --workdir <workdir> \
    output-coordinates \
    --selected-components 1 2 3 4 \
    --output-fl <path/to/output>
```

## Manhattan Plots
Manhattan plots show you the p-values of the SNPs across the chromosome.  Spatial correlation can indicate structural variations such as inversions.  Manhattan plots are generated using the p-values of all of the SNPs, not just SNPs that pass the significance testing above.


You can use the a script in the `utils` directory to make the Manhattan plot:

```
$ python utils/manhattan_plot.py \
    --input-tsv <workdir>/analysis/snp_pc_1_${i}_logreg_assoc_tests.tsv \
    --plot-fl <workdir>/figures/manhattan_pc_1_${i}.png
```

Inversions will be indicated by a step function-like pattern in the Manhattan plot.  Different karyotypes of the same inversion may be captured by separate PCs, so you may see the inversion present in more than one plot.
