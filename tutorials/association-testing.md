# Associating Testing

In our 2018 JBCB paper [Adjusted Likelihood-Ratio Test for Variants with Unknown Genotypes](https://www.worldscientific.com/doi/10.1142/S0219720018400206), we described an association test that is adjusted for missing genotypes to avoid false positives.  In this tutorial, we describe how to use Asaph to run the association tests.

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

## Association Tests
The `snp_assocation_tests` script implements several variations on our adjusted association test.  The tests can be run with our recommended defaults like so: 

```
$ bin/snp_association_tests --workdir <workdir> \
    --populations <populations file>
```

Each SNP is tested against the population labels.  The p-values will be written to tab-separated value (TSV) files under `<workdir>/statistics`.  The file will contain three columns (the chromosome, position, and p-value) with one row per SNP.  The tests from each PC will be written to its own file.

The association tests take a LONG, LONG time to run.  You are best running them in a `screen` or `tmux` session.

## Significance Testing of Scripts
Once Asaph calculates p-values for SNPs, either against a principal component or a population, we need to perform significant testing to identify which SNPs have a significant association.  The `sig_test_snps.py` can be used to do this.  Since we are testing multiple SNPs, we need to apply a correction for multiple hypothesis testing.  I liked to use the Bonferroni correction, with one caveat.  Instead of using the number of tests, I used the number of samples.  In our work, the number of samples is often far smaller than the number of SNPs and determines the true degrees of freedom.  The formula is as follows:

```
alpha_corrected = alpha / (n_samples - 1)
```

where `alpha` is your desired threshold (e.g., 0.05, 0.01).

Run the script like so:

```
$ python utils/sig_test_snps.py --input <workdir>/analysis/snp_pc_1_logreg_assoc_tests.tsv \
    --output <workdir>/analysis/snp_pc_1_logreg_assoc_tests.sig.tsv \
    --significance 0.000333
```

