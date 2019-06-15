# Asaph
Asaph is a toolbox for developing and evaluating algorithms for population genetics analysis of SNPs from insect species.  Asaph is developed primarily by [RJ Nowling](http://rnowling.github.io/).

Build status: [![build status](https://travis-ci.org/rnowling/asaph.svg?branch=master)](https://travis-ci.org/rnowling/asaph)

## Getting Asaph
To download Asaph, all you need to do is clone the git repository:

    $ git clone https://github.com/rnowling/asaph.git

Asaph requires that [numpy and scipy](http://www.numpy.org/), [matplotlib](http://matplotlib.org/), [scikit-learn](http://scikit-learn.org/stable/), and [seaborn](https://seaborn.pydata.org/index.html) are installed.

## Running Tests
Asaph has a test suite implemented using the [Bats](https://github.com/sstephenson/bats) framework.  You can run the test suite like so:

    bats tests/*.bats

## Using Asaph

Tutorials for using Asaph for analyses described in our papers are located in the `tutorials` directory.

1. [Detecting Inversions with PCA](tutorials/detecting-inversions-with-pca.md)
2. [Association Testing](tutorials/association-testing.md)

## Citing Asaph
Please cite the following papers if you are using Asaph in your work:

* RJ Nowling and SJ Emrich. [Detecting Chromosomal Inversions from Dense SNPs by Combining PCA and Association Tests](/publications/ACMBCB_2018.pdf). *Proceedings of the 9th ACM International Conference on Bioinformatics, Computational Biology and Health Informatics* (ACM-BCB 2018), August 2018.
* RJ Nowling and SJ Emrich. [Adjusted Likelihood-Ratio Test for Variants with Unknown Genotypes](https://www.worldscientific.com/doi/10.1142/S0219720018400206). *Journal of Bioinformatics and Computational Biology* (JBCB), 16(5) 2018.


