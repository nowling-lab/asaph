# Asaph
[![Build Status](https://travis-ci.com/rnowling/asaph.svg?branch=main)](https://travis-ci.com/rnowling/asaph) [![DOI](https://zenodo.org/badge/42882932.svg)](https://zenodo.org/badge/latestdoi/42882932)

Asaph enables analysis of large inversions from SNP data.  Asaph employs streaming and dimensionality reduction to be fast and maintain a small memory footprint.  Asaph is developed primarily by [RJ Nowling](http://rnowling.github.io/) and students at [MSOE](https://msoe.edu).

## Documentation
Tutorials describing how to use Asaph and associated resources  are located in the [`tutorials`](tutorials/README.md) directory.

1. [Installing Asaph](tutorials/installing-asaph.md)
1. [Preparing, Importing, and PCA of SNPs](tutorials/pca.md)
1. [Detecting and Localizing Inversions](tutorials/localizing-inversions.md)
1. [Genotyping Inversions](tutorials/genotyping-inversions.md)

Other tutorials related to Asaph:

* [Testing Changes to Asaph](tutorials/testing-asaph.md) (mostly for developers)
* [Preparing the Inversion Benchmark Data Set](inversion-benchmark/README.md) (optional, only if you don't have your own data)

## Citing Asaph
Asaph and associated methods have been described in the following papers:

* RJ Nowling, F Fallas-Moya, et al. [Fast, low-memory detection and localization of large, polymorphic inversions from SNPs](https://peerj.com/articles/12831/). *PeerJ* 2022.
* RJ Nowling, KR Manke, and SJ Emrich. [Detecting inversions with PCA in the presence of population structure](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0240429). *PLOS One* 2020.
* RJ Nowling and SJ Emrich. [Detecting Chromosomal Inversions from Dense SNPs by Combining PCA and Association Tests](/publications/ACMBCB_2018.pdf). *Proceedings of the 9th ACM International Conference on Bioinformatics, Computational Biology and Health Informatics* (ACM-BCB 2018), August 2018.
* RJ Nowling and SJ Emrich. [Adjusted Likelihood-Ratio Test for Variants with Unknown Genotypes](https://www.worldscientific.com/doi/10.1142/S0219720018400206). *Journal of Bioinformatics and Computational Biology* (JBCB), 16(5) 2018.


