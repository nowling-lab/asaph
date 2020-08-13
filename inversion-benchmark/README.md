# Inversion Benchmark Data Set

We composed a benchmark data set for evaluating SNP-based inversion detection methods.  We collected variant data from three insect data sets.  We provide links to the original data sources, a pipeline for processing the raw data into the final usable product for the benchmark, and labels of samples' inversion genotypes.

## Setup Instructions
You will need to download:

* The `dgrp2.vcf` VCF file from the [Drosophila Genetics Reference Panel v2](http://dgrp2.gnets.ncsu.edu/data.html)
* The data set for the Fontaine, et al. 2015 paper from the 16 Anopheles genomes project [DRYAD](https://datadryad.org/stash/dataset/doi:10.5061/dryad.f4114)
* The biallelic 3L, 2R, and 2L VCF files from the FTP site for the 1000 Anopheles Genomes project [phase 1 AR3 data release](ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/variation/main/vcf/).  These files are named: ag1000g.phase1.ar3.pass.biallelic.{chrom}.vcf.gz where chrom is 3L, 3R, or 2L.

Place this files in the directory `data/raw_data`.  The directory contains an empty file named `PLACE_INPUT_FILES_HERE`.

You can verify that you downloaded all of the files and put them in the correct directory by running:

```bash
$ snakemake check_inputs
```

If the command complains of missing input files, then you are either missing a file, the file is not named correctly, or the file is in the wrong place.

To run the pipeline, you will need to install:

* [VCFTools](https://vcftools.github.io/index.html)

## Preparing VCF Files
Once you've downloaded the data, you can process all of the VCF files with the following command:

```bash
$ snakemake prepare_vcfs
```
The Snakemake pipeline will filter the VCFs to retain only biallelic sites and the appropriate samples.  When complete, the three test cases will each be placed in their own directory:

* `data/negatives`
* `data/positives_single_pop`
* `data/positives_multiple_pops`

## Citing

If you use this benchmark and its, please cite the original papers from which the data are derived:

* **Anopheles 1000 Genomes**: Miles, A., Harding, N., Bottà, G. et al. [Genetic diversity of the African malaria vector Anopheles gambiae.](https://doi.org/10.1038/nature24995) Nature 552, 96–100 (2017).
* **16 Anopheles Genomes**: Fontaine, M. C., Pease, J. B., Steele, A., et al. [Extensive introgression in a malaria vector species complex revealed by phylogenomics.](https://doi.org/10.1126/science.1258524) Science 347:6217. (2015).
* **Drosophila Genetics Reference Panel**: Mackay, T., Richards, S., Stone, E., et al. [The Drosophila melanogaster Genetic Reference Panel.](https://doi.org/10.1038/nature10811) Nature 482, 173–178 (2012).
* **Drosophila Genetics Reference Panel**: Huang, W., Massouras, A., Inoue, Y., et al. [Natural variation in genome architecture among 205 Drosophila melanogaster Genetic Reference Panel lines.](https://doi.org/10.1101/gr.171546.113) Genome Research 24:1193-1208 (2014).
