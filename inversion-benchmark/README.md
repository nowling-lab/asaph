# Inversion Benchmark Data Set

We composed a benchmark data set for evaluating SNP-based inversion detection methods.  We collected variant data from three insect data sets.  We provide links to the original data sources, a pipeline for processing the raw data into the final usable product for the benchmark, and labels of samples' inversion genotypes.

## Setup Instructions
You will need to download:

* The `dgrp2.vcf` VCF file from the [Drosophila Genetics Reference Panel v2](http://dgrp2.gnets.ncsu.edu/data.html)
* The data set for the Fontaine, et al. 2015 paper from the 16 Anopheles genomes project [DRYAD](https://datadryad.org/stash/dataset/doi:10.5061/dryad.f4114)
* The biallelic 3L, 2R, and 2L VCF files from the FTP site for the 1000 Anopheles Genomes project phase 1 AR3 data release (ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/variation/main/vcf/).  These files are named: ag1000g.phase1.ar3.pass.biallelic.{chrom}.vcf.gz where chrom is 3L, 2R, or 2L.

Place this files in the directory `data/raw_data`.  The directory contains an empty file named `PLACE_INPUT_FILES_HERE`.

You can verify that you downloaded all of the files and put them in the correct directory by running:

```bash
$ snakemake check_inputs
```

If the command complains of missing input files, then you are either missing a file, the file is not named correctly, or the file is in the wrong place.

To run the pipeline, you will need to install:

* [VCFTools](https://vcftools.github.io/index.html)
* [Plink v1.9](https://www.cog-genomics.org/plink/1.9/)

## Running Pipeline to Generate Output Files
Once you've downloaded the data, you can process the individual data sets with the following commands:

```bash
$ snakemake prepare_ag16g
$ snakemake prepare_ag1000g
$ snakemake prepare_dgrp2
```

Each task is assigned one thread.  If you want to run multiple tasks concurrently, use Snakemake's `--cores` flag:

```bash
$ snakemake --cores 4 prepare_ag16g
$ snakemake --cores 4 prepare_ag1000g
$ snakemake --cores 4 prepare_dgrp2
```

## Output File Formats
The pipeline will convert the variants into three file formats:

* VCF (potentially gzipped): This file format can be read by Asaph
* Plink bed: This file format can be read by [pcadapt](https://bcm-uga.github.io/pcadapt/index.html)
* inveRsion: This text file format can be read by [inveRsion](https://bioconductor.org/packages/release/bioc/html/inveRsion.html)

## Description of Test Cases
The data constitute three test cases:

* Negatives: these data do not have an inversions
  * DGRP2 3L: Samples with no inversions (we removed any samples with inversions) and from a single population
  * Ag16 3L: 34 _Anopheles gambiae_ and _coluzzii_ samples from four geographic areas
  * Ag1000 3L: 150 _Anopheles gambiae_ and _coluzzii_ samples from a single geographic area (Burkina Faso)
* Positives from single population
  * DGRP2 2L: 198 samples with a single inversion In(2L)t
  * DGRP2 2R: 198 samples with a single inversion In(2R)NS
  * DGRP2 3R: 198 samples with three overlapping and mutually-exclusive inversions In(3R)Mo, In(3R)p, and In(3R)k
  * Ag1000 2L _An. gambiae_: 81 samples from a single geographic area (Burkina Faso) with a single inversion 2La
  * Ag1000 2L _An. coluzzii_: 69 samples from a single geographic area (Burkina Faso) and only one sample has a different genotype for 2La
  * Ag1000 2R _An. gambiae_: 81 from a single geographic area (Burkina Faso) with the 2Rb inversion
  * Ag1000 2R _An. coluzzii_: 69 from a single geographic area (Burkina Faso) with the 2Rbc inversion system and 2Rd inversion (no labels provided)
* Positives from multiple populations
  * Ag1000 2L: 150 _Anopheles gambiae_ and _coluzzii_ samples from a single geographic area (Burkina Faso) with a single inversion 2La
  * Ag16 2L: 34 _Anopheles gambiae_ and _coluzzii_ samples from four geographic areas with a single inversion 2La
  * Ag1000 2R: 150 _Anopheles gambiae_ and _coluzzii_ samples from a single geographic area (Burkina Faso) with a combination of the 2Rb and 2Rbc inversions
  
Inversion genotype labels are provided under the `sample_labels` directory.  Note that we do not have labels for the 2Rd inversion.

## Citing

If you use this data set, please cite the original papers from which the data are derived:

* **Anopheles 1000 Genomes**: Miles, A., Harding, N., Bottà, G. et al. [Genetic diversity of the African malaria vector Anopheles gambiae.](https://doi.org/10.1038/nature24995) Nature 552, 96–100 (2017).
* **16 Anopheles Genomes**: Fontaine, M. C., Pease, J. B., Steele, A., et al. [Extensive introgression in a malaria vector species complex revealed by phylogenomics.](https://doi.org/10.1126/science.1258524) Science 347:6217. (2015).
* **Drosophila Genetics Reference Panel**: Mackay, T., Richards, S., Stone, E., et al. [The Drosophila melanogaster Genetic Reference Panel.](https://doi.org/10.1038/nature10811) Nature 482, 173–178 (2012).
* **Drosophila Genetics Reference Panel**: Huang, W., Massouras, A., Inoue, Y., et al. [Natural variation in genome architecture among 205 Drosophila melanogaster Genetic Reference Panel lines.](https://doi.org/10.1101/gr.171546.113) Genome Research 24:1193-1208 (2014).

and our paper describing the data set:

* Nowling, R.J., Manke, K.R., and Emrich, S.J. [Detection of Inversions with PCA in the Presence of Population Structure.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0240429) PLOS One (2020).
