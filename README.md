# Asaph
Software for ranking SNPs using Random Forests

Build status: ![build status](https://travis-ci.org/rnowling/asaph.svg?branch=dev)

## Importing data
To create an Asaph project, we first need to import the data.  A minimal command for importing biallelic SNPs from a VCF file would look like so:

    bin/import --vcf <path/to/vcf> \
               --populations <path/to/populations_file> \
               --workdir <path/to/workdir> \
               --features-type categories

Asaph currently supports encoding SNPs as features in two ways: categories and counts.  In the categorical scheme, each genotype (e.g., A/T, A/A, T/T) is represented as a binary variable. In the counts scheme, each allele (e.g., A, T) is represented as an integer giving the number of copies (e.g., 0, 1, 2) of each allele the sample has.

A file listing the sample ids for each population must also be specified.  The populations file contains one group per line, with the first column indicating the population name, and the sample names separated by commas like so:

    Population 1,Sample1,Sample2
    Population 2,Sample3,Sample4

The sample ids have to match the ids in the VCF file.

The work directory will be created and contain the resulting Asaph data structures such as a feature matrix, feature labels, and sample labels.

To improve run times, reduce memory usage, and stabilize rankings, Asaph supports dictionary encoding of the feature matrix.  With dictionary encoding, all SNPs with the same genotypes across all samples are replaced by a single instance.  When the rankings are generated, the original SNPs are given the variable importance score of the remaining instance.  In practice, we've found that the compression can reduce run times from weeks to hours.  Compression can be enabled with the `--compress` flag:

    bin/import --vcf <path/to/vcf> \
               --populations <path/to/populations_file> \
               --workdir <path/to/workdir> \
               --features-type categories \
               --compress

## Random Forests
Asaph's main purpose is to support calculation of variable importances scores and ranking of SNPs using Random Forests.  Once data is imported, Random Forest models can be trained with the command:

    bin/random_forests train \
                       --workdir <path/to/workdir> \
                       --trees <number of trees>


Generally, you will want to sweep over the trees parameter, so you'll run the above command with a range of values for the `trees` parameter.  Asaph actually trains two Random Forests each time, to use in checking convergence.  You can check the convergence of the SNP rankings using the `analyzing-rankings` mode:

    bin/random_forests analyze-rankings \
                       --workdir <path/to/workdir>

Each The `analyze-rankings` mode will generate two plots, comparisons of the number of SNPs used and the agreement of the top 5%, 10%, 25%, and 50% of ranked SNPs between each pair of models.  The plots are written in PDF and PNG formats and stored in `<workdir>/figures`. If the rankings do not demonstrate convergence, run the training command with a larger number of trees.  Once the rankings have converged, you can output the rankings to a text file:

    bin/random_forests output-rankings \
                       --workdir <path/to/workdir> \
                       --ranks-file <path/for/output> \
                       --trees <select model with this many trees>

## Running Tests
Asaph has a test suite implemented using the [Bats](https://github.com/sstephenson/bats) framework.  You can run the test suite like so:

    bats tests/*.bats