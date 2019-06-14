# Asaph
Software for population genetics analysis of SNPs from insect species

Build status: [![build status](https://travis-ci.org/rnowling/asaph.svg?branch=master)](https://travis-ci.org/rnowling/asaph)

## Getting Asaph
To download Asaph, all you need to do is clone the git repository:

    $ git clone https://github.com/rnowling/asaph.git

Asaph requires that [numpy and scipy](http://www.numpy.org/), [matplotlib](http://matplotlib.org/), [scikit-learn](http://scikit-learn.org/stable/), and [seaborn](https://seaborn.pydata.org/index.html) are installed.

## Importing data
To create an Asaph project, we first need to import the data.  A minimal command for importing biallelic SNPs from a VCF file would look like so:

    bin/import --vcf <path/to/vcf> \
               --populations <path/to/populations_file> \
               --workdir <path/to/workdir> \
               --features-type categories

Asaph currently supports encoding SNPs as features in two ways: categories and counts.  In the categorical scheme, each genotype (e.g., A/T, A/A, T/T) is represented as a binary variable. In the counts scheme, each allele (e.g., A, T) is represented as an integer giving the number of copies (e.g., 0, 1, 2) of each allele the sample has.  We recommend using the categories encoding for Random Forests and the counts encoding for Logistic Regression ensembles.

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
               
## Single-SNP Association Tests
Asaph can perform single-SNP association tests using Likelihood-Ratio Tests on Logistic Regression models. Once data is imported:

    bin/snp_association_tests --workdir <path/to/workdir>
    
The p-value for each SNP is written to a tab-separated value file under 

    <workdir>/statistics/snp_lrtests_gt.tsv
    
OR

    <workdir>/statistics/snp_lrtests_pop.tsv

## SNP Rankings with Random Forests Variable Importance Scores
Asaph's original purpose, which it has since outgrown, was to support calculation of variable importances scores and ranking of SNPs using Random Forests.  Once data is imported, Random Forest models can be trained with the command:

    bin/random_forests --workdir <path/to/workdir> \
                       train \
                       --trees <number of trees> \
                       --populations <populations file>


Generally, you will want to sweep over the trees parameter, so you'll run the above command with a range of values for the `trees` parameter.  Asaph actually trains two Random Forests each time, to use in checking convergence.  You can check the convergence of the SNP rankings using the `analyzing-rankings` mode:

    bin/random_forests --workdir <path/to/workdir> \
                       analyze-rankings
                       

The `analyze-rankings` mode will generate two plots, comparisons of the number of SNPs used and the agreement of the top 0.01%, 0.1%, 1%, and 10% of ranked SNPs between each pair of models.  The plots are written in PDF and PNG formats and stored in `<workdir>/figures`. If the rankings do not demonstrate convergence, run the training command with a larger number of trees.  Once the rankings have converged, you can output the rankings to a text file:

    bin/random_forests --workdir <path/to/workdir> \
                       output-rankings \
                       --trees <select model with this many trees>

The rankings will be output to a text file in the `<workdir>/rankings` directory.

## SNP Ranking From Logistic Regression Weights (Ridge and Lasso)
Asaph can also be used for training ensembles of Logistic Regression models.  By training an ensemble and averaging over the feature weights, we can ensure that the rankings of the SNPs are consistent.  The LR workflow follows the RF workflow.  Once data is imported, you can train a LR model like so:

    bin/logistic_regression --workdir <path/to/workdir> \
                            train \
                            --populations <populations file>
                            --n-models <number of models>

Convergence of the SNP rankings can be evaluated like with the command:

    bin/logistic_regression --workdir <path/to/workdir> \
                            analyze-rankings

The `analyze-rankings` mode will generate two plots, comparisons of the number of SNPs used and the agreement of the top 0.01%, 0.1%, 1%, and 10% of ranked SNPs between each pair of models.  The plots are written in PDF and PNG formats and stored in `<workdir>/figures`. If the rankings do not demonstrate convergence, run the training command with a larger number of models.  Once the rankings have converged, you can output the rankings to a text file:

    bin/logistic_gression --workdir <path/to/workdir> \
                          output-rankings \
                          --n-models <select ensemble with this many models>

The rankings will be output to a text file in the `<workdir>/rankings` directory.

By default, LR models are trained with Stochastic Gradient Descent (SGD) and a L2 penalty.  Asaph also supports the average Stochastic Gradient Descent (ASGD) and Stochastic Average Gradient Descent (SAG) optimization algorithms. SGD additionally supports the elastic-net penalty. You can select different methods by using the `--methods` flag. If you do so, you'll need to use `--methods` each time you invoke `train`, `analyze-rankings`, and `output-rankings`.

You can also enable bagging, where the dataset is bootstrapped before each model is trained, with the `--bagging` flag for the `train` function. Bagging is disabled by default since we have found little impact from its usage.

## Running Tests
Asaph has a test suite implemented using the [Bats](https://github.com/sstephenson/bats) framework.  You can run the test suite like so:

    bats tests/*.bats
