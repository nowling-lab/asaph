# Asaph
Software for ranking SNPs using Random Forests

## Importing data
To create an Asaph project, we first need to import the data.  You can use the import mode for this like so:

    python driver.py --mode import \
                     --vcf <path/to/vcf> \
                     --groups <path/to/groups_file> \
                     --workdir <path/to/workdir> 

For input, we need a VCF file containing biallelic SNPs and a file specifying which populations (groups) the samples belong to.  The groups file contains one group per line, with the first column indicating the group name, and the sample names separated by commas like so:

    Group 1,Sample1,Sample2
    Group 2,Sample3,Sample4

The work directory will be created and contain the resulting Asaph data structures such as a feature matrix, feature labels, and sample labels.

Asaph supports two more features which can be useful: imputing unknown genotypes and compressing the feature matrix using dictionary encoding.  In cases where genotypes for an individual SNP are known for all but a few samples, Asaph supports imputing the missing genotypes using the most common value found in the sample's class.  The `--impute-unknown` flag can be used to enable imputation.  Imputation is only performed if the frequency of the majority genotype is above a user-provided threshold to prevent inaccurate results.

Secondly, to improve run times, reduce memory usage, and stabilize rankings, Asaph supports dictionary encoding of the feature matrix.  With dictionary encoding, all SNPs with the same genotypes across all samples are replaced by a single instance.  When the rankings are generated, the original SNPs are given the variable importance score of the remaining instance.  In practice, we've found that the compression can reduce run times from weeks to hours.  Compression can be enabled with the `--compress` flag.

## Random Forests

Asaph's main purpose is to support calculation of variable importances scores and ranking of SNPs using Random Forests.  Once data is imported, Random Forest models can be trained with the command:

    python driver.py --mode train \
    		     --workdir <path/to/workdir> \
		     --trees <number of trees>


Generally, you will want to sweep over the trees parameter, so you'll run the above command with a range of values for the `trees` parameter.  Asaph actually trains two Random Forests each time, to use in checking convergence.  You can check the convergence of the SNP rankings using the `analyzing-rankings` mode:

    python driver.py --mode  analyze-rankings \
    		     --workdir <path/to/workdir>

Each The `analyze-rankings` mode will generate two plots, comparisons of the number of SNPs used and the agreement of the top 5%, 10%, 25%, and 50% of ranked SNPs between each pair of models.  The plots are written in PDF and PNG formats and stored in `<workdir>/figures`. If the rankings do not demonstrate convergence, run the training command with a larger number of trees.  Once the rankings have converged, you can output the rankings to a text file:

    python driver.py --mode output-rankings \
    		     --workdir <path/to/workdir> \
		     --ranks-file <path/for/output> \
		     --trees <select model with this many trees>
