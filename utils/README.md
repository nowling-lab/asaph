# Utility Scripts
This directory contains several small but useful scripts I created for analyzing the output of Asaph.

* `sig_test_snps.py`: Perform significance testing of SNPs based on the p-values estimated by the association tests.
* `manhattan_plot.py`: Plot the p-values of the SNPs along a chromosome arm
* `split_by_chromosome.py`: Select SNP data from an association test file from a single chromosome. Supports renaming the chromosome in the output.

## Significance Testing of Scripts
Once Asaph calculates p-values for SNPs, either against a principal component or a population, we need to perform significant testing to identify which SNPs have a significant association.  The `sig_test_snps.py` can be used to do this.  Since we are testing multiple SNPs, we need to apply a correction for multiple hypothesis testing.  I liked to use the Bonferroni correction, with one caveat.  Instead of using the number of tests, I used the number of samples.  In our work, the number of samples is often far smaller than the number of SNPs and determines the true degrees of freedom.  The formula is as follows:

```
alpha_corrected = alpha / (n_samples - 1)
```

where `alpha` is your desired threshold (e.g., 0.05, 0.01).

Run the script like so:

```
$ python sig_test_snps.py --input <workdir>/analysis/snp_pc_1_logreg_assoc_tests.tsv \
                          --output <workdir>/analysis/snp_pc_1_logreg_assoc_tests.sig.tsv \
                          --significance 0.000333
```

## Manhattan Plots
Manhattan plots show you the p-values of the SNPs across the chromosome.  Spatial correlation can indicate structural variations such as inversions.  Manhattan plots are generated using the p-values of all of the SNPs, not just SNPs that pass the significance testing above.

If you did association testing using SNPs from across the genome, you will need to first divide the SNPs by chromosome:

```
$ for i in 2L 2R X 3R 3L;
  do
      python split_by_chromsome.py --input <workdir>/analysis/snp_pc_1_logreg_assoc_tests.tsv \
                                   --output <workdir>/analysis/snp_pc_1_${i}_logreg_assoc_tests.tsv \
                                   --select-id $i \
                                   --output-id $i
  done
```

You can use the included script like so:

```
$ for i in 2L 2R X 3R 3L;
  do
      python manhattan_plot.py --input-tsv <workdir>/analysis/snp_pc_1_${i}_logreg_assoc_tests.tsv \
                               --plot-fl <workdir>/figures/manhattan_pc_1_${i}.png
  done
```
