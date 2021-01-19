This repository contains the scripts used to generate the results and plots for the paper ['High trait variability in optimal polygenic prediction strategy within multiple-ancestry cohorts'](https://doi.org/10.1101/2021.01.15.426781).

Scripts should be run in order. Please note the dependencies down below. 

In the sim_04_run.R script, you can vary 
```
iter <- 1 # vary from 1 to 50
beta_cor <- 0.5 # takes values in 0.5 0.6 0.7 0.8 0.9
prop_afr <- 0.1 # takes values in 0.1 0.2 0.3 0.4 0.5
```

In the ukb_03 - ukb_08 scripts, you can vary
```
pheno <- "height"
prop_min <- 0.1 # one of 0 (White only), 1 (Black only) or 0.1
f <- 1 # outer CV fold, vary from 1 to 5
```
In ukb_05b_lasso_cv.R, there is also the inner CV-fold variable `i`.

Please feel free to send any questions to brieuc dot lehmann at bdi dot ox dot ac dot uk.

# Dependencies

## Simulations

Simulated DNA sequences were generated using the [stdpopsim](https://stdpopsim.readthedocs.io/en/latest/) toolkit.

The simulated DNA sequences were processed using the python tree sequence toolkit [tskit](https://tskit.readthedocs.io/en/latest/).

## UK Biobank

The data is available to researchers - see [UK Biobank](https://www.ukbiobank.ac.uk) for details on how to apply for access. This analysis was performed under Application 12788.

Genotypes were preprocessed using [plink](https://www.cog-genomics.org/plink/) and [plink2](https://www.cog-genomics.org/plink/2.0/).

Phenotypes were preprocessed using [fmrib-unpack](https://pypi.org/project/fmrib-unpack/).

The LASSO is fit using a modified version of [snpnet](https://github.com/junyangq/snpnet), which must be installed using:
```
devtools::install_github("brieuclehmann/snpnet")
```
See [snpnet](https://github.com/junyangq/snpnet) for further installation requirements (zstdcat, pgenlibr, cindex).
