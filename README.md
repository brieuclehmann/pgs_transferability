This repository contains the scripts used to generate the results and plots for the preprint ['Optimal strategies for learning multi-ancestry polygenic scores vary across traits'](https://doi.org/10.1101/2021.01.15.426781).

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
