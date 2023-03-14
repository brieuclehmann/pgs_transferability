This repository contains the scripts used to generate the results and plots for the preprint ['Optimal strategies for learning multi-ancestry polygenic scores vary across traits'](https://doi.org/10.1101/2021.01.15.426781).

# Installation

To run these scripts, you will need R version 4.0.0 or later, widely available on Unix-like, Windows and Mac families of operating systems. If you have R version 4.0.0 or above on Windows, you will also need to install Rtools. The analyses for the manuscript were performed on a CentOS Linux 7 system. The dependencies for the analysis scripts are listed below. Installation for the main snpnet package and all dependencies may take up to 30 minutes, depending on the internet connection and computer.

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
See [snpnet](https://github.com/junyangq/snpnet) for further installation requirements (zstdcat, pgenlibr, cindex, plink2).

# Demo

We have included a toy dataset (a subset of the simulated dataset) to try out the code on. See `scripts/demo.R` for some example code. The entire script will run the code on a dual-ancestry training set and produce a data frame of predictions for a test set of individuals from both ancestries. It should take roughly 2 minutes to run.

# Workflow to reproduce manuscript results

All scripts to reproduce the manuscript results are contained in the `scripts` directory. We made extensive use of the [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system to perform the analyses on a research computing cluster. See the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) and [this useful github repository](https://github.com/Snakemake-Profiles/slurm) for help on setting up cluster profiles with snakemake.

The main snakefiles are `Snakefile_simulations.smk` and `Snakefile_ukb.smk`. The relevant conda environments can be found in the `envs` directory. took several weeks to run using several hundreds of CPUs. You may want to run on a particular subset of the input arguments that you are interested in!


## Simulations
To reproduce the simulation study, first run the `scripts/simulations/01_set_sim_params.R` in R, and then run the simulations snakefile in a shell as follows:

```
snakemake --snakefile Snakefile_simulations.smk --use-conda simul_wrapper
```

Finally, run `scripts/simulations/05_combine_output.R` to combine all the simulation output ahead of plotting.

## UK Biobank

To reproduce the UK Biobank analyses, first run the `scripts/simulations/00_set_sim_params.R` in R, and then run the UKB analyses snakefiles in a shell as follows:

```
snakemake --snakefile Snakefile_simulations.smk --use-conda ukb_wrapper
snakemake --snakefile Snakefile_simulations.smk --use-conda ss_wrapper
```

## Plots

The scripts in the `scripts/plots` directory contain all the plotting code to reproduce the manuscript figures, based on the output of both the simulations and UK Biobank (`ukb`) workflows.
