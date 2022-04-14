# Script to generate test set predictions

### Set up parameters ---

.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(snpnet)
library(purrr)
library(foreach)
library(readr)
library(dplyr)
set.seed(1)

beta_file <- snakemake@output[[1]]
dir.create(dirname(beta_file), recursive = TRUE, showWarnings = FALSE)
pheno <- snakemake@wildcards[["pheno"]]
min_ancestry <- snakemake@wildcards[["min_ancestry"]]
f <- as.integer(snakemake@wildcards[["fold"]])
type <- as.character(snakemake@wildcards[["type"]])
this_frac <- as.double(snakemake@wildcards[["frac"]])

pfile <- "data/all_vars.tsv"
### Combine chromosome predictions ###
outdir <- file.path(
  "output", "ukb",
  paste0("pheno~", pheno),
  paste0("min_ancestry~", min_ancestry),
  paste0("fold~", f),
  paste0("type~", type),
  paste0("frac~", format(this_frac, nsmall = 1))
)

chroms <- c(1:22, "X")

beta_df <- foreach(chrom = chroms, .combine = rbind) %do% {

  gfile <- file.path(
    "data", "genotypes",
    paste0("ancestry~", min_ancestry),
    paste0("chrom~", chrom)
  )
  outfile <- file.path(outdir, paste0("chrom~", chrom, ".RDS"))

  mod <- readRDS(outfile)
  lambda_ind <- which.max(mod$metric.val)
  covars <- mod$configs$covariates

  beta_fit <- mod$beta[[lambda_ind]][-(1:length(covars))]
  beta_nonzero <- beta_fit[beta_fit != 0]
  

  tibble(varname = names(beta_nonzero),
         beta = beta_nonzero,
         chrom = chrom)

}

### Save output
#beta_file <- file.path(outdir, "beta.tsv")
write_tsv(beta_df, beta_file)
