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

pred_file <- snakemake@output[[1]]
dir.create(dirname(pred_file), recursive = TRUE, showWarnings = FALSE)
variants <- snakemake@wildcards[["v"]]
pheno <- snakemake@wildcards[["pheno"]]
min_ancestry <- snakemake@wildcards[["min_ancestry"]]
f <- as.integer(snakemake@wildcards[["fold"]])
type <- as.character(snakemake@wildcards[["type"]])
this_frac <- as.double(snakemake@wildcards[["frac"]])

pfile <- "data/all_vars.tsv"
### Combine chromosome predictions ###
outdir <- file.path(
    "output", "ukb",
    paste0("v~", variants),
    paste0("pheno~", pheno),
    paste0("min_ancestry~", min_ancestry),
    paste0("fold~", f),
    paste0("type~", type),
    paste0("frac~", format(this_frac, nsmall = 1))
)

chroms <- c(1:22, "X")

out_df <- foreach(chrom = chroms, .combine = rbind) %do% {

  gfile <- file.path(
    "data", "genotypes",
    paste0("v~", variants),
    paste0("min_ancestry~", min_ancestry),
    paste0("chrom~", chrom)
  )
  outfile <- file.path(outdir, paste0("chrom~", chrom, ".RDS"))

  mod <- readRDS(outfile)
  lambda_ind <- which.max(mod$metric.val)
  covars <- mod$configs$covariates
  mod$beta[[lambda_ind]][covars] <- 0
  mod$a0[[lambda_ind]] <- 0

  pred <- predict_snpnet(mod,
                         new_genotype_file = gfile,
                         new_phenotype_file = pfile,
                         phenotype = pheno,
                         covariate_names = covars,
                         idx = lambda_ind)

  tibble(ID = rownames(pred$prediction$train),
         pred = pred$prediction$train[, 1],
         chrom = chrom)

}

pred_df <- out_df %>%
  group_by(ID) %>%
  summarise(pred = sum(pred))

### Save output
write_tsv(pred_df, pred_file)
