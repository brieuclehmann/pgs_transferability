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
variants <- snakemake@wildcards[["v"]]
pheno <- snakemake@wildcards[["pheno"]]
prop_min <- snakemake@wildcards[["prop_min"]]
min_ancestry <- snakemake@wildcards[["min_ancestry"]]
f <- as.integer(snakemake@wildcards[["fold"]])
pow <- as.double(snakemake@wildcards[["pow"]])

pfile <- "data/all_vars.tsv"

### Combine chromosome predictions ###
outdir <- file.path(
  "output", "ukb",
  paste0("v~", variants),
  paste0("pheno~", pheno),
  paste0("min_ancestry~", min_ancestry),
  paste0("prop_min~", prop_min),
  paste0("fold~", f),
  paste0("pow~", format(pow, nsmall = 1))
)

chroms <- c(1:22, "X")

pan_ukbb_variant_file <- "data/full_variant_qc_metrics.txt"
pan_ukbb_df <- read_tsv(
  pan_ukbb_variant_file,
  col_types = cols_only(
    chrom = "c", pos = "i", ref = "c", alt = "c", rsid = "c", varid = "c",
    high_quality = "l", info = "d",
    af_AFR = "d", af_AMR = "d", af_CSA = "d", af_EAS = "d", af_EUR = "d",
    af_MID = "d"
  )
)

this_df <- pan_ukbb_df %>%
  filter(high_quality) %>%
  rename(maf = paste0("af_", min_ancestry)) %>%
  mutate(
    maf_maj = pmin(af_EUR, 1 - af_EUR),
    maf_min = pmin(maf, 1 - maf)
  ) %>%
  filter(maf_min >= 0.01 | maf_maj >= 0.01) %>%
  mutate(varid2 = paste0(chrom, ":", pos, "_", ref, "_", alt, "_", alt))

beta_df <- foreach(chrom = chroms, .combine = rbind) %do% {

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

  beta_fit <- mod$beta[[lambda_ind]][- (seq_along(covars))]
  beta_nonzero <- beta_fit[beta_fit != 0]

  tibble(varname = names(beta_nonzero),
         beta = beta_nonzero,
         chrom = chrom)
}

beta_df <- beta_df %>%
  left_join(this_df, by = c("varname" = "varid2", "chrom")) %>%
  select(varname, chrom, beta, maf_min, maf_maj)

### Save output
#beta_file <- file.path(outdir, "beta.tsv")
write_tsv(beta_df, beta_file)
