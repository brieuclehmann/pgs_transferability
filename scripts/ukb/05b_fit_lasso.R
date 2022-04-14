# This script fits the LASSO using a validation set to select lambda.

### Set up parameters ---
.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(snpnet)
library(readr)
library(dplyr)
set.seed(1)

outfile <- snakemake@output[[1]]
outdir <- dirname(outfile)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
gfile <- tools::file_path_sans_ext(snakemake@input[[1]])
kfile <- snakemake@input[[2]]

pheno <- snakemake@wildcards[["pheno"]]
min_ancestry <- snakemake@wildcards[["min_ancestry"]]
chrom <- snakemake@wildcards[["chrom"]]
frac <- snakemake@wildcards[["frac"]]
type <- snakemake@wildcards[["type"]]
f <- as.integer(snakemake@wildcards[["fold"]])
ncores <- as.integer(Sys.getenv("NSLOTS"))

pfile  <- "data/all_vars.tsv"
covars <- c("sex", "age", paste0("PC", 1:10), paste0("PC", 1:10, "_sex"))

keep_ids <- read.table(kfile)[, 1]
ids <- readIDsFromPsam(paste0(gfile, ".psam"))
id_df <- data.frame(ID = ids, stringsAsFactors = F) %>%
  mutate(sort_order = 1:n())


results.dir <- file.path("temp", substr(outfile, 1, nchar(outfile) - 4))
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)

glmnetPlus_flag <- TRUE
if (type == "maj" & frac == "0.0") glmnetPlus_flag <- FALSE

configs <- list(results.dir = results.dir,
                stopping.lag = 2,
                keep = kfile,
                nCores = ncores,
                mem = 12e3 * ncores,
                use.glmnetPlus = glmnetPlus_flag,
                #glmnet.thresh = 1e-6,
                verbose = TRUE,
                save = TRUE)

### Set up train-validation split ---
phe <- readPheMaster(
  phenotype.file = pfile, psam.ids = ids,
  family = NULL, covariates = c(covars, "pop"),
  phenotype = pheno,
  status = NULL, split.col = NULL, configs = configs
)
ntrain <- nrow(phe)

if (length(unique(phe$sex)) == 1) {
  covars <- c("age", paste0("PC", 1:10))
}

train_ids <- phe %>%
  filter(pop %in% c("EUR", min_ancestry)) %>%
  group_by(pop) %>%
  sample_n(round(0.8 * n())) %>%
  pull(FID)

new_pfile <- file.path(results.dir, "pheno.tsv")
phe %>%
  dplyr::mutate(
    split = if_else(FID %in% train_ids, "train", "val")
  ) %>%
  data.table::fwrite(new_pfile, sep = "\t")

if (length(unique(phe$Sex)) == 1)
  covars <- c("Age", paste0("PC", 1:10))

rm(phe)
gc()

### Fit lasso ---

if (!file.exists(outfile)) {
  intermediate_files <- list.files(file.path(results.dir, "results"),
                                   pattern = "output*") %>%
    tools::file_path_sans_ext()
  prev_iters <- as.integer(gsub("output_iter_", "", intermediate_files))
  if (length(prev_iters) > 0) {
    configs$prev_iter <- max(prev_iters)
  }

  system.time(
    mod <- snpnet(gfile,
                  new_pfile,
                  pheno,
                  covariates = covars,
                  split.col = "split",
                  configs = configs)
  )
  saveRDS(mod, outfile)
  snpnet:::cleanUpIntermediateFiles(mod$configs)
}
