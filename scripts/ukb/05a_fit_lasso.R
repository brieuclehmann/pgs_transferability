# This script fits the LASSO using a validation set to select lambda.

### Set up parameters ---

library(snpnet)
library(readr)
library(dplyr)
set.seed(1)

outfile <- snakemake@output[[1]]
outdir <- dirname(outfile)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
gfile <- tools::file_path_sans_ext(snakemake@input[[1]])
kfile <- snakemake@input[[2]]
results.dir <- tools::file_path_sans_ext(sub("output/ukb", "temp", outfile))
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)

pheno <- snakemake@wildcards[["pheno"]]
prop_min <- snakemake@wildcards[["prop_min"]]
min_ancestry <- snakemake@wildcards[["min_ancestry"]]
pow <- as.double(snakemake@wildcards[["pow"]])
ncores <- as.integer(Sys.getenv("NSLOTS"))

pfile  <- "data/all_vars.tsv"
covars <- c("sex", "age", paste0("PC", 1:10), paste0("PC", 1:10, "_sex"))

keep_ids <- read.table(kfile)[, 1]
ids <- readIDsFromPsam(paste0(gfile, ".psam"))
id_df <- data.frame(ID = ids, stringsAsFactors = F) %>%
  mutate(sort_order = 1:n())

glmnetPlus_flag <- TRUE
if (!(prop_min %in% c("0.0", "-1.0")) & pheno != "NC1111") {
  glmnetPlus_flag <- FALSE
}
glmnet.thresh <- 1e-07
if (pheno == "NC1111" & pow >= 1.2) {
  glmnet.thresh <- 1e-06
}

configs <- list(results.dir = results.dir,
                stopping.lag = 2,
                keep = kfile,
                nCores = ncores,
                mem = 12e3 * ncores,
                use.glmnetPlus = glmnetPlus_flag,
                glmnet.thresh = glmnet.thresh,
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

weights <- phe %>%
  group_by(pop) %>%
  transmute(
    ID = paste(FID, FID, sep = "_"),
    weight = ntrain / n()
  ) %>%
  left_join(id_df, by = "ID") %>%
    ungroup() %>%
    mutate(weight = ntrain * (weight^pow / sum(weight^pow))) %>%
    arrange(sort_order) %>%
    pull(weight)

if (prop_min %in% c("0.0", "1.0") && !all.equal(weights, rep(1, ntrain))) {
  stop("Check weights")
}

new_pfile <- file.path(results.dir, "pheno.tsv")
phe %>%
  dplyr::mutate(
    split = if_else(FID %in% train_ids, "train", "val"),
    weights = weights
  ) %>%
  data.table::fwrite(new_pfile, sep = "\t")

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

  start_time <- proc.time()
  system.time(
    mod <- snpnet(gfile,
                  new_pfile,
                  pheno,
                  weights = weights,
                  covariates = covars,
                  split.col = "split",
                  configs = configs)
  )
  diff_time <- proc.time() - start_time
  mod$time <- diff_time

  saveRDS(mod, outfile)
  snpnet:::cleanUpIntermediateFiles(mod$configs)
}
