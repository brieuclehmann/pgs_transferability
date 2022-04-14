# This script fits the LASSO using a validation set to select lambda.

### Set up parameters ---
.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(snpnet)
library(readr)
library(dplyr)
set.seed(1)


pheno <- snakemake@wildcards[["pheno"]]
prop_min <- snakemake@wildcards[["prop_min"]]
min_ancestry <- snakemake@wildcards[["min_ancestry"]]
chrom <- snakemake@wildcards[["chrom"]]
f <- as.integer(snakemake@wildcards[["fold"]])
pow <- as.double(snakemake@wildcards[["pow"]])
ncores <- as.integer(Sys.getenv("NSLOTS"))

kfile  <- file.path("data", "train_ids",
                    paste0("pheno~", pheno),
                    paste0("min_ancestry~", min_ancestry),
                    paste0("prop_min~", prop_min),
                    paste0("fold~", f, ".txt"))
pfile  <- "data/all_vars.tsv"
gfile  <- file.path("data", "genotypes",
                    paste0("ancestry~", min_ancestry),
                    paste0("chrom~", chrom))
covars <- c("sex", "age", paste0("PC", 1:10), paste0("PC", 1:10, "_sex"))

keep_ids <- read.table(kfile)[, 1]
ids <- readIDsFromPsam(paste0(gfile, ".psam"))
id_df <- data.frame(ID = ids, stringsAsFactors = F) %>%
  mutate(sort_order = 1:n())

results.dir <- file.path(
  "temp",
  paste0("pheno~", pheno),
  paste0("min_ancestry~", min_ancestry),
  paste0("prop_min~", prop_min),
  paste0("fold~", f),
  paste0("pow~", format(pow, nsmall = 1)),
  paste0("chrom~", chrom)
)
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)

configs <- list(results.dir = results.dir,
                stopping.lag = 2,
                keep = kfile,
                nCores = ncores,
                mem = 12e3 * ncores,
                use.glmnetPlus = TRUE,
                verbose = TRUE,
                save = TRUE)

### Set up train-validation split ---
phe <- readPheMaster(
  phenotype.file = pfile, psam.ids = ids,
  family = NULL, covariates = c(covars, "pop"),
  phenotype = pheno,
  status = NULL, split.col = NULL, configs = configs
) %>%
  dplyr::mutate_at(all_of(c(pheno, covars)), scale)
ntrain <- nrow(phe)

train_ids <- phe %>%
  filter(pop %in% c("EUR", min_ancestry)) %>%
  group_by(pop) %>%
  sample_n(round(0.8 * n())) %>%
  pull(FID)

 #ntrain <- length(train_ids)
weights <- phe %>%
  group_by(pop) %>%
  transmute(
    ID = paste(FID, FID, sep = "_"),
    weight = ntrain / n(),
    split = if_else(FID %in% train_ids, "train", "val")
  ) %>%
  left_join(id_df, by = "ID") %>%
    ungroup() %>%
   # group_by(split) %>%
    mutate(weight = ntrain * (weight^pow / sum(weight^pow))) %>%
    arrange(sort_order) %>%
    pull(weight)

if (prop_min %in% c("0.0", "1.0") & any(weights != 1)) {
  stop("Check weights")
}

new_pfile <- file.path(results.dir, "pheno.tsv")
phe %>%
  dplyr::mutate(
    split = if_else(FID %in% train_ids, "train", "val"),
    weights = weights
  ) %>%
  data.table::fwrite(new_pfile, sep = "\t")

if (length(unique(phe$Sex)) == 1)
  covars <- c("Age", paste0("PC", 1:10))


########################
### Initialise LASSO ###
########################
phe_train <- phe[phe$FID %in% train_ids, ]
features_train <- cbind(
  as.matrix(phe_train[, covars, with = F]),
  rnorm(nrow(phe_train))
)
weights_train <- weights[phe$FID %in% train_ids]
response_train <- as.matrix(phe_train[, pheno, with = F])

glmmod <- glmnetPlus::glmnet(
  features_train, response_train, weights = weights_train, standardize = F,
  lambda = 100 * max(abs(response_train)),
  penalty.factor = c(rep(0, length(covars)), length(covars) + 1)
)
glmmod2 <- stats::glm(
  stats::as.formula(paste(pheno, " ~ ", paste(c(1, covars), collapse = " + "))),
  data = phe_train, weights = weights_train
)
residual <- response_train - predict(glmmod, features_train)
rownames(residual) <- phe_train$ID
colnames(residual) <- c("0")

residual2 <- matrix(glmmod2$residual)
rownames(residual2) <- phe_train$ID
colnames(residual2) <- c("0")

configs <- snpnet:::setupConfigs(
  configs, gfile, new_pfile, pheno, covars,
  alpha = 1, nlambda = 100, split.col = NULL, p.factor = NULL,
  status.col = NULL, mem = NULL
)
vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd = paste0(configs[["zstdcat.path"]], " ", paste0(gfile, ".pvar.zst"))), "CHROM" = "#CHROM"), VAR_ID = paste(ID, ALT, sep = "_"))$VAR_ID
stats <- snpnet:::computeStats(gfile, phe_train$ID, configs = configs)
prod.full <- snpnet:::computeProduct(
  residual, gfile, vars, stats, configs, weights_train,
  iter = 0
) / nrow(phe_train)
prod.full2 <- snpnet:::computeProduct(
  residual2, gfile, vars, stats, configs, weights_train,
  iter = 0
) / nrow(phe_train)
score <- abs(prod.full[, 1])
score2 <- abs(prod.full2[, 1])

nobs <- nrow(phe_train)
nvars <- length(vars) - length(stats[["excludeSNP"]]) - length(covars)
lambda.min.ratio <- ifelse(nobs < nvars, 0.01, 1e-04)
configs[["lambda.min.ratio"]] <- lambda.min.ratio
full.lams <- snpnet:::computeLambdas(
  score, configs[["nlambda"]], configs[["lambda.min.ratio"]]
)

rm(phe, stats, vars, prod.full, score, weights_train, phe_train)
gc()

#################
### Fit lasso ###
#################

outdir <- file.path(
  "output", "ukb",
  paste0("pheno~", pheno),
  paste0("min_ancestry~", min_ancestry),
  paste0("prop_min~", prop_min),
  paste0("fold~", f),
  paste0("pow~", format(pow, nsmall = 1))
)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, paste0("chrom~", chrom, ".RDS"))

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
      weights = weights,
      covariates = covars,
      split.col = "split",
      configs = configs,
      full.lams = full.lams)
    )
  saveRDS(mod, outfile)
  snpnet:::cleanUpIntermediateFiles(mod$configs)
}
