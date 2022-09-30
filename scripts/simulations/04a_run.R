# LOAD PACKAGES AND DATA ------------------------------------------------------

# old_path <- .libPaths()
# new_path <- c("/well/holmes/users/rxa753/R/3.6/skylake", old_path)
# .libPaths(new_path)
.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(tibble)
library(dplyr)
library(tidyr)
library(mvtnorm)
library(pROC)
library(tidyselect)
source("scripts/simulations/utils.R")

tskit <- reticulate::import("tskit")
ts <- tskit$load("data/ooa.trees")
gfile <- "data/ooa"

maf_df <- readr::read_delim("data/maf.txt",
    delim = " ",
    col_names = c("YRI", "CEU", "CHB")
) %>%
    mutate_all(function(x) pmin(x, 1 - x))
maf_df$SNP <- 1:nrow(maf_df)

# SPECIFY SIMULATION HYPERPARAMETERS ------------------------------------------
maf_threshold <- 1e-2
n_test <- 2e4L
prop_wt <- 0.1
n_sim <- 10
n_eur <- 18000

iter <- as.double(snakemake@wildcards[["iter"]])
beta_cor <- as.double(snakemake@wildcards[["beta_cor"]])
prop_afr <- as.double(snakemake@wildcards[["prop_afr"]])
pow <- as.double(snakemake@wildcards[["pow"]])
r <- as.double(snakemake@wildcards[["rep"]])
h2 <- as.double(snakemake@wildcards[["h2"]])
n_causal <- as.double(snakemake@wildcards[["n_causal"]])

n_cores <- as.integer(Sys.getenv("NSLOTS"))
mem <- 12e3 * n_cores

n_afr <- (prop_afr / (1 - prop_afr)) * n_eur

n_pop <- 3
pop <- c("YRI", "CEU", "CHB")
pow_range <- seq(0, 1, 0.1)

include <- with(maf_df, pmax(YRI, CEU, CHB) >= maf_threshold)
ind_snp <- maf_df$SNP[include]
nSNP <- length(ind_snp)

beta_cov <- matrix(beta_cor, n_pop, n_pop)
diag(beta_cov) <- 1

# GENERATE DATA ----------------------------------------------------------------
set.seed(iter)

ind_sub <- get_subjects(ts)
names(ind_sub) <- pop

true_beta <- simul_effects(nSNP, n_causal, beta_cov)
colnames(true_beta) <- pop
true_beta_df <- as_tibble(true_beta)

snp_risk <- mapply(function(x, y) get_snp_risk(ts, y, ind_snp, x),
    ind_sub, true_beta_df,
    SIMPLIFY = FALSE
)
snp_risk_var <- sapply(snp_risk, var)
env_var <- mean(snp_risk_var) / h2 - mean(snp_risk_var)
pheno <- lapply(snp_risk, function(x) x + rnorm(length(x), sd = sqrt(env_var)))

sam_df <- readr::read_tsv("data/ooa.psam",
    col_names = c("FID", "IID"),
    skip = 1
)
pheno_df <- tibble(
    IID = sam_df$IID,
    FID = sam_df$IID,
    pheno = unlist(pheno),
    pop = rep(pop, each = length(ind_sub[[1]])),
    split = NA_character_,
    value = 1
) %>%
    pivot_wider(names_from = pop, values_from = value, values_fill = 0)

# RUN SIMUL --------------------------------------------------------------------
pred_df <- tibble(
    IID = character(),
    pop = character(),
    pheno = double(),
    pred = double()
)

if (pow == Inf) {
    n_train <- round(c(n_afr, 0, 0)) # n_eur <- 0
    weights <- rep(1, n_afr)
    pow <- "inf"
} else {
    n_train <- round(c(n_afr, n_eur, 0))
    IPW <- rep(1 / c(prop_wt, 1 - prop_wt, 0), n_train)

    # normalised weights
    weights <- sum(n_train) * (IPW^pow) / sum(IPW^pow)
}

set.seed(r)

ind_train <- mapply(
    function(x, y) sort(sample(x, y)),
    ind_sub, n_train
)
ind_val <- mapply(
    function(x, y) sort(sample(x, y)),
    ind_train, round(n_train * 0.2)
)
ind_test <- mapply(
    function(x, y) sort(sample(setdiff(x, y), n_test)),
    ind_sub, ind_train,
    SIMPLIFY = FALSE
)

pheno_df$split <- NA_character_
pheno_df$split[unlist(ind_train)] <- "train"
pheno_df$split[unlist(ind_val)] <- "val"

pfile <- tempfile()
readr::write_tsv(filter(pheno_df, split %in% c("train", "val")), pfile)
all_weights <- double(nrow(pheno_df))
all_weights[unlist(ind_train)] <- weights

start_time <- proc.time()
  system.time(
      fit <- snpnet::snpnet(
          gfile, pfile, "pheno",
          covariates = "YRI",
          weights = weights,
          split.col = "split",
          mem = mem
      )
  )
diff_time <- proc.time() - start_time

pheno_df$split[unlist(ind_test)] <- "test"
pfile <- tempfile()
readr::write_tsv(filter(pheno_df, split == "test"), pfile)

lambda_ind <- which.max(fit$metric.val)
pred <- snpnet::predict_snpnet(
    fit,
    new_genotype_file = gfile, new_phenotype_file = pfile,
    phenotype = "pheno", split_col = "split", split_name = "test",
    covariate_names = "YRI", idx = lambda_ind
)

this_pred_df <- pheno_df %>%
    filter(split == "test") %>%
    mutate(
        IID_FID = rownames(pred$prediction[[1]]),
        pred = pred$prediction[[1]][, 1],
        check = paste(IID, FID, sep = "_") == IID_FID
    ) %>%
    pivot_longer(
        cols = c(YRI, CEU, CHB), names_to = "pop", values_to = "keep"
    ) %>%
    filter(keep == 1) %>%
    select(IID, pop, pheno, pred, check)

if (!all(this_pred_df$check)) stop("Error combining predictions")

pred_df <- bind_rows(pred_df, select(this_pred_df, -check))

# SAVE OUTPUT ------------------------------------------------------------------
out_dir <- file.path(
    "output", "simulations",
    paste0("n_causal~", format(n_causal, nsmall = 1)),
    paste0("h2~", format(h2, nsmall = 1)),
    paste0("beta_cor~", format(beta_cor, nsmall = 1)),
    paste0("prop_afr~", format(prop_afr, nsmall = 1)),
    paste0("pow~", format(pow, nsmall = 1)),
    paste0("iter~", format(iter, nsmall = 1))
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pred_file <- file.path(out_dir, paste0("rep~", format(r, nsmall = 1), ".csv"))
readr::write_csv(pred_df, pred_file)

time_file <- file.path(out_dir, paste0("time~", format(r, nsmall = 1), ".txt"))
write(diff_time, time_file)