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
ts <- tskit$load("data/ooa_demo.trees")
gfile <- "data/ooa_demo"


# SPECIFY HYPERPARAMETERS ------------------------------------------
n_test <- 100
prop_wt <- 0.1
n_sim <- 1
n_eur <- 1000

iter <- 1
beta_cor <- 0.8
prop_afr <- 0.2
pow <- 0.2
r <- 1
h2 <- 0.4
n_causal <- 10

n_cores <- 1
mem <- 12e3 * n_cores

n_afr <- (prop_afr / (1 - prop_afr)) * n_eur

n_pop <- 2
pop <- c("YRI", "CEU")
pow_range <- seq(0, 1, 0.1)

ind_snp <- 1:ts$num_sites
nSNP <- ts$num_sites

beta_cov <- matrix(beta_cor, n_pop, n_pop)
diag(beta_cov) <- 1

# GENERATE DATA ----------------------------------------------------------------

ind_sub <- list()
ind_sub$YRI <- 1:2000
ind_sub$CEU <- 2001:4000

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

sam_df <- readr::read_tsv("data/ooa_demo.psam",
    col_names = c("FID", "IID"),
    skip = 1
)
pheno_df <- tibble(
    IID = sam_df$IID,
    FID = 0,
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

n_train <- round(c(n_afr, n_eur))
IPW <- rep(1 / c(prop_wt, 1 - prop_wt), n_train)

# normalised weights
weights <- sum(n_train) * (IPW^pow) / sum(IPW^pow)

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
        check = paste(FID, IID, sep = "_") == IID_FID
    ) %>%
    pivot_longer(
        cols = c(YRI, CEU), names_to = "pop", values_to = "keep"
    ) %>%
    filter(keep == 1) %>%
    select(IID, pop, pheno, pred, check)

if (!all(this_pred_df$check)) stop("Error combining predictions")

pred_df <- bind_rows(pred_df, select(this_pred_df, -check))

# SAVE OUTPUT ------------------------------------------------------------------
out_dir <- file.path(
    "output", "simulations"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pred_file <- file.path(out_dir, paste0("demo_pred.csv"))
readr::write_csv(pred_df, pred_file)

time_file <- file.path(out_dir, paste0("demo_time.txt"))
write(diff_time, time_file)