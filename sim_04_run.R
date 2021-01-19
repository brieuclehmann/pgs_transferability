# LOAD PACKAGES AND DATA ------------------------------------------------------

library(tibble)
library(dplyr)
library(tidyr)
library(mvtnorm)
library(pROC)
library(tidyselect)
source("sim_utils.R")

tskit <- reticulate::import("tskit")
ts <- tskit$load("data/ooa_chr20_10.trees")

maf_df   <- readr::read_csv("data/maf_chr20_10.csv", col_types=readr::cols())

# SPECIFY SIMULATION HYPERPARAMETERS ------------------------------------------
maf_threshold <- 1e-3
n_test <- 2e3L
n_causal <- 100
prop_wt <- 0.1
h2 <- 0.5
n_sim <- 20
n_eur <- 1800

iter <- 1 # vary from 1 to 50
beta_cor <- 0.5 # takes values in 0.5 0.6 0.7 0.8 0.9
prop_afr <- 0.1 # takes values in 0.1 0.2 0.3 0.4 0.5

n_afr <- (prop_afr / (1 - prop_afr)) * n_eur

n_pop <- 3
pop <- c("AFR", "EUR", "EAS")
pow_range <- seq(0, 1, 0.1)

include <- with(maf_df, pmin(mafEUR, 1 - mafEUR, 
                             mafAFR, 1 - mafAFR, 
                             mafEAS, 1 - mafEAS) > maf_threshold)
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
                   ind_sub, true_beta_df, SIMPLIFY = FALSE)
snp_risk_var <- sapply(snp_risk, var)
env_var <- mean(snp_risk_var)/h2 - mean(snp_risk_var)
pheno <- lapply(snp_risk, function(x) x + rnorm(length(x), sd = sqrt(env_var)))

ind_test <- mapply(function(x, y) sort(sample(x, y)), 
                   ind_sub, n_test, SIMPLIFY = FALSE)
Gtest <- lapply(ind_test, function(x) get_genotype(ts, x, ind_snp))

# Add subpopulation-intercept to genotype matrix
is_afr <- as.integer(pop == "AFR")
Gtest <- mapply(function(x, y) cbind(rep(y, n_test), x), Gtest, is_afr)

snp_risk_test  <- mapply(function(x, y, z) x[y %in% z], 
                         snp_risk, ind_sub, ind_test, SIMPLIFY = FALSE)
pheno_test  <- mapply(function(x, y, z) x[y %in% z], 
                      pheno, ind_sub, ind_test, SIMPLIFY = FALSE)

# RUN SIMUL --------------------------------------------------------------------
pred_df <- tibble(r = integer(),
                  pow = double(),
                  pop = character(),
                  n = character(),
                  pred = double())

n_train <- round(c(n_afr, n_eur, 0))
IPW <- rep(1 / c(prop_wt, 1 - prop_wt, 0), n_train)

for (pow in pow_range) {
  
  weights <- n_train * (IPW^pow) / sum(IPW^pow) # normalised weights

  for (r in seq_len(n_sim)) {
    set.seed(r)

    ind_train <- mapply(function(x, y, z) sort(sample(setdiff(x, y), z)),
                        ind_sub, ind_test, n_train, SIMPLIFY = FALSE)
    Gtrain <- get_genotype(ts, unlist(ind_train), ind_snp)
    Gtrain <- cbind(rep(is_afr, n_train), Gtrain)
    
    pheno_train <- mapply(function(x, y, z) x[y %in% z], 
                          pheno, ind_sub, ind_train)
    pheno_train <- unlist(pheno_train)
    
    fit <- glmnet::cv.glmnet(Gtrain, pheno_train, nfolds = 10, 
                             weights = weights, 
                             penalty.factor = c(0, rep(1, nSNP)),
                             standardize = FALSE)

    pred <- lapply(Gtest, function(x) drop(predict(fit, x, s = "lambda.min")))

    this_pred_df <- as_tibble(pred) %>% 
      pivot_longer(all_of(pop), values_to = "pred", names_to = "pop") %>%
      mutate(r = r, pow = pow) %>%
      rownames_to_column(var = "n")
    
    pred_df <- bind_rows(pred_df, this_pred_df)
  }
}

# SAVE OUTPUT ------------------------------------------------------------------
dir.create("temp", showWarnings = FALSE)
this_sim <- paste0("beta_cor", beta_cor, "-prop_afr", prop_afr)

pred_df$simID <- this_sim
pred_df$iter <- iter
pred_file <- paste0("temp/", this_sim, iter, "_pred.tsv")
readr::write_tsv(pred_df, pred_file)

risk_df <- as_tibble(snp_risk_test) %>%
  pivot_longer(all_of(pop), values_to = "risk", names_to = "pop") %>%
  rownames_to_column(var = "n")
pheno_df <- as_tibble(pheno_test) %>%
  pivot_longer(all_of(pop), values_to = "pheno", names_to = "pop") %>%
  rownames_to_column(var = "n") %>%
  inner_join(risk_df, by = c("n", "pop"))
pheno_file <- paste0("temp/", this_sim, iter, "_pheno.tsv")
readr::write_tsv(pheno_df, pheno_file)
