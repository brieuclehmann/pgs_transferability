.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")

library(tibble)
library(dplyr)
source("scripts/simulations/utils.R")

sim_df <- readr::read_csv("data/sim_params.csv") %>%
    distinct(iter, beta_cor, h2, n_causal)

# Load and preprocess data
tskit <- reticulate::import("tskit")
ts <- tskit$load("data/ooa.trees")

maf_df <- readr::read_delim("data/maf.txt",
    delim = " ",
    col_names = c("YRI", "CEU", "CHB")
) %>%
    mutate_all(function(x) pmin(x, 1 - x))
maf_df$SNP <- 1:nrow(maf_df)

n_pop <- 3
pop <- c("YRI", "CEU", "CHB")
maf_threshold <- 1e-2
n_causal <- 100
h2 <- 0.3

include <- with(maf_df, pmax(YRI, CEU, CHB) >= maf_threshold)
ind_snp <- maf_df$SNP[include]
nSNP <- length(ind_snp)

# Get SNP heritability
h2_df <- tibble(
    iter = integer(),
    beta_cor = double(),
    pop = character(),
    n_causal = n_causal,
    h2 = double(),
    h2_pop = double()
)

for (i in seq(nrow(sim_df))) {
    iter <- sim_df$iter[i]
    beta_cor <- sim_df$beta_cor[i]
    h2 <- sim_df$h2[i]
    n_causal <- sim_df$n_causal[i]

    set.seed(iter)

    beta_cov <- matrix(beta_cor, n_pop, n_pop)
    diag(beta_cov) <- 1

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
    pheno <- lapply(
        snp_risk, function(x) x + rnorm(length(x), sd = sqrt(env_var))
    )

    this_df <- tibble(
        iter = iter,
        beta_cor = beta_cor,
        pop = pop,
        n_causal = n_causal,
        h2 = h2,
        h2_pop = snp_risk_var / (snp_risk_var + env_var),
    )
    h2_df <- bind_rows(h2_df, this_df)

}

readr::write_csv(h2_df, "output/simulations/h2.csv")