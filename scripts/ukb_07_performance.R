# Script to evaluate predictive accuracy of PGS

### Set up parameters ---

.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(readr)
library(dplyr)
set.seed(1)

this_pheno <- pheno <- snakemake@wildcards[["pheno"]]
this_min_ancestry <- min_ancestry <- snakemake@wildcards[["min_ancestry"]]
pfile <- "data/all_vars.tsv"
nfolds <- 5

compute_r2 <- function(y, pred) {
    ss_tot <- sum((y - mean(y))^2)
    ss_res <- sum((y - pred)^2)

    1 - (ss_res / ss_tot)
}

outdir <- file.path(
    "output", "ukb",
    paste0("pheno~", pheno),
    paste0("min_ancestry~", min_ancestry)
)

covars <- c(
    "pop", "age", "sex", "age_sex", "age2", "age2_sex",
    paste0("PC", 1:10), paste0("PC", 1:10, "_sex")
)
pheno_df <- read_tsv("data/all_vars.tsv") %>%
    rename(y = all_of(pheno)) %>%
    select(eid, y, all_of(covars)) %>%
    mutate(ID = paste0(eid, "_", eid)) %>%
    filter(!is.na(y)) %>%
    group_by(pop) %>%
        mutate(
            fold = sample(cut(seq_len(n()), breaks = nfolds, labels = FALSE))
        )

pred_files <- list.files(outdir, pattern = "pred.tsv", recursive = TRUE)

covar_formula <- y ~ age + sex *
    (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
full_formula <- y ~ age + pred + sex *
    (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)

fit_lm_on_bootstrap <- function(split) {
    lm(y ~ age + sex *
        (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10), split)
}

out_df <- tibble(pop = character())
boot_df <- tibble(pop = character())
for (pred_file in pred_files) {
    prop_min <- gsub(".*prop_min~(.*?)/.*", "\\1", pred_file)
    f <- gsub(".*fold~(.*?)/.*", "\\1", pred_file)
    pow <- gsub(".*pow~(.*?)/.*", "\\1", pred_file)

    pred_df <- read_tsv(file.path(outdir, pred_file)) %>%
        inner_join(pheno_df, by = "ID") %>%
        filter(fold == f)

    boot_df <- pred_df %>%
        group_by(pop) %>%
        group_modify(
            ~ bootstraps(., times = 100, apparent = TRUE) %>%
                mutate(
                    model = map(splits, ~ lm(covar_formula, data = .)),
                    model2 = map(splits, ~ lm(full_formula, data = .)),
                    full_glance = map(model2, glance),
                    covar_glance = map(model, glance)
                ) %>%
                unnest(covar_glance) %>%
                select(id, covar = r.squared, full_glance) %>%
                unnest(full_glance) %>%
                select(id, covar, full = r.squared)
        ) %>%
        mutate(
            bootstrap = id != "Apparent",
            partial = full - covar) %>%
        group_by(pop, bootstrap) %>%
            summarise(
                mean = mean(partial),
                upper = quantile(partial, 0.95),
                lower = quantile(partial, 0.05)
            ) %>%
        mutate(prop_min = prop_min, pow = pow, fold = f) %>%
        bind_rows(boot_df)

    out_df <- pred_df %>%
        group_by(pop) %>%
        mutate(
            resid_full = residuals(lm(y ~ age + pred + sex *
                (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10))),
            resid_covar = residuals(lm(y ~ age + sex *
                (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10))),
            diff_resid = resid_covar^2 - resid_full^2 / var(y)
        ) %>%
        summarise(
            partial = mean(diff_resid) / var(y),
            var_resid = var(diff_resid / var(y)),
            partial_upper = quantile(diff_resid, 0.95) / var(y),
            partial_lower = quantile(diff_resid, 0.05) / var(y)
            )
    
   out_df <- pred_df %>%
       group_by(pop) %>%
       summarise(
           r2 = compute_r2(y, pred),
           covar = summary(lm(y ~ age + sex *
               (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)))$r.squared,
           full = summary(lm(y ~ age + pred + sex *
               (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)))$r.squared,
           partial = full - covar,
           bias = mean(pred) - mean(y),
           mse = mean((pred - y)^2),
           cor = cor(pred, y)
       ) %>%
       mutate(prop_min = prop_min, pow = pow, fold = f) %>%
       bind_rows(out_df)
}

out_df$pheno <- pheno
out_df$min_ancestry <- min_ancestry
scores_file <- file.path(outdir, "scores.tsv")
write_tsv(out_df, scores_file)

boot_df$pheno <- pheno
boot_df$min_ancestry <- min_ancestry
boot_scores_file <- file.path(outdir, "boot_scores.tsv")
write_tsv(boot_df, boot_scores_file)