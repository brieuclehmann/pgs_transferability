# Script to evaluate predictive accuracy of PGS

### Set up parameters ---

.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(readr)
library(dplyr)
set.seed(1)

scores_file <- snakemake@output[[1]]
dir.create(dirname(scores_file), recursive = TRUE, showWarnings = FALSE)
this_pheno <- pheno <- snakemake@wildcards[["pheno"]]
this_min_ancestry <- min_ancestry <- snakemake@wildcards[["min_ancestry"]]
variants <- snakemake@wildcards[["v"]]
pfile <- "data/all_vars.tsv"
nfolds <- 5

compute_r2 <- function(y, pred) {
    ss_tot <- sum((y - mean(y))^2)
    ss_res <- sum((y - pred)^2)

    1 - (ss_res / ss_tot)
}

outdir <- file.path(
    "output", "ukb",
    paste0("v~", variants),
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
    filter(!is.na(y))

fold_df <- pheno_df %>%
    select(eid, ID, pop) %>%
    mutate(fold = NA_integer_)
for (this_f in seq(nfolds)) {
    eur_kfile <- paste0(
        "data/train_ids/pheno~", pheno,
        "/min_ancestry~", min_ancestry,
        "/prop_min~0.0/fold~", this_f, ".txt"
    )
    maj_train_df <- read_delim(eur_kfile, delim = " ", col_names = c("eid", "fid"))
    fold_df <- fold_df %>%
        mutate(fold = if_else(!(eid %in% maj_train_df$eid) & pop == "EUR", this_f, fold))
        

    if (pheno %in% c("A50", "A30040", "NC1111", "A30070")) {
        min_ancestries <- c("AFR", "CSA", "EAS", "AMR", "MID")
    } else {
        min_ancestries <- min_ancestry
    }

    for (this_ancestry in min_ancestries) {
        min_kfile <- paste0(
            "data/train_ids/pheno~", pheno,
            "/min_ancestry~", this_ancestry,
            "/prop_min~1.0/fold~", this_f, ".txt"
        )
        min_train_df <- read_delim(min_kfile, delim = " ", col_names = c("eid", "fid"))
        fold_df <- fold_df %>%
            mutate(fold = if_else(!(eid %in% min_train_df$eid) & pop == this_ancestry, this_f, fold))
    }
}

pheno_df <- pheno_df %>%
    inner_join(fold_df, by = c("eid", "ID", "pop"))

if (pheno == "N81") {
    pheno_df <- pheno_df %>%
        filter(sex == 0)
}

pred_files <- list.files(outdir, pattern = "pred.tsv", recursive = TRUE)


if (pheno %in% c("NC1226", "NC1111", "I48", "K57", "N81")) {
    fam <- "binomial"
} else {
    fam <- "gaussian"
}

out_df <- tibble(pop = character())
out_sex_df <- tibble(pop = character(), sex = double())
for (pred_file in pred_files) {
    prop_min <- gsub(".*prop_min~(.*?)/.*", "\\1", pred_file)
    f <- gsub(".*fold~(.*?)/.*", "\\1", pred_file)
    pow <- gsub(".*pow~(.*?)/.*", "\\1", pred_file)

    pred_df <- read_tsv(file.path(outdir, pred_file)) %>%
        inner_join(pheno_df, by = "ID") %>%
        filter(fold == f)

        out_df <- pred_df %>%
            group_by(pop) %>%
            mutate(
                resid_covar = residuals(lm(y ~ age + sex *
                    (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10))),
                fit_covar = fitted.values(lm(y ~ age + sex *
                    (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10))),
                resid_pgs = residuals(lm(pred ~ age + sex *
                    (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10))),
                resid_full = residuals(lm(y ~ age + pred + sex *
                    (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)))
            ) %>%
            summarise(
                r2 = compute_r2(y, pred),
                covar_r2 = 1 - sum(resid_covar^2) / sum((y - mean(y))^2),
                covar = summary(glm(y ~ age + sex *
                    (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10), family = fam))$deviance,
                full_r2 = 1 - sum(resid_full^2) / sum((y - mean(y))^2),
                full = summary(glm(y ~ age + pred + sex *
                    (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10), family = fam))$deviance,
                inc_r2 = full_r2 - covar_r2,
                partial_r2 = 1 - (sum(resid_full^2) / sum(resid_covar^2)),
                pseudo_r2 = 1 - exp((full - covar) / n()),
                partial_cor = cor(resid_covar, resid_pgs),
                bias = mean(pred) - mean(y),
                VE = 1 - (mean((resid_pgs - resid_covar)^2) / var(resid_covar)),
                cor = cor(pred, y),
                implied_var = var(pred),
                implied_var_cov = var(fit_covar)
            ) %>%
            mutate(prop_min = prop_min, pow = pow, fold = f) %>%
            bind_rows(out_df)

}

out_df$pheno <- pheno
out_df$min_ancestry <- min_ancestry
write_tsv(out_df, scores_file)
