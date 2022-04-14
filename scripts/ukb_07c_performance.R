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

match_dirs <- paste0(outdir, "/fold~", 1:5)
pred_files <- do.call(c, lapply(match_dirs, function(x) list.files(x, pattern = "pred.tsv$", recursive = TRUE, full.names = TRUE)))

covar_formula <- y ~ age + sex *
    (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
full_formula <- y ~ age + pred + sex *
    (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)

fit_lm_on_bootstrap <- function(split) {
    lm(y ~ age + sex *
        (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10), split)
}

if (pheno %in% c("NC1226", "NC1111", "I48", "K57", "N81")) {
    fam <- "binomial"
} else {
    fam <- "gaussian"
}

out_df <- tibble(pop = character())
boot_df <- tibble(pop = character())
for (pred_file in pred_files) {
    type <- gsub(".*type~(.*?)/.*", "\\1", pred_file)
    f <- gsub(".*fold~(.*?)/.*", "\\1", pred_file)
    frac <- gsub(".*frac~(.*?)/.*", "\\1", pred_file)

    pred_df <- read_tsv(pred_file) %>%
        inner_join(pheno_df, by = "ID") %>%
        filter(fold == f)

        out_df <- pred_df %>%
            group_by(pop) %>%
            mutate(
                resid_covar = residuals(lm(y ~ age + sex *
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
                cor = cor(pred, y)
            ) %>%
            mutate(type = type, frac = frac, fold = f) %>%
            bind_rows(out_df)
   
}

out_df$pheno <- pheno
out_df$min_ancestry <- min_ancestry
#scores_file <- file.path(outdir, "scores.tsv")
write_tsv(out_df, scores_file)