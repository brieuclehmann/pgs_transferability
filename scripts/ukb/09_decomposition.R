# Script to generate test set predictions

### Set up parameters ---

.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(snpnet)
library(purrr)
library(foreach)
library(readr)
library(dplyr)
set.seed(1)

decomp_file <- snakemake@output[[1]]
dir.create(dirname(decomp_file), recursive = TRUE, showWarnings = FALSE)
variants <- snakemake@wildcards[["v"]]
pheno <- snakemake@wildcards[["pheno"]]
prop_min <- snakemake@wildcards[["prop_min"]]
min_ancestry <- snakemake@wildcards[["min_ancestry"]]
f <- as.integer(snakemake@wildcards[["fold"]])
pow <- as.double(snakemake@wildcards[["pow"]])
nfolds <- 5

pfile <- "data/all_vars.tsv"
### Combine chromosome predictions ###
outdir <- file.path(
  "output", "ukb",
  paste0("v~", variants),
  paste0("pheno~", pheno),
  paste0("min_ancestry~", min_ancestry),
  paste0("prop_min~", prop_min),
  paste0("fold~", f),
  paste0("pow~", format(pow, nsmall = 1))
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

beta_file <- file.path(outdir, "beta.tsv")

chroms <- c(1:22, "X")
maf_cuts <- c(0, 0.01, 0.05, 0.2, 0.5)

pred_df <- tibble()
for (maf_pop in c("min", "maj")) {

    if (maf_pop == "min") {
        beta_df <- read_tsv(beta_file, col_types="ccddd") %>%
            mutate(maf_grp = cut(maf_min, maf_cuts))
    } else {
        beta_df <- read_tsv(beta_file, col_types="ccddd") %>%
            mutate(maf_grp = cut(maf_maj, maf_cuts))
    }

    for (grp in levels(beta_df$maf_grp)) {

        out_df <- foreach(this_chrom = chroms, .combine = rbind) %do% {
            gfile <- file.path(
                "data", "genotypes",
                paste0("v~", variants),
                paste0("min_ancestry~", min_ancestry),
                paste0("chrom~", this_chrom)
            )

            snp_null <- beta_df %>%
                filter(chrom == this_chrom & maf_grp != grp) %>%
                pull(varname)

            outfile <- file.path(outdir, paste0("chrom~", this_chrom, ".RDS"))

            mod <- readRDS(outfile)
            lambda_ind <- which.max(mod$metric.val)
            covars <- mod$configs$covariates
            mod$beta[[lambda_ind]][covars] <- 0
            mod$beta[[lambda_ind]][snp_null] <- 0
            mod$a0[[lambda_ind]] <- 0

            pred <- predict_snpnet(mod,
                new_genotype_file = gfile,
                new_phenotype_file = pfile,
                phenotype = pheno,
                covariate_names = covars,
                idx = lambda_ind
            )

            tibble(
                ID = rownames(pred$prediction$train),
                pred = pred$prediction$train[, 1],
                chrom = this_chrom
            )
        }

        this_pred_df <- out_df %>%
            group_by(ID) %>%
            summarise(pred = sum(pred)) %>%
            mutate(maf_pop = maf_pop, maf_grp = grp)

        pred_df <- bind_rows(this_pred_df, pred_df)
    }
}

### Save predictions
pred_file <- file.path(outdir, "pred_split.tsv")
write_tsv(pred_df, pred_file)


### Compute performance

if (pheno %in% c("NC1226", "NC1111", "I48", "K57", "N81")) {
    fam <- "binomial"
} else {
    fam <- "gaussian"
}

compute_r2 <- function(y, pred) {
    ss_tot <- sum((y - mean(y))^2)
    ss_res <- sum((y - pred)^2)

    1 - (ss_res / ss_tot)
}

out_df <- pred_df %>%
    inner_join(pheno_df, by = ("ID")) %>%
    filter(fold == f) %>%
    group_by(pop, maf_pop, maf_grp) %>%
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
            (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
            family = fam))$deviance,
        full_r2 = 1 - sum(resid_full^2) / sum((y - mean(y))^2),
        full = summary(glm(y ~ age + pred + sex *
            (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
            family = fam))$deviance,
        inc_r2 = full_r2 - covar_r2,
        partial_r2 = 1 - (sum(resid_full^2) / sum(resid_covar^2)),
        pseudo_r2 = 1 - exp((full - covar) / n()),
        partial_cor = cor(resid_covar, resid_pgs),
        bias = mean(pred) - mean(y),
        VE = 1 - (mean((resid_pgs - resid_covar)^2) / var(resid_covar)),
        cor = cor(pred, y)
    )

write_tsv(out_df, decomp_file)