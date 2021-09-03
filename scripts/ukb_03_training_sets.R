.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(Matching)
library(readr)
library(dplyr)
set.seed(1) # for reproducibility

# Set parameters
maj_ancestry <- "EUR"
nfolds <- 5

pheno <- as.character(snakemake@wildcards[["pheno"]])
prop_min <- as.double(snakemake@wildcards[["prop_min"]])
min_ancestry <- as.character(snakemake@wildcards[["min_ancestry"]])
f <- as.integer(snakemake@wildcards[["fold"]])

print(min_ancestry)

covars <- c(
    "pop", "age", "sex", "age_sex", "age2", "age2_sex",
    paste0("PC", 1:10), paste0("PC", 1:10, "_sex")
)
pheno_df <- read_tsv("data/all_vars.tsv") %>%
    rename(y = all_of(pheno)) %>%
    select(eid, y, all_of(covars)) %>%
    filter(!is.na(y))

out_dir <- file.path(
    "data", "train_ids",
    paste0("pheno~", pheno),
    paste0("min_ancestry~", min_ancestry),
    paste0("prop_min~", format(prop_min, nsmall = 1))
)
out_file <- file.path(out_dir, paste0("fold~", f, ".txt"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fold_df <- pheno_df %>%
    group_by(pop) %>%
    mutate(
        fold = sample(cut(seq_len(n()), breaks = nfolds, labels = FALSE))
        )

# Construct training sets
if (prop_min == 0) {

    train_id <- fold_df %>%
        filter(pop == maj_ancestry & fold != f) %>%
        pull(eid)

} else if (prop_min == 1) {

    train_id <- fold_df %>%
        filter(pop == min_ancestry & fold != f) %>%
        pull(eid)

} else {
    # Train-test split
    train_ratio <- round((1 - prop_min) / prop_min)

    match_df <- fold_df %>%
        filter(fold != f & pop %in% c(min_ancestry, maj_ancestry)) %>%
        mutate(Tr = (pop == min_ancestry)) %>%
        dplyr::select(eid, Tr, sex, age)

    match_out <- match_df %>%
        with(Match(
            Tr = Tr, X = cbind(sex, age), exact = c(TRUE, FALSE),
            M = train_ratio, Z = eid, ties = FALSE, replace = FALSE
        ))

    idx <- unique(c(match_out$index.treated, match_out$index.control))
    train_id <- match_df$eid[idx]
}

write(rep(train_id, each = 2), out_file, ncol = 2)
