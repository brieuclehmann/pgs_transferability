.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(Matching)
library(readr)
library(dplyr)
set.seed(1) # for reproducibility

# Set parameters
maj_ancestry <- "EUR"
prop_min <- 0.1
nfolds <- 5
nfracs <- 5

out_file <- snakemake@output[[1]]
type <- as.character(snakemake@wildcards[["type"]])
pheno <- as.character(snakemake@wildcards[["pheno"]])
this_frac <- as.double(snakemake@wildcards[["frac"]])
min_ancestry <- as.character(snakemake@wildcards[["min_ancestry"]])
f <- as.integer(snakemake@wildcards[["fold"]])

# frac_type <- ifelse(
#     frac_type == "frac_maj",
#     paste0("frac_", maj_ancestry), 
#     paste0("frac_", min_ancestry)
# )
print(min_ancestry)

covars <- c(
    "pop", "age", "sex", "age_sex", "age2", "age2_sex",
    paste0("PC", 1:10), paste0("PC", 1:10, "_sex")
)
pheno_df <- read_tsv("data/all_vars.tsv") %>%
    rename(y = all_of(pheno)) %>%
    dplyr::select(eid, y, all_of(covars)) %>%
    filter(!is.na(y))

unique_sex <- unique(pheno_df$sex[pheno_df$y == 1])
if (length(unique_sex) == 1) {
    pheno_df <- pheno_df %>%
        filter(sex == unique_sex)
}

fold_df <- pheno_df %>%
    group_by(pop) %>%
    mutate(
        fold = sample(cut(seq_len(n()), breaks = nfolds, labels = FALSE))
    ) %>%
    group_by(pop, fold) %>%
    mutate(
        frac = sample(cut(seq_len(n()), breaks = nfracs, labels = FALSE))
    ) %>%
    ungroup()

# Construct training sets
train_ratio <- round((1 - prop_min) / prop_min)

frac_filter <- ifelse(
    type == "maj",
    maj_ancestry,
    min_ancestry
)
match_df <- fold_df %>%
    filter(fold != f & pop %in% c(min_ancestry, maj_ancestry)) %>%
    mutate(
        Tr = (pop == min_ancestry),
        frac = if_else(pop != frac_filter, 0L, frac)
    ) %>%
    dplyr::select(eid, frac, Tr, sex, age)

match_out <- match_df %>%
    with(Match(
            Tr = Tr, X = cbind(sex, age), exact = c(TRUE, FALSE),
            M = train_ratio, Z = eid, ties = FALSE, replace = FALSE
    ))

idx <- unique(c(match_out$index.treated, match_out$index.control))

frac_ind <- round(this_frac * nfracs)
frac_df <- match_df %>%
    filter(eid %in% match_df$eid[idx] & frac <= frac_ind)


# Save output
out_dir <- dirname(out_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write(rep(frac_df$eid, each = 2), out_file, ncol = 2)
