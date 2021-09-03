### This script constructs the single-ancestry training sets and the dual-ancestry training set with 10% Black or Black British and 90% White individuals.

library(Matching)
library(readr)
library(dplyr)
set.seed(1) # for reproducibility

# Set parameters
min_coding <- 4 # Black or Black British
maj_coding <- 1 # White
nfolds <- 5

pheno <- "height"
prop_min <- 0.1 # one of 0 (White only), 1 (Black only) or 0.1
f <- 1 # outer CV fold, vary from 1 to 5

pheno_df <- read_tsv("data/all_vars.tsv") %>%
  rename(y = all_of(pheno)) %>%
  filter(!is.na(y) & parent_coding %in% c(min_coding, maj_coding)) %>%
  mutate(Tr = (parent_coding == min_coding))

dir.create(file.path("data", pheno, paste0("fold", f)),
           showWarnings = FALSE, recursive = TRUE)

# Construct training sets
if (prop_min == 0) {

  sample_df <- pheno_df %>%
    filter(parent_coding == maj_coding) %>%
    mutate(fold = sample(cut(seq_len(n()),
                             breaks = nfolds, labels = FALSE))) %>%
    dplyr::select(eid, fold)

  outfile <- file.path("data", pheno,
                       paste0("fold", f),
                       paste0("prop", prop_min, ".txt"))

  sample_df %>%
    filter(fold != f) %>%
    pull(eid) %>%
    rep(each = 2) %>%
    write(outfile, ncol = 2)
  
} else {
  # Train-test split
  train_ratio <- round((1 - prop_min) / prop_min)

  min_df <- pheno_df %>%
    filter(parent_coding == min_coding) %>%
    mutate(fold = sample(cut(seq_len(n()),
                             breaks = nfolds, labels = FALSE))) %>%
    dplyr::select(eid, Tr, Sex, Age, fold)

  maj_df <- pheno_df %>%
    filter(parent_coding == maj_coding) %>%
    mutate(fold = NA_integer_) %>%
    dplyr::select(eid, Tr, Sex, Age, fold)

  match_df <- min_df %>%
    filter(fold != f) %>%
    bind_rows(maj_df)

  if (prop_min == 1) {
    train_id <- min_df %>%
      filter(fold != f) %>%
      pull(eid)
  } else {
    match_out <- match_df %>%
      with(Match(Tr = Tr, X = cbind(Sex, Age), exact = c(TRUE, FALSE),
                 M = train_ratio, Z = eid, ties = FALSE, replace = FALSE))

    idx <- unique(c(match_out$index.treated, match_out$index.control))
    train_id <- match_df$eid[idx]
  }
  outfile <- file.path("data", pheno,
                       paste0("fold", f),
                       paste0("prop", prop_min, ".txt"))
  write(rep(train_id, each = 2), outfile, ncol = 2)
}
