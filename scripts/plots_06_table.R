library(readr)
library(dplyr)
library(xtable)


count_df <- read_tsv("output/counts.tsv", show_col_types = FALSE) %>%
  mutate(trait = phenos[pheno],
         is_binary = pheno %in% binary_pheno_codes,
         is_full = pheno %in% full_pheno_codes,
         n_cases = if_else(is_binary, n_cases, NA_real_)) %>%
  arrange(-is_full, is_binary, trait) %>%
  filter(pop %in% c("EUR", "AFR") | pheno %in% full_pheno_codes) %>%
  select(ancestry = pop, trait, n_total, n_cases)

xtable(count_df, digits = 0)
