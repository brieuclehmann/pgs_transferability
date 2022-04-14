library(readr)
library(dplyr)
library(xtable)

ukbb_df <- read_csv("data/Pan-UK Biobank phenotype manifest - phenotype_manifest.csv")

pan_quant_pheno_codes <- c(50, 21001, 30100, 30040, 30090, 
                       30300, 30130, 30070, 30210, 30120)
pan_binary_pheno_codes <-c("I48", "K57", "N81")
nc_codes <- c(1111, 1226)

h2_df <- ukbb_df %>% 
  filter(phenocode %in% c(pan_quant_pheno_codes, pan_binary_pheno_codes) |
           (phenocode == 20002 & coding %in% nc_codes)) %>%
  mutate(description = if_else(phenocode == 20002, 
                               coding_description, 
                               description),
         pheno = if_else(phenocode == 20002, 
                         paste0("NC", coding), 
                         phenocode),
         pheno = if_else(pheno %in% pan_quant_pheno_codes,
                         paste0("A", pheno),
                         pheno)) %>%
  select_at(vars(pheno, description,
                 starts_with("saige"))) %>%
  rename_with(~ gsub("saige_heritability_", "", .x), 
              starts_with("saige")) %>%
  tidyr::pivot_longer(-c("pheno", "description"),
                      names_to = "pop", values_to = "h2")

count_df <- read_tsv("output/counts.tsv", show_col_types = FALSE) %>%
  mutate(trait = phenos[pheno],
         is_binary = pheno %in% binary_pheno_codes,
         is_full = pheno %in% full_pheno_codes,
         n_cases = if_else(is_binary, n_cases, NA_real_)) %>%
  arrange(-is_full, is_binary, trait) %>%
  filter(pop %in% c("EUR", "AFR") | pheno %in% full_pheno_codes) %>%
  left_join(h2_df, by = c("pheno", "pop")) %>%
  select(ancestry = pop, trait, n_total, n_cases, h2)

print(xtable(count_df, digits = c(rep(0, 5), 4)), include.rownames = FALSE)
