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


afr_implied_var_df <- basic_df %>% 
  filter(pop %in% c("AFR", "EUR") & prop_min %in% c(0,1)) %>% 
  select(pop, trait, `Training set`, implied_var) %>% 
  tidyr::pivot_wider(names_from = c("Training set", "pop"), values_from = "implied_var", values_fn = mean) %>%
  rename(g_min_min = 2, g_min_eur = 3, g_maj_min = 4, g_maj_eur = 5) %>%
  mutate(min_ancestry = "AFR", g_min_min = g_min_min / g_maj_eur, 
         g_min_eur = g_min_eur / g_maj_eur, g_maj_min = g_maj_min / g_maj_eur)

full_score_df <- tibble()
for (code in full_pheno_codes) {
  for (anc in all_ancestries) {
    this_file <- paste0("output/ukb/v~imputed/pheno~", code, "/min_ancestry~", anc, "/scores.tsv")
    if (file.exists(this_file)) {
      this_df <- read_tsv(this_file, show_col_types = FALSE) %>%
        mutate(prop_min = as.double(prop_min),
               pow = as.double(pow))
      full_score_df <- bind_rows(full_score_df, this_df)
    }
  }
}
full_score_df <- full_score_df %>%
  mutate(trait = phenos[pheno],
         `Training set` = case_when(prop_min == 0 ~ "European ancestry only",
                                    prop_min == 0.1 ~ "Multi-ancestry",
                                    TRUE ~ "Minority ancestry only")) %>%
  mutate(`Training set` = factor(`Training set`, 
                                 levels = c("European ancestry only",
                                            "Multi-ancestry", 
                                            "Minority ancestry only")))

full_score_df %>%
  filter(pop %in% c(min_ancestry, "EUR") & prop_min %in% c(0,1)) %>% 
  mutate(min_maj = if_else(min_ancestry == pop, "min", "maj")) %>%
  select(pop, min_ancestry, trait, min_maj, `Training set`, implied_var) %>% 
  tidyr::pivot_wider(names_from = c("Training set", "min_maj"), values_from = "implied_var", values_fn = mean)
  
g_eur_df <- full_score_df %>% 
  filter(pop == "EUR" & prop_min %in% c(0,1)) %>% 
  select(min_ancestry, `Training set`, trait, implied_var) %>%
  pivot_wider(names_from = "Training set", values_from = "implied_var", values_fn = mean) %>% 
  group_by(trait) %>% 
  mutate(g_maj_eur = max(`European ancestry only`, na.rm = TRUE),
         g_min_eur = `Minority ancestry only`) %>%
  select(trait, min_ancestry, g_min_eur, g_maj_eur)

g_min_df <- full_score_df %>% 
  filter((pop == min_ancestry & prop_min == 1) | (prop_min == 0 & pop != "EUR")) %>% 
  select(pop, `Training set`, trait, implied_var) %>%
  pivot_wider(names_from = "Training set", values_from = "implied_var", values_fn = mean) %>%
  rename(g_min_min = 3, g_maj_min = 4)

var_df <- g_eur_df %>%
  left_join(g_min_df, by = c("trait", "min_ancestry" = "pop")) %>%
  mutate(g_min_min = g_min_min / g_maj_eur, g_min_eur = g_min_eur / g_maj_eur, g_maj_min = g_maj_min / g_maj_eur) %>%
  bind_rows(afr_implied_var_df) %>%
  arrange(trait, min_ancestry) %>%
  select(trait, min_ancestry, g_maj_min, g_min_min, g_min_eur)

print(xtable(var_df, digits = 3), include.rownames = FALSE)
  
