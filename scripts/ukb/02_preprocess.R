# Script to preprocess phenotypes

library(readr)
library(dplyr)
library(tidyr)

trait_df <- read_csv("data/selected_traits.csv") %>%
  mutate(
    phenotype = gsub("_irnt", "", phenotype),
    ukb_code = paste0(phenotype, "-0.0")
  )
traits <- trait_df$phenotype
icd10_codes <- traits[grep("^[A-Za-z].*$", traits)]
noncancer_codes <- traits[grep("20002", traits)]
dvt_code <- "6152_5"
nopain_code <- "6159_100"
traits <- traits[!traits %in% c(icd10_codes, noncancer_codes, dvt_code, nopain_code)]


ukb_df <- read_tsv("data/all_ukb_vars.tsv")

### Blood and body traits ###

phys_df <- ukb_df %>%
  select(eid, paste0(traits, "-0.0"))

names(phys_df)[-1] <- paste0("A", gsub("-0.0", "", names(phys_df[-1])))
### DVT and no pain ###
dvt_id <- "6152"
nopain_id <- "6159"
dvt_df <- ukb_df %>%
  select(eid, starts_with(dvt_id), starts_with(nopain_id)) %>%
  transmute(eid, dvt = `6152-0.0` == 5, nopain = `6159-0.0` == -7)


### Self-reported non-cancer disease ###
noncancer_id <- gsub("20002_", "", noncancer_codes)
noncancer_df <- ukb_df %>%
  select(eid, starts_with("20002")) %>%
  tidyr::pivot_longer(-eid) %>%
  filter(!is.na(value)) %>%
  select(eid, value) %>%
  filter(value %in% noncancer_id) %>%
  mutate(disease = TRUE) %>%
  distinct()

full_noncancer_df <- ukb_df %>%
  select(eid) %>%
  full_join(tibble(value = as.double(noncancer_id)), by = character()) %>%
  left_join(noncancer_df, by = c("eid", "value")) %>%
  mutate(disease = tidyr::replace_na(disease, FALSE)) %>%
  tidyr::pivot_wider(names_from = value, values_from = disease) 

names(full_noncancer_df)[-1] <- paste0("NC", names(full_noncancer_df[-1]))
### ICD10 codes ###

icd10_names <- trait_df %>%
  filter(phenotype %in% icd10_codes) %>%
  select(code = phenotype, description) %>%
  mutate(description = substr(description, 29, 100))

icd10_df <- ukb_df %>%
  select(eid, starts_with("41270")) %>%
  tidyr::pivot_longer(-eid) %>%
  filter(!is.na(value)) %>%
  mutate(code = substr(value, 1, 3)) %>%
  select(eid, code) %>%
  filter(code %in% icd10_codes) %>%
  mutate(disease = TRUE) %>%
  distinct()

full_icd10_df <- ukb_df %>%
  select(eid) %>%
  full_join(tibble(code = icd10_codes), by = character()) %>%
  left_join(icd10_df, by = c("eid", "code")) %>%
  mutate(disease = tidyr::replace_na(disease, FALSE)) %>%
  tidyr::pivot_wider(names_from = code, values_from = disease)


### Inferred genetic ancestry ---
bridge <- read_delim("data/links/ukb12788bridge31063.txt",
  delim = " ",
  col_names = c("eid", "s")
)

pan_ukbb_ancestry_file <- "data/all_pops_non_eur_pruned_within_pop_pc_covs.tsv"
ancestry_df <- read_tsv(pan_ukbb_ancestry_file) %>%
  filter(!related)
ancestry_df <- ancestry_df %>%
  left_join(bridge, by = "s") %>%
  select(-c(s, related, paste0("PC", 11:20))) %>%
  mutate_at(vars(starts_with("PC")),
    .funs = list(sex = ~ . * ancestry_df$sex)
  )
  

pheno_df <- phys_df %>%
  left_join(dvt_df, by = "eid") %>%
  left_join(full_noncancer_df, by = "eid") %>%
  left_join(full_icd10_df, by = "eid") %>%
  right_join(ancestry_df, by = "eid") %>%
  mutate(FID = eid, IID = eid) %>%
  filter(!is.na(eid)) # remove 15 individuals not found in bridge file

write_tsv(pheno_df, "data/all_vars.tsv")