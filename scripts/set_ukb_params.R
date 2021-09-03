### Set UKB analysis parameters ###
library(dplyr)
library(readr)

# Basic analysis
covars <- c(
    "pop", "age", "sex", "age_sex", "age2", "age2_sex",
    paste0("PC", 1:10), paste0("PC", 1:10, "_sex")
)
pheno_codes <- read_tsv("data/all_vars.tsv") %>%
    dplyr::select(-c(all_of(covars), eid)) %>%
    colnames()

params <- expand.grid(
    pheno = pheno_codes,
    min_ancestry = "AFR",
    prop_min = c(0, 0.1, 1),
    fold = seq(5)
)

params$prop_min <- format(params$prop_min, digits = 2)

write_csv(params, "data/ukb_basic_params.csv")

# Enhanced analysis