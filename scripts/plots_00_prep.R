###################
## Load packages ##
###################

library(readr)
library(dplyr)
library(ggplot2)
library(purrr)
library(patchwork)

theme_set(theme_minimal())

scale_colour_discrete <- function(...) scale_colour_viridis_d(..., end = 0.9, option = "D")
scale_fill_discrete <- function(...) scale_fill_viridis_d(..., end = 0.9)

###################
## Get UKB codes ##
###################

quant_pheno_codes <- c("A50", "A21001", "A30100", "A30040", "A30090", 
                       "A30300", "A30130", "A30070", "A30210", "A30120")
quant_phenos <- c("height", "BMI", "MPV", "MCV", "platelet",
                  "reticulocyte", "monocyte", "erythrocyte", "eosinophill", "lymphocyte")

binary_pheno_codes <- c("NC1226", "NC1111", "I48", "K57", "N81")
binary_phenos <- c("hypothyroidism", "asthma", "AFib", "DDI", "FGP")

pheno_codes <- c(quant_pheno_codes, binary_pheno_codes)
phenos <- c(quant_phenos, binary_phenos)
names(phenos) <- pheno_codes

full_ind <- c(1, 4, 12, 15)
full_phenos <- phenos[full_ind]
full_pheno_codes <- pheno_codes[full_ind]

all_ancestries <- c("AFR", "CSA", "EAS", "AMR", "MID")

ntrain_df <- read_tsv("output/ntrain.tsv", show_col_types = FALSE) %>%
  mutate(n_min = 0.1 * ntrain,
         n_eur = 0.9 * ntrain)
