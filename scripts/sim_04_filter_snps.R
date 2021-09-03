# Script to filter SNPs according to MAF

library(dplyr)
library(readr)

maf_thresh <- 0.01

maf_df <- read_csv("data/maf.csv") %>%
    filter(pmin(mafAFR, 1 - mafAFR) >= maf_thresh |
        pmin(mafEUR, 1 - mafEUR) >= maf_thresh |
        pmin(mafEAS, 1 - mafEAS) >= maf_thresh)