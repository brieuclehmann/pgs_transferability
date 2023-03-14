library(dplyr)

# Read in available fields from UKBB application
# NB: this will be different for each application
ukb_fields <- scan("data/ukb_fields.txt")

# Filter out parent illness phenotypes and Number of self-reported illnesses
ukb_fields <- ukb_fields[!ukb_fields %in% c("20107", "20110", "135")]

gcor_file <- "data/genetic_correlation_matrices/geno_correlation.r2"
gcor <- readr::read_table2(gcor_file) %>%
    mutate(
        p1 = gsub(".*files/(.+).ldsc.*", "\\1", p1),
        p2 = gsub(".*files/(.+).ldsc.*", "\\1", p2)
    )

h2_df <- readr::read_tsv("data/ukb31063_h2_topline.02Oct2019.tsv") %>%
    mutate(phenotype2 = stringr::str_extract(phenotype, "[0-9]*"))

phe <- h2_df %>%
    filter((phenotype2 %in% ukb_fields | source == "icd10") &
        !is.na(confidence) & confidence %in% c("medium", "high") &
        (h2_observed > 0.05 | h2_liability > 0.05) &
        phenotype %in% gcor$p1) %>%
    arrange(desc(h2_observed)) %>%
    pull(phenotype)


i <- 1
while (i <= length(phe)) {
    id <- phe[i]
    corr_phe1 <- gcor %>%
        filter(p1 == id & abs(rg) >= 0.5) %>%
        pull(p2) %>%
        as.character()

    corr_phe2 <- gcor %>%
        filter(p2 == id & abs(rg) >= 0.5) %>%
        pull(p1) %>%
        as.character()

    phe <- phe[!phe %in% c(corr_phe1, corr_phe2)]
    i <- i + 1
}

trait_df <- h2_df %>%
    filter(phenotype2 %in% phe | phenotype %in% phe &
        confidence %in% c("medium", "high")) %>%
    select(phenotype, description, h2_observed, h2_liability) %>%
    arrange(desc(h2_observed))

readr::write_csv(trait_df, "data/selected_traits.csv")