library(tibble)
library(dplyr)
library(readr)
library(tidyr)

all_phenos <- c(
    "A21001", "A30100", "A30040", "A30090", "A30300", "A30130",
    "A30070", "A30210", "A30120", "A50", "NC1226", "NC1111", "I48", "K57", "N81"
)
min_ancestries <- c("AFR", "EAS", "MID", "CSA", "AMR")
full_phenos <- c("A50", "A30040", "NC1111", "N81")

fracs <- seq(0, 1, 0.2)
nfolds <- 5

ntrain_df <- tibble(
    pheno = character(),
    min_ancestry = character(),
    fold = integer(),
    ntrain = integer()
)

pheno_df <- read_tsv("data/all_vars.tsv") %>%
    select(pop, sex, all_of(all_phenos)) %>%
    pivot_longer(-c(1, 2), names_to = "pheno") %>%
    mutate(value = if_else(sex == 1 & pheno == "N81", NA_real_, value))


count_df <- pheno_df %>%
    group_by(pop, pheno) %>%
    summarise(
        n_total = sum(!is.na(value)),
        n_cases = sum(value, na.rm = TRUE),
        .groups = "drop"
    )

for (pheno in all_phenos) {
    for (fold in seq_len(nfolds)) {
        for (min_ancestry in min_ancestries) {
            kfile <- file.path(
                "data", "train_ids",
                paste0("pheno~", pheno),
                paste0("min_ancestry~", min_ancestry),
                paste0("fold~", fold),
                paste0("type~min"),
                paste0("frac~1.0.txt")
            )
            if (file.exists(kfile)) {
                ntrain <- length(read.table(kfile)[, 1])
                ntrain_df <- ntrain_df %>%
                    add_row(pheno, min_ancestry, fold, ntrain)
            }
        }
    }
}

write_tsv(count_df, "output/counts.tsv")
write_tsv(ntrain_df, "output/ntrain.tsv")
