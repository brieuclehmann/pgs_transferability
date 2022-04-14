.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(readr)
library(dplyr)
library(tidyr)
library(bigsnpr)

fst_file <- snakemake@output[[1]]
pheno <- snakemake@wildcards[["pheno"]]
min_ancestry <- snakemake@wildcards[["min_ancestry"]]

pan_ukbb_variant_file <- "data/full_variant_qc_metrics.txt"
min_ukbb_df <- read_tsv(
    pan_ukbb_variant_file,
    col_types = cols_only(
        varid = "c",
        af_AFR = "d",
        an_AFR = "d",
        af_CSA = "d",
        an_CSA = "d",
        af_EAS = "d",
        an_EAS = "d",
        af_AMR = "d",
        an_AMR = "d",
        af_MID = "d",
        an_MID = "d"
    )
) %>%
    pivot_longer(-varid, names_sep = "_", names_to = c(".value", "pop")) %>%
        filter(pop == min_ancestry) %>%
        mutate(N = an / 2) %>%
        select(-c(pop, an))

eur_ukbb_df <- read_tsv(
    pan_ukbb_variant_file,
    col_types = cols_only(
        varid = "c",
        af_EUR = "d",
        an_EUR = "d"
    )
) %>%
    transmute(varid = varid, N = an_EUR / 2, af = af_EUR)

outdir <- file.path(
    "output", "ukb",
    paste0("pheno~", pheno),
    paste0("min_ancestry~", min_ancestry)
)
beta_files <- list.files(outdir, pattern = "beta.tsv", recursive = TRUE)

out_df <- tibble()
for (betafile in beta_files) {
    prop_min <- gsub(".*prop_min~(.*?)/.*", "\\1", betafile)
    f <- gsub(".*fold~(.*?)/.*", "\\1", betafile)
    pow <- gsub(".*pow~(.*?)/.*", "\\1", betafile)

    beta_df <- read_tsv(file.path(outdir, betafile)) %>%
        mutate(varid = sub("_[^_]+$", "", varname))

    eur_df <- eur_ukbb_df %>%
        filter(varid %in% beta_df$varid)

    min_df <- min_ukbb_df %>%
        filter(varid %in% beta_df$varid)

    list_af_df <- list(min_df, eur_df)
    fst <- snp_fst(list_af_df, overall = TRUE)

    out_df <- bind_rows(
        out_df,
        tibble(
            pheno = pheno, min_ancestry = min_ancestry,
            prop_min = prop_min, fold = f, pow = pow, fst = fst
        )
    )
}

write_tsv(out_df, fst_file)