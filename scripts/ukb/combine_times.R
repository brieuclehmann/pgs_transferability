

### Set up parameters ---

.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(readr)
library(dplyr)
set.seed(1)

out_file <- snakemake@output[[1]]
this_pheno <- pheno <- snakemake@wildcards[["pheno"]]
this_min_ancestry <- min_ancestry <- snakemake@wildcards[["min_ancestry"]]
variants <- snakemake@wildcards[["v"]]


outdir <- file.path(
    "output", "ukb",
    paste0("v~", variants),
    paste0("pheno~", pheno),
    paste0("min_ancestry~", min_ancestry)
)

time_files <- list.files(outdir, pattern = "runtimes.tsv", recursive = TRUE, )

time_df = tibble()
for (time_file in time_files) {
    prop_min <- gsub(".*prop_min~(.*?)/.*", "\\1", time_file)
    f <- gsub(".*fold~(.*?)/.*", "\\1", time_file)
    pow <- gsub(".*pow~(.*?)/.*", "\\1", time_file)

    this_df <- read_tsv(paste0(outdir, "/", time_file))
    this_df$pow <- pow
    this_df$prop_min <- prop_min
    this_df$fold <- f

    time_df <- rbind(time_df, this_df)
}

time_df$vtype <- variants
time_df$pheno <- pheno
time_df$min_ancestry <- min_ancestry

write_tsv(time_df, out_file)
