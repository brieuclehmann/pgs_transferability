# Script to calculate MAF
.libPaths("/gpfs2/well/holmes/users/rxa753/pgs_transferability/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu")
print(.libPaths())

library(readr)
library(dplyr)

ancestry <- snakemake@wildcards[["ancestry"]]
this_chrom <- snakemake@wildcards[["chrom"]]

pan_ukbb_variant_file <- "data/full_variant_qc_metrics.txt"
pan_ukbb_df <- read_tsv(
  pan_ukbb_variant_file,
  col_types = cols_only(
    chrom = "c", pos = "i", rsid = "c", varid = "c", high_quality = "l",
    af_AFR = "d", af_AMR = "d", af_CSA = "d", af_EAS = "d", af_EUR = "d",
    af_MID = "d"
  )
)

this_df <- pan_ukbb_df %>%
  filter(high_quality & chrom == this_chrom) %>%
  rename(maf = paste0("af_", ancestry)) %>%
    mutate(
      maf_maj = pmin(af_EUR, 1 - af_EUR),
      maf_min = pmin(maf, 1 - maf)
    ) %>%
  filter(maf_min >= 0.01 | maf_maj >= 0.01) %>%
    mutate(varid2 = paste0(chrom, ":", pos, "_", ref, "_", alt))

out_rsids <- this_df$varid2
out_dir <- file.path("data", "snp_ids", paste0("ancestry~", ancestry))
out_file <- file.path(out_dir, paste0("chrom~", this_chrom, ".txt"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write(out_rsids, out_file, ncol = 1)