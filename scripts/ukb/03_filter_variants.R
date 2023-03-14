# Script to calculate MAF and filter very rare SNPs

library(readr)
library(dplyr)

variants <- snakemake@wildcards[["v"]]
ancestry <- snakemake@wildcards[["min_ancestry"]]
this_chrom <- snakemake@wildcards[["chrom"]]
out_file <- snakemake@output[[1]]

pan_ukbb_variant_file <- "data/full_variant_qc_metrics.txt"
pan_ukbb_df <- read_tsv(
  pan_ukbb_variant_file,
  col_types = cols_only(
    chrom = "c", pos = "i", ref = "c", alt = "c", rsid = "c", varid = "c",
    high_quality = "l", info = "d", af_AFR = "d", af_AMR = "d", af_CSA = "d",
    af_EAS = "d", af_EUR = "d", af_MID = "d"
  )
)
this_df <- pan_ukbb_df %>%
  filter(high_quality & chrom == this_chrom)


if (variants == "tagged") {
  # QC filter for tagged SNPs

  # Downloaded from http://kunertgraf.com/data/files/snp_maf_comparison.csv
  af_file <- "data/snp_maf_comparison.csv"
  af_df <- read_csv(af_file) %>%
    filter(oob == 1)

  batches_qc <- c(
    paste0("Batch_b", formatC(1:95, width = 3, flag = "0"), "_qc"),
    paste0("UKBiLEVEAX_b", 1:11, "_qc")
  )
  call_ukbb_file <- "/well/ukbb-wtchg/v2/qc/ukb_snp_qc.txt"
  call_ukbb_df <- read_delim(call_ukbb_file, delim = " ")
  qc_all <- call_ukbb_df %>%
    select(all_of(batches_qc)) %>%
    rowMeans()
  # 'array == 2' genotyped on both arrays
  call_ukbb_df <- call_ukbb_df %>%
    select(-batches_qc) %>%
    mutate(qc_all = qc_all) %>%
    filter(qc_all == 1 & array == 2 & !rs_id %in% af_df$rsid) %>%
    select(
      rsid = rs_id,
      chrom = chromosome,
      pos = position,
      ref = allele1_ref,
      alt = allele2_alt
    ) %>%
    mutate(chrom = if_else(chrom == 23, "X", as.character(chrom)))

  this_df <- this_df %>%
    filter(info == 1) %>%
    semi_join(call_ukbb_df, by = c("chrom", "pos", "ref", "alt"))
}

# MAF filter
this_df <- this_df %>%
  rename(maf = paste0("af_", ancestry)) %>%
  mutate(
    maf_maj = pmin(af_EUR, 1 - af_EUR),
    maf_min = pmin(maf, 1 - maf)
  ) %>%
  filter(maf_min >= 0.01 | maf_maj >= 0.01) %>%
  mutate(varid2 = paste0(chrom, ":", pos, "_", ref, "_", alt))

# Save output
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write(this_df$varid2, out_file, ncol = 1)