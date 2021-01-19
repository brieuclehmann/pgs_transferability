# Script to compute minor allele frequency for simulated tree sequence

maf <- function(x, indSub = NULL, minor = TRUE){

  # Filter subjects
  if (!is.null(indSub))
    x <- filter_sub(x, indSub)

  # Define allele frequency function using call to python
  np <- reticulate::import("numpy")
  maf_def <- "maf = lambda x: ([g.genotypes.mean() for g in x.variants()])"
  reticulate::py_run_string(maf_def)
  get_maf <- reticulate::py$maf

  freq <- np$array(get_maf(x))

  # Return (minor) allele frequency
  if (minor)
    freq <- pmin(freq, 1-freq)

  freq
}

tskit <- reticulate::import("tskit")
ts <- tskit$load("data/ooa_chr20_10.trees")

mafAFR <- maf(ts, indSub = 1:200000, minor = FALSE)
mafEUR <- maf(ts, indSub = 200001:400000, minor = FALSE)
mafEAS <- maf(ts, indSub = 400001:600000, minor = FALSE)
mafTOT <- maf(ts, minor = FALSE)

maf_df <- tibble::as_tibble(cbind(mafTOT, mafAFR, mafEUR, mafEAS))
maf_df <- tibble::rownames_to_column(maf_df, "SNP")
maf_df$SNP <- as.integer(maf_df$SNP)

readr::write_csv(maf_df, "data/maf_chr20_10.csv")
