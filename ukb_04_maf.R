# Script to calculate MAF

library(dplyr)
library(readr)

ncores <- 2
pfile <- "data/all_vars.tsv"
gfile <- "data/ukb_cal_pgen/ukb_cal_all"

pheno_df <- read_tsv(pfile)
train_black_id <- pheno_df %>%
  filter(parent_meaning == "Black or Black British") %>%
  pull(eid) %>%
  sort()
train_black_id <- sapply(train_black_id, function(x) paste(x, x, sep = "_"))

train_white_id <- pheno_df %>%
  filter(parent_meaning == "White") %>%
  pull(eid) %>%
  sort()
train_white_id <- sapply(train_white_id, function(x) paste(x, x, sep = "_"))
ids <- snpnet::readIDsFromPsam(paste0(gfile, '.psam'))
# Note that first 130 IDs are 'missing'

vars <- mutate(rename(data.table::fread(cmd=paste0('zstdcat ',
                                                   paste0(gfile, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID

pvar <- pgenlibr::NewPvar(paste0(gfile, ".pvar.zst"))

pgen_black <- pgenlibr::NewPgen(paste0(gfile, ".pgen"), pvar = pvar,
                                sample_subset = match(train_black_id, ids))

pgen_white <- pgenlibr::NewPgen(paste0(gfile, ".pgen"), pvar = pvar,
                                sample_subset = match(train_white_id, ids))

beta_df <- tibble(snp = character(),
                  maf_black = double(),
                  maf_white = double())

for (snp in split(vars, ceiling(seq_along(vars)/100))) {
  X_black <- pgenlibr::ReadList(pgen_black, match(snp, vars))
  maf_black <- colMeans(X_black, na.rm = T)
        
  X_white <- pgenlibr::ReadList(pgen_white, match(snp, vars))
  maf_white <- colMeans(X_white, na.rm = T)
        
  beta_df <- beta_df %>%
    add_row(snp, maf_black, maf_white)
}

betafile <- "output/beta_all.tsv"
write_tsv(beta_df, betafile)
