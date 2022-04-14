### Set up parameters ---

library(snpnet)
library(purrr)
library(foreach)
library(readr)
library(dplyr)
devtools::load_all()
set.seed(1)

pheno  <- Sys.getenv("pheno")
min_coding <- as.integer(Sys.getenv("min_coding"))
this_prop   <- as.double(Sys.getenv("prop"))
this_fold   <- as.integer(Sys.getenv("fold"))
ncores <- as.integer(Sys.getenv("NSLOTS"))

kfile  <- file.path("data", paste0("min", min_coding), pheno,
                    paste0("fold", this_fold),
                    paste0("prop", this_prop, ".txt"))
pfile  <- "data/all_vars.tsv"
gfile  <- "data/ukb_cal_pgen/ukb_cal_all"

train_id <- read.table(kfile)[, 1]
pheno_df <- read_tsv(pfile) %>%
  mutate(ID = paste(FID, IID, sep = "_"), train = (eid %in% train_id)) %>%
  rename(y = all_of(pheno))

### Prepare test set ---

kfile_test <- file.path("data", paste0("min", min_coding), pheno,
                        paste0("fold", this_fold),
                        paste0("prop", this_prop, "_test.txt"))
pheno_df %>%
  filter(!train) %>%
  pull(eid) %>%
  rep(each = 2) %>%
  write(kfile_test, ncol = 2)

gfile_test <- file.path("data",paste0("min", min_coding), pheno,
                        paste0("fold", this_fold),
                        paste0("prop", this_prop, "_ukb_cal_test"))

cmd_plink2 <- paste("plink2",
                    "--pfile", gfile,
                    "--threads", ncores,
                    "--keep", kfile_test,
                    "--make-pgen vzs",
                    "--out", gfile_test)
system(cmd_plink2, intern = F)

### Evaluate fit on test set ---
beta_file <- file.path("output", paste0("min", min_coding), "beta_all.tsv")

out_files <- list.files(file.path("output", paste0("min", min_coding), pheno,
                                  paste0("fold", this_fold),
                                  paste0("prop", this_prop)),
                        pattern = "^pow.*.RDS$")
pow_range <- gsub(".RDS", "", gsub("pow", "", out_files))
pow_range <- pow_range[!is.na(pow_range)]

out_df <- foreach(this_pow = pow_range, .combine = rbind) %do% {

  outfile <- file.path("output", paste0("min", min_coding), pheno,
                       paste0("fold", this_fold),
                       paste0("prop", this_prop),
                       paste0("pow", this_pow, ".RDS"))
  mod <- readRDS(outfile)

  if (this_prop %in% c(-1, 0)) {
    lambda_ind <- which.max(mod$metric.val)
    snpnet.object <- mod
  } else {
    lambda_ind <- which(mod$full.lams == mod$lambda.min)
    snpnet.object <- mod$snpnet.object
  }

  configs <- snpnet.object$configs
  covars <- configs$covariates
  a0 <- snpnet.object$a0
  beta <- snpnet.object$beta[[lambda_ind]]
  stats <- snpnet.object$stats
  feature_names <- names(beta[beta != 0])
  feature_names <- setdiff(feature_names, covars)

  ids <- readIDsFromPsam(paste0(gfile_test, '.psam'))
  configs_temp <- list(zstdcat.path = "zstdcat")
  test_ids <- readPheMaster(pfile, ids, NULL, covars, pheno, NULL, NULL, configs_temp)$ID

  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0('zstdcat ', paste0(gfile_test, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID 
  pvar <- pgenlibr::NewPvar(paste0(gfile_test, '.pvar.zst'))
  pgen <- pgenlibr::NewPgen(paste0(gfile_test, '.pgen'), pvar = pvar, sample_subset = match(test_ids, ids))
  pgenlibr::ClosePvar(pvar)

  maf_cuts <- c(0, 0.01, 0.05, 0.2, 0.5)
  beta_df <- read_tsv(beta_file) %>%
    filter(snp %in% feature_names) %>%
    mutate(maf_black_grp = cut(0.5 * pmin(maf_black, 2 - maf_black), maf_cuts),
           maf_white_grp = cut(0.5 * pmin(maf_white, 2 - maf_white), maf_cuts))

  pred_split_df <- tibble(ID = character(), pow = double(), 
                          ethnicity = character(), 
                          grp = character(), pred = double())
  for (grp in levels(beta_df$maf_black_grp)) {
    snps <- beta_df %>%
      filter(maf_black_grp == grp) %>%
      pull(snp)
    
    if (length(snps) > 0) {
        features <- snpnet:::prepareFeatures(pgen, vars, snps, stats)
        features_single <- as.matrix(features[, snps, with = F])
      } else {
        features_single <- matrix(0, length(test_ids), 0)
    }

    pred_single <- features_single %*% beta[snps]
    pred_split_df <- pred_split_df %>%
      add_row(ID = test_ids, pow = this_pow, ethnicity = "Black",
              grp = grp, pred = pred_single[ ,1])
   
    snps <- beta_df %>%
      filter(maf_white_grp == grp) %>%
      pull(snp)
    if (length(snps) > 0) {
        features <- snpnet:::prepareFeatures(pgen, vars, snps, stats)
        features_single <- as.matrix(features[, snps, with = F])
      } else {
        features_single <- matrix(0, length(test_ids), 0)
    }

    pred_single <- features_single %*% beta[snps]
    pred_split_df <- pred_split_df %>%
      add_row(ID = test_ids, pow = this_pow, ethnicity = "White", 
              grp = grp, pred = pred_single[ ,1])
  }

  temp_df <- pred_split_df %>%
    inner_join(pheno_df, by = "ID") %>%
    select(eid, pow, ethnicity, grp, pred)

  temp_df
}

### Save output and remove temporary test files

pred_file <- file.path("output", paste0("min", min_coding), pheno,
                      paste0("fold", this_fold),
                      paste0("prop", this_prop),
                      "pred_split.tsv")
write_tsv(out_df, pred_file)

file.remove(list.files(dirname(gfile_test), basename(gfile_test),
            full.names = TRUE))
file.remove(kfile_test)
