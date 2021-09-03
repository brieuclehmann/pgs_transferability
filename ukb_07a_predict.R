# Script to generate test set predictions

### Set up parameters ---

library(snpnet)
library(purrr)
library(foreach)
library(readr)
library(dplyr)
set.seed(1)

pheno  <- "height"
min_coding <- 4
prop   <- 0.1
fold   <- 1
ncores <- 2

kfile  <- file.path("data", pheno,
                    paste0("fold", fold),
                    paste0("prop", prop, ".txt"))
pfile  <- "data/all_vars.tsv"
gfile  <- "data/ukb_cal_pgen/ukb_cal_all"

train_id <- read.table(kfile)[, 1]
pheno_df <- read_tsv(pfile) %>%
  mutate(ID = paste(FID, IID, sep = "_"), train = (eid %in% train_id)) %>%
  rename(y = all_of(pheno))

### Prepare test set ---

kfile_test <- file.path("data", pheno,
                        paste0("fold", fold),
                        paste0("prop", prop, "_test.txt"))
pheno_df %>%
  filter(!train) %>%
  pull(eid) %>%
  rep(each = 2) %>%
  write(kfile_test, ncol = 2)

gfile_test <- file.path("data", pheno,
                        paste0("fold", fold),
                        paste0("prop", prop, "_ukb_cal_test"))

cmd_plink2 <- paste("plink2",
                    "--pfile", gfile,
                    "--threads", ncores,
                    "--keep", kfile_test,
                    "--make-pgen vzs",
                    "--out", gfile_test)
system(cmd_plink2, intern = F)

### Evaluate fit on test set ---

out_files <- list.files(file.path("output", pheno,
                                  paste0("fold", fold),
                                  paste0("prop", prop)),
                        pattern = "^pow.*.RDS$")
pow_range <- gsub(".RDS", "", gsub("pow", "", out_files))
pow_range <- pow_range[!is.na(pow_range)]

out_df <- foreach(pow = pow_range, .combine = rbind) %do% {

  outfile <- file.path("output", pheno,
                       paste0("fold", fold),
                       paste0("prop", prop),
                       paste0("pow", pow, ".RDS"))
  mod <- readRDS(outfile)

  if (prop == 0) {
    lambda_ind <- which.max(mod$metric.val)
    snpnet.object <- mod
  } else {
    lambda_ind <- which(mod$full.lams == mod$lambda.min)
    snpnet.object <- mod$snpnet.object
  }
  
  covars <- snpnet.object$configs$covariates
  pred <- predict_snpnet(snpnet.object,
                         new_genotype_file = gfile_test,
                         new_phenotype_file = pfile,
                         phenotype = pheno,
                         covariate_names = covars,
                         idx = lambda_ind)

  this_pred_df <- tibble(ID = rownames(pred$prediction$train),
                         pred = pred$prediction$train[, 1]) %>%
    inner_join(pheno_df, by = "ID") %>%
    mutate(pred_bin = 1 / (1 + exp(-pred)), type = "lasso")
  
  this_pred_df %>%
    select(eid, type, pred) %>%
    mutate(pow = pow)

}

### Save output and remove temporary test files

pred_file <- file.path("output", pheno,
                      paste0("fold", fold),
                      paste0("prop", prop),
                      "pred.tsv")
write_tsv(pred_df, pred_file)

file.remove(list.files(dirname(gfile_test), basename(gfile_test),
            full.names = TRUE))
file.remove(kfile_test)
