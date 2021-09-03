# This script fits the LASSO on the inner folds in order to select the regularisation parameter lambda.

### Set up parameters ---

# Install version of snpnet: devtools::install_github("brieuclehmann/snpnet")
library(snpnet)
library(readr)
library(dplyr)
set.seed(1)

ncores <- 2 # Number of available cores
mem <- ncores * 10e3
min_coding <- 4 # Black or Black British
pheno <- "height"
prop <- 0.1
pow <- 0
fold <- 1L # outer fold
i <- 1L # inner fold

kfile  <- file.path("data", pheno,
                    paste0("fold", fold),
                    paste0("prop", prop, ".txt"))
pfile <- "data/all_vars.tsv"
gfile <- "data/ukb_cal_pgen/ukb_cal_all"

results.dir <- file.path("temp", pheno,
                         paste0("fold", fold),
                         paste0("prop", prop),
                         paste0("pow", format(pow, nsmall = 1)))
configs <- readRDS(file.path(results.dir, "configs.RDS"))
full.lams <- scan(file.path(results.dir, "full.lams.txt"))
covars <- configs$covariates

cv.phenotype.file <- paste0(configs[["cv.full.prefix"]], ".tsv")
weights <- read_tsv(cv.phenotype.file) %>%
  pull(weights)

cv_configs <- configs
cv_configs$results.dir <- paste0(configs$results.dir, "/fold", i)
cv_configs$nCores <- ncores
cv_configs$mem <- mem
cv_configs$stopping.lag <- 4
cv_configs$gcount.full.prefix <- NULL
cv_configs$verbose <- TRUE
cv_configs$save <- TRUE
      
### Fit lasso for ith fold ---
outfile <- file.path(cv_configs$results.dir, "metric.val.txt")

if (!file.exists(outfile)) {
  intermediate_files <- list.files(file.path(results.dir, "results"),
                                   pattern = "output*") %>%
    tools::file_path_sans_ext()
  prev_iters <- as.integer(gsub("output_iter_", "", intermediate_files))
  if (length(prev_iters) > 0) {
    cv_configs$prev_iter <- max(prev_iters)
  }

  out <- snpnet(gfile,
                cv.phenotype.file,
                pheno,
                covariates = covars,
                weights = weights,
                full.lams = full.lams,
                split.col = paste0("fold", i),
                configs = cv_configs)
  write(out$metric.val, file = outfile, ncolumns = 1)

  snpnet:::cleanUpIntermediateFiles(out$configs)
}
