# This script fits the final lasso for the regularisation parameter chosen in 05c_lasso_cv.R

### Set up parameters ---

# Install version of snpnet: devtools::install_github("brieuclehmann/snpnet")
library(snpnet)
library(readr)
library(dplyr)
set.seed(1)

ncores <- 2
min_coding <- 4 # Black or Black British
pheno  <- "height"
prop   <- 0.1
pow    <- 0
fold   <- 1
nfolds <- 5

kfile  <- file.path("data", pheno,
                    paste0("fold", fold),
                    paste0("prop", prop, ".txt"))
pfile  <- "data/all_vars.tsv"
gfile  <- "data/ukb_cal_pgen/ukb_cal_all"

results.dir <- file.path("temp", pheno,
                         paste0("fold", fold),
                         paste0("prop", prop),
                         paste0("pow", format(pow, nsmall = 1)))
configs   <- readRDS(file.path(results.dir, "configs.RDS"))
full.lams <- scan(file.path(results.dir, "full.lams.txt"))
covars <- configs$covariates

cv.phenotype.file <- paste0(configs[["cv.full.prefix"]], ".tsv")
weights <- read_tsv(cv.phenotype.file) %>%
  pull(weights)

### Gather cross-validation results ---

inner_folds <- seq_len(nfolds)
cvout <- inner_folds %>%
  sapply(function(x) file.path(results.dir, paste0("fold", x),
                               "metric.val.txt")) %>%
  sapply(scan)

cvm <- apply(cvout, 1, mean)
cvsd <- apply(cvout, 1, sd)
lambda.min <- full.lams[which.max(cvm)]
lambda_na <- apply(cvout, 1, function(x) !all(is.na(x)))

### Fit lasso on full dataset ---

fit.lams <- full.lams[lambda_na]
configs$nCores <- ncores
configs$verbose <- TRUE
configs$save <- TRUE

outdir <- file.path("output", pheno, paste0("fold", fold), paste0("prop", prop))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, paste0("pow", format(pow, nsmall = 1), ".RDS"))

#if (!file.exists(outfile)) {
  intermediate_files <- list.files(file.path(results.dir, "results"),
                                   pattern = "output*") %>%
    tools::file_path_sans_ext()
  prev_iters <- as.integer(gsub("output_iter_", "", intermediate_files))
  if (length(prev_iters) > 0) {
    configs$prev_iter <- max(prev_iters)
  }

  snpnet.object <- snpnet(gfile,
                          pfile,
                          pheno,
                          covariates = covars,
                          weights = weights,
                          full.lams = fit.lams,
                          configs = configs,
                          mem = 8e3 * ncores)

  max_cv_ind <- which.max(cvm)
  if (max_cv_ind == length(cvm) || 
      is.na(cvm[min(length(cvm), max_cv_ind + 1)]))
    warning("Cross-validation may have stopped early. Consider increasing   stopping.lag.")

  ### Save output ---
  mod <- list(snpnet.object = snpnet.object,
              cvm = cvm,
              cvsd = cvsd,
              cvout = cvout,
              lambda.min = lambda.min,
              full.lams = full.lams,
              fit.lams = fit.lams)

  saveRDS(mod, outfile)
  snpnet:::cleanUpIntermediateFiles(mod$configs)
#}
