# This script fits the LASSO using a validation set to select lambda.

### Set up parameters ---

library(snpnet)
library(readr)
library(dplyr)
set.seed(1)

min_coding <- 4 # Black or Black British
ncores <- 2
pheno  <- "height"
fold   <- 1
prop   <- 0
pow    <- 0

kfile  <- file.path("data", pheno, 
                    paste0("fold", fold),
                    paste0("prop", prop, ".txt"))
pfile  <- "data/all_vars.tsv"
gfile  <- "data/ukb_cal_pgen/ukb_cal_all"
covars <- c("Sex", "Age", paste0("PC", 1:10), paste0("PC", 1:10, "_Sex"))

keep_ids <- read.table(kfile)[, 1]
ntrain <- length(keep_ids)

ids <- readIDsFromPsam(paste0(gfile, ".psam"))
id_df <- data.frame(ID = ids, stringsAsFactors = F) %>%
  mutate(sort_order = 1:n())

results.dir <- file.path("temp", pheno,
                         paste0("fold", fold),
                         paste0("prop", prop),
                         paste0("pow", format(pow, nsmall = 1)))
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)

configs <- list(results.dir = results.dir,
                stopping.lag = 2,
                keep = kfile,
                nCores = ncores,
                mem = 10e3 * ncores,
                use.glmnetPlus = TRUE,
                verbose = TRUE,
                save = TRUE)

### Set up train-validation split ---
phe <- readPheMaster(phenotype.file = pfile, psam.ids = ids,
                     family = NULL, covariates = covars, phenotype = pheno,
                     status = NULL, split.col = NULL, configs = configs)
train_ids <- sample(phe$ID, round(0.8 * nrow(phe)))

new_pfile <- file.path(results.dir, "pheno.tsv")
phe %>%
  dplyr::mutate(split = if_else(ID %in% train_ids, "train", "val")) %>%
  data.table::fwrite(new_pfile, sep = "\t")

if (length(unique(phe$Sex)) == 1)
  covars <- c("Age", paste0("PC", 1:10))

rm(phe)
gc()

### Fit lasso ---

outdir <- file.path("output",
                    pheno, paste0("fold", fold), 
                    paste0("prop", prop))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, paste0("pow", format(pow, nsmall = 1), ".RDS"))

if (!file.exists(outfile)) {
  intermediate_files <- list.files(file.path(results.dir, "results"),
                                   pattern = "output*") %>%
    tools::file_path_sans_ext()
  prev_iters <- as.integer(gsub("output_iter_", "", intermediate_files))
  if (length(prev_iters) > 0) {
    configs$prev_iter <- max(prev_iters)
  }

  system.time(
    mod <- snpnet(gfile,
                  new_pfile,
                  pheno,
                  covariates = covars,
                  split.col = "split",
                  configs = configs)
  )
  saveRDS(mod, outfile)
  snpnet:::cleanUpIntermediateFiles(mod$configs)
}
