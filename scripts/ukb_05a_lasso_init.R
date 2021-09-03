# This script initialises the lasso cross-validation for snpnet.

### Set up parameters ---

# Install version of snpnet: devtools::install_github("brieuclehmann/snpnet")
library(snpnet)
library(readr)
library(dplyr)
set.seed(1)

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
covars <- c("Sex", "Age", paste0("PC", 1:10), paste0("PC", 1:10, "_Sex"))

keep_ids <- read.table(kfile)[, 1]
ntrain <- length(keep_ids)

ids <- readIDsFromPsam(paste0(gfile, ".psam"))
id_df <- data.frame(ID = ids, stringsAsFactors = F) %>%
  mutate(sort_order = 1:n())

weights <- read_tsv(pfile) %>%
  filter(eid %in% keep_ids) %>%
  group_by(parent_coding) %>%
  transmute(ID = paste(eid, eid, sep = "_"),
            weight = ntrain / n()) %>%
  left_join(id_df, by = "ID") %>%
  ungroup() %>%
  mutate(weight = ntrain * (weight^pow / sum(weight^pow))) %>%
  arrange(sort_order) %>%
  pull(weight)

configs <- list(results.dir = file.path("temp", pheno,
                                        paste0("fold", fold),
                                        paste0("prop", prop),
                                        paste0("pow", format(pow, nsmall = 1))),
                cv.dir = "cv",
                keep = kfile)
configs$cv.full.prefix <- file.path(configs$results.dir, configs$cv.dir,
                                    "snpnet.cvpheno")

### Prepare phenotype data.table ---
  
phe <- readPheMaster(phenotype.file = pfile, psam.ids = ids,
                     family = NULL, covariates = covars,
                     phenotype = pheno, status = NULL,
                     split.col = NULL, configs = configs)
  
folds <- cut(seq(1, nrow(phe)), breaks = nfolds, labels = FALSE)
foldids <- sample(folds, length(folds))

dir.create(configs$cv.full.prefix, showWarnings = FALSE, recursive = TRUE)
cv.phenotype.file <- paste0(configs[["cv.full.prefix"]], ".tsv")
phe %>%
  dplyr::mutate(fold = foldids, tmp = "val") %>%
  tidyr::pivot_wider(names_prefix = "fold", names_from = fold,
                     values_from = tmp, values_fill = list(tmp = "train")) %>%
  dplyr::mutate(weights = weights) %>%
  data.table::fwrite(cv.phenotype.file, sep = "\t")

### Pre-compute lambda path ---

if (length(unique(phe$Sex)) == 1)
  covars <- c("Age", paste0("PC", 1:10))

rm(phe)
gc()

configs$covariates <- covars
full.lams <- snpnet(genotype.pfile = gfile,
                    phenotype.file = pfile,
                    phenotype = pheno,
                    covariates = covars,
                    weights = weights,
                    configs = configs,
                    lambda_only = TRUE)

dir.create(configs$results.dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(configs, file.path(configs$results.dir, "configs.RDS"))

full.lams.file <- file.path(configs$results.dir, "full.lams.txt")
write(full.lams, file = full.lams.file, ncolumns = 1)
