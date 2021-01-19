# Script to evaluate predictive accuracy of PGS

### Set up parameters ---

library(readr)
library(dplyr)
set.seed(1)

pheno  <- "height"
min_coding <- 4
prop   <- 0.1
fold   <- 1

pfile  <- "data/all_vars.tsv"
kfile_min <- file.path("data", pheno,
                       paste0("fold", fold),
                       paste0("prop1.txt"))
min_id <- read.table(kfile_min)[,1]

pheno_df <- read_tsv(pfile) %>%
  select(eid, parent_meaning, Age, Sex, paste0("PC", 1:10), 
         paste0("PC", 1:10, "_Sex"), y = all_of(pheno))

pred_file <- file.path("output", pheno,
                      paste0("fold", fold),
                      paste0("prop", prop),
                      "pred.tsv")
pred_df <- read_tsv(pred_file) %>%
  filter(type == "lasso" & !(eid %in% min_id)) %>%
  select(-type) %>%
  inner_join(pheno_df, by = "eid")

covar_formula <- y ~ Age + Sex *
                   (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
full_formula <- y ~ Age + pred + Sex *
                   (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)

out_df <- pred_df %>%
  group_by(parent_meaning, pow) %>%
  summarise(r2 = compute_r2(y, pred),
            covar = summary(lm(y ~ Age + Sex *
                   (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)))$r.squared,
            full = summary(lm(y ~ Age + pred + Sex *
                   (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)))$r.squared,
            partial = full - covar,
            bias = mean(pred) - mean(y),
            mse = mean((pred - y)^2),
            cor = cor(pred, y))

scores_file <- file.path("output", pheno,
                      paste0("fold", fold),
                      paste0("prop", prop),
                      "scores.tsv")
write_tsv(out_df, scores_file)
