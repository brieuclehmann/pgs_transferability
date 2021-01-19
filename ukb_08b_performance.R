# Script to evaluate predictive accuracy of PGS based on subsets of SNPs

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

pred_split_file <- file.path("output", pheno,
                             paste0("fold", fold),
                             paste0("prop", prop),
                             "pred_split.tsv")
pred_split_df <- read_tsv(pred_split_file) %>%
   filter(!(eid %in% min_id)) %>%
   rename(part = pred)

out_df <- pred_df %>%
  group_by(parent_meaning, pow) %>%
  summarise(r2 = compute_r2(y, pred),
            covar = summary(lm(y ~ Age + Sex *
                   (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)))$r.squared,
            full = summary(lm(y ~ Age + pred + Sex *
                   (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)))$r.squared,
            semi_partial = full - covar,
            partial = semi_partial / (1 - covar),
            bias = mean(pred) - mean(y),
            mse = mean((pred - y)^2),
            cor = cor(pred, y))

part_df <- pred_df %>%
  inner_join(pred_split_df, by = c("eid", "pow")) %>%
  group_by(parent_meaning, pow, ethnicity, grp) %>%
  summarise(full = summary(lm(y ~ Age + part + Sex *
                   (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)))$r.squared) %>%
  inner_join(select(out_df, parent_meaning, pow, covar, partial),
             by = c("parent_meaning", "pow")) %>%
  mutate(semi_partial_split = full - covar,
         partial_split = semi_partial_split / (1 - covar))

scores_file <- file.path("output", pheno,
                      paste0("fold", fold),
                      paste0("prop", prop),
                      "scores_split.tsv")
write_tsv(part_df, scores_file)
