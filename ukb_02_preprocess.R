# Script to preprocess phenotypes

library(readr)
library(dplyr)
library(tidyr)

### Blood measures ---

# 30000 - WBC
# 30010 - RBC
# 30020 - haemoglobin
# 30030 - haematocrit
# 30040 - MCV
# 30050 - MCH
# 30080 - platelet
# 30120 - lymphocyte
# 30130 - monocyte
# 30710 - C reactive protein

blood_df <- read_tsv("data/blood_phenotypes_raw.tsv")
names(blood_df) <- c("eid", "WBC", "RBC", "haemoglobin", "haematocrit",
                     "MCV", "MCH", "platelet", "lymphocyte", 
                     "monocyte", "CRP")

### Body measures ---

body_df <- read_tsv("data/body_phenotypes_raw.tsv")

names(body_df) <- c("eid", "waist_circum", "height",
                    "sys_man1", "sys_man2", "dia_man1", "dia_man2",
                    "dia_auto1", "dia_auto2", "sys_auto1",
                    "sys_auto2", "BMI")

body_df <- body_df %>%
  mutate(systolic = rowMeans(select(body_df, starts_with("sys")),
                             na.rm = TRUE),
         diastolic = rowMeans(select(body_df, starts_with("dia")),
                             na.rm = TRUE)) %>%
  select(eid, waist_circum, height, BMI,
         systolic, diastolic)

### Covariates ---

covariate_df <- read_tsv("data/covariates_raw.tsv")
n_pc <- covariate_df %>%
  select(starts_with("22009")) %>%
  ncol
codings_df  <- read_tsv("Codings/coding1001.tsv")

names(covariate_df) <- c("eid", "Sex", "YearOfBirth", "MonthOfBirth",
                         "coding", "Batch", paste0("PC", seq_len(n_pc)))

covariate_df <- covariate_df %>%
  mutate_at(vars(starts_with("PC")),
            .funs = list(Sex = ~ . * covariate_df$Sex)) %>%
  inner_join(codings_df, by = "coding") %>%
  rename(main_coding = coding, main_meaning = meaning) %>%
  mutate(parent_coding = if_else(parent_id != 0, parent_id, node_id),
         Batch = if_else(Batch < 0 | Batch == 1000, "BiLEVE", "Axiom"),
         Age = 2020 - YearOfBirth - (MonthOfBirth - 1) / 12)  %>%
  select(eid, Sex, Age, Batch, starts_with("PC"),
         main_coding, main_meaning, parent_coding) %>%
  inner_join(select(codings_df, coding, meaning),
             by = c("parent_coding" = "coding")) %>%
  rename(parent_meaning =  meaning)

### QC on samples

fam <- read.table("data/ukb12788_cal_chr1_v2_s488264.fam")
qc <- read_delim("ukb/v2/qc/ukb_sqc_v2.txt", delim = " ")
qc$eid <- fam$V1

qc_keep <- qc %>%
  filter(used.in.pca.calculation == 1) %>%
  filter(putative.sex.chromosome.aneuploidy == 0) %>%
  pull(eid)

pheno_df <- covariate_df %>%
  filter(eid %in% qc_keep) %>%
  left_join(blood_df,  by = "eid") %>%
  left_join(body_df, by = "eid") %>%
  mutate(FID = eid, IID = eid)

write_tsv(pheno_df, "data/all_vars.tsv")
 