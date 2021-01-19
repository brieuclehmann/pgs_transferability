# Scripts to construct training sets with varying numbers of White and Black/Black British individuals

pheno <- "height"
min_coding <- 4 # Black or Black British
maj_coding <- 1 # White
nfolds <- 5
nfracs <- 5
frac_range <- seq(0, 1, 0.2)

prop_min <- 0.1
f <- 1
frac_type <- "frac_min" # one of frac_min or frac_maj to vary fraction of minority (Black or Black British) or majority (White)

library(Matching)
library(readr)
library(dplyr)
set.seed(1) # for reproducibility

pheno_df <- read_tsv("data/all_vars.tsv") %>%
  rename(y = all_of(pheno)) %>%
  filter(!is.na(y) & parent_coding %in% c(min_coding, maj_coding)) %>%
  mutate(Tr = (parent_coding == min_coding))

# Train-test split
train_ratio <- round((1 - prop_min) / prop_min)

min_df <- pheno_df %>%
  filter(parent_coding == min_coding) %>%
  mutate(fold = sample(cut(seq_len(n()), breaks = nfolds, labels = FALSE))) %>%
  group_by(fold) %>%
  mutate(frac = sample(cut(seq_len(n()), breaks = nfracs, labels = FALSE))) %>%
  dplyr::select(eid, Tr, Sex, Age, fold, frac)

maj_df <- pheno_df %>%
  filter(parent_coding == maj_coding) %>%
  mutate(fold = NA_integer_, frac = NA_integer_) %>%
  dplyr::select(eid, Tr, Sex, Age, fold, frac)

match_df <- min_df %>%
  filter(fold != f) %>%
  bind_rows(maj_df)

match_out <- match_df %>%
  with(Match(Tr = Tr, X = cbind(Sex, Age), exact = c(TRUE, FALSE),
             M = train_ratio, Z = eid, ties = FALSE, replace = FALSE))

idx <- unique(c(match_out$index.treated, match_out$index.control))
train_id <- match_df$eid[idx]

maj_train_df <- maj_df %>%
  filter(eid %in% train_id) %>%
  mutate(frac = sample(cut(seq_len(n()), breaks = nfracs, labels = FALSE)))

for (frac_ind in seq_along(frac_range)) {
  this_frac <- frac_range[frac_ind]
  outdir <- file.path("data",
                       pheno, paste0("fold", f),
                       paste0("prop", prop_min))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(outdir, 
                       paste0(frac_type, format(this_frac, nsmall = 1), ".txt"))

  if (frac_type == "frac_min") {

    min_df %>%
      filter(eid %in% train_id & frac < frac_ind) %>%
      bind_rows(maj_train_df) %>%
      pull(eid) %>%
      rep(each = 2) %>%
      write(outfile, ncol = 2)

  } else if (frac_type == "frac_maj") {

    maj_train_df %>%
      filter(frac < frac_ind) %>%
      bind_rows(min_df) %>%
      filter(eid %in% train_id) %>%
      pull(eid) %>%
      rep(each = 2) %>%
      write(outfile, ncol = 2) 

  }

}
