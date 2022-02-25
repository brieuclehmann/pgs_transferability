library(readr)
library(dplyr)
library(ggplot2)
library(purrr)

pheno_codes <- c("A50", "A21001", "A30100", "A30040", "A30090")
phenos <- c("height", "BMI", "MPV", "MCV", "platelet")
names(phenos) <- pheno_codes

score_files <- paste0("output/ukb/pheno~", pheno_codes, "/min_ancestry~AFR/scores.tsv")
score_df <- map_dfr(score_files, read_tsv) %>%
  mutate(trait = phenos[pheno])

beta_files <- paste0("output/ukb/pheno~", pheno_codes, "/min_ancestry~AFR/beta.tsv")
beta_df <- map_dfr(beta_files, read_tsv) %>%
  mutate(trait = phenos[pheno])

height_df <- read_tsv("output/ukb/pheno~A50/min_ancestry~AFR/scores.tsv")
boot_df <- read_tsv("output/ukb/pheno~A50/min_ancestry~AFR/boot_scores.tsv")

phe <- "height"
min_maj_df <- score_df %>%
  filter(pop %in% c("AFR", "EUR")  & prop_min %in% c(0, 1)) %>%
  select(-pow) %>%
  full_join(tibble(pow = seq(0, 1, 0.2)), by = character())

p2 <- score_df %>%
  filter(pop %in% c("AFR", "EUR") & prop_min == 0.1) %>%
  bind_rows(min_maj_df) %>%
  group_by(pow, prop_min, trait, pop) %>%
  summarise(mean = mean(partial), 
            max = max(partial),
            min = min(partial), .groups = "drop") %>%
  #  left_join(train_name_df, by = "prop") %>%
  #  mutate(Training = fct_relevel(Training, train_name_df$Training)) %>%
  #  mutate(Training = paste0(Training, "\n")) %>%
  ggplot(aes(pow, mean)) + 
  geom_line(aes(color = factor(prop_min), linetype = factor(prop_min))) +
  facet_wrap(c("pop", "trait"), nrow = 2) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.2)) +
  #  scale_linetype_manual(values = c(2, 2, 1)) +
  #  scale_color_manual(values = gg_color_hue(4)[c(4,1,3)]) +
  xlab(expression(gamma)) +
  ylab(expression(r^2)) 

fold1_df <- boot_df %>%
  filter(bootstrap & pop %in% c("AFR", "EUR") & fold == 1) %>%
  select(pow, prop_min, pop, mean, upper, lower) %>%
  mutate(type = "bootstrap")

min_maj_df2 <- fold1_df %>%
  filter(pop %in% c("AFR", "EUR")  & prop_min %in% c(0, 1)) %>%
  select(-pow) %>%
  full_join(tibble(pow = seq(0, 1, 0.2)), by = character())

score_df %>%
  filter(pop %in% c("AFR", "EUR") & prop_min == 0.1) %>%
  bind_rows(min_maj_df) %>%
  filter(trait == "height") %>%
  group_by(pow, prop_min, pop) %>%
  summarise(mean = mean(partial), 
            upper = max(partial),
            lower = min(partial), .groups = "drop") %>%
  mutate(type = "cv") %>%
  bind_rows(fold1_df) %>%
  bind_rows(min_maj_df2) %>%
  ggplot(aes(pow, mean)) + 
  geom_line(aes(color = factor(prop_min), linetype = factor(prop_min))) +
  facet_wrap(c("pop", "type"), nrow = 2) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.2)) +
  #  scale_linetype_manual(values = c(2, 2, 1)) +
  #  scale_color_manual(values = gg_color_hue(4)[c(4,1,3)]) +
  xlab(expression(gamma)) +
  ylab(expression(r^2))
  


