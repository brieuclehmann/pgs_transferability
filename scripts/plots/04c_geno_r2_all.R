#####################
## Load full files ##
#####################

full_score_df <- tibble()
for (code in full_pheno_codes) {
  for (anc in all_ancestries) {
    this_file <- paste0("output/ukb/v~tagged/pheno~", code, "/min_ancestry~", anc, "/scores.tsv")
    if (file.exists(this_file)) {
      this_df <- read_tsv(this_file, show_col_types = FALSE)
      full_score_df <- bind_rows(full_score_df, this_df)
    }
  }
}
full_score_df <- full_score_df %>%
  filter(pop == min_ancestry | (prop_min == 0 & min_ancestry == "AFR")) %>%
  mutate(trait = phenos[pheno],
         `Training set` = case_when(prop_min == 0 ~ "European ancestry only",
                                    prop_min == 0.1 ~ "Multi-ancestry",
                                    TRUE ~ "Minority ancestry only")) %>%
  mutate(`Training set` = factor(`Training set`, 
                                 levels = c("European ancestry only",
                                            "Multi-ancestry", 
                                            "Minority ancestry only")))

full_min_maj_df <- full_score_df %>%
  filter(prop_min %in% c(0, 1)) %>%
  select(-pow) %>%
  full_join(tibble(pow = seq(0, 1, 0.2)), by = character()) %>%
  filter(!(trait == "FGP" & pop %in% c("AMR", "EAS", "MID")))

##########################
## Prep plotting output ##
##########################

plot_df <- full_score_df %>%
  filter(prop_min == 0.1) %>%
  bind_rows(full_min_maj_df) %>%
 # left_join(eur_df, by = c("pheno", "fold")) %>%
  mutate(r2 = if_else(trait %in% quant_phenos, partial_r2, pseudo_r2)) %>%
  #  mutate(partial_r2 = partial_eur - partial_r2) %>%
  group_by(pow, prop_min, trait, `Test set` = pop, `Training set`) %>%
  summarise(mean = mean(r2), 
            max = max(r2),
            min = min(r2), .groups = "drop") %>%
  mutate(trait = factor(trait, trait_order)) %>%
  filter(`Test set` != "EUR")

p1 <- plot_df %>%
  ggplot(aes(pow, mean)) + 
  geom_line(aes(color = `Training set`),
            position = position_dodge(width = 0*0.1)) +
  geom_point(size = 0.3, aes(color = `Training set`)) +
  # geom_errorbar(aes(ymin = min, ymax = max, color = `Training set`),
  #               position = position_dodge(width = 0.1)) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = `Training set`), alpha = 0.2) +
  facet_grid(trait ~ `Test set`, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.2)) +
  xlab(expression(gamma)) +
  ylab(expression(r^2)) +
  scale_colour_discrete(name = "Polygenic score",
                        labels = c(expression("PGS"["EUR"]),
                                   expression("PGS"["dual"]), 
                                   expression("PGS"["min"]))) +
  scale_fill_discrete(name = "Polygenic score",
                      labels = c(expression("PGS"["EUR"]),
                                 expression("PGS"["dual"]), 
                                 expression("PGS"["min"]))) +
  theme(legend.position = "bottom")

dir.create("plots", showWarnings = FALSE)
ggsave("plots/SI_fig1_geno_r2.pdf", width = 8, height = 6)
