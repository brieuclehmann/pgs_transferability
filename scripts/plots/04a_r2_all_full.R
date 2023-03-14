######################
## Load basic files ##
######################

basic_files <- paste0("output/ukb/v~imputed/pheno~", full_pheno_codes, "/min_ancestry~AFR/scores.tsv")
basic_df <- map_dfr(basic_files, read_tsv, show_col_types = FALSE) %>%
  mutate(trait = phenos[pheno],
         `Training set` = case_when(prop_min == 0 ~ "European ancestry only",
                                    prop_min == 0.1 ~ "Multi-ancestry (reduced)",
                                    prop_min == -1 ~ "Multi-ancestry (full)",
                                    TRUE ~ "African ancestry only")) %>%
  mutate(`Training set` = factor(`Training set`, 
                                 levels = c("European ancestry only",
                                            "Multi-ancestry (reduced)",
                                            "Multi-ancestry (full)", 
                                            "African ancestry only")))

min_maj_df <- basic_df %>%
  filter(pop %in% c("AFR", "EUR")  & prop_min %in% c(0, 1)) %>%
  select(-pow) %>%
  full_join(tibble(pow = seq(0, 1, 0.2)), by = character())

eur_df <- min_maj_df %>%
  filter(pop == "EUR" & prop_min == 0) %>%
  mutate(r2_eur = if_else(trait %in% quant_phenos, partial_r2, pseudo_r2)) %>%
  distinct(pheno, fold, r2_eur)

pheno_order <- eur_df %>%
  group_by(pheno) %>%
  summarise(r2 = mean(r2_eur)) %>%
  arrange(desc(r2)) %>%
  pull(pheno)
trait_order <- phenos[pheno_order]

##########################
## Prep plotting output ##
##########################

plot_df <- basic_df2 %>%
  filter(pop %in% c("AFR", "EUR") & prop_min %in% c(0.1, -1) & pow <= 1.0) %>%
  bind_rows(min_maj_df) %>%
  left_join(eur_df, by = c("pheno", "fold")) %>%
  mutate(r2 = if_else(trait %in% quant_phenos, partial_r2, pseudo_r2)) %>%
  group_by(pow, prop_min, trait, `Test set` = pop, `Training set`) %>%
  summarise(mean = mean(r2), 
            max = max(r2),
            min = min(r2), .groups = "drop") %>%
  mutate(trait = factor(trait, trait_order))

p1 <- ggplot(filter(plot_df, `Test set` == "AFR"), aes(pow, mean)) + 
  geom_line(aes(color = `Training set`, linetype = `Test set`)) +
  geom_point(size = 0.3, aes(color = `Training set`)) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = `Training set`), alpha = 0.2) +
  geom_line(data = filter(plot_df, `Test set` == "EUR" & prop_min == 0), 
            aes(pow, mean, linetype = `Test set`),
            inherit.aes = FALSE) +
  geom_ribbon(data = filter(plot_df, `Test set` == "EUR" & prop_min == 0),
              aes(pow, mean, ymin = min, ymax = max), alpha = 0.2,
              inherit.aes = FALSE) +
  facet_wrap(c("trait"), nrow = 3, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.2)) +
  scale_colour_discrete(name = "Polygenic score",
                        labels = c(expression("PGS"["EUR"]),
                                   expression("PGS"["dual (reduced)"]), 
                                   expression("PGS"["dual (full)"]),
                                   expression("PGS"["AFR"]))) +
  scale_fill_discrete(name = "Polygenic score",
                      labels = c(expression("PGS"["EUR"]),
                                 expression("PGS"["dual (reduced)"]), 
                                 expression("PGS"["dual (full)"]),
                                 expression("PGS"["AFR"]))) +
  xlab(expression(gamma)) +
  ylab(expression(r^2)) +
  theme(legend.position = "bottom")

dir.create("plots", showWarnings = FALSE)
ggsave("plots/SI_fig4_afr_r2_full.pdf", width = 8, height = 6)
