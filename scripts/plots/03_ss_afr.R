######################
## Load basic files ##
######################

basic_files <- paste0("output/ukb/pheno~", 
                      pheno_codes, 
                      "/min_ancestry~AFR/samplesize_scores.tsv")
basic_df <- map_dfr(basic_files, read_tsv, show_col_types = FALSE) %>%
  mutate(trait = phenos[pheno],
         r2 = if_else(trait %in% quant_phenos, partial_r2, pseudo_r2),
         trait = factor(trait, trait_order)) %>%
  left_join(ntrain_df, by = c("pheno", "fold", "min_ancestry")) %>%
  mutate(n_min = if_else(type == "min", frac * n_min, n_min),
         n_eur = if_else(type == "maj", frac * n_eur, n_eur))

out_df <- basic_df %>% 
  filter(pop %in% c("AFR", "EUR")) %>%
  select(trait, pop, frac, type, fold, r2, n_min, n_eur) %>%
  group_by(trait, pop, frac, type) %>%
  summarise(out_up = max(r2), 
            out_low = min(r2), 
            out = mean(r2), 
            n_min = round(mean(n_min)),
            n_eur = round(mean(n_eur)),
            .groups = "drop")
  # inner_join(filter(ntrain_df, min_coding == 4), by = c("pheno", "prop")) %>%
  # mutate(nblack = ntrain * 0.1,
  #        nwhite = (ntrain - nblack) * frac,
  #        prop_white = nwhite / (nwhite + nblack)) %>%
p1 <- out_df %>%
  filter(type == "min") %>%
  ggplot(aes(n_min, out)) +
  geom_line(aes(color = trait, group = trait)) +
  geom_ribbon(aes(ymin = out_low, ymax = out_up, fill = trait), alpha = 0.2) +
  geom_point(aes(color = trait, group = trait), size = 0.5) +
  facet_wrap(c("pop"), nrow = 2) +
  ylab(expression(r^2)) +
  xlab("Number of AFR in training set") +
#  scale_color_viridis_d(option = "A") +
#  scale_fill_viridis_d(option = "A") +
  theme(legend.position = "bottom") +
  

p2 <- out_df %>%
  filter(type == "maj") %>%
  ggplot(aes(n_eur, out)) +
  geom_line(aes(color = trait, group = trait)) +
  geom_ribbon(aes(ymin = out_low, ymax = out_up, fill = trait), alpha = 0.2) +
  geom_point(aes(color = trait, group = trait), size = 0.5) +
  facet_wrap(c("pop"), nrow = 2) +
  ylab(expression(r^2)) +
  xlab("Number of EUR in training set") +
  theme(legend.position = "bottom")

pout <- p1 + p2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom') &
  ggsci::scale_color_d3(palette = "category20") &
  ggsci::scale_fill_d3(palette = "category20")

dir.create("plots", showWarnings = FALSE)
ggsave("plots/fig2_afr_ss.pdf", width = 8, height = 6)
