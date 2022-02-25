#######################
## Load output files ##
#######################
decomp_all <- "output/ukb/decomp.tsv"

if (file.exists(decomp_all)) {
  decomp_df <- read_tsv(decomp_all, show_col_types = FALSE) %>%
    mutate(`Training set` = case_when(prop_min == 0 ~ "European ancestry only",
                                      prop_min == 0.1 ~ "Multi-ancestry",
                                      TRUE ~ "Minority ancestry only"))
} else {
  decomp_files <- list.files("output/ukb", "*decomp.tsv", 
                             full.names = TRUE, 
                             recursive = TRUE)
  
  decomp_df <- tibble()
  for (decomp_file in decomp_files) {
    pheno <- gsub(".*pheno~(.*?)/.*", "\\1", decomp_file)
    prop_min <- gsub(".*prop_min~(.*?)/.*", "\\1", decomp_file)
    min_ancestry <- gsub(".*min_ancestry~(.*?)/.*", "\\1", decomp_file)
    f <- gsub(".*fold~(.*?)/.*", "\\1", decomp_file)
    pow <- gsub(".*pow~(.*?)/.*", "\\1", decomp_file)
    
    this_df <- readr::read_tsv(decomp_file, show_col_types = FALSE) %>%
      mutate(pheno = pheno,
             prop_min = prop_min,
             min_ancestry = min_ancestry,
             fold = f,
             pow = pow,
             pow_aux = if_else(prop_min == "0.0", -0.25, as.double(pow))) %>%
      mutate(pow_aux = if_else(prop_min == "1.0", 1.25, pow_aux),
             trait = phenos[pheno]) %>%
      select(pheno, trait, min_ancestry, prop_min, pow_aux, pop, 
             maf_pop, maf_grp, partial_r2, pseudo_r2)
    decomp_df <- bind_rows(decomp_df, this_df)
  }
  
  write_tsv(decomp_df, decomp_all)
}

full_score_df <- tibble()
for (code in full_pheno_codes) {
  for (anc in all_ancestries) {
    this_file <- paste0("output/ukb/pheno~", code, "/min_ancestry~", anc, "/scores.tsv")
    if (file.exists(this_file)) {
      this_df <- read_tsv(this_file, show_col_types = FALSE)
      full_score_df <- bind_rows(full_score_df, this_df)
    }
  }
}
full_score_df <- full_score_df %>%
  mutate(min_ancestry = if_else(prop_min == 0 & pop != "EUR", pop, min_ancestry),
         trait = phenos[pheno],
         pow_aux = if_else(prop_min == 0, -0.25, as.double(pow)),
         `Training set` = case_when(prop_min == 0 ~ "European ancestry only",
                                    prop_min == 0.1 ~ "Multi-ancestry",
                                    TRUE ~ "Minority ancestry only")) %>%
  mutate(pow_aux = if_else(prop_min == 1, 1.25, pow_aux),
         `Training set` = factor(`Training set`, 
                                 levels = c("European ancestry only",
                                            "Multi-ancestry", 
                                            "Minority ancestry only")),
         r2 = if_else(trait %in% quant_phenos, partial_r2, pseudo_r2)) %>%
  group_by(trait, min_ancestry, pop, pow_aux) %>%
  summarise(base_r2 = mean(r2), .groups = "drop")

eur_df <- full_score_df %>% 
  filter(pow_aux == -0.25 & pop == "EUR" & min_ancestry == "AFR") %>%
  select(-min_ancestry) %>%
  full_join(tibble(min_ancestry = c("AFR", "AMR", "EAS", "CSA", "MID")),
            by = character())
  
full_score_df <- full_score_df %>%
  filter(!(pow_aux == -0.25 & pop == "EUR")) %>%
  bind_rows(eur_df)
#######################
## Preprocess output ##
#######################

out_df <- decomp_df %>%
  filter(pop == "EUR" | pop == min_ancestry) %>%
  mutate(r2 = if_else(trait %in% quant_phenos, partial_r2, pseudo_r2),
         trait = factor(trait, trait_order)) %>%
  group_by(trait, min_ancestry, pop, pow_aux, maf_pop, maf_grp) %>%
  summarise(mean = mean(r2), 
            max = max(r2),
            min = min(r2), .groups = "drop") %>%
  left_join(full_score_df, 
            by = c("trait", "min_ancestry", "pop", "pow_aux")) %>%
  mutate(trait = factor(trait, trait_order))

plot_maf <- function(this_maf_pop, this_min_ancestry, this_pop, this_trait) {
  colours <- ifelse(this_maf_pop == "min", "Blues", "Greens")
  maf_label <- ifelse(this_maf_pop == "min", 
                      paste0("MAF (", this_min_ancestry, ")"),
                      "MAF (EUR)")
  
  out_df %>%
    filter(maf_pop == this_maf_pop & pop == this_pop & 
             min_ancestry == this_min_ancestry & trait == this_trait) %>%
    ggplot(aes(pow_aux, mean)) + 
    geom_bar(aes(fill = maf_grp), stat = "identity") +
    geom_point(aes(pow_aux, base_r2)) +
    #  geom_point(aes(y = partial)) +
    #    facet_wrap("trait", scales = "free_y", nrow = 3) +
    scale_fill_brewer(direction = 1, palette = colours,
                      limits=c("a", "b", unique(out_df$maf_grp)), 
                      breaks = unique(out_df$maf_grp)) +
    labs(fill = maf_label) +
    ylab(expression(r^2))  +
    #    xlab("Training set") +
    scale_x_continuous(breaks = c(-0.25, 0, 0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.25), 
                       minor_breaks = c(-0.25, 0, 0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.25),
                       labels = c("EUR-only", 
                                  c(0, 0.2, 0.4, "\n\u03B3\nMulti-ancestry\n", 0.6, 0.8, 1), 
                                  paste0(this_min_ancestry, "-only"))) +
    theme(legend.position="bottom",
          axis.title.y = element_text(angle = 0, vjust=0.5),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_line(color = c(rep("gray92", 4), 
                                                      NA, 
                                                      rep("gray92", 4))))
}


plot_maf_wrapper <- function(min_ancestry, pheno1, pheno2) {
  p1 <- plot_maf("min", min_ancestry, min_ancestry, pheno1)
  p2 <- plot_maf("maj", min_ancestry, min_ancestry, pheno1)
  
  ylims <- c(0, max(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,
                    ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range))
  
  p1 <- p1 + ylim(ylims) + theme(legend.position = "none") + ggtitle(pheno1)
  p2 <- p2 + ylim(ylims) + theme(legend.position = "none")
  
  p3 <- plot_maf("min", min_ancestry, "EUR", pheno1)
  p4 <- plot_maf("maj", min_ancestry, "EUR", pheno1)
  ylims <- c(0, max(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range,
                    ggplot_build(p4)$layout$panel_scales_y[[1]]$range$range))
  
  p3 <- p3 + ylim(ylims) + theme(legend.position = "none")
  p4 <- p4 + ylim(ylims) + theme(legend.position = "none")
  
  q1 <- plot_maf("min", min_ancestry, min_ancestry, pheno2) + ggtitle(pheno2)
  q2 <- plot_maf("maj", min_ancestry, min_ancestry, pheno2)
  ylims <- c(0, max(ggplot_build(q1)$layout$panel_scales_y[[1]]$range$range,
                    ggplot_build(q2)$layout$panel_scales_y[[1]]$range$range))
  
  q1 <- q1 + ylim(ylims) + theme(legend.position = "none")
  q2 <- q2 + ylim(ylims) + theme(legend.position = "none")
  
  q3 <- plot_maf("min", min_ancestry, "EUR", pheno2)
  q4 <- plot_maf("maj", min_ancestry, "EUR", pheno2)
  ylims <- c(0, max(ggplot_build(q3)$layout$panel_scales_y[[1]]$range$range,
                    ggplot_build(q4)$layout$panel_scales_y[[1]]$range$range))
  
  q3 <- q3 + ylim(ylims)
  q4 <- q4 + ylim(ylims)
  
  # for (p in c("p1", "p2", "p3", "p4", "q1", "q2", "q3", "q4")) {
  #   gt = ggplot_gtable(ggplot_build(get(p)))
  #   gt$widths[7] <- gt$widths[7] * 0.5
  #   assign(p, wrap_ggplot_grob(gt))
  # }
  
  p_top <- (p1 + q1) / (p2 + q2) + 
    plot_annotation(title = paste0(min_ancestry, " test set"), 
                    theme = theme(plot.title = element_text(hjust = 0.5)))
  p_bottom <- (p3 + q3) / (p4 + q4) / guide_area() + 
    plot_layout(guides = 'collect', heights = c(2, 2, 1)) +
    plot_annotation(title = "EUR test set", 
                    theme = theme(plot.title = element_text(hjust = 0.5)))
  
  out_file <- paste0("plots/SI_maf_upper_", min_ancestry, 
                     "_", pheno1, "-", pheno2, ".pdf")
  ggsave(out_file, p_top, device = cairo_pdf, width = 8, height = 4)
  
  out_file <- paste0("plots/SI_maf_lower_", min_ancestry, 
                     "_", pheno1, "-", pheno2, ".pdf")
  ggsave(out_file, p_bottom, device = cairo_pdf, width = 8, height = 4)
}

for (min_ancestry in c("EAS", "CSA", "MID", "AMR", "AFR")) {
  for (pheno_pair in list(c("height", "MCV"), c("asthma", "FGP"))) {
    plot_maf_wrapper(min_ancestry, pheno_pair[[1]], pheno_pair[[2]])
  }
}
