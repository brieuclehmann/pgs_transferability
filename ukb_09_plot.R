library(foreach)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(theme_minimal(base_size = 12))

### LOAD DATA
phenos <- c("systolic", "CRP", "height", "BMI", "platelet",
            "MCV", "WBC", "monocyte", "diastolic", "RBC", 
            "haematocrit", "MCH", "lymphocyte", "haemoglobin", "waist_circum")
props <- c(0, 0.1, 1)
folds <- seq_len(5)
opts <- expand.grid(phenos, props, folds)

ntrain_df <- tibble(pheno = character(), prop = double(),
                    fold = integer(), ntrain = integer())
for (pheno in phenos) {
    for (prop in props) {
        for (fold in seq_len(nfolds)) {
            kfile  <- file.path("data", pheno,
                                paste0("fold", fold),
                                paste0("prop", prop, ".txt"))
            if (file.exists(kfile)) {
                ntrain <- length(read.table(kfile)[ ,1])
                ntrain_df <- ntrain_df %>%
                add_row(pheno, prop, fold, ntrain)
            }
        }
    }
}

power_min_df <- foreach(phe = opts$Var1, prop = opts$Var2, fold = opts$Var3,
                        .combine = bind_rows) %do% {
                          out <- NULL
                          out_file <- file.path("output", phe, 
                                                paste0("fold", fold), 
                                                paste0("prop", prop), "scores_frac_min.tsv")
                          
                          if (file.exists(out_file)) {
                            out <- read_tsv(out_file) %>%
                              mutate(prop = prop, pheno = phe, fold = fold)
                          }
                        }

power_maj_df <- foreach(phe = opts$Var1, prop = opts$Var2, fold = opts$Var3,
                        .combine = bind_rows) %do% {
                          out <- NULL
                          out_file <- file.path("output", phe, 
                                                paste0("fold", fold), 
                                                paste0("prop", prop), "scores_frac_maj.tsv")
                          
                          if (file.exists(out_file)) {
                            out <- read_tsv(out_file) %>%
                              mutate(prop = prop, pheno = phe, fold = fold)
                          }
                        }

out_df <- foreach(phe = opts$Var1, prop = opts$Var2, fold = opts$Var3,
                  .combine = bind_rows) %do% {
                    out <- NULL
                    out_file <- file.path("output", phe, paste0("fold", fold), 
                                          paste0("prop", prop), "scores.tsv")
                    
                    if (file.exists(out_file)) {
                      out <- read_tsv(out_file) %>%
                        mutate(prop = prop, pheno = phe, fold = fold) 
                    }
                    
                    out
                  }

split_df <- foreach(phe = opts$Var1, prop = opts$Var2, fold = opts$Var3,
                    .combine = bind_rows) %do% {
                    out <- NULL
                    out_file <- file.path("output", phe, paste0("fold", fold), 
                                          paste0("prop", prop), "scores_split.tsv")
                    
                    if (file.exists(out_file)) {
                      out <- read_tsv(out_file) %>%
                        mutate(prop = prop, pheno = phe, fold = fold) %>%
                        filter(pow <= 1)
                    }
                  
                    out
                  }


# Get order of phenotypes

df <- out_df %>%
  filter(parent_meaning == "Black or Black British" & type == "aall") %>%
  group_by(pow, prop, pheno) %>%
  summarise(mean = mean(partial), 
            max = max(partial),
            min = min(partial), .groups = "drop")

order_df <- df %>%
  filter(pow == 0) %>%
  group_by(pheno) %>%
  filter(mean == max(mean)) %>%
  arrange(desc(prop), desc(mean))

### FIGURE 1 ###

p2 <- power_df %>% 
  filter(prop == 0.1 &
         parent_meaning %in% c("White", "Black or Black British")) %>%
  select(pheno, parent_meaning, frac, fold, prop, partial) %>%
  group_by(parent_meaning, pheno, frac, prop) %>%
  summarise(out_up = max(partial), 
            out_low = min(partial), 
            out = mean(partial), .groups = "drop") %>%
  inner_join(filter(ntrain_df, min_coding == 4), by = c("pheno", "prop")) %>%
  mutate(nwhite = ntrain * 0.9,
         nblack = (ntrain - nwhite) * frac,
         prop_white = nwhite / (nwhite + nblack)) %>%
  mutate(pheno = factor(pheno, levels = order_df$pheno)) %>% 
  rename(Phenotype = pheno) %>%
  ggplot(aes(nblack, out)) +
  geom_line(aes(color = Phenotype, group = Phenotype)) +
  geom_point(aes(color = Phenotype, group = Phenotype), size = 0.5) +
  facet_wrap("parent_meaning", nrow = 1) +
  xlab("Number of Black individuals in training set") +
  ylab(expression(r^2)) +
  theme(legend.position = "top")

p1 <- power_maj_df %>% 
  filter(parent_meaning %in% c("White", "Black or Black British") & 
         prop == 0.1) %>%
  select(pheno, parent_meaning, frac, fold, prop, partial) %>%
  group_by(parent_meaning, pheno, frac, prop) %>%
  summarise(out_up = max(partial), 
            out_low = min(partial), 
            out = mean(partial), .groups = "drop") %>%
  inner_join(filter(ntrain_df, min_coding == 4), by = c("pheno", "prop")) %>%
  mutate(nblack = ntrain * 0.1,
         nwhite = (ntrain - nblack) * frac,
         prop_white = nwhite / (nwhite + nblack)) %>%
  mutate(pheno = factor(pheno, levels = order_df$pheno)) %>% 
  rename(Phenotype = pheno) %>%
  ggplot(aes(nwhite, out)) +
  geom_line(aes(color = Phenotype, group = Phenotype)) +
  geom_point(aes(color = Phenotype, group = Phenotype), size = 0.5) +
  facet_wrap("parent_meaning", nrow = 1) +
  xlab("Number of White individuals in training set") +
  ylab(expression(r^2))+
  theme(legend.position = "top")

(p1 / p2) +
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') &
  theme(legend.position='bottom', plot.tag = element_text(size = 12),
        axis.title.y = element_text(angle = 0, vjust=0.5)) 

ggsave("plots/fig1.png", width = 8, height = 6)


### Figure 3 ###

phe <- "MCV"
min_maj_df <- out_df %>%
  filter(parent_meaning %in% c("Black or Black British", "White") & 
         pheno == phe & prop %in% c(0, 1) & pow == 0) %>%
  select(-pow) %>%
  full_join(tibble(pow = seq(0, 1, 0.2)), by = character())

train_name_df <- ntrain_df %>%
  filter(pheno == phe) %>%
  inner_join(tibble(prop = c(1, 0.1, 0), 
                    train_lab = c("prop_black = 100%",
                                  "prop_black = 10%\nprop_white = 90%",
                                  "prop_white = 100%")), by = "prop") %>%
  group_by(prop, train_lab) %>%
  summarise(ntrain = round(mean(ntrain))) %>%
  mutate(n_lab = paste0("n = ", ntrain),
         Training = paste(n_lab, train_lab, sep = "\n")) %>%
  select(prop, Training)

p1 <- out_df %>%
  filter(parent_meaning %in% c("Black or Black British", "White") & 
           type == "aall" & pheno == phe & prop == 0.1) %>%
  bind_rows(min_maj_df) %>%
  group_by(pow, prop, pheno, parent_meaning) %>%
  summarise(mean = mean(partial), 
            max = max(partial),
            min = min(partial), .groups = "drop") %>%
  left_join(train_name_df, by = "prop") %>%
  mutate(Training = fct_relevel(Training, train_name_df$Training)) %>%
  mutate(Training = paste0(Training, "\n")) %>%
  ggplot(aes(pow, mean)) + 
  geom_line(aes(color = Training, linetype = Training)) +
  facet_wrap(c("parent_meaning"), nrow = 1) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.2)) +
  scale_linetype_manual(values = c(2, 2, 1)) +
  scale_color_manual(values = gg_color_hue(4)[c(4,1,3)]) +
  xlab(expression(gamma)) +
  ylab(expression(r^2)) +
  ggtitle(phe)


### height

phe <- "height"
min_maj_df <- out_df %>%
  filter(parent_meaning %in% c("Black or Black British", "White") & 
         pheno == phe & prop %in% c(0, 1) & pow == 0) %>%
  select(-pow) %>%
  full_join(tibble(pow = seq(0, 1, 0.2)), by = character())

train_name_df <- ntrain_df %>% 
  filter(pheno == phe) %>%
  inner_join(tibble(prop = c(1, 0.1, 0), 
                    train_lab = c("prop_black = 100%",
                                  "prop_black = 10%\nprop_white = 90%",
                                  "prop_white = 100%")), by = "prop") %>%
  group_by(prop, train_lab) %>%
  summarise(ntrain = round(mean(ntrain))) %>%
  mutate(n_lab = paste0("n = ", ntrain),
         Training = paste(n_lab, train_lab, sep = "\n")) %>%
  select(prop, Training)

p2 <- out_df %>%
  filter(parent_meaning %in% c("Black or Black British", "White") & 
           type == "aall" & pheno == phe & prop == 0.1) %>%
  bind_rows(min_maj_df) %>%
  group_by(pow, prop, pheno, parent_meaning) %>%
  summarise(mean = mean(partial), 
            max = max(partial),
            min = min(partial), .groups = "drop") %>%
  left_join(train_name_df, by = "prop") %>%
  mutate(Training = fct_relevel(Training, train_name_df$Training)) %>%
  mutate(Training = paste0(Training, "\n")) %>%
  ggplot(aes(pow, mean)) + 
  geom_line(aes(color = Training, linetype = Training)) +
  facet_wrap(c("parent_meaning"), nrow = 1) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.2)) +
  scale_linetype_manual(values = c(2, 2, 1)) +
  scale_color_manual(values = gg_color_hue(4)[c(4,1,3)]) +
  xlab(expression(gamma)) +
  ylab(expression(r^2)) +
  ggtitle(phe)

(p1 / p2) + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12),
        axis.title.y = element_text(angle = 0, vjust=0.5)) 

ggsave("plots/fig3.png", width = 8, height = 5)



### Figure 4 ###

train_name_df <- tibble(prop = c(1, 0.1, 0),
                        Training = c("prop_black = 100%\nn \u2248 5k\n",
                                     "prop_black = 10%\nprop_white = 90%,\nn \u2248 53k\n",
                                     "prop_white = 100%\nn \u2248 300k\n"))

rewt_df <- out_df %>%
  filter(parent_meaning == "Black or Black British" & prop == 0.1) %>%
  group_by(pheno, pow) %>%
  mutate(mean = mean(partial)) %>%
  group_by(pheno) %>%
  filter(mean == max(mean)) %>%
  group_by(pheno, pow) %>%
  summarise(mean = mean(partial), 
            max = max(partial),
            min = min(partial), .groups = "drop") %>%
  mutate(Training = "prop_black = 10%\nprop_white = 90%,\nn \u2248 53k (reweighted)\n")

### Summary plot
plot_df <- df %>%
  filter(pow == 0) %>%
  bind_rows(rewt_df) %>%
  mutate(pheno = factor(pheno, levels = order_df$pheno)) 

xlabels <- plot_df %>%
  filter(Training == "prop_black = 10%\nprop_white = 90%,\nn \u2248 53k (reweighted)\n") %>%
  mutate(label = paste0(pheno, " (\u03b3 = ", pow, ")")) %>%
  arrange(pheno) %>%
  pull(label)

plot_df %>%
  mutate(Training = factor(Training, levels = unique(Training)[c(3,2,4,1)])) %>%
  ggplot(aes(pheno, mean)) +
  geom_pointrange(aes(ymin = min, ymax = max, colour = Training), 
                  fatten = 1, position = position_dodge2(width = 0.5)) +
  geom_vline(xintercept = 5.5, linetype = 2, size = 0.5, color =  "gray") +
  theme(axis.text.x = element_text(angle = -30, hjust = 0),
        axis.title.y = element_text(angle = 0, vjust=0.5)) +
  scale_x_discrete(labels = xlabels) +
  xlab("") +
  ylab(expression(r^2))

ggsave("plots/fig4.png", width = 8, height = 4)


### Figure 5 ###

phe <- "MCV"
train_name_df <- ntrain_df %>% 
  filter(pheno == phe) %>%
  inner_join(tibble(prop = c(1, 0.1, 0), 
                    train_lab = c("prop_black:\n100%",
                                  "prop_black: 10%\nprop_white: 90%",
                                  "prop_white:\n100%")), by = "prop") %>%
  mutate(n_lab = paste0("n: ", ntrain),
         Training = paste(n_lab, train_lab, sep = "\n"))

sum_df <- split_df %>%
  mutate(semi_partial = partial * (1 - covar)) %>%
  filter(parent_meaning %in% c("Black or Black British", "White") & 
         pheno == phe) %>%
  filter(pow == 0 | prop == 0.1) %>%
  group_by(pow, ethnicity, grp, prop, parent_meaning) %>%
  summarise(partial = mean(semi_partial),
            mean = mean(semi_partial_split), 
            max = max(partial_split),
            min = min(partial_split), .groups = "drop") %>%
  left_join(train_name_df, by = "prop") %>%
  mutate(Training = paste0(Training, "\n"))

p1 <- sum_df %>% 
  filter(ethnicity == "Black" & parent_meaning == "Black or Black British") %>%
  mutate(Training = factor(Training, levels = sort(unique(Training))[c(1,3,2)])) %>%
  ggplot(aes(factor(pow), mean)) + 
  geom_bar(aes(fill = grp), stat = "identity") +
  geom_point(aes(y = partial)) +
  facet_grid( ~ Training, scales = "free_x", space = "free",
             labeller = labeller(parent_meaning = label_wrap_gen(width = 14))) +
  scale_fill_brewer(direction = 1, palette = "Blues",
                    limits=c("a", "b", unique(sum_df$grp)), 
                    breaks=unique(sum_df$grp)) +
  labs(fill = "MAF (Black)") +
  ylab(expression(r^2))  +
  theme(legend.position="bottom", axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(angle = 0, vjust=0.5),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Mean corpuscular volume")

p2 <- sum_df %>% 
  filter(ethnicity == "White" & parent_meaning == "Black or Black British") %>%
  mutate(Training = factor(Training, levels = sort(unique(Training))[c(1,3,2)])) %>%
  ggplot(aes(factor(pow), mean)) + 
  geom_bar(aes(fill = grp), stat = "identity") +
  geom_point(aes(y = partial)) +
  facet_grid( ~ Training, scales = "free_x", space = "free") +
  scale_fill_brewer(direction = 1, palette = "Greens",
                    limits=c("a", "b", unique(sum_df$grp)), 
                    breaks=unique(sum_df$grp)) +
  labs(fill = "MAF (White)") +
  xlab(expression(gamma)) +
  ylab(expression(r^2))  +
  theme(legend.position="bottom",
        strip.text.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust=0.5))

ylims <- c(0, max(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,
                  ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range))

p1 <- p1 + ylim(ylims) + theme(legend.position = "none")
p2 <- p2 + ylim(ylims) + theme(legend.position = "none")

p3 <- sum_df %>% 
  filter(ethnicity == "Black" & parent_meaning == "White") %>%
  mutate(Training = factor(Training, levels = sort(unique(Training))[c(1,3,2)])) %>%
  ggplot(aes(factor(pow), mean)) + 
  geom_bar(aes(fill = grp), stat = "identity") +
  geom_point(aes(y = partial)) +
  facet_grid( ~ Training, scales = "free_x", space = "free",
              labeller = labeller(parent_meaning = label_wrap_gen(width = 14))) +
  scale_fill_brewer(direction = 1, palette = "Blues",
                    limits=c("a", "b", unique(sum_df$grp)), 
                    breaks=unique(sum_df$grp)) +
  labs(fill = "MAF (Black)") +
  ylab(expression(r^2))  +
  theme(legend.position="bottom", axis.title.x=element_blank(),
        strip.text.x = element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(angle = 0, vjust=0.5))

p4 <- sum_df %>% 
  filter(ethnicity == "White" & parent_meaning == "White") %>%
  mutate(Training = factor(Training, levels = sort(unique(Training))[c(1,3,2)])) %>%
  ggplot(aes(factor(pow), mean)) + 
  geom_bar(aes(fill = grp), stat = "identity") +
  geom_point(aes(y = partial)) +
  facet_grid( ~ Training, scales = "free_x", space = "free") +
  scale_fill_brewer(direction = 1, palette = "Greens",
                    limits=c("a", "b", unique(sum_df$grp)), 
                    breaks=unique(sum_df$grp)) +
  labs(fill = "MAF (White)") +
  xlab(expression(gamma)) +
  ylab(expression(r^2))  +
  theme(legend.position="bottom",
        strip.text.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust=0.5))

ylims <- c(0, max(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range,
                  ggplot_build(p4)$layout$panel_scales_y[[1]]$range$range))

p3 <- p3 + ylim(ylims) + theme(legend.position = "none")
p4 <- p4 + ylim(ylims) + theme(legend.position = "none")

# HEIGHT
phe <- "height"
train_name_df <- ntrain_df %>% 
  filter(pheno == phe) %>%
  inner_join(tibble(prop = c(1, 0.1, 0), 
                    train_lab = c("prop_black:\n100%",
                                  "prop_black: 10%\nprop_white: 90%",
                                  "prop_white:\n100%")), by = "prop") %>%
  mutate(n_lab = paste0("n: ", ntrain),
         Training = paste(n_lab, train_lab, sep = "\n"))


sum_df <- split_df %>%
  mutate(semi_partial = partial * (1 - covar)) %>%
  filter(parent_meaning %in% c("Black or Black British", "White") & 
           type == "aall" & pheno == phe) %>%
  filter(pow == 0 | prop == 0.1) %>%
  group_by(pow, ethnicity, grp, prop, parent_meaning) %>%
  summarise(partial = mean(semi_partial),
            mean = mean(semi_partial_split), 
            max = max(partial_split),
            min = min(partial_split), .groups = "drop") %>%
  left_join(train_name_df, by = "prop") %>%
  mutate(Training = paste0(Training, "\n"))

q1 <- sum_df %>% 
  filter(ethnicity == "Black" & parent_meaning == "Black or Black British") %>%
  mutate(Training = factor(Training, levels = sort(unique(Training))[c(1,3,2)])) %>%
  ggplot(aes(factor(pow), mean)) + 
  geom_bar(aes(fill = grp), stat = "identity") +
  geom_point(aes(y = partial)) +
  facet_grid( ~ Training, scales = "free_x", space = "free",
              labeller = labeller(parent_meaning = label_wrap_gen(width = 14))) +
  scale_fill_brewer(direction = 1, palette = "Blues",
                    limits=c("a", "b", unique(sum_df$grp)), 
                    breaks=unique(sum_df$grp)) +
  labs(fill = "MAF (Black)") +
  ylab(expression(r^2))  +
  theme(legend.position="bottom", axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(angle = 0, vjust=0.5),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Height")

q2 <- sum_df %>% 
  filter(ethnicity == "White" & parent_meaning == "Black or Black British") %>%
  mutate(Training = factor(Training, levels = sort(unique(Training))[c(1,3,2)])) %>%
  ggplot(aes(factor(pow), mean)) + 
  geom_bar(aes(fill = grp), stat = "identity") +
  geom_point(aes(y = partial)) +
  facet_grid( ~ Training, scales = "free_x", space = "free") +
  scale_fill_brewer(direction = 1, palette = "Greens",
                    limits=c("a", "b", unique(sum_df$grp)), 
                    breaks=unique(sum_df$grp)) +
  labs(fill = "MAF (White)") +
  xlab(expression(gamma)) +
  ylab(expression(r^2))  +
  theme(legend.position="bottom",
        strip.text.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust=0.5))

ylims <- c(0, max(ggplot_build(q1)$layout$panel_scales_y[[1]]$range$range,
                  ggplot_build(q2)$layout$panel_scales_y[[1]]$range$range))

q1 <- q1 + ylim(ylims) + theme(legend.position = "none")
q2 <- q2 + ylim(ylims) + theme(legend.position = "none")

q3 <- sum_df %>% 
  filter(ethnicity == "Black" & parent_meaning == "White") %>%
  mutate(Training = factor(Training, levels = sort(unique(Training))[c(1,3,2)])) %>%
  ggplot(aes(factor(pow), mean)) + 
  geom_bar(aes(fill = grp), stat = "identity") +
  geom_point(aes(y = partial)) +
  facet_grid( ~ Training, scales = "free_x", space = "free",
              labeller = labeller(parent_meaning = label_wrap_gen(width = 14))) +
  scale_fill_brewer(direction = 1, palette = "Blues",
                    limits=c("a", "b", unique(sum_df$grp)), 
                    breaks=unique(sum_df$grp)) +
  labs(fill = "MAF (Black)") +
  ylab(expression(r^2))  +
  theme(legend.position="bottom", axis.title.x=element_blank(),
        strip.text.x = element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(angle = 0, vjust=0.5))

q4 <- sum_df %>% 
  filter(ethnicity == "White" & parent_meaning == "White") %>%
  mutate(Training = factor(Training, levels = sort(unique(Training))[c(1,3,2)])) %>%
  ggplot(aes(factor(pow), mean)) + 
  geom_bar(aes(fill = grp), stat = "identity") +
  geom_point(aes(y = partial)) +
  facet_grid( ~ Training, scales = "free_x", space = "free") +
  scale_fill_brewer(direction = 1, palette = "Greens",
                    limits=c("a", "b", unique(sum_df$grp)), 
                    breaks=unique(sum_df$grp)) +
  labs(fill = "MAF (White)") +
  xlab(expression(gamma)) +
  ylab(expression(r^2))  +
  theme(legend.position="bottom",
        strip.text.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust=0.5))

ylims <- c(0, max(ggplot_build(q3)$layout$panel_scales_y[[1]]$range$range,
                  ggplot_build(q4)$layout$panel_scales_y[[1]]$range$range))
q3 <- q3 + ylim(ylims)
q4 <- q4 + ylim(ylims)

for (p in c("p1", "p2", "p3", "p4", "q1", "q2", "q3", "q4")) {
  gt = ggplot_gtable(ggplot_build(get(p)))
  gt$widths[7] <- gt$widths[7] * 0.5
  assign(p, wrap_ggplot_grob(gt))
}

p_top <- (p1 + q1) / (p2 + q2) + 
  plot_annotation(title = "Black or Black British test set", 
                  theme = theme(plot.title = element_text(hjust = 0.5)))
p_bottom <- (p3 + q3) / (p4 + q4) / guide_area() + 
  plot_layout(guides = 'collect', heights = c(2, 2, 1)) +
  plot_annotation(title = "White test set", 
                  theme = theme(plot.title = element_text(hjust = 0.5)))

out_file <- paste0("plots/fig5upper.png")
ggsave(out_file, p_top, width = 8, height = 4)

out_file <- paste0("plots/fig5lower.png")
ggsave(out_file, p_bottom, width = 8, height = 4)


### Figure 6 ###

get_labs <- function(x, n = 20) {
  brks <- ggplot2:::breaks(range(x), "width", n = n)

  round(brks[-length(brks)] + diff(brks)/2, 2)
}

pheno <- "MCV"
fold <- 1
prop <- 0

pred_path <- file.path("output", pheno,
                       paste0("fold", fold), paste0("prop", prop),
                       paste0("pred.tsv"))


kfile <- file.path("data", pheno,
                   paste0("fold", fold),
                   paste0("prop0.txt"))
train_id <- read.table(kfile)[, 1]

pfile <- "data/all_vars.tsv"
pheno_df <- read_tsv(pfile)

this_df <- pheno_df %>%
 rename(y = all_of(pheno))
  
train_mean <- this_df %>%
  filter(eid %in% train_id) %>%
  pull(y) %>%
  mean()

pred_df <- read_tsv(pred_path) %>%
  filter(type == "lasso" & pow == 0) %>%
  inner_join(this_df, by = "eid") %>%
  mutate(bias = pred - y,
         abs_bias = abs(bias))
  
out_df <- pred_df %>%
  mutate(covar = lm(y ~ Age + Sex *
                 (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10))$fitted.values,
         full = lm(y ~ Age + pred + Sex *
                 (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10))$fitted.values,
         PRS = y - full,
         Covariate = y - covar) %>%
  select(eid, parent_meaning, PC1, y, covar, full, pred, PRS, Covariate)

p1 <- out_df %>%
  mutate(Residual = PRS) %>%
  ggplot(aes(PC1, Residual)) + 
  geom_point(aes(color = parent_meaning), size = 0.1) +
  guides(colour = guide_legend(override.aes = list(size = 2),
                               title = "Self-reported ethnicity")) +
  xlab(NULL) +
  ggtitle(pheno)

ylims <- out_df %>% 
  group_by(cut_interval(PC1, 20)) %>% 
  summarise(lower = boxplot.stats(PRS)$stats[1], 
            upper = boxplot.stats(PRS)$stats[5]) %>%
  summarise(min(lower), max(upper)) %>%
  as.double()

p2 <- out_df %>%
  mutate(y_centered = y - train_mean) %>%
  pivot_longer(c(PRS, Covariate), 
               names_to = "Type", values_to = "Residual") %>%
  mutate(Type = factor(Type, levels = c("PRS", "Covariate"))) %>%
  ggplot(aes(PC1, Residual)) +
  geom_boxplot(aes(group = interaction(Type, cut_interval(PC1, 20, 
      labels = get_labs(PC1, 20))), fill = Type), outlier.shape = NA) +
    geom_hline(yintercept = 0, alpha = 0.5) +
    coord_cartesian(ylim = ylims*1.05)

ggsave("fig6upper", p1 / p2, width = 10, height = 5)


# Height

pheno <- "height"
fold <- 1
prop <- 0

pred_path <- file.path("output", pheno,
                       paste0("fold", fold), paste0("prop", prop),
                       paste0("pred.tsv"))


kfile <- file.path("data", pheno,
                   paste0("fold", fold),
                   paste0("prop0.txt"))
train_id <- read.table(kfile)[, 1]

pfile <- "data/all_vars.tsv"
pheno_df <- read_tsv(pfile)

this_df <- pheno_df %>%
 rename(y = all_of(pheno))
  
train_mean <- this_df %>%
  filter(eid %in% train_id) %>%
  pull(y) %>%
  mean()

pred_df <- read_tsv(pred_path) %>%
  filter(type == "lasso" & pow == 0) %>%
  inner_join(this_df, by = "eid") %>%
  mutate(bias = pred - y,
         abs_bias = abs(bias))
  
out_df <- pred_df %>%
  mutate(covar = lm(y ~ Age + Sex *
                 (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10))$fitted.values,
         full = lm(y ~ Age + pred + Sex *
                 (PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10))$fitted.values,
         PRS = y - full,
         Covariate = y - covar) %>%
  select(eid, parent_meaning, PC1, y, covar, full, pred, PRS, Covariate)

p1 <- out_df %>%
  mutate(Residual = PRS) %>%
  ggplot(aes(PC1, Residual)) + 
  geom_point(aes(color = parent_meaning), size = 0.1) +
  guides(colour = guide_legend(override.aes = list(size = 2),
                               title = "Self-reported ethnicity")) +
  xlab(NULL) +
  ggtitle(pheno)

ylims <- out_df %>% 
  group_by(cut_interval(PC1, 20)) %>% 
  summarise(lower = boxplot.stats(PRS)$stats[1], 
            upper = boxplot.stats(PRS)$stats[5]) %>%
  summarise(min(lower), max(upper)) %>%
  as.double()

p2 <- out_df %>%
  mutate(y_centered = y - train_mean) %>%
  pivot_longer(c(PRS, Covariate), 
               names_to = "Type", values_to = "Residual") %>%
  mutate(Type = factor(Type, levels = c("PRS", "Covariate"))) %>%
  ggplot(aes(PC1, Residual)) +
  geom_boxplot(aes(group = interaction(Type, cut_interval(PC1, 20, 
      labels = get_labs(PC1, 20))), fill = Type), outlier.shape = NA) +
    geom_hline(yintercept = 0, alpha = 0.5) +
    coord_cartesian(ylim = ylims*1.05)

ggsave("fig6lower", p1 / p2, width = 10, height = 5)