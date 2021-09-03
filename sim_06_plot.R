library(driftlasso)
library(tidyverse)
library(patchwork)
theme_set(theme_minimal(base_size = 12))

n_eur <- 1800
beta_cor_range <- c(0.5, 0.6, 0.7, 0.8, 0.9)
prop_afr_range <- c(0.1, 0.2, 0.3, 0.4, 0.5)

r2_df <- read_csv("output/beta_cor0.5-prop_afr0.1.csv")
r2_list <- list()
for (beta_cor in beta_cor_range) {
    for (prop_afr in prop_afr_range) {
        sim_id <- paste0("beta_cor", beta_cor, "-prop_afr", prop_afr)
        n_afr <- (prop_afr / (1 - prop_afr)) * n_eur
        n_train <- n_eur + n_afr

        this_df <- read_csv(paste0("output/", sim_id, ".csv")) %>%
          mutate(beta_cor = beta_cor,
                 n_train = n_train,
                 sim_id = sim_id)
        
        r2_list <- c(r2_list, this_df)
    }
}
r2_df <- bind_rows(r2_list)

pow_range <- unique(r2_df$pow)

pop_string <- c(AFR = "Test: AFR", EUR = "Test: EUR")
n_string <- r2_df %>%
  distinct(n_train) %>%
  mutate(label = paste0("n_afr = ", n_train - 1800, "\nn_eur = 1800")) %>%
  deframe()

afr_df <- r2_df %>%
  filter(pop == "AFR") %>%
  group_by(pop, beta_cor, n_train, pow) %>%
  summarise(r2 = mean(r2_y), .groups = "drop_last") 

eur_df <- r2_df %>%
  filter(pop == "EUR" & pow == 0) %>%
  group_by(pop, beta_cor, n_train, pow) %>%
  summarise(r2 = mean(r2_y), .groups = "drop_last") 

afr_df %>%
  ggplot(aes(pow, r2)) +
  geom_line(aes(group = beta_cor, color = factor(beta_cor), linetype = pop)) +
  geom_hline(data = eur_df, 
             aes(yintercept = r2, group = beta_cor, 
                 color = factor(beta_cor), linetype = pop)) +
  facet_wrap("n_train", nrow = 1,
             labeller = labeller(pop = pop_string, n_train = n_string)) + 
  scale_x_continuous(minor_breaks = seq(0, 1, 0.1), breaks = seq(0, 1, 0.2)) +
  scale_color_grey(name = "\u03c1", start = 0.8, end = 0.2) + 
  xlab(expression(gamma)) +
  ylab(expression(r^2)) +
  theme(legend.position = "bottom",
        axis.title.y = element_text(angle = 0, vjust=0.5))

dir.create("plots", showWarnings=FALSE)
ggsave("plots/sim-r2vntrain.png", width = 8, height = 4)
