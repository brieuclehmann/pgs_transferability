library(dplyr)
library(ggplot2)

n_eur <- 18000
r2_df <- readr::read_csv("output/simulations/full.csv") %>%
    filter(prop_afr != 0) %>%
    mutate(
        n_train = round(n_eur / (1 - prop_afr)),
        h2_gap = h2 - r2
    )

n_afr <- round(unique(n_eur * r2_df$prop_afr / (1 - r2_df$prop_afr)))
pop_string <- c(AFR = "Test: YRI", EUR = "Test: CEU")
n_string <- paste0("n_yri = ", n_afr, "\nn_ceu = 18000")
names(n_string) <- n_afr + n_eur

afr_df <- r2_df %>%
    filter(pop == "YRI") %>%
    group_by(pop, beta_cor, n_train, pow) %>%
    summarise(r2 = mean(h2_gap), .groups = "drop_last")

eur_df <- r2_df %>%
    filter(pop == "CEU" & pow == 0) %>%
    group_by(pop, beta_cor, n_train, pow) %>%
    summarise(r2 = mean(h2_gap), .groups = "drop_last")

p1 <- afr_df %>%
    ggplot(aes(pow, r2)) +
    geom_line(aes(group = beta_cor, color = factor(beta_cor), linetype = pop)) +
    geom_hline(
        data = eur_df,
        aes(
            yintercept = r2, group = beta_cor,
            color = factor(beta_cor), linetype = pop
        )
    ) +
    facet_wrap("n_train",
        nrow = 1,
        labeller = labeller(pop = pop_string, n_train = n_string)
    ) +
    scale_x_continuous(minor_breaks = seq(0, 1, 0.1), breaks = seq(0, 1, 0.2)) +
    scale_color_grey(name = "\u03c1", start = 0.8, end = 0.2) +
    xlab(expression(gamma)) +
    ylab("Predictive gap") +
    #ylab(expression(r^2)) +
    theme(
        legend.position = "bottom",
        axis.title.y = element_text(angle = 0, vjust = 0.5)
    )

dir.create("plots", showWarnings = FALSE)
ggsave("plots/sim-gapVntrain.png", type = "cairo-png", width = 8, height = 4)