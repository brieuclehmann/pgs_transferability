source("scripts/plots_00_prep.R")

n_eur <- 18000
pop_map <- c("CEU" = "EUR", "YRI" = "AFR", "CHB" = "EAS")
r2_df <- readr::read_csv("output/simulations/full.csv", 
                         show_col_types = FALSE) %>%
    filter(prop_afr != 0) %>%
    mutate(
        n_train = round(n_eur / (1 - prop_afr)),
        h2_gap = h2 - r2,
        beta_cor = factor(beta_cor),
        pop = pop_map[pop]
    )

n_afr <- round(unique(n_eur * r2_df$prop_afr / (1 - r2_df$prop_afr)))
pop_string <- c(AFR = "Test: AFR", EUR = "Test: EUR")
n_string <- paste0("No. AFR = ", n_afr, "\nNo. EUR = 18000")
n_train_values <- c(20000, 22500, 25714, 30000, 36000)
n_string <- expression(atop(paste0("n" [AFR], " = ", n_afr), "No. EUR = 18000"))
label_ntrain <- sapply(n_train_values - 18000, 
                       function(x) bquote(atop(n [AFR] == .(x), n [EUR] == 18000)))

names(n_string) <- n_afr + n_eur

afr_df <- r2_df %>%
  filter(pop == "AFR" & pow != Inf) %>%
  group_by(pop, beta_cor, n_train, pow) %>%
  summarise(r2 = mean(h2_gap), .groups = "drop_last") %>%
  mutate(`Training set` = "Multi-ancestry",
         n_train_label = factor(n_train,
                                levels = n_train_values,
                                labels = label_ntrain))

eur_df <- r2_df %>%
    filter(pop == "EUR" & pow == 0) %>%
    group_by(pop, beta_cor, n_train, pow) %>%
    summarise(r2 = mean(h2_gap), .groups = "drop_last") %>%
    mutate(`Training set` = "Multi-ancestry",
           n_train_label = factor(n_train,
                                  levels = n_train_values,
                                  labels = label_ntrain))

min_df <- r2_df %>%
    filter(pop == "AFR" & pow == Inf) %>%
    group_by(pop, beta_cor, n_train, pow) %>%
    summarise(r2 = mean(h2_gap), .groups = "drop_last") %>%
    mutate(pow = 0,
           `Training set` = "AFR only",
           n_train_label = factor(n_train,
                                  levels = n_train_values,
                                  labels = label_ntrain))

p1 <- afr_df %>%
    ggplot(aes(pow, r2)) +
    geom_line(aes(group = beta_cor, 
                  alpha = beta_cor, 
                  linetype = pop, 
                  color = `Training set`)) +
    geom_hline(
        data = eur_df,
        aes(
            yintercept = r2, group = beta_cor,
            alpha = beta_cor, linetype = pop, color = `Training set`
        )
    ) +
    geom_hline(
        data = min_df,
        aes(
            yintercept = r2, group = beta_cor,
            alpha = beta_cor, linetype = pop, color = `Training set`
        )
    ) +
    facet_wrap("n_train_label",
        nrow = 1,
        labeller = labeller(pop = pop_string, n_train_label = label_parsed)
    ) +
    scale_x_continuous(minor_breaks = seq(0, 1, 0.1), breaks = seq(0, 1, 0.2)) +
 #   scale_color_grey(name = "\u03c1", start = 0.8, end = 0.2) +
    scale_linetype_manual(name = "Test set", values = c(1,2,3)) +
    xlab(expression(gamma)) +
    ylab("Predictive gap") +
    #ylab(expression(r^2)) +
    theme(
        legend.position = "bottom",
        axis.title.y = element_text(vjust = 0.5),
        legend.box="vertical",
        legend.margin = margin(),
        legend.box.margin = margin()
    ) +
    scale_colour_viridis_d(begin = 0.45, end = 0.9, direction = -1,
                           labels = c(expression("PGS"["AFR"]), 
                                      expression("PGS"["dual"])),
                           name = "Polygenic score") +
    scale_alpha_discrete(name = "\u03c1") 

dir.create("plots", showWarnings = FALSE)
ggsave("plots/fig2_sim.pdf", device = cairo_pdf, width = 8, height = 5)
