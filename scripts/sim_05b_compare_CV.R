# Combine output for a simulation

library(foreach)
library(dplyr)
library(ggplot2)
source("sim_utils.R")

beta_cor <- 0.7
prop_afr <- 0.1

pow_range <- seq(0, 1, 0.1)
iter_range <- as.numeric(seq(5))
rep_range <- as.numeric(seq(10))

for (prefix in c("cv", "")) {
    out_dir <- file.path("output", "simulations", prefix)

    params <- expand.grid(pow = pow_range, iter = iter_range, rep = rep_range)

    out_df <- foreach(n = seq(nrow(params)), .combine = rbind) %do% {
        pow <- params$pow[n]
        iter <- params$iter[n]
        rep <- params$rep[n]

        pred_dir <- file.path(
            out_dir,
            paste0("beta_cor~", format(beta_cor, nsmall = 1)),
            paste0("prop_afr~", format(prop_afr, nsmall = 1)),
            paste0("pow~", format(pow, nsmall = 1)),
            paste0("iter~", format(iter, nsmall = 1))
        )
        pred_file <- file.path(pred_dir, paste0("rep~", format(rep, nsmall = 1), ".csv"))
        pred_df <- readr::read_csv(pred_file)

        out <- pred_df %>%
            group_by(pop) %>%
            summarise(
                r2 = compute_r2(pheno, pred),
                mse = mean((pred - pheno)^2),
                cor = cor(pred, pheno),
                .groups = "drop"
            ) %>%
            mutate(pow = pow, iter = iter, rep = rep)

        out
    }

    out_file <- file.path(out_dir, "combined.csv")
    readr::write_csv(out_df, out_file)
}

cv_file <- file.path("output", "simulations", "cv", "combined.csv")
cv_df <- readr::read_csv(cv_file) %>%
    mutate(lambda = "cv")

val_file <- file.path("output", "simulations", "combined.csv") 
val_df <- readr::read_csv(val_file) %>%
    mutate(lambda = "val")

combi_df <- bind_rows(cv_df, val_df) %>%
    group_by(lambda, pop, pow) %>%
    summarise(r2 = mean(r2), .groups = "drop") %>%
    filter(pop != "CHB")

p1 <- ggplot(combi_df, aes(pow, r2, color = lambda)) +
    geom_line() +
    facet_wrap("pop")

ggsave("plots/cv.png", p1, type = "cairo-png", width = 10)

combi_iter_df <- bind_rows(cv_df, val_df) %>%
    group_by(lambda, pop, pow, iter) %>%
    arrange(r2) %>%
    summarise(upper = nth(r2, 9), lower = nth(r2, 2), r2 = mean(r2),  
              .groups = "drop") %>%
    filter(pop != "CHB")

p2 <- ggplot(combi_iter_df, aes(pow, r2, color = lambda, fill = lambda)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    facet_wrap(c("iter", "pop"), nrow = 5, 
               labeller = labeller(.multi_line=FALSE))
ggsave("plots/cv_iter.png", p2, type = "cairo-png", height = 10)