library(dplyr)

sim_df <- readr::read_csv("data/sim_params.csv")

compute_r2 <- function(y, pred) {
    ss_tot <- sum((y - mean(y))^2)
    ss_res <- sum((y - pred)^2)

    1 - (ss_res / ss_tot)
}

extract_results <- function(beta_cor, prop_afr, pow, iter, rep) {

    if (is.infinite(pow)) {
        pow_txt <- "inf"
    } else {
        pow_txt <- format(pow, nsmall = 1)
    }
    out_dir <- file.path(
        "output", "simulations",
        paste0("beta_cor~", format(beta_cor, nsmall = 1)),
        paste0("prop_afr~", format(prop_afr, nsmall = 1)),
        paste0("pow~", pow_txt),
        paste0("iter~", format(iter, nsmall = 1))
    )

    out_file <- file.path(
        out_dir,
        paste0("rep~", format(rep, nsmall = 1), ".csv")
    )

    readr::read_csv(out_file) %>%
        group_by(pop) %>%
        summarise(
            bias = mean(pred) - mean(pheno),
            r2 = compute_r2(pheno, pred),
            mse = mean((pred - pheno)^2),
            cor = cor(pred, pheno),
            var = var(pheno),
            .groups = "drop"
        ) %>%
        mutate(
            beta_cor = beta_cor,
            prop_afr = prop_afr,
            pow = pow,
            iter = iter,
            rep = rep
        )
}

h2_df <- readr::read_csv("output/simulations/h2.csv")

out_df <- sim_df %>%
    rowwise() %>%
    summarise(
        out = extract_results(beta_cor, prop_afr, pow, iter, rep),
        .groups = "drop") %>%
    pull(out) %>%
    left_join(h2_df, by = c("pop", "beta_cor", "iter"))

readr::write_csv(out_df, "output/simulations/full.csv")