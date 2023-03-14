# Set simulation parameters

params <- expand.grid(
    prop_afr = seq(0, 0.5, 0.1),
    pow = seq(0, 1.5, 0.1),
    iter = seq(5),
    rep = seq(5)
)
params_afr <- expand.grid(
    prop_afr = seq(0.1, 0.5, 0.1),
    pow = Inf,
    iter = seq(5),
    rep = seq(5)
)

params <- rbind(params, params_afr)

extra_params <- rbind(
    data.frame(n_causal = 100, h2 = c(0.1, 0.3, 0.5), beta_cor = 0.8),
    data.frame(n_causal = c(10, 100, 1000), h2 = 0.3, beta_cor = 0.8),
    data.frame(n_causal = 100, h2 = 0.3, beta_cor = seq(0.5, 1, 0.1))
)

all_params <- merge(extra_params, params)

all_params$pow <- format(all_params$pow, digits = 2)
all_params$prop_afr <- format(all_params$prop_afr, digits = 2)

readr::write_csv(all_params, "data/sim_params.csv")