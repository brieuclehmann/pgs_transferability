# Set simulation parameters

params <- expand.grid(
    beta_cor = seq(0.5, 0.9, 0.1),
    prop_afr = seq(0, 0.5, 0.1),
    pow = seq(0, 1, 0.1),
    iter = seq(5),
    rep = seq(5)
)
params_afr <- expand.grid(
    beta_cor = seq(0.5, 0.9, 0.1),
    prop_afr = seq(0.1, 0.5, 0.1),
    pow = Inf,
    iter = seq(5),
    rep = seq(5)
)

params <- rbind(params, params_afr)

params$pow <- format(params$pow, digits = 2)
params$prop_afr <- format(params$prop_afr, digits = 2)

readr::write_csv(params, "data/sim_params.csv")