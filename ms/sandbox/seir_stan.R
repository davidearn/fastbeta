source("seir.R"); stopifnot(is.list(data))
options(mc.cores = parallel::detectCores())
fit <- rstan::stan("seir_stan.stan", data = data, init = "0")
