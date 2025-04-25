source("seir.R"); stopifnot(is.list(data))
options(mc.cores = parallel::detectCores())
fit <- rstan::stan("seir.stan", data = data, init = "0")
