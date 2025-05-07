source("seir.R"); stopifnot(is.list(data))
options(mc.cores = parallel::detectCores())

set.seed(0L)
fit <- rstan::stan("seir_stan.stan", data = data, init = "0")
sit <- rstan::summary(fit, probs = c(0.025, 0.5, 0.975))
pos <- grep("^log_trans[[][1-9][0-9]*[]]$", names(fit))
stopifnot(length(pos) == length(times))

matplot(times, cbind(beta(times), exp(sit$summary[pos, c("50%", "2.5%", "97.5%")]), deparse.level = 0L),
        lty = c(1, 2, 3, 3), lwd = c(2, 2, 1, 1), col = c(1, 2, 2, 2),
        type = "l", xlab = "Time", ylab = "Transmission rate")
