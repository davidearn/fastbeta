source("seir.R"); stopifnot(is.list(data))

dll <- "seir_tmb"
cpp <- paste0(dll, ".cpp")
lib <- paste0(dll, .Platform[["dynlib.ext"]])
TMB::compile(cpp)
dyn.load(lib)
parameters <- list(log_sd = 0,
                   log_size = 0,
                   b0 = double(data$R0),
                   b1 = double(data$R1))
obj <- TMB::MakeADFun(data = data, parameters = parameters,
                      random = "b1", hessian = TRUE, DLL = dll)
opt <- optim(obj$par, obj$fn, obj$gr, method = obj$method)
sdr <- TMB::sdreport(obj)
sdr. <- summary(sdr)
sdr.. <- unname(sdr.)[rownames(sdr.) == "log_trans", , drop = FALSE]

wald <-
function (value, se, level) {
	h <- 0.5 * (1 - level)
	p <- c(h, 1 - h)
	q <- qnorm(p)
	ans <- value + rep(q, each = length(se)) * se
	dim(ans) <- c(length(se), 2L)
	ans
}

pe <- exp(sdr..[, 1L])
ci <- exp(wald(sdr..[, 1L], sdr..[, 2L], level = 0.95))

matplot(times, cbind(beta(times), pe, ci, deparse.level = 0L),
        lty = c(1, 2, 3, 3), lwd = c(2, 2, 1, 1), col = c(1, 2, 2, 2),
        type = "l", xlab = "Time", ylab = "Transmission rate")

options(mc.cores = parallel::detectCores())

set.seed(0L)
fit <- tmbstan::tmbstan(obj, init = "0")
sit <- rstan::summary(fit, probs = c(0.025, 0.5, 0.975))
pos <- grep("^log_trans[[][1-9][0-9]*[]]$", names(fit))
stopifnot(length(pos) == length(times))

matplot(times, cbind(beta(times), exp(sit$summary[pos, c("50%", "2.5%", "97.5%")]), deparse.level = 0L),
        lty = c(1, 2, 3, 3), lwd = c(2, 2, 1, 1), col = c(1, 2, 2, 2),
        type = "l", xlab = "Time", ylab = "Transmission rate")
