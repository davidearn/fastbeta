if (!requireNamespace("adaptivetau"))
	q("no")

library(fastbeta)
options(warn = 2L, error = if (interactive()) utils::recover)

utils::data(seir.ts02, package = "fastbeta")
a <- attributes(seir.ts02)
m <- a[["m"]]
n <- a[["n"]]
p <- m + n + 2L

series <- cbind(seir.ts02[, c("Z.obs", "B")], mu = a[["mu"]](0))
colnames(series) <- c("Z.obs", "B", "mu")

args <- c(list(series = series),
          a[c("sigma", "gamma", "delta", "init", "m", "n", "prob", "delay")])
X <- do.call(fastbeta, args)
str(X)

stopifnot(exprs = {
	is.double(X)
	stats::is.mts(X)
	identical(dim(X), c(nrow(seir.ts02), p + 1L))
	identical(dimnames(X), list(NULL, rep.int(c("S", "E", "I", "R", "beta"), c(1L, m, n, 1L, 1L))))
	identical(stats::tsp(X), stats::tsp(seir.ts02))
	!anyNA(X[-length(X)])
	min(0, X, na.rm = TRUE) >= 0
})

if (grDevices::dev.interactive(TRUE))
	plot(X)

proc.time()
