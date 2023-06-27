bootstrap <-
function (series, constants)
{
	stopifnot(exprs = {
		is.mts(series)
		is.double(series)
		ncol(series) == 3L
		min(0, series, na.rm = TRUE) >= 0
		is.double(constants)
		length(constants) == 4L
		is.finite(constants)
		all(constants >= 0)
	})

}
