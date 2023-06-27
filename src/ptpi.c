#include <math.h> /* fabs */
#include <Rinternals.h>

static void ptpi0(double *Z, double *B, double *mu,
                  double start, int a, int b, double tol, int itermax,
                  double *value, double *delta, int *iter)
{
	int k;
	double halfmu, end = start;
	*iter = 0;
	while (*iter < itermax) {
		halfmu = 0.5 * mu[a];
		for (k = a + 1; k <= b; ++k) {
			end  = (1.0 - halfmu) * end + B[k] - Z[k];
			end /= 1.0 + (halfmu = 0.5 * mu[k]);
		}
		++(*iter);
		*delta = (end - start) / start;
		if (ISNAN(*delta) || fabs(*delta) < tol)
			break;
		start = end;
	}
	halfmu = 0.5 * mu[a];
	for (k = a; k > 0; --k) {
		start  = (1.0 + halfmu) * start - B[k] + Z[k];
		start /= 1.0 - (halfmu = 0.5 * mu[k - 1]);
	}
	*value = start;
	return;
}

static void ptpi1(double *Z, double *B, double *mu,
                  double start, int a, int b, double tol, int itermax,
                  double *value, double *delta, int *iter, double *X)
{
	int k;
	double halfmu;
	*iter = 0;
	X[a] = start;
	while (*iter < itermax) {
		halfmu = 0.5 * mu[a];
		for (k = a + 1; k <= b; ++k) {
			X[k]  = (1.0 - halfmu) * X[k - 1] + B[k] - Z[k];
			X[k] /= 1.0 + (halfmu = 0.5 * mu[k]);
		}
		++(*iter);
		*delta = (X[b] - X[a]) / X[a];
		if (ISNAN(*delta) || fabs(*delta) < tol)
			break;
		X[a] = X[b];
	}
	start = X[a];
	halfmu = 0.5 * mu[a];
	for (k = a; k > 0; --k) {
		start  = (1.0 + halfmu) * start - B[k] + Z[k];
		start /= 1.0 - (halfmu = 0.5 * mu[k - 1]);
	}
	*value = start;
	return;
}

SEXP R_ptpi(SEXP series, SEXP start, SEXP a, SEXP b, SEXP tol, SEXP itermax,
            SEXP complete) {
	int n = INTEGER(getAttrib(series, R_DimSymbol))[0] - 1,
		a_ = INTEGER(a)[0], b_ = INTEGER(b)[0],
		itermax_ = INTEGER(itermax)[0];
	double start_ = REAL(start)[0], tol_ = REAL(tol)[0],
		*s0 = REAL(series), *s1 = s0 + n + 1, *s2 = s1 + n + 1;
	SEXP res = PROTECT(allocVector(VECSXP, 4)),
		nms = PROTECT(allocVector(STRSXP, 4)),
		value = PROTECT(allocVector(REALSXP, 1)),
		delta = PROTECT(allocVector(REALSXP, 1)),
		iter = PROTECT(allocVector(INTSXP, 1));

	SET_STRING_ELT(nms, 0, mkChar("value"));
	SET_STRING_ELT(nms, 1, mkChar("delta"));
	SET_STRING_ELT(nms, 2, mkChar("iter"));
	SET_STRING_ELT(nms, 3, mkChar("X"));
	setAttrib(res, R_NamesSymbol, nms);

	SET_VECTOR_ELT(res, 0, value);
	SET_VECTOR_ELT(res, 1, delta);
	SET_VECTOR_ELT(res, 2, iter);

	if (LOGICAL(complete)[0]) {
		SEXP X = PROTECT(allocMatrix(REALSXP, b_ - a_ + 1, itermax_ + 1));
		SET_VECTOR_ELT(res, 3, X);
		ptpi1(s0, s1, s2, start_, a_, b_, tol_, itermax_,
		      REAL(value), REAL(delta), INTEGER(iter), REAL(X) - a_);
		UNPROTECT(1);
	} else {
		ptpi0(s0, s1, s2, start_, a_, b_, tol_, itermax_,
		      REAL(value), REAL(delta), INTEGER(iter));
	}

	UNPROTECT(5);
	return res;
}
