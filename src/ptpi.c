#include <math.h> /* sqrt */
#include <string.h> /* memcpy */
#include <Rinternals.h>
#include <R_ext/RS.h>

static
void ptpi(const double *series, int lengthOut,
          double sigma, double gamma, double delta,
          int m, int n, const double *init,
          int a, int b, double tol, int iterMax, int backcalc,
          double *value, double *diff, int *iter,
          double *x)
{
	memcpy(value, init, (m + n + 2) * sizeof(double));
	*diff = 0.0;
	*iter = (lengthOut <= 0) ? iterMax : 0;

	if (lengthOut <= 0)
		return;

	const
	double *Z = series, *B = Z + lengthOut, *mu = B + lengthOut;

	double *S = (x) ? x - a : R_Calloc((size_t) 2 * (m + n + 2), double),
		halfsigma = 0.5 * sigma * (double) m,
		halfgamma = 0.5 * gamma * (double) n,
		halfdelta = 0.5 * delta * (double) 1,
		work[4];

	int d[2];
	d[0] = b - a;
	d[1] = m + n + 2;

	if (x) {

	while (*iter < iterMax) {

	x = S;
	for (int k = 0; k < d[1]; ++k) {
		x[a] = value[k];
		x += d[0];
	}

	for (int s = a, t = a + 1; t < b; ++s, ++t) {
		x = S + d[0];
		work[0] = 1.0 - 0.5 * mu[s];
		work[1] = 1.0 + 0.5 * mu[t];
		work[2] = Z[t];
		for (int i = 0; i < m; ++i) {
		x[t] = ((work[0] - halfsigma) * x[s] + work[2])/(work[1] + halfsigma);
		work[2] = halfsigma * (x[s] + x[t]);
		x += d[0];
		}
		for (int j = 0; j < n; ++j) {
		x[t] = ((work[0] - halfgamma) * x[s] + work[2])/(work[1] + halfgamma);
		work[2] = halfgamma * (x[s] + x[t]);
		x += d[0];
		}
		x[t] = ((work[0] - halfdelta) * x[s] + work[2])/(work[1] + halfdelta);
		work[2] = halfdelta * (x[s] + x[t]) - Z[t] + B[t];
		x += d[0];
		S[t] = ( work[0]              * S[s] + work[2])/ work[1]             ;
	}

	*iter += 1;

	x = S;
	work[0] = work[1] = 0.0;
	for (int k = 0; k < d[1]; ++k) {
		work[0] += (x[b - 1] - value[k]) * (x[b - 1] - value[k]);
		work[1] +=             value[k]  *             value[k] ;
		value[k] = x[b - 1];
		x += d[0];
	}

	*diff = sqrt(work[0] / work[1]);
	if (ISNAN(*diff) || *diff < tol)
		break;

	S += (R_xlen_t) d[0] * d[1];

	}

	}

	else {

	while (*iter < iterMax) {

	x = S;
	for (int k = 0; k < d[1]; ++k) {
		x[0] = value[k];
		x += 2;
	}

	for (int s = a, t = a + 1; t < b; ++s, ++t) {
		x = S + 2;
		work[0] = 1.0 - 0.5 * mu[s];
		work[1] = 1.0 + 0.5 * mu[t];
		work[2] = Z[t];
		for (int i = 0; i < m; ++i) {
		x[1] = ((work[0] - halfsigma) * x[0] + work[2])/(work[1] + halfsigma);
		work[2] = halfsigma * (x[0] + x[1]);
		x[0] = x[1]; x += 2;
		}
		for (int j = 0; j < n; ++j) {
		x[1] = ((work[0] - halfgamma) * x[0] + work[2])/(work[1] + halfgamma);
		work[2] = halfgamma * (x[0] + x[1]);
		x[0] = x[1]; x += 2;
		}
		x[1] = ((work[0] - halfdelta) * x[0] + work[2])/(work[1] + halfdelta);
		work[2] = halfdelta * (x[0] + x[1]) - Z[t] + B[t];
		x[0] = x[1]; x += 2;
		S[1] = ( work[0]              * S[0] + work[2])/ work[1]             ;
		S[0] = S[1];
	}

	*iter += 1;

	x = S;
	work[0] = work[1] = 0.0;
	for (int k = 0; k < d[1]; ++k) {
		work[0] += (x[1] - value[k]) * (x[1] - value[k]);
		work[1] +=         value[k]  *         value[k] ;
		value[k] = x[1];
		x += 2;
	}

	*diff = sqrt(work[0] / work[1]);
	if (ISNAN(*diff) || *diff < tol)
		break;

	}

	R_Free(S);

	}

	if (backcalc) {

	for (int s = a - 1, t = a; s >= 0; --s, --t) {
		work[0] = 1.0 - 0.5 * mu[s];
		work[1] = 1.0 + 0.5 * mu[t];
		work[2] = Z[t];

		for (int i = 0; i < m; ++i) {
		work[3] = value[1 + i    ];
		value[1 + i    ] = ((work[1] + halfsigma) * value[1 + i    ] - work[2])/(work[0] - halfsigma);
		work[2] = halfsigma * (work[3] + value[1 + i    ]);
		}
		for (int j = 0; j < n; ++j) {
		work[3] = value[1 + m + j];
		value[1 + m + j] = ((work[1] + halfgamma) * value[1 + m + j] - work[2])/(work[0] - halfgamma);
		work[2] = halfgamma * (work[3] + value[1 + m + j]);
		}
		work[3] = value[1 + m + n];
		value[1 + m + n] = ((work[1] + halfdelta) * value[1 + m + n] - work[2])/(work[0] - halfdelta);
		work[2] = halfdelta * (work[3] + value[1 + m + n]) - Z[t] + B[t];
		value[0        ] = ( work[1]              * value[0        ] - work[2])/ work[0]             ;
	}

	}

	return;
}

SEXP R_ptpi(SEXP s_series, SEXP s_sigma, SEXP s_gamma, SEXP s_delta,
            SEXP s_m, SEXP s_n, SEXP s_init,
            SEXP s_a, SEXP s_b, SEXP s_tol, SEXP s_iterMax,
            SEXP s_backcalc, SEXP s_complete)
{
	int m = INTEGER(s_m)[0], n = INTEGER(s_n)[0],
		lengthOut = INTEGER(getAttrib(s_series, R_DimSymbol))[0],
		a = INTEGER(s_a)[0], b = INTEGER(s_b)[0],
		iterMax = INTEGER(s_iterMax)[0],
		complete = LOGICAL(s_complete)[0];

	SEXP ans = PROTECT(allocVector(VECSXP, 4)),
		nms = PROTECT(allocVector(STRSXP, 4)),
		value = PROTECT(allocVector(REALSXP, m + n + 2)),
		diff = PROTECT(allocVector(REALSXP, 1)),
		iter = PROTECT(allocVector(INTSXP, 1)),
		x = R_NilValue;

	int d[3];
	if (complete) {
		d[0] = b - a;
		d[1] = m + n + 2;
		d[2] = iterMax;
		x = allocVector(REALSXP, (R_xlen_t) d[0] * d[1] * d[2]);
	}

	PROTECT(x);

	SET_STRING_ELT(nms, 0, mkChar("value"));
	SET_STRING_ELT(nms, 1, mkChar("diff"));
	SET_STRING_ELT(nms, 2, mkChar("iter"));
	SET_STRING_ELT(nms, 3, mkChar("x"));
	setAttrib(ans, R_NamesSymbol, nms);

	SET_VECTOR_ELT(ans, 0, value);
	SET_VECTOR_ELT(ans, 1, diff);
	SET_VECTOR_ELT(ans, 2, iter);

	ptpi(REAL(s_series), lengthOut,
	     REAL(s_sigma)[0], REAL(s_gamma)[0], REAL(s_delta)[0],
	     m, n, REAL(s_init),
	     a, b, REAL(s_tol)[0], iterMax, LOGICAL(s_backcalc)[0],
	     REAL(value), REAL(diff), INTEGER(iter), (complete) ? REAL(x) : NULL);

	if (complete) {
		d[2] = INTEGER(iter)[0];

		SEXP dim = PROTECT(allocVector(INTSXP, 3));
		memcpy(INTEGER(dim), &d, 3 * sizeof(int));

		if (d[2] >= iterMax) {
		setAttrib(x, R_DimSymbol, dim);
		SET_VECTOR_ELT(ans, 3, x);
		}
		else {
		SEXP y = PROTECT(allocVector(REALSXP, (R_xlen_t) d[0] * d[1] * d[2]));
		memcpy(REAL(y), REAL(x), (size_t) d[0] * d[1] * d[2] * sizeof(double));
		setAttrib(y, R_DimSymbol, dim);
		SET_VECTOR_ELT(ans, 3, y);
		UNPROTECT(1);
		}

		UNPROTECT(1);
	}

	UNPROTECT(6);
	return ans;
}
