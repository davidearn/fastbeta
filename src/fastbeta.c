#include <Rinternals.h>

static
void fastbeta(const double *series, int lengthOut,
              double sigma, double gamma, double delta,
              int m, int n, const double *init,
              double *x)
{
	if (lengthOut <= 0)
		return;

	const
	double *Z  = series, *B  = Z + lengthOut, *mu = B + lengthOut;

	double *S = x,
		halfsigma = 0.5 * sigma * (double) m,
		halfgamma = 0.5 * gamma * (double) n,
		halfdelta = 0.5 * delta * (double) 1,
		work[4];

	int d[2];
	d[0] = lengthOut;
	d[1] = m + n + 2; /* not counting last column storing 'beta' */

	for (int k = 0; k < d[1]; ++k) {
		x[0] = init[k];
		x += d[0];
	}
	x[d[0] - 1] = NA_REAL;

	for (int s = 0, t = 1; t < d[0]; ++s, ++t) {
		x = S + d[0];

		work[0] = 1.0 - 0.5 * mu[s];
		work[1] = 1.0 + 0.5 * mu[t];
		work[2] = Z[t];
		work[3] = 0.0;

		for (int i = 0; i < m; ++i) {
		x[t] = ((work[0] - halfsigma) * x[s] + work[2])/(work[1] + halfsigma);
		work[2] = halfsigma * (x[s] + x[t]);
		x += d[0];
		}
		for (int j = 0; j < n; ++j) {
		x[t] = ((work[0] - halfgamma) * x[s] + work[2])/(work[1] + halfgamma);
		work[2] = halfgamma * (x[s] + x[t]);
		work[3] += x[s];
		x += d[0];
		}
		x[t] = ((work[0] - halfdelta) * x[s] + work[2])/(work[1] + halfdelta);
		work[2] = halfdelta * (x[s] + x[t]) - Z[t] + B[t];
		x += d[0];
		S[t] = ( work[0]              * S[s] + work[2])/ work[1]             ;
		x[s] = (Z[s] + Z[t]) / (2.0 * S[s] * work[3]);
	}

	return;
}

SEXP R_fastbeta(SEXP s_series,
                SEXP s_sigma, SEXP s_gamma, SEXP s_delta,
                SEXP s_m, SEXP s_n, SEXP s_init)
{
	int m = INTEGER(s_m)[0], n = INTEGER(s_n)[0],
		lengthOut = INTEGER(Rf_getAttrib(s_series, R_DimSymbol))[0];
	SEXP x = PROTECT(Rf_allocMatrix(REALSXP, lengthOut, m + n + 3));

	fastbeta(REAL(s_series), lengthOut,
	         REAL(s_sigma)[0], REAL(s_gamma)[0], REAL(s_delta)[0],
	         m, n, REAL(s_init),
	         REAL(x));

	double *px = REAL(x);
	for (R_xlen_t k = 0, end = XLENGTH(x) - 1; k < end; ++k) {
		if (!ISNAN(px[k]) && px[k] < 0.0) {
			Rf_warning("entry [%d, %d] of result is negative",
			           (int) (k % lengthOut) + 1,
			           (int) (k / lengthOut) + 1);
			break;
		}
	}

	UNPROTECT(1);
	return x;
}
