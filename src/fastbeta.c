#include <Rinternals.h>

static
void fastbeta(double *series,
              double sigma, double gamma, double delta,
              double *init, int m, int n, int lengthOut,
              double *x)
{
	if (lengthOut <= 0)
		return;

	double *x_ = x;
	for (int k = 0, end = m + n + 2; k < end; ++k) {
		*x = *(init++);
		x += lengthOut;
	}
	x[lengthOut - 1] = NA_REAL;
	x = x_;

	if (lengthOut <= 1)
		return;

	double work[4],
		*Z    = series,
		*B    = Z + lengthOut,
		*mu   = B + lengthOut,
		*S    = x,
		halfmu    = 0.5 * mu[0] * (double) 1,
		halfsigma = 0.5 * sigma * (double) m,
		halfgamma = 0.5 * gamma * (double) n,
		halfdelta = 0.5 * delta * (double) 1;

	x = x_ = x_ + lengthOut;

	for (int s = 0, t = 1; t < lengthOut; ++s, ++t) {
		work[0] = 1.0 -  halfmu               ;
		work[1] = 1.0 + (halfmu = 0.5 * mu[t]);
		work[2] = Z[t];
		work[3] = 0.0;

		for (int i = 0; i < m; ++i) {
		x[t] = ((work[0] - halfsigma) * x[s] + work[2])/(work[1] + halfsigma);
		work[2] = halfsigma * (x[s] + x[t]);
		x += lengthOut;
		}
		for (int j = 0; j < n; ++j) {
		x[t] = ((work[0] - halfgamma) * x[s] + work[2])/(work[1] + halfgamma);
		work[2] = halfgamma * (x[s] + x[t]);
		work[3] += x[s];
		x += lengthOut;
		}
		x[t] = ((work[0] - halfdelta) * x[s] + work[2])/(work[1] + halfdelta);
		work[2] = halfdelta * (x[s] + x[t]);
		x += lengthOut;
		S[t] = ( work[0]              * S[s] + work[2])/ work[1]             ;
		x[s] = (Z[s] + Z[t]) / (2.0 * S[s] * work[3]);

		x = x_;
	}

	return;
}

SEXP R_fastbeta(SEXP s_series,
                SEXP s_sigma, SEXP s_gamma, SEXP s_delta,
                SEXP s_init, SEXP s_m, SEXP s_n)
{
	int m = INTEGER(s_m)[0], n = INTEGER(s_n)[0],
		lengthOut = INTEGER(getAttrib(s_series, R_DimSymbol))[0];
	SEXP x = allocMatrix(REALSXP, lengthOut, m + n + 3);

	fastbeta(REAL(s_series),
	         REAL(s_sigma)[0], REAL(s_gamma)[0], REAL(s_delta)[0],
	         REAL(s_init), m, n, lengthOut,
	         REAL(x));

	double *px = REAL(x);
	for (R_xlen_t k = 0, end = XLENGTH(x) - 1; k < end; ++k) {
		if (ISNAN(px[k]) || px[k] < 0.0) {
			warning("entry [%d, %d] of result is negative or NA",
			        (int) (k % lengthOut) + 1, (int) (k / lengthOut) + 1);
			break;
		}
	}

	return x;
}
