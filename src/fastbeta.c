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
	for (int i = 0; i < m + n + 2; ++i) {
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

	for (int i = 0, j = 1; j < lengthOut; ++i, ++j) {
		work[0] = 1.0 -  halfmu;
		work[1] = 1.0 + (halfmu = 0.5 * mu[j]);
		work[2] = Z[j];
		work[3] = 0.0;

		for (int k = 0; k < m; ++k) {
		x[j] = ((work[0] - halfsigma) * x[i] + work[2])/(work[1] + halfsigma);
		work[2] = halfsigma * (x[i] + x[j]);
		x += lengthOut;
		}
		for (int k = 0; k < n; ++k) {
		x[j] = ((work[0] - halfgamma) * x[i] + work[2])/(work[1] + halfgamma);
		work[2] = halfgamma * (x[i] + x[j]);
		work[3] += x[i];
		x += lengthOut;
		}
		x[j] = ((work[0] - halfdelta) * x[i] + work[2])/(work[1] + halfdelta);
		work[2] = halfdelta * (x[i] + x[j]);
		x += lengthOut;
		S[j] = ( work[0]              * S[i] + work[2])/ work[1]             ;
		x[i] = (Z[i] + Z[j]) / (2.0 * S[i] * work[3]);

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
	for (R_xlen_t k = 0, l = XLENGTH(x) - 1; k < l; ++k) {
		if (ISNAN(px[k]) || px[k] < 0.0) {
			warning("entry [%d, %d] of result is negative or NA",
			        (int) (k % lengthOut) + 1, (int) (k / lengthOut) + 1);
			break;
		}
	}

	return x;
}
