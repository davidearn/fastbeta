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
	for (int i = 0; i < m + n + 2) {
		*x_ = *(x++);
		x_ += lengthOut;
	}
	x_[lengthOut - 1] = NA_REAL;

	if (lengthOut <= 1)
		return;

	double tmp0, tmp1,
		*Z    = series,
		*B    = Z + lengthOut,
		*mu   = B + lengthOut,
		*S    = x,
		*E    = S + lengthOut,
		*I    = E + lengthOut * m,
		*R    = I + lengthOut * n,
		*beta = R + lengthOut,
		hmu    = 0.5 * mu[0] * (double) 1,
		hsigma = 0.5 * sigma * (double) m,
		hgamma = 0.5 * gamma * (double) n,
		hdelta = 0.5 * delta * (double) 1;
	char name[] = { 'S', 'E', 'I', 'R' }, warn[] = { 1, 1, 1, 1 };

	for (i = 0, j = 1; i < n; ++i, ++j) {
		tmp0 = 1.0 -  halfmu;
		tmp1 = 1.0 + (halfmu = 0.5 * mu[j]);

		I[j]  = (tmp0 - halfgamma) * I[i]                             + Z[j];
		I[j] /= tmp1 + halfgamma;

		R[j]  = (tmp0 - halfdelta) * R[i] + halfgamma * (I[i] + I[j]);
		R[j] /= tmp1 + halfdelta;

		S[j]  =  tmp0              * S[i] + halfdelta * (R[i] + R[j]) - Z[j] + B[j];
		S[j] /= tmp1;

		r_ = r;
		for (k = 0; k < 2; ++k) {
			if (warn[k] && (ISNAN(*r) || *r < 0.0)) {
				warning("%c[%d] is NA or negative", name[k], k);
				warn[k] = 0;
			}
			r += n + 1;
		}
		r = r_ + 1;

		beta[i] = 0.5 * (Z[i] + Z[j]) / (S[i] * I[i]);
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
	return x;
}
