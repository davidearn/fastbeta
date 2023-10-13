#include <Rinternals.h>

static
void fastbeta(double *s, double *c, int n, double *r)
{
	double
		*Z = s, *B = Z + n + 1, *mu = B + n + 1,
		*S = r, *I = S + n + 1, * R = I + n + 1, *beta = R + n + 1;
	S[0] = c[2]; I[0] = c[3]; R[0] = c[4]; beta[n] = NA_REAL;

	if (n == 0)
		return;

	int i, j, k;
	char name[] = { 'S', 'I', 'R' }, warn[] = { 1, 1, 1 };
	double halfmu = 0.5 * mu[0], halfgamma = 0.5 * c[0], halfdelta = 0.5 * c[1],
		tmp, *r_;

	for (i = 0, j = 1; i < n; ++i, ++j) {
		tmp = 1.0 - halfmu;
		I[j] = (tmp - halfgamma) * I[i]                             + Z[j];
		R[j] = (tmp - halfdelta) * R[i] + halfgamma * (I[i] + I[j]);
		S[j] =  tmp              * S[i] + halfdelta * (R[i] + R[j]) - Z[j] + B[j];
		r_ = r;
		for (k = 0; k < 2; ++k) {
			if (warn[k] && (ISNAN(*r) || *r < 0.0)) {
				warning("%c[%d] is NA or negative", name[k], k);
				warn[k] = 0;
			}
			r += n + 1;
		}
		r = r_ + 1;
		tmp = 1.0 + (halfmu = 0.5 * mu[j]);
		S[j] /= tmp;
		I[j] /= tmp + halfgamma;
		R[j] /= tmp + halfdelta;
		beta[j] = 0.5 * (Z[i] + Z[j]) / (S[i] * I[i]);
	}
	return;
}

SEXP R_fastbeta(SEXP series, SEXP constants)
{
	int n = INTEGER(getAttrib(series, R_DimSymbol))[0] - 1;
	SEXP res = allocMatrix(REALSXP, n + 1, 4);
	if (n >= 0)
		fastbeta(REAL(series), REAL(constants), n, REAL(res));
	return res;
}
