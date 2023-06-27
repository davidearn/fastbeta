#include <Rinternals.h>

static void fastbeta(double *Z, double *B, double *mu, double gamma,
                     double *S, double *I, double *beta, int n)
{
	int k, Sw = 1, Iw = 1;
	double tmp, Z_ = Z[0], S_ = S[0], I_ = I[0],
		halfmu = 0.5 * mu[0], halfgamma = 0.5 * gamma;
	++Z; ++B; ++mu; ++S; ++I;
	for (k = 0; k < n; ++k) {
		tmp = 1.0 - halfmu;
		S[k] = tmp * S_ + B[k] - Z[k];
		I[k] = (tmp - halfgamma) * I_ + Z[k];
		if (Sw && (ISNAN(S[k]) || S[k] < 0.0)) {
			warning("S[k] is NA or negative, k=%d", k + 1);
			Sw = 0;
		}
		if (Iw && (ISNAN(I[k]) && I[k] < 0.0)) {
			warning("I[k] is NA or negative, k=%d", k + 1);
			Iw = 0;
		}
		tmp = 1.0 + (halfmu = 0.5 * mu[k]);
		S[k] /= tmp;
		I[k] /= tmp + halfgamma;
		beta[k] = 0.5 * (Z_ + Z[k]) / (S_ * I_);
		Z_ = Z[k];
		S_ = S[k];
		I_ = I[k];
	}
	return;
}

SEXP R_fastbeta(SEXP data, SEXP gamma, SEXP S0, SEXP I0)
{
	SEXP dim = PROTECT(getAttrib(data, R_DimSymbol));
	int n = INTEGER(dim)[0] - 1;
	UNPROTECT(1);

	SEXP res = PROTECT(allocMatrix(REALSXP, n + 1, 3));
	if (n >= 0) {
		double *r0 = REAL(res), *r1 = r0 + n + 1, *r2 = r1 + n + 1;
		r0[0] = REAL(S0)[0];
		r1[0] = REAL(I0)[0];
		r2[n] = NA_REAL;
		if (n >= 1) {
			double *d0 = REAL(data), *d1 = d0 + n + 1, *d2 = d1 + n + 1;
			fastbeta(d0, d1, d2, REAL(gamma)[0], r0, r1, r2, n);
		}
	}
	UNPROTECT(1);

	return res;
}
