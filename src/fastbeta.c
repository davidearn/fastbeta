#include <Rinternals.h>

static void fastbeta(double *Z, double *B, double *mu, double gamma,
                     double *beta, double *S, double *I, int n)
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

SEXP R_fastbeta(SEXP series, SEXP constants)
{
	int n = INTEGER(getAttrib(series, R_DimSymbol))[0] - 1;
	SEXP res = PROTECT(allocMatrix(REALSXP, n + 1, 3));
	if (n >= 0) {
		double *r0 = REAL(res), *r1 = r0 + n + 1, *r2 = r1 + n + 1,
			*cc = REAL(constants);
		r0[0] = NA_REAL;
		r1[0] = cc[1];
		r2[0] = cc[2];
		if (n >= 1) {
			double *s0 = REAL(series), *s1 = s0 + n + 1, *s2 = s1 + n + 1;
			fastbeta(s0, s1, s2, cc[0], r0, r1, r2, n);
		}
	}
	UNPROTECT(1);

	return res;
}
