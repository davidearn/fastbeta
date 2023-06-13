#include <math.h>
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

SEXP R_fastbeta(SEXP Z, SEXP B, SEXP mu, SEXP gamma, SEXP S0, SEXP I0)
{
    int n = LENGTH(Z) - 1;
    SEXP res = PROTECT(allocVector(VECSXP, 3)),
	nms = PROTECT(allocVector(STRSXP, 3)),
	S = PROTECT(allocVector(REALSXP, n + 1)),
	I = PROTECT(allocVector(REALSXP, n + 1)),
	beta = PROTECT(allocVector(REALSXP, n + 1));

    SET_STRING_ELT(nms, 0, mkChar("S"));
    SET_STRING_ELT(nms, 1, mkChar("I"));
    SET_STRING_ELT(nms, 2, mkChar("beta"));
    setAttrib(res, R_NamesSymbol, nms);

    SET_VECTOR_ELT(res, 0, S);
    SET_VECTOR_ELT(res, 1, I);
    SET_VECTOR_ELT(res, 2, beta);

    if (n >= 0) {
	REAL(S)[0] = REAL(S0)[0];
	REAL(I)[0] = REAL(I0)[0];
	REAL(beta)[n] = NA_REAL;
	if (n >= 1)
	    fastbeta(REAL(Z), REAL(B), REAL(mu), REAL(gamma)[0],
		     REAL(S), REAL(I), REAL(beta), n);
    }
    
    UNPROTECT(5);
    return res;
}

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
    		  double *value, double *delta, int *iter, double *S)
{
    int k;
    double halfmu;
    *iter = 0;
    S[a] = start;
    while (*iter < itermax) {
	halfmu = 0.5 * mu[a];
	for (k = a + 1; k <= b; ++k) {
	    S[k]  = (1.0 - halfmu) * S[k - 1] + B[k] - Z[k];
	    S[k] /= 1.0 + (halfmu = 0.5 * mu[k]);
	}
	++(*iter);
	*delta = (S[b] - S[a]) / S[a];
	if (ISNAN(*delta) || fabs(*delta) < tol)
	    break;
	S[a] = S[b];
    }
    start = S[a];
    halfmu = 0.5 * mu[a];
    for (k = a; k > 0; --k) {
	start  = (1.0 + halfmu) * start - B[k] + Z[k];
	start /= 1.0 - (halfmu = 0.5 * mu[k - 1]);
    }
    *value = start;
    return;
}

SEXP R_ptpi(SEXP Z, SEXP B, SEXP mu,
	    SEXP start, SEXP a, SEXP b, SEXP tol, SEXP itermax,
	    SEXP complete) {
    int a_ = INTEGER(a)[0], b_ = INTEGER(b)[0],
	itermax_ = INTEGER(itermax)[0];
    double start_ = REAL(start)[0], tol_=  REAL(tol)[0];
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
	ptpi1(REAL(Z), REAL(B), REAL(mu),
	      start_, a_, b_, tol_, itermax_,
	      REAL(value), REAL(delta), INTEGER(iter), REAL(X) - a_);
	UNPROTECT(1);
    } else {
	ptpi0(REAL(Z), REAL(B), REAL(mu),
	      start_, a_, b_, tol_, itermax_,
	      REAL(value), REAL(delta), INTEGER(iter));
    }
    
    UNPROTECT(5);
    return res;    
}
