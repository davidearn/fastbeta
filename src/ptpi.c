#include <math.h> /* sqrt */
#include <string.h> /* memcpy */
#include <Rinternals.h>

static
void ptpi0(double *s, double *c, int n, int a, int b, double tol, int itermax,
           double *value, double *delta, int *iter)
{
	int i, j;
	double
		*Z = s, *B = Z + n + 1, *mu = B + n + 1,
		S_a_ = c[0], I_a_ = c[1], R_a_ = c[2],
		S_i_, S_j_, I_i_, I_j_, R_i_, R_j_,
		halfmu, halfgamma = 0.5 * c[3], halfdelta = 0.5 * c[4],
		tmp0, tmp1;

	*iter = 0;
	while (*iter < itermax) {
	S_i_ = S_a_;
	I_i_ = I_a_;
	R_i_ = R_a_;
	halfmu = 0.5 * mu[a];
	for (j = a + 1; j <= b; ++j) {
		tmp0 = 1.0 -  halfmu;
		tmp1 = 1.0 + (halfmu = 0.5 * mu[j]);

		I_j_  = (tmp0 - halfgamma) * I_i_                             + Z[j];
		I_j_ /= tmp1 + halfgamma;

		R_j_  = (tmp0 - halfdelta) * R_i_ + halfgamma * (I_i_ + I_j_);
		R_j_ /= tmp1 + halfdelta;

		S_j_  =  tmp0              * S_i_ + halfdelta * (R_i_ + R_j_) - Z[j] + B[j];
		S_j_ /= tmp1;

		S_i_ = S_j_;
		I_i_ = I_j_;
		R_i_ = R_j_;
	}
	*delta = sqrt(
		((S_i_ - S_a_) * (S_i_ - S_a_) +
		 (I_i_ - I_a_) * (I_i_ - I_a_) +
		 (R_i_ - R_a_) * (R_i_ - R_a_)) /
		(S_a_ * S_a_ + I_a_ * I_a_ + R_a_ * R_a_));
	S_a_ = S_i_;
	I_a_ = I_i_;
	R_a_ = R_i_;
	++(*iter);
	if (ISNAN(*delta) || *delta < tol)
		break;
	}
	halfmu = 0.5 * mu[a];
	for (i = a - 1, j = a; i >= 0; --i, --j) {
		tmp0 = 1.0 +  halfmu;
		tmp1 = 1.0 - (halfmu = 0.5 * mu[i]);

		I_i_  = (tmp0 + halfgamma) * I_j_                             - Z[j];
		I_i_ /= tmp1 - halfgamma;

		R_i_  = (tmp0 + halfdelta) * R_j_ - halfgamma * (I_i_ + I_j_);
		R_i_ /= tmp1 - halfdelta;

		S_i_  =  tmp0              * S_j_ - halfdelta * (R_i_ + R_j_) + Z[j] - B[j];
		S_i_ /= tmp1;

		S_j_ = S_i_;
		I_j_ = I_i_;
		R_j_ = R_i_;
	}
	value[0] = S_j_;
	value[1] = I_j_;
	value[2] = R_j_;
	return;
}

static
void ptpi1(double *s, double *c, int n, int a, int b, double tol, int itermax,
           double *value, double *delta, int *iter, double *x)
{
	x -= a;

	int i, j;
	R_xlen_t ldx = (b - a + 1) * 3;
	double
		*Z = s, *B = Z + n + 1, *mu = B + n + 1,
		*S = x, *I = S + b - a + 1, * R = I + b - a + 1,
		S_a_ = c[0], I_a_ = c[1], R_a_ = c[2],
		halfmu, halfgamma = 0.5 * c[3], halfdelta = 0.5 * c[4],
		tmp0, tmp1;

	*iter = 0;
	while (*iter < itermax) {
	S[a] = S_a_;
	I[a] = I_a_;
	R[a] = R_a_;
	halfmu = 0.5 * mu[a];
	for (i = a, j = a + 1; i < b; ++i, ++j) {
		tmp0 = 1.0 -  halfmu;
		tmp1 = 1.0 + (halfmu = 0.5 * mu[j]);

		I[j]  = (tmp0 - halfgamma) * I[i]                             + Z[j];
		I[j] /= tmp1 + halfgamma;

		R[j]  = (tmp0 - halfdelta) * R[i] + halfgamma * (I[i] + I[j]);
		R[j] /= tmp1 + halfdelta;

		S[j]  =  tmp0              * S[i] + halfdelta * (R[i] + R[j]) - Z[j] + B[j];
		S[j] /= tmp1;
	}
	*delta = sqrt(
		((S[b] - S_a_) * (S[b] - S_a_) +
		 (I[b] - I_a_) * (I[b] - I_a_) +
		 (R[b] - R_a_) * (R[b] - R_a_)) /
		(S_a_ * S_a_ + I_a_ * I_a_ + R_a_ * R_a_));
	S_a_ = S[b];
	I_a_ = I[b];
	R_a_ = R[b];
	++(*iter);
	if (ISNAN(*delta) || *delta < tol)
		break;
	S += ldx;
	I += ldx;
	R += ldx;
	}
	double S_i_, S_j_ = S_a_, I_i_, I_j_ = I_a_, R_i_, R_j_ = R_a_;
	halfmu = 0.5 * mu[a];
	for (i = a - 1, j = a; i >= 0; --i, --j) {
		tmp0 = 1.0 +  halfmu;
		tmp1 = 1.0 - (halfmu = 0.5 * mu[i]);

		I_i_  = (tmp0 + halfgamma) * I_j_                             - Z[j];
		I_i_ /= tmp1 - halfgamma;

		R_i_  = (tmp0 + halfdelta) * R_j_ - halfgamma * (I_i_ + I_j_);
		R_i_ /= tmp1 - halfdelta;

		S_i_  =  tmp0              * S_j_ - halfdelta * (R_i_ + R_j_) + Z[j] - B[j];
		S_i_ /= tmp1;

		S_j_ = S_i_;
		I_j_ = I_i_;
		R_j_ = R_i_;
	}
	value[0] = S_j_;
	value[1] = I_j_;
	value[2] = R_j_;
	return;
}

SEXP R_ptpi(SEXP series, SEXP constants, SEXP a, SEXP b,
            SEXP tol, SEXP itermax, SEXP complete)
{
	int n_ = INTEGER(getAttrib(series, R_DimSymbol))[0] - 1,
		a_ = INTEGER(a)[0], b_ = INTEGER(b)[0],
		itermax_ = INTEGER(itermax)[0];
	double tol_ = REAL(tol)[0];
	SEXP res = PROTECT(allocVector(VECSXP, 4)),
		nms = PROTECT(allocVector(STRSXP, 4)),
		value = PROTECT(allocVector(REALSXP, 3)),
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
		int d[3]; d[0] = b_ - a_ + 1; d[1] = 3; d[2] = itermax_;
		SEXP x = PROTECT(allocVector(REALSXP, (R_xlen_t) d[0] * d[1] * d[2])),
			dim = PROTECT(allocVector(INTSXP, 3));
		ptpi1(REAL(series), REAL(constants), n_, a_, b_, tol_, itermax_,
		      REAL(value), REAL(delta), INTEGER(iter), REAL(x));
		d[2] = INTEGER(iter)[0];
		memcpy(INTEGER(dim), &d, 3 * sizeof(int));
		if (d[2] < itermax_) {
			SEXP y = allocVector(REALSXP, (R_xlen_t) d[0] * d[1] * d[2]);
			memcpy(REAL(y), REAL(x), (size_t) d[0] * d[1] * d[2] * sizeof(double));
			UNPROTECT(1);
			PROTECT(x = y);
		}
		setAttrib(x, R_DimSymbol, dim);
		SET_VECTOR_ELT(res, 3, x);
		UNPROTECT(2);
	}
	else {
		ptpi0(REAL(series), REAL(constants), n_, a_, b_, tol_, itermax_,
		      REAL(value), REAL(delta), INTEGER(iter));
	}
	UNPROTECT(5);
	return res;
}
