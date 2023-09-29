#include <math.h> /* exp */
#include <string.h> /* memset */
#include <Rinternals.h>

static SEXP F, DF, betaCall, betaArg, nuCall, nuArg, muCall, muArg, s;
static double *pF, *pDF, *pbetaArg, betaVal, *pnuArg, nuVal, *pmuArg, muVal,
	gammaVal, deltaVal, lastTimeDot, lastTimeJac, tmp, *pt, *px;

#define INIT_CALL(_V_)                        \
do {                                          \
	_V_ ## Call = allocVector(LANGSXP, 2);    \
	R_PreserveObject(_V_ ## Call);            \
	_V_ ## Arg = allocVector(REALSXP, 1);     \
	R_PreserveObject(_V_ ## Arg);             \
	SETCAR(_V_ ## Call, _V_);                 \
	SETCADR(_V_ ## Call, _V_ ## Arg);         \
	p ## _V_ ## Arg = REAL(_V_ ## Arg);       \
} while (0)

#define TRY_EVAL_CALL(_V_)                                                \
do {                                                                      \
	s = eval(_V_ ## Call, R_GlobalEnv);                                   \
	if (TYPEOF(s) != REALSXP)                                             \
		error("'%s' did not evaluate to type \"%s\"", #_V_, "double");    \
	if (LENGTH(s) != 1)                                                   \
		error("'%s' did not evaluate to length %d", #_V_, 1)              \
	_V_ ## Val = REAL(s)[0];                                              \
	if (!R_FINITE(_V_ ## Val) || _V_ ## Val < 0.0)                        \
		error("'%s' returned a nonfinite or negative value", #_V_);       \
} while (0)

#define MAYBE_EVAL_CALL(_T1_, _T0_)              \
do {                                             \
	if (_T1_ != _T0_) {                          \
		*pbetaArg = *pnuArg = *pmuArg = _T1_;    \
		TRY_EVAL_CALL(beta);                     \
		TRY_EVAL_CALL(nu);                       \
		TRY_EVAL_CALL(mu);                       \
		 _T0_ = _T1_;                            \
	}                                            \
} while (0)

SEXP R_adsir_initialize(SEXP beta, SEXP nu, SEXP mu, SEXP gamma, SEXP delta)
{
	F = allocVector(REALSXP, 7);
	R_PreserveObject(F);
	pF = REAL(F);
	memset(pF, 0, 7 * sizeof(double));

	DF = allocMatrix(REALSXP, 5, 7);
	R_PreserveObject(DF);
	pDF = REAL(DF);
	memset(pDF, 0, 35 * sizeof(double));

	INIT_CALL(beta);
	INIT_CALL(nu);
	INIT_CALL(mu);
	gammaVal = REAL(gamma)[0];
	deltaVal = REAL(delta)[0];

	lastTimeDot = lastTimeJac = -1.0;
	return R_NilValue;
}

SEXP R_adsir_finalize(void)
{
	R_ReleaseObject(F);
	R_ReleaseObject(DF);
	R_ReleaseObject(betaCall);
	R_ReleaseObject(betaArg);
	R_ReleaseObject(nuCall);
	R_ReleaseObject(nuArg);
	R_ReleaseObject(muCall);
	R_ReleaseObject(muArg);
	return R_NilValue;
}

SEXP R_adsir_dot(SEXP t, SEXP x)
{
	pt = REAL(t);
	px = REAL(x);
	MAYBE_EVAL_CALL(*pt, lastTimeDot);
	pF[0] = nuVal;
	pF[1] = betaVal * px[0] * px[1];
	pF[2] = (px[1] > 1.0) ? gammaVal * px[1] : 0.0;
	pF[3] = muVal * px[0];
	pF[4] = (px[1] > 1.0) ? muVal * px[1] : 0.0;
	pF[5] = muVal * px[2];
	pF[6] = deltaVal * px[2];
	return F;
}

SEXP R_adsir_jac(SEXP t, SEXP x)
{
	pt = REAL(t);
	px = REAL(x);
	MAYBE_EVAL_CALL(*pt, lastTimeJac);
	pDF[ 5] = betaVal * px[1];
	pDF[ 6] = betaVal * px[0];
	pDF[11] = gammaVal;
	pDF[15] = muVal;
	pDF[21] = muVal;
	pDF[27] = muVal;
	pDF[32] = deltaVal;
	return DF;
}

SEXP R_desir_initialize(SEXP beta, SEXP nu, SEXP mu, SEXP gamma, SEXP delta)
{
	INIT_CALL(beta);
	INIT_CALL(nu);
	INIT_CALL(mu);
	gammaVal = REAL(gamma)[0];
	deltaVal = REAL(delta)[0];

	lastTimeDot = lastTimeJac = -1.0;
	return R_NilValue;
}

SEXP R_desir_finalize(void)
{
	R_ReleaseObject(betaCall);
	R_ReleaseObject(betaArg);
	R_ReleaseObject(nuCall);
	R_ReleaseObject(nuArg);
	R_ReleaseObject(muCall);
	R_ReleaseObject(muArg);
	return R_NilValue;
}

/* vignette("compiledCode", package = "deSolve") */
void R_desir_dot(int *neq, double *t, double *y, double *ydot,
                 double *yout, int *ip)
{
	MAYBE_EVAL_CALL(*t, lastTimeDot);
	tmp = y[1];
	y[1] = exp(tmp);
	ydot[0] = nuVal - betaVal * y[0] * y[1] + deltaVal * y[2] - muVal * y[0];
	ydot[1] = betaVal * y[0] - gammaVal - muVal;
	ydot[2] = gammaVal * y[1] - deltaVal * y[2] - muVal * y[2];
	ydot[3] = nuVal;
	ydot[4] = betaVal * y[0] * y[1];
	y[1] = tmp;
	return;
}

/* vignette("compiledCode", package = "deSolve") */
void R_desir_jac(int *neq, double *t, double *y, int *ml,
                 int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{
	MAYBE_EVAL_CALL(*t, lastTimeJac);
	tmp = y[1];
	y[1] = exp(tmp);
	pd[ 0] = -betaVal * y[1] - muVal;
	pd[ 1] = betaVal;
	pd[ 4] = betaVal * y[1];
	pd[ 5] = -betaVal * y[0] * y[1];
	pd[ 7] = gammaVal * y[1];
	pd[ 9] = betaVal * y[0] * y[1];
	pd[10] = deltaVal;
	pd[12] = -deltaVal - muVal;
	y[1] = tmp;
	return;
}
