#include <math.h>
#include <Rinternals.h>

static SEXP F, D, betaCall, betaArg, nuCall, nuArg, muCall, muArg, s;
static double *px, *pF, *pJ, *pbetaArg, betaVal, *pnuArg, nuVal, *pmuArg, muVal,
    gammaVal, lastTimeRate, lastTimeJaco, tmp;

#define INIT_CALL(_V_)				\
    do {					\
	_V_ ## Call = allocVector(LANGSXP, 2);	\
	R_PreserveObject(_V_ ## Call);		\
	_V_ ## Arg = allocVector(REALSXP, 1);	\
	R_PreserveObject(_V_ ## Arg);		\
	SETCAR(_V_ ## Call, _V_);		\
	SETCADR(_V_ ## Call, _V_ ## Arg);	\
	p # _V_ ## Arg = REAL(_V_ ## Arg);	\
    } while (0)

#define TRY_EVAL_CALL(_V_, _VALUE_)
do {
    s = eval(_V_ ## Call, R_GlobalEnv);
    if (TYPEOF(s) != REALSXP)
	error("'%s' did not evaluate to type \"double\"", #_V_);
    if (LENGTH(s) != 1)
	error("'%s' did not evaluate to length 1", #_V_);
    _V_ ## Val = *REAL(s);
} while (0)

#define MAYBE_EVAL_CALL(_T1_, _T0_)				\
    do {							\
	if (_T1_ != _T0_) {					\
	    *pbetaArg = *pnuArg = *pmuArg = _T1_;		\
	    TRY_EVAL_CALL(beta);				\
	    TRY_EVAL_CALL(nu);					\
	    TRY_EVAL_CALL(mu);					\
	    _T0_ = _T1_;					\
	}							\
    } while (0)

SEXP R_adsir_initialize(SEXP beta, SEXP nu, SEXP mu, SEXP gamma)
{
    F = allocVector(REALSXP, 6);
    R_PreserveObject(F);
    pF = REAL(F);
    memset(pF, 0, 6 * sizeof(double));
    
    J = allocMatrix(REALSXP, 4, 6);
    R_PreserveObject(J);
    pJ = REAL(J);
    memset(pJ, 0, 24 * sizeof(double));

    INIT_CALL(beta);
    INIT_CALL(nu);
    INIT_CALL(mu);
    gammaVal = *REAL(gamma);
    
    lastTimeRate = lastTimeJaco = -1.0;
    return R_NilValue;
}

SEXP R_adsir_finalize(void)
{
    R_ReleaseObject(F);
    R_ReleaseObject(J);
    R_ReleaseObject(betaCall);
    R_ReleaseObject(betaArg);
    R_ReleaseObject(nuCall);
    R_ReleaseObject(nuArg);
    R_ReleaseObject(muCall);
    R_ReleaseObject(muArg);
    return R_NilValue;
}

SEXP R_adsir_rate(SEXP t, SEXP x)
{
    MAYBE_EVAL_CALL(*REAL(t), lastTimeRate);
    px = REAL(x);
    pF[0] = nuVal;
    pF[1] = betaVal * px[0] * px[1];
    pF[2] = (px[1] > 1.0) ? gammaVal * px[1] : 0.0;
    pF[3] = muVal * px[0];
    pF[4] = (px[1] > 1.0) ? muVal * px[1] : 0.0;
    pF[5] = muVal * px[2];
    return F;
}

SEXP R_adsir_jaco(SEXP t, SEXP x)
{
    MAYBE_EVAL_CALL(*REAL(t), lastTimeJaco);
    px = REAL(x);
    pJ[ 0] = nuVal;
    pJ[ 4] = betaVal * px[1];
    pJ[ 5] = betaVal * px[0];
    pJ[ 9] = gammaVal;
    pJ[12] = pJ[17] = pJ[22] = muVal;
    return J;
}

SEXP R_desir_initialize(SEXP beta, SEXP nu, SEXP mu, SEXP gamma)
{
    INIT_CALL(beta);
    INIT_CALL(nu);
    INIT_CALL(mu);
    gammaVal = REAL(gamma)[0];

    lastTimeRate = lastTimeJaco = -1.0;
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
void R_desir_rate(int *neq, double *t, double *y, double *ydot,
		  double *yout, int *ip)
{
    MAYBE_EVAL_CALL(*t, lastTimeRate);
    tmp = y[1];
    y[1] = exp(tmp);
#if 0
    ydot[0] = nuVal - betaVal * y[0] * y[1] - muVal * y[0];
    ydot[1] = betaVal * y[0] - gammaVal - muVal;
    ydot[2] = gammaVal * y[1] - muVal * y[2];
    ydot[3] = betaVal * y[0] * y[1];
#else
    ydot[3] = betaVal * y[0] * y[1];
    ydot[0] = nuVal - ydot[3] - muVal * y[0];
    ydot[1] = betaVal * y[0] - gammaVal - muVal;
    ydot[2] = gammaVal * y[1] - muVal * y[2];
#endif
    y[1] = tmp;
    return;
}

/* vignette("compiledCode", package = "deSolve") */
void R_desir_jaco(int *neq, double *t, double *y, int *ml,
		  int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{
    MAYBE_EVAL_CALL(*t, lastTimeJaco);
    tmp = y[1];
    y[1] = exp(tmp);
#if 0
    pd [0] = -betaVal * y[1] - muVal;
    pd [1] = betaVal;
    pd [3] = betaVal * y[1];
    pd [4] = -betaVal * y[0] * y[1];
    pd [6] = gammaVal * y[1];
    pd [7] = betaVal * y[0] * y[1];
    pd[10] = -muVal;
#else
    pd [1] = betaVal;
    pd [3] = betaVal * y[1];
    pd [0] = -pd[3] - muVal;
    pd [7] = pd[3] * y[0];
    pd [4] = -pd[7];
    pd [6] = gammaVal * y[1];
    pd[10] = -muVal;
#endif
    y[1] = tmp;
    return;
}
