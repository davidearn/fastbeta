#include <math.h>
#include <Rinternals.h>

static SEXP F, D, betaCall, betaArg, nuCall, nuArg, muCall, muArg;
static double *px, *pF, *pJ, *pbetaArg, betaVal, *pnuArg, nuVal, *pmuArg, muVal,
    gammaVal, lastTime, thisTime, tmp;

#define INIT_CALL(_V_)				\
    do {					\
	_V_ # Call = allocVector(LANGSXP, 2);	\
	R_PreserveObject(_V_ # Call);		\
	_V_ # Arg = allocVector(REALSXP, 1);	\
	R_PreserveObject(_V_ # Arg);		\
	SETCAR(_V_ # Call, _V_);		\
	SETCADR(_V_ # Call, _V_ # Arg);		\
	p # _V_ # Arg = REAL(_V_ # Arg);	\
    } while (0)

#define MAYBE_EVAL_CALL(_T_)					\
    do {							\
	thisTime = (_T_);					\
	if (thisTime != lastTime) {				\
	    *pbetaArg = *pnuArg = *pmuArg = thisTime;		\
	    betaVal = *REAL(eval(betaCall, R_GlobalEnv));	\
	    nuVal = *REAL(eval(  nuCall, R_GlobalEnv));		\
	    muVal = *REAL(eval(  muCall, R_GlobalEnv));		\
	    lastTime = thisTime;				\
	}							\
    } while (0)

void R_sir_ad_initialize(SEXP beta, SEXP nu, SEXP mu, SEXP gamma)
{
    F = allocVector(REALSXP, 6);
    R_PreserveObject(F);
    pF = REAL(F);
    memset(pF, 0, 6 * sizeof(double));
    
    J = allocMatrix(REALSXP, 4, 6);
    R_PreserveObject(J);
    pJ = REAL(J);
    memset(pJ, 0, 4 * 6 * sizeof(double));

    INIT_CALL(beta);
    INIT_CALL(nu);
    INIT_CALL(mu);
    gammaVal = *REAL(gamma);
    lastTime = -1.0;
    return;
}

SEXP R_sir_ad_rate(SEXP t, SEXP x)
{
    MAYBE_EVAL_CALL(*REAL(t));
    px = REAL(x);
    pF[0] = nuVal;
    pF[1] = betaVal * px[0] * px[1];
    pF[2] = (px[1] > 1.0) ? gammaVal * px[1] : 0.0;
    pF[3] = muVal * px[0];
    pF[4] = (px[1] > 1.0) ? muVal * px[1] : 0.0;
    pF[5] = muVal * px[2];
    return F;
}

SEXP R_sir_ad_jaco(SEXP t, SEXP x)
{
    MAYBE_EVAL_CALL(*REAL(t));
    px = REAL(x);
    pJ[ 0] = nuVal;
    pJ[ 4] = betaVal * px[1];
    pJ[ 5] = betaVal * px[0];
    pJ[ 9] = gammaVal;
    pJ[12] = pJ[17] = pJ[22] = muVal;
    return J;
}

void R_sir_de_initialize(SEXP beta, SEXP nu, SEXP mu, SEXP gamma)
{
    INIT_CALL(beta);
    INIT_CALL(nu);
    INIT_CALL(mu);
    gammaVal = REAL(gamma)[0];
    return;
}

void R_sir_de_rate(int *neq, double *t, double *y, double *ydot,
		   double *yout, int *ip)
{
    MAYBE_EVAL_CALL(*t);
    tmp = y[1];
    y[1] = exp(y[1]);
    ydot[0] = nuVal - (ydot[3] = betaVal * y[0] * y[1]) - muVal * y[0];
    ydot[1] = betaVal * y[0] - gammaVal - muVal;
    ydot[2] = gammaVal * y[1] - muVal * y[2];
    y[1] = tmp;
    return;
}

void R_sir_de_jaco(int *neq, double *t, double *y, int *ml,
		   int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{
    MAYBE_EVAL_CALL(*t);
    tmp = y[1];
    y[1] = exp(y[1]);
    pd[ 0] = -betaVal * y[1] - muVal;
    pd[ 1] = betaVal;
    pd[ 3] = betaVal * y[1];
    pd[ 4] = -betaVal * y[0] * y[1];
    pd[ 6] = gammaVal * y[1];
    pd[ 7] = betaVal * y[0] * y[1];
    pd[10] = -muVal;
    y[1] = tmp;
    return;
}
