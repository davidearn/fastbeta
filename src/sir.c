#include <math.h>
#include <Rinternals.h>

static SEXP F, J, betaCall, betaArg, nuCall, nuArg, muCall, muArg;
static double *pF, *pJ, *pbetaArg, *pnuArg, *pmuArg, gamma0;

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

    betaCall = allocVector(LANGSXP, 2);
    R_PreserveObject(betaCall);
    betaArg = allocVector(REALSXP, 1);
    R_PreserveObject(betaArg);
    SETCAR(betaCall, beta);
    SETCADR(betaCall, betaArg);
    pbetaArg = REAL(betaArg);
    
    nuCall = allocVector(LANGSXP, 2);
    R_PreserveObject(nuCall);
    nuArg = allocVector(REALSXP, 1);
    R_PreserveObject(nuArg);
    SETCAR(nuCall, nu);
    SETCADR(nuCall, nuArg);
    pnuArg = REAL(nuArg);
    
    muCall = allocVector(LANGSXP, 2);
    R_PreserveObject(muCall);
    muArg = allocVector(REALSXP, 1);
    R_PreserveObject(muArg);
    SETCAR(muCall, mu);
    SETCADR(muCall, muArg);
    pmuArg = REAL(muArg);

    gamma0 = REAL(gamma)[0];

    return;
}

SEXP R_sir_ad_rate(SEXP t, SEXP x)
{

}

SEXP R_sir_ad_jaco(SEXP t, SEXP x)
{
    
}

void R_sir_de_initialize(SEXP beta, SEXP nu, SEXP mu, SEXP gamma)
{

}

void R_sir_de_rate(int *neq, double *t, double *y, double *ydot,
		   double *yout, int *ip)
{

}

void R_sir_de_jaco(int *neq, double *t, double *y, int *ml,
		   int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{
    
}
