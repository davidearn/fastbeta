#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

SEXP R_fastbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_ptpi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallMethods[] =
{
    { "R_fastbeta", (DL_FUNC) &R_fastbeta, 6 },
    { "R_ptpi", (DL_FUNC) &R_ptpi, 9 },
    { NULL, NULL, 0 }
};

void attribute_visible R_init_fastbeta(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallMethods, NULL, NULL);
}
