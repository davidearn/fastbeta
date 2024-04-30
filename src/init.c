#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

SEXP R_fastbeta(SEXP, SEXP);
SEXP R_ptpi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_adseir_initialize(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_adseir_finalize(void);
SEXP R_adseir_dot(SEXP, SEXP);
SEXP R_adseir_jac(SEXP, SEXP);

SEXP R_deseir_initialize(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_deseir_finalize(void);
void R_deseir_dot(int *, double *, double *, double *, double *, int *);
void R_deseir_jac(int *, double *, double *, int *, int *, double *, int *,
                  double *, int *);

static const R_CallMethodDef CallMethods[] =
{
	{ "R_fastbeta",          (DL_FUNC) &R_fastbeta,          2 },
	{ "R_ptpi",              (DL_FUNC) &R_ptpi,              8 },

	{ "R_adseir_initialize", (DL_FUNC) &R_adseir_initialize, 8 },
	{ "R_adseir_finalize",   (DL_FUNC) &R_adseir_finalize,   0 },
	{ "R_adseir_dot",        (DL_FUNC) &R_adseir_dot,        2 },
	{ "R_adseir_jac",        (DL_FUNC) &R_adseir_jac,        2 },

	{ "R_deseir_initialize", (DL_FUNC) &R_deseir_initialize, 8 },
	{ "R_deseir_finalize",   (DL_FUNC) &R_deseir_finalize,   0 },

	{ NULL, NULL, 0 }
};

static R_NativePrimitiveArgType
	type0[] = { INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP },
	type1[] = { INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP };

static const R_CMethodDef CMethods[] =
{
	{ "R_deseir_dot",        (DL_FUNC) &R_deseir_dot,        6, type0 },
	{ "R_deseir_jac",        (DL_FUNC) &R_deseir_jac,        9, type1 },

	{ NULL, NULL, 0, NULL }
};

void attribute_visible R_init_fastbeta(DllInfo *info)
{
	R_registerRoutines(info, CMethods, CallMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
#if 0
	/* No, because deSolve tests is.loaded(<character string>) : */
	R_forceSymbols(info, TRUE);
#endif
	return;
}
