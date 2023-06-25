#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

SEXP R_fastbeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_ptpi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_adsir_initialize(SEXP, SEXP, SEXP, SEXP);
SEXP R_adsir_finalize(void);
SEXP R_adsir_dot(SEXP, SEXP);
SEXP R_adsir_jac(SEXP, SEXP);

SEXP R_desir_initialize(SEXP, SEXP, SEXP, SEXP);
SEXP R_desir_finalize(void);
void R_desir_dot(int *, double *, double *, double *, double *, int *);
void R_desir_jac(int *, double *, double *, int *, int *, double *, int *,
                 double *, int *);

static const R_CallMethodDef CallMethods[] =
{
	{ "R_fastbeta",         (DL_FUNC) &R_fastbeta,         6 },
	{ "R_ptpi",             (DL_FUNC) &R_ptpi,             9 },

	{ "R_adsir_initialize", (DL_FUNC) &R_adsir_initialize, 4 },
	{ "R_adsir_finalize",   (DL_FUNC) &R_adsir_finalize,   0 },
	{ "R_adsir_dot",        (DL_FUNC) &R_adsir_dot,        2 },
	{ "R_adsir_jac",        (DL_FUNC) &R_adsir_jac,        2 },

	{ "R_desir_initialize", (DL_FUNC) &R_desir_initialize, 4 },
	{ "R_desir_finalize",   (DL_FUNC) &R_desir_finalize,   0 },

	{ NULL, NULL, 0 }
};

static R_NativePrimitiveArgType
	type0[] = { INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP },
	type1[] = { INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP };

static const R_CMethodDef CMethods[] =
{
	{ "R_desir_dot",        (DL_FUNC) &R_desir_dot,        6, type0 },
	{ "R_desir_jac",        (DL_FUNC) &R_desir_jac,        9, type1 },

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
