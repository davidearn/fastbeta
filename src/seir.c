#include <math.h> /* exp */
#include <string.h> /* memset */
#include <Rinternals.h>
#include <R_ext/RS.h>

static SEXP FF, DF, betaCall, betaArg, nuCall, nuArg, muCall, muArg, s;
static double *pFF, *pDF, *pbetaArg, betaVal, *pnuArg, nuVal, *pmuArg, muVal,
	sigmaVal, gammaVal, deltaVal, lastTimeDot, lastTimeJac, *pt, *px,
	tmp, *work0, *work1, swork1, *work2, swork2;
static int m, n, p, nrow, ncol, ok;

#define INIT_CALL(_V_) \
do { \
	_V_ ## Call = allocVector(LANGSXP, 2); \
	R_PreserveObject(_V_ ## Call); \
	_V_ ## Arg = allocVector(REALSXP, 1); \
	R_PreserveObject(_V_ ## Arg); \
	SETCAR(_V_ ## Call, s_ ## _V_); \
	SETCADR(_V_ ## Call, _V_ ## Arg); \
	p ## _V_ ## Arg = REAL(_V_ ## Arg); \
} while (0)

#define EVAL_CALL(_V_) \
do { \
	s = eval(_V_ ## Call, R_GlobalEnv); \
	if (TYPEOF(s) != REALSXP) \
		error("'%s' did not evaluate to type \"%s\"", #_V_, "double"); \
	if (LENGTH(s) != 1) \
		error("'%s' did not evaluate to length %d", #_V_, 1); \
	_V_ ## Val = REAL(s)[0]; \
	if (!R_FINITE(_V_ ## Val) || _V_ ## Val < 0.0) \
		error("'%s' returned a nonfinite or negative value", #_V_); \
} while (0)

#define SETUP(_MODE_, _T0_, _T1_, _X_) \
do { \
	if (_T0_ != _T1_) { \
		*pbetaArg = *pnuArg = *pmuArg = _T1_; \
		EVAL_CALL(beta); \
		EVAL_CALL(nu); \
		EVAL_CALL(mu); \
		_T0_ = _T1_; \
		R_ ## _MODE_ ## _setup(_X_); \
	} \
} while (0)

SEXP R_adseir_initialize(SEXP s_beta, SEXP s_nu, SEXP s_mu,
                         SEXP s_sigma, SEXP s_gamma, SEXP s_delta,
                         SEXP s_m, SEXP s_n)
{
	m = INTEGER(s_m)[0];
	n = INTEGER(s_n)[0];
	p = m + n + 2;

	nrow = p + 2;
	ncol = p + p + 1;

	FF = allocVector(REALSXP, ncol);
	R_PreserveObject(FF);
	pFF = REAL(FF);
	memset(pFF, 0, LENGTH(FF) * sizeof(double));

	DF = allocMatrix(REALSXP, nrow, ncol);
	R_PreserveObject(DF);
	pDF = REAL(DF);
	memset(pDF, 0, LENGTH(DF) * sizeof(double));

	INIT_CALL(beta);
	INIT_CALL(nu);
	INIT_CALL(mu);
	sigmaVal = REAL(s_sigma)[0] * (double) m;
	gammaVal = REAL(s_gamma)[0] * (double) n;
	deltaVal = REAL(s_delta)[0] * (double) 1;

	pDF += (R_xlen_t) nrow * (2 + p) + 1;
	for (int i = 0; i < m; ++i) {
		*(pDF++) = sigmaVal;
		pDF += nrow;
	}
	for (int j = 0; j < n; ++j) {
		*(pDF++) = gammaVal;
		pDF += nrow;
	}
	*pDF = deltaVal;

	lastTimeDot = lastTimeJac = -1.0;
	return R_NilValue;
}

SEXP R_adseir_finalize(void)
{
	R_ReleaseObject(FF);
	R_ReleaseObject(DF);
	R_ReleaseObject(betaCall);
	R_ReleaseObject(betaArg);
	R_ReleaseObject(nuCall);
	R_ReleaseObject(nuArg);
	R_ReleaseObject(muCall);
	R_ReleaseObject(muArg);
	return R_NilValue;
}

static
void R_adseir_setup(const double *px)
{
	px++;

	swork1 = swork2 = 0.0;
	for (int i = 0; i < m; ++i)
		swork1 += *(px++);
	for (int j = 0; j < n; ++j)
		swork2 += *(px++);
	ok = swork1 + swork2 > 1.0;

	return;
}

SEXP R_adseir_dot(SEXP s_t, SEXP s_x)
{
	pt = REAL(s_t);
	px = REAL(s_x);
	SETUP(adseir, lastTimeDot, *pt, px);

	pFF = REAL(FF);
	*(pFF++) = betaVal * px[0] * swork2;
	*(pFF++) = nuVal;
	if (ok)
		for (int k = 0; k < p; ++k)
			*(pFF++) = muVal * px[k];
	else {
		memset(pFF, 0, p * sizeof(double));
		pFF[0    ] = muVal * px[0    ];
		pFF[p - 1] = muVal * px[p - 1];
		pFF += p;
	}
	for (int i = 0; i < m; ++i)
		*(pFF++) = sigmaVal * px[1 + i];
	if (ok)
		for (int j = 0; j < n; ++j)
			*(pFF++) = gammaVal * px[1 + m + j];
	else {
		memset(pFF, 0, n * sizeof(double));
		pFF += n;
	}
	*(pFF++) = deltaVal * px[1 + m + n];

	return FF;
}

SEXP R_adseir_jac(SEXP t, SEXP x)
{
	pt = REAL(t);
	px = REAL(x);
	SETUP(adseir, lastTimeDot, *pt, px);

	tmp = betaVal * px[0];

	pDF = REAL(DF);
	*(pDF++) = betaVal * swork2;
	pDF += m;
	for (int j = 0; j < n; ++j)
		*(pDF++) = tmp;

	pDF += 3 + nrow;
	for (int k = 0; k < p; ++k) {
		*(pDF++) = muVal;
		pDF += nrow;
	}

	return DF;
}

SEXP R_deseir_initialize(SEXP s_beta, SEXP s_nu, SEXP s_mu,
                         SEXP s_sigma, SEXP s_gamma, SEXP s_delta,
                         SEXP s_m, SEXP s_n)
{
	m = INTEGER(s_m)[0];
	n = INTEGER(s_n)[0];
	p = m + n + 2;

	nrow = p + 2;
	ncol = p + 2;

	INIT_CALL(beta);
	INIT_CALL(nu);
	INIT_CALL(mu);
	sigmaVal = REAL(s_sigma)[0] * (double) m;
	gammaVal = REAL(s_gamma)[0] * (double) n;
	deltaVal = REAL(s_delta)[0] * (double) 1;

	work0 = R_Calloc((size_t) m + n + n + n, double);
	work1 = work0 + m + n;
	work2 = work1     + n;

	lastTimeDot = lastTimeJac = -1.0;
	return R_NilValue;
}

SEXP R_deseir_finalize(void)
{
	R_ReleaseObject(betaCall);
	R_ReleaseObject(betaArg);
	R_ReleaseObject(nuCall);
	R_ReleaseObject(nuArg);
	R_ReleaseObject(muCall);
	R_ReleaseObject(muArg);
	R_Free(work0);
	return R_NilValue;
}

static
void R_deseir_setup(const double *px)
{
	px++;
	const double *py = px + 1;

	tmp = *px;

	swork1 = swork2 = 0.0;
	for (int i = 0; i < m; ++i)
		work0[i    ] = exp(*(px++) - *(py++));
	for (int j = 0; j < n; ++j) {
		swork1 += (work0[j] = exp(*px      ));
		swork2 += (work1[j] = exp(*px - tmp));
		work0[m + j] = exp(*(px++) - *(py++));
	}

	return;
}

/* vignette("compiledCode", package = "deSolve") */
void R_deseir_dot(const int *neq, const double *t, const double *y,
                  double *ydot, double *yout, const int *ip)
{
	SETUP(deseir, lastTimeDot, *t, y);

	*(ydot++) = nuVal - betaVal * y[0] * swork1 + deltaVal * y[p - 1] - muVal * y[0];
	if (m == 0)
	*(ydot++) = betaVal * y[0] * swork2 - gammaVal - muVal;
	else {
	*(ydot++) = betaVal * y[0] * swork2 - sigmaVal - muVal;
	for (int i = 0; i < m - 1; ++i)
	*(ydot++) = sigmaVal * work0[i    ] - sigmaVal - muVal;
	*(ydot++) = sigmaVal * work0[m - 1] - gammaVal - muVal;
	}
	for (int j = 0; j < n - 1; ++j)
	*(ydot++) = gammaVal * work0[m + j] - gammaVal - muVal;
	*(ydot++) = gammaVal * work0[m + n] - deltaVal - muVal;
	*(ydot++) = betaVal * y[0] * swork1;
	*(ydot++) = nuVal;

	return;
}

/* vignette("compiledCode", package = "deSolve") */
void R_deseir_jac(const int *neq, const double *t, const double *y,
                  const int *ml, const int *mu,
                  double *pd, const int *nrowpd, double *yout, const int *ip)
{
	SETUP(deseir, lastTimeDot, *t, y);

	pd[0] = -(pd[p] = betaVal * swork1) - muVal;
	pd[1] =           betaVal * swork2         ;
	pd += nrow + 2;
	for (int i = 0; i < m; ++i) {
		*(pd + nrow) = -(*pd = sigmaVal * work0[i    ]);
		*pd += nrow + 1;
	}
	for (int j = 0; j < n; ++j) {
		*(pd + nrow) = -(*pd = gammaVal * work0[m + j]);
		*pd += nrow + 1;
	}
	pd -= p;
	pd[0] = deltaVal * exp(y[p - 1]);
	pd -= (size_t) nrow * n;
	for (int j = 0; j < n; ++j) {
		pd[0] = -(pd[p] = betaVal * y[0] * work1[j]);
		pd[1] =           betaVal * y[0] * work2[j] ;
		pd += nrow;
	}
	pd -= (size_t) nrow * (m + n);
	pd[1] = (m == 0) ? 0.0 : -betaVal * swork2;

	return;
}
