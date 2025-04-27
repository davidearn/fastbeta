#include <TMB.hpp>

template<class Type>
vector<Type> rungekutta(vector<Type> (*f)(Type,
                                          const vector<Type> &,
                                          const vector<Type> &,
                                          const vector<int> &),
                        matrix<Type> &F,
                        Type t,
                        const vector<Type> &y,
                        const vector<Type> &dtheta,
                        const vector<int> &itheta,
                        const vector<Type> &a,
                        const vector<Type> &b,
                        const vector<Type> &c,
                        Type h)
{
	Eigen::Index i, j, k, m = F.rows(), n = F.cols(), p = 0;
	vector<Type> ans(m);
	for (j = 0; j < n; ++j) {
		for (i = 0; i < m; ++i)
			ans(i) = y(i);
		for (k = 0; k < j; ++k, ++p)
		for (i = 0; i < m; ++i)
			ans(i) += a(p) * h * F(i, k);
		F.col(j) = f(t + c(j) * h, ans, dtheta, itheta);
	}
	for (i = 0; i < m; ++i)
		ans(i) = y(i);
	for (j = 0; j < n; ++j)
	for (i = 0; i < m; ++i)
		ans(i) += b(j) * h * F(i, j);
	return ans;
}

template<class Type>
vector<Type> dot(Type t,
                 const vector<Type> &y,
                 const vector<Type> &dtheta,
                 const vector<int> &itheta)
{
	int m = itheta(0);
	int n = itheta(1);
	Type s1 = (y.segment(1 + m, n)       ).exp().sum();
	Type s2 = (y.segment(1 + m, n) - y(1)).exp().sum();
	vector<Type> dydt(1 + m + n + 1 + 1);
	dydt(0) = dtheta(1) + dtheta(3 + m + n) * exp(y(1 + m + n)) -
		(dtheta(0) * s1 + dtheta(2)) * y(0);
	dydt(1) = dtheta(0) * s2 * y(0) - (dtheta(3) + dtheta(2));
	for (int k = 0; k < m + n; ++k)
	dydt(2 + k) = dtheta(3 + k) * exp(y(1 + k) - y(2 + k)) -
		(dtheta(4 + k) + dtheta(2));
	dydt(2 + m + n) = dtheta(0) * s1 * y(0);
	return dydt;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
	/* .... Data .................................................... */

	/* Number of time points */
	DATA_INTEGER(T);

	/* Time series */
	DATA_IVECTOR(incidence);
	DATA_VECTOR(birth);
	DATA_VECTOR(death);

	/* Constants */
	DATA_SCALAR(intercept);
	DATA_SCALAR(sigma);
	DATA_SCALAR(gamma);
	DATA_SCALAR(delta);

	/* Number of E[i], I[j] compartments */
	DATA_INTEGER(m);
	DATA_INTEGER(n);

	/* State (S, E, I, R) at first time point */
	DATA_VECTOR(init);

	/* Nullity, rank */
#if 0
	DATA_INTEGER(R0);
#endif
	DATA_INTEGER(R1);

	/* Model matrices */
	DATA_MATRIX(X0);
	DATA_MATRIX(X1);

	/* Runge-Kutta coefficients */
	DATA_VECTOR(rka);
	DATA_VECTOR(rkb);
	DATA_VECTOR(rkc);


	/* .... Transformed data ........................................ */
	vector<Type> dtheta(3 + m + n + 1);
	for (int i = 0; i < m; ++i)
		dtheta(3 + i + 0 + 0) = m * sigma;
	for (int j = 0; j < n; ++j)
		dtheta(3 + m + j + 0) = n * gamma;
	for (int k = 0; k < 1; ++k)
		dtheta(3 + m + n + k) = 1 * delta;
	vector<int> itheta(2);
	itheta(0) = m;
	itheta(1) = n;


	/* .... Parameters .............................................. */

	PARAMETER(log_sd);
	PARAMETER(log_size);
	PARAMETER_VECTOR(b0);
	PARAMETER_VECTOR(b1);


	/* .... Transformed parameters .................................. */

	/* Spline */
	vector<Type> log_trans = intercept + X0 * b0 + X1 * b1;

	/* State (S, log(E), log(I), log(R), cumulative incidence) */
	matrix<Type> state(1 + m + n + 1 + 1, T);

	state(0, 0) = init(0);
	for (int i = 0; i < m; ++i)
		state(1 + i + 0 + 0, 0) = log(init(1 + i + 0 + 0));
	for (int j = 0; j < n; ++j)
		state(1 + m + j + 0, 0) = log(init(1 + m + j + 0));
	for (int k = 0; k < 1; ++k)
		state(1 + m + n + k, 0) = log(init(1 + m + n + k));
	state(1 + m + n + 1, 0) = Type(0.0);

	matrix<Type> F(state.rows(), rkb.size());
	vector<Type> y = state.col(0);
	for (int t = 1; t < T; ++t) {
		dtheta(0) = exp(log_trans(t - 1));
		dtheta(1) = birth(t - 1);
		dtheta(2) = death(t - 1);
		y = rungekutta(dot, F, Type(0.0), y, dtheta, itheta,
		               rka, rkb, rkc, Type(1.0));
		state.col(t) = y;
	}

	REPORT(log_trans);
	REPORT(state);
	ADREPORT(log_trans);
	ADREPORT(state);


	/* .... Model ................................................... */
	Type ans = Type(0.0);

	Type sd = exp(log_sd);
	for (int j = 0; j < R1; ++j)
		ans -= dnorm(b1(j), Type(0.0), sd, 1);

#if 0
	Type size = exp(log_size);
#endif
	for (int t = 1; t < T; ++t) {
		Type log_mu = log(state(1 + m + n + 1, t - 0) -
		                  state(1 + m + n + 1, t - 1));
		ans -= dnbinom_robust(Type(incidence(t)),
		                      log_mu, Type(2.0) * log_mu - log_size, 1);
	}


	return ans;
}
