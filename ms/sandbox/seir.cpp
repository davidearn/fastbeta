#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
	/* .... Data .................................................... */

	/* Number of time points */
	DATA_INTEGER(T);

	/* Time series */
	DATA_IVECTOR(incidence);
	DATA_IVECTOR(birth);
	DATA_VECTOR(death);

	/* Constants */
	DATA_SCALAR(sigma);
	DATA_SCALAR(gamma);
	DATA_SCALAR(delta);

	/* Number of E[i], I[j] compartments */
	DATA_INTEGER(m);
	DATA_INTEGER(n);

	/* State (S, E, I, R) at first time point */
	DATA_VECTOR(init);

	/* Nullity, rank */
	DATA_INTEGER(R0);
	DATA_INTEGER(R1);

	/* Model matrices */
	DATA_MATRIX(X0);
	DATA_MATRIX(X1);


	/* .... Transformed data ........................................ */

	vector<Type> a(m + n + 1);
	for (int i = 0; i < m; ++i)
		a(i + 0 + 0) = m * sigma;
	for (int j = 0; j < n; ++j)
		a(m + j + 0) = n * gamma;
	for (int k = 0; k < 1; ++k)
		a(m + n + k) = 1 * delta;


	/* .... Parameters .............................................. */

	PARAMETER(log_sd);
	PARAMETER(log_size);
	PARAMETER_VECTOR(b0);
	PARAMETER_VECTOR(b1);


	/* .... Transformed parameters .................................. */

	/* Random effect standard deviation */
	Type sd = exp(log_sd);

	/* Negative binomial dispersion parameter */
	Type size = exp(log_size);

	/* Spline */
	vector<Type> trans = gamma/init(0) * exp(X0 * b0 + X1 * b1);

	/* State (S, log(E), log(I), log(R), cumulative incidence) */
	matrix<Type> state(T, 1 + m + n + 1 + 1);
	state(0, 0) = init(0);
	for (int i = 0; i < m; ++i)
		state(0, 1 + i + 0 + 0) = log(init(1 + i + 0 + 0));
	for (int j = 0; j < n; ++j)
		state(0, 1 + m + j + 0) = log(init(1 + m + j + 0));
	for (int k = 0; k < 1; ++k)
		state(0, 1 + m + n + k) = log(init(1 + m + n + k));
	state(0, 1 + m + n + 1) = Type(0.0);

	/* RK step goes here */


	/* .... Model ................................................... */
	Type ans = Type(0.0);
	for (int j = 0; j < R1; ++j)
		ans -= dnorm(b1(j), Type(0.0), sd, 1);
	for (int t = 1; t < T; ++t) {
		Type log_mu = log(state(t - 0, 1 + m + n + 1) -
		                  state(t - 1, 1 + m + n + 1));
		Type log_var_minus_mu = Type(2.0) * log_mu - log_size;
		ans -= dnbinom_robust(incidence(t), log_mu, log_var_minus_mu, 1);
	}


	return ans;
}
