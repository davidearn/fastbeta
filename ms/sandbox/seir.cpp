#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
	/* Incidence */
	DATA_IVECTOR(series1);
	/* Births */
	DATA_IVECTOR(series2);
	/* Natural mortality */
	DATA_VECTOR(series3);
	/* Rates out of E, I, R */
	DATA_SCALAR(sigma);
	DATA_SCALAR(gamma);
	DATA_SCALAR(delta);
	/* Number of E[i], I[j] subcompartments */
	DATA_INTEGER(m);
	DATA_INTEGER(n);
	/* Initial state (S, E, I, R, cumulative incidence) */
	DATA_VECTOR(init);

	/* Number of time points */
	DATA_INTEGER(T);
	/* Number of knots */
	DATA_INTEGER(K);
	/* Model matrix */
	DATA_MATRIX(X);
	/* Penalty matrix */
	DATA_MATRIX(S);
	/* Hyperparameter: negative binomial dispersion */
	PARAMETER(disp);

}
