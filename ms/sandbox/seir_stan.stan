functions {
	vector dot(real t, vector y,
	           real log_beta, data real nu, data real mu,
	           data vector dtheta, data array[] int itheta)
	{
		int m = itheta[1];
		int n = itheta[2];
		real log_lambda = log_beta + log_sum_exp(y[(1+m+1):(1+m+n)]);
		vector[1+m+n+1+1] dydt;
		dydt[1] = nu * exp(-y[1]) +
			dtheta[m+n+1] * exp(y[1+m+n+1] - y[1]) -
			(exp(log_lambda) + mu);
		dydt[2] = exp(log_lambda + y[1] - y[2]) - (dtheta[1] + mu);
		dydt[(2+1):(2+m+n)] =
			dtheta[(0+1):(0+m+n)] .* exp(y[(1+1):(1+m+n)] - y[(2+1):(2+m+n)]) -
			(dtheta[(1+1):(1+m+n)]);
		dydt[2+m+n+1] = exp(log_lambda + y[1]);
		return dydt;
	}
}

data {
	/* Number of time points */
	int<lower=1> T;

	/* Time series */
	array[T] int<lower=0> incidence;
	array[T] real<lower=0.0> birth;
	array[T] real<lower=0.0> death;

	/* Constants */
	real intercept;
	real<lower=0.0> sigma;
	real<lower=0.0> gamma;
	real<lower=0.0> delta;

	/* Number of E[i], I[j] compartments */
	int<lower=0> m;
	int<lower=1> n;

	/* State (log S, log E, log I, log R, cumulative incidence) */
	vector[1+m+n+1+1] init;

	/* Nullity, rank of spline penalty matrix */
	int<lower=1> R0;
	int<lower=1> R1;

	/* Model matrices */
	matrix[T, R0] X0;
	matrix[T, R1] X1;
}

transformed data {
	vector[m+n+1] dtheta;
	for (i in 1:m)
		dtheta[i+0+0] = m * sigma;
	for (j in 1:n)
		dtheta[m+j+0] = n * gamma;
	for (k in 1:1)
		dtheta[m+n+k] = 1 * delta;
	array[2] int itheta;
	itheta[1] = m;
	itheta[2] = n;
	array[1] real times = { 1.0 };
}

parameters {
	real log_sd;
	real log_size;
	vector[R0] b0;
	vector[R1] b1;
}

transformed parameters {
	/* Spline */
	vector[T] log_trans = intercept + X0 * b0 + X1 * b1;

	/* State (log S, log E, log I, log R, cumulative incidence) */
	matrix[1+m+n+1+1, 1+T] state;

	state[:, 1] = init;
	for (t in 1:T)
		state[:, 1+t] =
		ode_rk45(dot, state[:, t], 0.0, times,
		         log_trans[t], birth[t], death[t],
		         dtheta, itheta)[1];
}

model {
	b1 ~ normal(0.0, exp(log_sd));
	incidence ~ neg_binomial_2(state[1+m+n+1+1, (1+1):(1+T)] -
	                           state[1+m+n+1+1, (0+1):(0+T)],
	                           exp(log_size));
}

/*
generated quantities {

}
*/
