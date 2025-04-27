functions {
	vector dot(real t, vector y,
	           real beta, data real nu, data real mu,
	           data vector dtheta, data array[] int itheta)
	{
		int m = itheta[1];
		int n = itheta[2];
		real s1 = sum(exp(y[(1+m+1):(1+m+n)]       ));
		real s2 = sum(exp(y[(1+m+1):(1+m+n)] - y[2]));
		vector[1+m+n+1+1] dydt;
		dydt[1] = nu + dtheta[m+n+1] * exp(y[1+m+n+1]) -
			(beta * s1 + mu) * y[1];
		dydt[2] = beta * s2 * y[1] - (dtheta[1] + mu);
		for (k in 1:(m + n))
		dydt[2+k] = dtheta[0+k] * exp(y[1+k] - y[2+k]) -
			(dtheta[1+k] + mu);
		dydt[2+m+n+1] = beta * s1 * y[1];
		return dydt;
	}
}

data {
	/* Number of time points */
	int<lower=2> T;

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

	/* State (S, E, I, R) at first time point */
	vector<lower=0.0>[1+m+n+1] init;

	/* Nullity, rank */
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
	vector[T] trans = exp(intercept + X0 * b0 + X1 * b1);

	/* State (S, log(E), log(I), log(R), cumulative incidence) */
	matrix[1+m+n+1+1, T] state;

	state[:, 1] =
	append_row(append_row(init[1], log(init[(1+1):(1+m+n+1)])), 0.0);

	for (t in 2:T)
		state[:, t] =
		ode_rk45(dot, state[:, t-1], 0.0, times,
		         trans[t-1], birth[t-1], death[t-1],
		         dtheta, itheta)[1];
}

model {
	real sd = exp(log_sd);
	for (j in 1:R1)
		b1[j] ~ normal(0.0, sd);

	real size = exp(log_size);
	for (t in 2:T) {
		real mu = state[1+m+n+1+1, t] - state[1+m+n+1+1, t-1];
		incidence[t] ~ neg_binomial_2(mu, size);
	}
}

/*
generated quantities {

}
*/
