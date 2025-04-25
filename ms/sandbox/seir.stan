functions {
	vector dot(real t, vector y,
	           real beta, data real nu, data real mu,
	           data vector a, data int m, data int n)
	{
		vector[1+m+n+1+1] dydt;
		real s1 = sum(exp(y[(1+m+1):(1+m+n)]         ));
		real s2 = sum(exp(y[(1+m+1):(1+m+n)] - y[1+1]));
		dydt[1] = nu + a[m+n+1] * exp(y[1+m+n+1]) - (beta * s1 + mu) * y[1];
		dydt[2] = beta * s2 * y[1] - (a[1] + mu);
		dydt[(1+2):(1+m+n+1)] = a[1:(m+n)] .* exp(y[(1+1):(1+m+n)] - y[(1+2):(1+m+n+1)]) - (a[2:(m+n+1)] + mu);
		dydt[1+m+n+1+1] = beta * s1 * y[1];
		return dydt;
	}
}

data {
	/* Number of time points */
	int<lower=2> T;
	/* Time series */
	array[T] int<lower=0> incidence;
	array[T] int<lower=0> birth;
	array[T] real<lower=0.0> death;
	/* Constants */
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
	/* Model matrix */
	matrix[T, R0] X0;
	matrix[T, R1] X1;
}

transformed data {
	vector<lower=0.0>[m+n+1] a;
	for (i in 1:m)
		a[i+0+0] = m * sigma;
	for (j in 1:n)
		a[m+j+0] = n * gamma;
	for (k in 1:1)
		a[m+n+k] = 1 * delta;
	array[1] real step = { 1.0 };
}

parameters {
	/* Coefficients */
	vector[R0] b0;
	vector[R1] b1;
	/* Random effect standard deviation */
	real<lower=0.0> sd;
	/* Negative binomial dispersion parameter */
	real<lower=0.0> dispersion;
}

transformed parameters {
	/* Spline */
	vector<lower=0.0>[T] beta = gamma/init[1] * exp(X0 * b0 + X1 * b1);
	/* State (S, log(E), log(I), log(R), cumulative incidence) */
	matrix[T, 1+m+n+1+1] state;
	state[1, :] = append_row(append_row(init[1], log(init[(1+1):(1+m+n+1)])), 0.0)';
	for (t in 2:T)
	state[t, :] = ode_rk45(dot, state[t-1, :]', 0.0, step,
	                       beta[t-1], birth[t-1], death[t-1], a,
	                       m, n)[1]';
}

model {
	for (j in 1:R1)
		b1[j] ~ normal(0.0, sd);
	for (t in 2:T) {
		real location = state[t, 1+m+n+1+1] - state[t-1, 1+m+n+1+1];
		incidence[t] ~ neg_binomial_2(location, dispersion);
	}
}

/*
generated quantities {

}
*/
