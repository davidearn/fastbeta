functions {
	vector dot(real t, vector y,
	           real beta, data real nu, data real mu,
	           data array[] real a, data int m, data int n) {
		vector[1+m+n+1+1] dydt;
		real s1 = sum(exp(y[(1+m+1):(1+m+n)]       ));
		real s2 = sum(exp(y[(1+m+1):(1+m+n)] - y[2]));
		dydt[1] = nu + a[m+n+1] * exp(y[1+m+n+1]) - (beta * s1 + mu) * y[1];
		dydt[2] = beta * s2 * y[1] - (a[1] + mu);
		dydt[3:(1+m+n+1)] = a[1:(m+n)] * exp(y[2:(1+m+n)] - y[3:(1+m+n+1)]) - (a[2:(m+n+1)] + mu);
		dydt[1+m+n+1+1] = beta * s1 * y[1];
		return dydt;
	}
}

data {
	/* Number of time points */
	int<lower=2> T;
	/* Incidence */
	array[T] int<lower=0> series1;
	/* Births */
	array[T] int<lower=0> series2;
	/* Natural mortality */
	array[T] real<lower=0.0> series3;
	/* Rates out of E, I, R */
	real<lower=0.0> sigma;
	real<lower=0.0> gamma;
	real<lower=0.0> delta;
	/* Number of E[i], I[j] subcompartments */
	int<lower=0> m;
	int<lower=1> n;
	/* Initial state (S, E, I, R, cumulative incidence) */
	vector[1+m+n+1+1] init;
	/* Nullity, rank */
	int<lower=1> R0;
	int<lower=1> R1;
	/* Model matrix */
	matrix[T, R0] X0;
	matrix[T, R1] X1;
}

transformed data {
	/* rep(c(sigma, gamma, delta), c(m, n, 1)) */
	array[m+n+1] real<lower=0.0> a;
	for (i in 1:m) {
		a[i+0+0] = sigma;
	}
	for (j in 1:n) {
		a[m+j+0] = gamma;
	}
	for (k in 1:1) {
		a[m+n+k] = delta;
	}
	/* Unit time step */
	array[1] real step = { 1.0 };
}

parameters {
	/* Coefficients */
	vector[R0] b0;
	vector[R1] b1;
	/* Random effect standard deviation */
	real<lower=0.0> sd;
	/* Negative binomial dispersion */
	real<lower=0.0> dispersion;
}

transformed parameters {
	/* Transmission */
	vector[T] real<lower=0.0> beta = exp(X0 * b0 + X1 * b1);
	/* List of state vectors */
	array[T] array[1] vector mu;
	mu[1][1] = init;
	for (t in 2:T) {
		mu[t] = ode_rk45(dot, mu[t-1], 0.0, step,
		                 beta[t-1], series2[t-1], series3[t-1], a, m, n);
	}
}

model {
	for (j in 1:R1)
		b1[j] ~ normal(0.0, sd);
	for (t in 2:T)
		series1[t] ~ neg_binomial_2(mu[t], dispersion);
}

/*
generated quantities {

}
*/
