// poisson_regression.stan
data {
  int<lower=0> N;            // Number of observations
  int<lower=0> y[N];         // Count outcome
  int<lower=1> K;            // Number of predictors
  matrix[N, K] X;            // Predictor matrix
}
parameters {
  vector[K] beta;            // Coefficients
}
model {
  beta ~ normal(0, 10);      // Prior on coefficients
  y ~ poisson_log(X * beta); // Poisson likelihood with log link
}
