data {
  int<lower=1> N; //Total sample size
  vector<lower=0>[N] X; //Observed X-values
  vector[N] Y; //Observed Y-values
}

parameters {
  //Mean parameter of X
  real<lower=0> mu;
  
  //Parameters for regression
  real beta_0;
  real beta_1;
  
  //Parameters for variance
  real<lower=0> sigma;
  real<lower=0> tau;
}

transformed parameters {
  vector[N] y_mean = rep_vector(beta_0, N) + beta_1*X;
}

model {
  beta_0 ~ normal(0, 100);
  beta_1 ~ normal(0, 100);
  sigma  ~ normal(0, 100);
  tau    ~ normal(0, 100);
  X      ~ lognormal(mu, tau);
  Y      ~ normal(y_mean, sigma);
}
