#include integrand.stan

data {
  int<lower=1> N;           //Total sample size
  vector<lower=0>[N] X_low; //Observed lower bound for X
  vector<lower=0>[N] X_up;  //Observed upper bound for X
  vector[N] Y;              //Observed Y-values
}


transformed data {
  //Normalize Y to help integration process
  real mu_Y = mean(Y);
  real sigma_Y = sd(Y);
  vector[N] y_centered = (Y - mu_Y)/sigma_Y;
  
  //Required for integration
  array[0] int x_i;
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

model {
  beta_0 ~ normal(0, 100);
  beta_1 ~ normal(0, 100);
  sigma  ~ normal(0, 100);
  tau    ~ normal(0, 100);
  
  //Add the likelihood for intervals
  for (i in 1:N){
    target += log(integrate_1d(
      integrand_pdf, X_low[i], X_up[i], 
        {sigma, mu, tau, beta_0, beta_1}, {y_centered[i]}, x_i,
        pow(machine_precision(), 0.25))
    );
  }
}

generated quantities {
  real beta_0_real = sigma_Y*beta_0 + mu_Y;
  real beta_1_real = sigma_Y*beta_1;
  real sigma_real  = sigma_Y*sigma;
}
