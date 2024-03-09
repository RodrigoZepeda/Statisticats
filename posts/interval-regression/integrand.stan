functions {
  // Lognormal probability density function
  real log_lognormal_pdf(real x, real mu_lognormal, real sigma_lognormal) {
    //1/(1 / (x * sigma * sqrt(2 * pi())) * exp(- 0.5*((log(x) - mu)/sigma)^2)
    
    //Values for logarithm. If logarithm is too small then we cancel
    real log_x = log(x);

    //This corresponds to log(1 / (x * sigma * sqrt(2 * pi())))
    real log_normalization = -(log_x + log(sigma_lognormal));
    
    //Exponent of the lognormal
    real exponent = -0.5*square((log_x - mu_lognormal) / sigma_lognormal);
    
    //For numerical stability
    return log_normalization + exponent;
  }
  
  // Normal probability density function
  real log_normal_pdf(real y, real mu, real sigma) {
    //1/(1 / (sqrt(2 * pi())*sigma) * exp(- 0.5*((x - mu)/sigma)^2)
    
    //Normal density
    real log_normalization = -log(sigma);
    real log_exponent      = -0.5 * square((y - mu) / sigma);

    return log_normalization + log_exponent;
  }
  
  real integrand_pdf(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    real mu  = theta[4] + theta[5]*x;
    return exp(log_lognormal_pdf(x, theta[2], theta[3]) + log_normal_pdf(x_r[1], mu, theta[1]));
  }
}
