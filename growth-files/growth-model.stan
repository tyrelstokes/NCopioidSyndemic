//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector<lower=0, upper =1>[N] y;
  vector[N] t;
  real t_cutoff;
}transformed data{
  vector[N] t_scale;
  real t_0;
  t_scale  = (t - min(t))/max(t);
  t_0 = (t_cutoff - min(t))/max(t);
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real k;
  real delta;
  real<lower=0,upper=1> L;
  real<lower=0,upper=1> x_0;
  
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[N] mu;
  real k_eff;
  for(i in 1:N){
    if(t_scale[i] >= t_0){
      k_eff = k + delta;
    }else{
      k_eff = k;
    }
    mu[i] = L/(1.0 + exp(k_eff*(t_scale[i]-x_0)));
  }
  
  
  y ~ normal(mu, sigma);
  
  L ~ normal(0.5,0.5);
  k ~ normal(0.0,0.25);
  delta ~ normal(0.0,0.1);
  x_0 ~ normal(0.5,0.6);
  
  sigma ~ normal(0.0,0.1);
}generated quantities{
  real p_delta;
  
  if(delta >= 0.0){
    p_delta = 1.0;
  }else{
    p_delta = 0.0;
  }
}

