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
  int<lower=0> n_time_points;
  array[n_time_points] potential_change_point;
  vector[N] y;
}transformed data{
  
int n_changes;

n_changes = 0;

for(i in 1:n_time_points){
  n_changes += potential_change_point[i];
}
  
array[n_changes] int changes;  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real alpha;
  vector[n_changes] beta;
  vector<lower =0,upper =1>theta;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  real<lower =0.0,upper = 1.0> theta_eff;
  theta_eff = 0.0;
  log_theta_eff;
  
  vector[N] lp;
  
  for(i in 1:n_changes){
    theta_eff += (1.0 - theta[i]);
  }
  
  log_theta_eff = log(theta_eff);
  target += log_sum_exp(log_theta_eff,normal_lpdf())
}

