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
  real eps;
  //vector<lower=0, upper =1>[N] y;
  vector[N] y;
  vector[N] t;
  real t_cutoff;
 // array[N] int state;
 // int n_s;
}transformed data{
  vector[N] t_scale;
  array[N] int hurdle;
  real t_0;
  t_scale  = (t - min(t))/(max(t)-min(t));
  t_0 = (t_cutoff - min(t))/(max(t)-min(t));
  
  
    
    for(i in 1:N){
    if(y[i] > eps){
      hurdle[i] = 1;
    }else{
      hurdle[i] = 0;
    }
  }
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
 // vector[n_s] k_0;
 //vector[n_s] k;
 // vector[n_s] delta_0;
 // vector[n_s] L_0;
  //vector[n_s] x_naught;
  
  real k;
  real delta;
  
 // real<lower =0, upper = 1.0> x_0;
  
//  real<lower=0> tau_k;
//  real<lower=0> tau_delta;
//  real<lower=0> tau_L;
//  real<lower=0> tau_x;
  
//  real<lower =0> alpha_k;
//  real<lower =0> alpha_delta;
 // real<lower = 0> alpha_l;
  //real<lower=0> alpha_x;
  
  //real<lower=0> sigma;
  real<lower=0> invphi;
 // real<lower=0> y0;
 
   real alpha;
  real b;
}transformed parameters{
  
//  vector[n_s] k;
//  vector[n_s] delta;
  //vector[n_s] L;
 // vector[n_s] x_0;
 // for(i in 1:n_s){
 //   k[i] = alpha_k + k_0[i]*tau_k;
 //   delta[i] = alpha_delta + delta_0[i]*tau_delta;
   // L[i] = inv_logit(alpha_l + L_0[i]*tau_L);
   // x_0[i] = inv_logit(alpha_x + x_naught[i]*tau_x);
//  }
  
  
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[N] mu;
  vector[N] beta;
  real k_eff;
  for(i in 1:N){
    if(t_scale[i] >= t_0){
      k_eff = k + delta;
    }else{
      k_eff = k;
    }
   // mu[i] = L[state[i]]/(1.0 + exp(k_eff*(t_scale[i]-x_0[state[i]])));
   mu[i] = y[1]*exp(k_eff*(t_scale[i]));
   beta[i] = invphi/mu[i];
  }
  
  
 // y ~ normal(mu, sigma)T[0.0,1.0];
 
 hurdle ~ bernoulli_logit(alpha + mu*b);
 
 for(i in 1:N){
   if(hurdle[i]==1){
 y[i] ~ gamma(invphi,beta[i])T[0.0,1.0];
   }
 }
  
  
  //y0 ~ normal(0.0,0.01);
  
 // L_0 ~ normal(0.0,1.0);
 // k_0 ~ normal(0.0,1.0);
 // delta_0 ~ normal(0.0,1.0);
  //x_naught ~ normal(0.0,1.0);
  
 // tau_k ~ normal(0.0,1.0);
 // tau_delta ~ normal(0.0,0.5);
 // tau_L ~ normal(0.0,0.5);
 // tau_x ~ normal(0.0,0.5);
  
  k ~ normal(7.0,2.0);
  
 // alpha_k ~ normal(7.0,2.0);
 // alpha_delta ~ normal(0.0,0.5);
 // alpha_l ~ normal(-1.0,0.5);
  //alpha_x ~ normal(0.0,0.5);
  
 
  
  //sigma ~ normal(0.0,0.02);
  invphi ~ normal(1.0,0.2);
  
  alpha ~ normal(-2.0, 1.0);
  beta ~ normal(1.0,1.0);
  
 // x_0 ~ normal(0.8,0.2);
  
}generated quantities{
  real p_delta;
  

  if(delta >= 0.0){
    p_delta = 1.0;
  }else{
    p_delta = 0.0;
  }
  

}
