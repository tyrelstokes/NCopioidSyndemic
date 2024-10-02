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
  //vector<lower=0, upper =1>[N] y;
  vector[N] y;
  vector[N] t;
  real t_cutoff;
  array[N] int state;
  int n_s;
  vector[n_s] y_0;
  real eps;
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
  vector[n_s] k_0;
 //vector[n_s] k;
  vector[n_s] delta_0;
 // vector[n_s] L_0;
  //vector[n_s] x_naught;
  
 // real<lower =0, upper = 1.0> x_0;
  
  real<lower=0> tau_k;
  real<lower=0> tau_delta;
//  real<lower=0> tau_L;
//  real<lower=0> tau_x;
  
  real<lower =0> alpha_k;
  real<lower =0> alpha_delta;
 // real<lower = 0> alpha_l;
  //real<lower=0> alpha_x;
  
  //real<lower=0> sigma;
  real<lower=0> invphi;
 // real<lower=0> y0;
  vector[n_s] alpha_0;
  real mu_alpha;
  real mu_b;
  vector[n_s] b_0;
  
  real<lower=0> tau_a;
  real<lower=0> tau_b;
  
}transformed parameters{
  
  vector[n_s] k;
  vector[n_s] delta;
  vector[n_s] alpha;
  vector[n_s] b;
  //vector[n_s] L;
 // vector[n_s] x_0;
  for(i in 1:n_s){
    k[i] = alpha_k + k_0[i]*tau_k;
    delta[i] = alpha_delta + delta_0[i]*tau_delta;
    alpha[i] = mu_alpha + alpha_0[i]*tau_a;
    b[i] = mu_b + b_0[i]*tau_b;
   // L[i] = inv_logit(alpha_l + L_0[i]*tau_L);
   // x_0[i] = inv_logit(alpha_x + x_naught[i]*tau_x);
  }
  
  
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[N] mu;
  vector[N] beta;
  vector[N] b_ln;
  real k_eff;
  for(i in 1:N){
    if(t_scale[i] >= t_0){
      k_eff = k[state[i]] + delta[state[i]];
    }else{
      k_eff = k[state[i]];
    }
   // mu[i] = L[state[i]]/(1.0 + exp(k_eff*(t_scale[i]-x_0[state[i]])));
   mu[i] = y_0[state[i]]*exp(k_eff*(t_scale[i]));
   beta[i] = invphi/mu[i];
   b_ln[i] = (alpha[state[i]] + mu[i]*b[state[i]]);
  }
  
  
  hurdle ~ bernoulli_logit(b_ln);
  
 // y ~ normal(mu, sigma)T[0.0,1.0];
 for(i in 1:N){
   if(hurdle[i] ==1){
 y[i] ~ gamma(invphi,beta[i])T[0.0,1.0];
   }
}
  
 // y0 ~ normal(0.0,0.01);
  
 // L_0 ~ normal(0.0,1.0);
  k_0 ~ normal(0.0,1.0);
  delta_0 ~ normal(0.0,1.0);
  //x_naught ~ normal(0.0,1.0);
  
  tau_k ~ normal(0.0,1.0);
  tau_delta ~ normal(0.0,0.5);
 // tau_L ~ normal(0.0,0.5);
 // tau_x ~ normal(0.0,0.5);
  
  //k ~ normal(alpha_k,0.25);
  
  alpha_k ~ normal(7.0,2.0);
  alpha_delta ~ normal(0.0,1.0);
 // alpha_l ~ normal(-1.0,0.5);
  //alpha_x ~ normal(0.0,0.5);
  
 
  
  //sigma ~ normal(0.0,0.02);
  invphi ~ normal(1.0,0.3);
  
  
  alpha_0 ~ normal(0.0,1.0);
  b_0 ~ normal(0.0,1.0);
  mu_alpha ~ normal(-2.0,0.5);
  mu_b ~ normal(2.0,0.5);
  
  tau_a ~ normal(0.0,0.5);
  tau_b ~ normal(0.0,0.5);


  
}generated quantities{
  vector[n_s] p_delta;
  
  for(i in 1:n_s){
  if(delta[i] >= 0.0){
    p_delta[i] = 1.0;
  }else{
    p_delta[i] = 0.0;
  }
  
  }
}
