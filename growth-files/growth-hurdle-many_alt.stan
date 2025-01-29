data {
  int<lower=0> N;
  vector[N] y;
  vector[N] t;
  real t_cutoff;
  real eps;
  array[N] int state;
  int n_s;
  int pred_ind;
} 

transformed data {
  array[N] int hurdle;
  vector[N] t_scale;
  real t_0;
  t_scale  = (t - min(t))/(max(t)-min(t));
  t_0 = (t_cutoff - min(t))/(max(t)-min(t));
  
  for(i in 1:N){
    if(y[i] > eps){
      hurdle[i] = 1;
    } else {
      hurdle[i] = 0;
    }
  }
}

parameters {
  vector[n_s] k_0;
  vector[n_s] delta_0;
  vector[n_s] L_0;
  vector[n_s] x_naught;
  
  real<lower=0> tau_k;
  real<lower=0> tau_delta;
  real<lower=0> tau_L;
  real<lower=0> tau_x;
  
  real<lower = 0> alpha_k;
  real<lower = 0> alpha_delta;
  real<lower = 0> alpha_l;
  real<lower = 0> alpha_x;
  
  real<lower = 0> invphi;
  real alpha;
  vector[n_s] beta;
}

transformed parameters {
  vector[n_s] k;
  vector[n_s] delta;
  vector[n_s] L;
  vector[n_s] x_0;
  for(i in 1:n_s){
    k[i] = alpha_k + k_0[i]*tau_k;
    delta[i] = alpha_delta + delta_0[i]*tau_delta;
    L[i] = inv_logit(alpha_l + L_0[i]*tau_L);
    x_0[i] = inv_logit(alpha_x + x_naught[i]*tau_x);
  }
}

model {
  vector[N] mu;
  vector[N] betas;
  vector[N] mu_ber;
  real k_eff;

  for(i in 1:N){
    if(t_scale[i] >= t_0){
      k_eff = k[state[i]] + delta[state[i]];
    } else {
      k_eff = k[state[i]];
    }
    mu[i] = L[state[i]]*inv_logit(k_eff*(t_scale[i]-x_0[state[i]]));
    betas[i] = invphi / mu[i];
    mu_ber[i] = alpha + t_scale[i] * beta[state[i]];
  }
  
  hurdle ~ bernoulli_logit(mu_ber);
  
  for(i in 1:N){
    if(hurdle[i] == 1){
      y[i] ~ gamma(invphi, betas[i])T[0.0, 1.0];
    }
  }
  
  beta ~ normal(0.0, 1.0);
  alpha ~ normal(0.0, 1.0);
  x_naught ~ normal(0.0, 1.0);
  
  tau_k ~ normal(0.0, 1.0);
  tau_delta ~ normal(0.0, 1.0);
  tau_L ~ normal(0.0, 0.5);
  tau_x ~ normal(0.0, 0.5);
  
  k ~ normal(alpha_k, 0.25);
  alpha_k ~ normal(5.0, 2.0);
  alpha_delta ~ normal(0.0, 0.5);
  alpha_l ~ normal(0.0, 0.5);
  alpha_x ~ normal(0.0, 0.5);
  
  invphi ~ normal(0.0, 1.0);
}

generated quantities {
  vector[n_s] p_delta;
  vector[N] log_lik;  // Log-likelihood for each observation
  vector[N] mu;
  vector[N] betas;
  vector[N] mu_ber;  // Declare as a vector to match the usage
  vector[N] pred;    // Declare pred as a vector, not inside the loop
  real k_eff;

  for(i in 1:N){
    if(t_scale[i] >= t_0){
      k_eff = k[state[i]] + delta[state[i]];
    } else {
      k_eff = k[state[i]];
    }
    mu[i] = L[state[i]] * inv_logit(k_eff * (t_scale[i] - x_0[state[i]]));
    betas[i] = invphi / mu[i];
    mu_ber[i] = alpha + t_scale[i] * beta[state[i]];
  }

  // Calculate log-likelihood for each observation
  for(i in 1:N){
    if(hurdle[i] == 1){
      log_lik[i] = gamma_lpdf(y[i] | invphi, betas[i]) - log_diff_exp(gamma_lcdf(1.0 | invphi, betas[i]), gamma_lcdf(0.0 | invphi, betas[i]));
    } else {
      log_lik[i] = bernoulli_logit_lpmf(0 | mu_ber[i]);
    }
  }

  // Calculate p_delta for each state
  for(i in 1:n_s){
    if(delta[i] >= 0.0){
      p_delta[i] = 1.0;
    } else {
      p_delta[i] = 0.0;
    }
  }

  // Predicted values for testing
  if(pred_ind == 1){
    real be;
    for(i in 1:N){
      if(t_scale[i] >= t_0){
        k_eff = k[state[i]] + delta[state[i]];
      } else {
        k_eff = k[state[i]];
      }
      mu[i] = L[state[i]] * inv_logit(k_eff * (t_scale[i] - x_0[state[i]]));
      be = invphi / mu[i];
      mu_ber[i] = alpha + t_scale[i] * beta[state[i]];
      pred[i] = inv_logit(mu_ber[i]) * mu[i];
    }
  }
}
