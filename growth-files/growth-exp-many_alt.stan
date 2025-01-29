//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
data {
  int<lower=0> N;
  vector[N] y;
  vector[N] t;
  real t_cutoff;
  array[N] int state;
  int n_s;
  vector[n_s] y_0;
  int pred_ind;
}
transformed data {
  vector[N] t_scale;
  real t_0;
  t_scale = (t - min(t)) / (max(t) - min(t));
  t_0 = (t_cutoff - min(t)) / (max(t) - min(t));
}

// Parameters
parameters {
  vector[n_s] k_0;
  vector[n_s] delta_0;

  real<lower=0> tau_k;
  real<lower=0> tau_delta;

  real<lower=0> alpha_k;
  real<lower=0> alpha_delta;

  real<lower=0> invphi;
}
transformed parameters {
  vector[n_s] k;
  vector[n_s] delta;

  for (i in 1:n_s) {
    k[i] = alpha_k + k_0[i] * tau_k;
    delta[i] = alpha_delta + delta_0[i] * tau_delta;
  }
}

// Model
model {
  vector[N] mu;
  vector[N] beta;
  real k_eff;

  for (i in 1:N) {
    if (t_scale[i] >= t_0) {
      k_eff = k[state[i]] + delta[state[i]];
    } else {
      k_eff = k[state[i]];
    }
    mu[i] = y_0[state[i]] * exp(k_eff * t_scale[i]);
    beta[i] = invphi / mu[i];
  }

  y ~ gamma(invphi, beta);

  k_0 ~ normal(0.0, 1.0);
  delta_0 ~ normal(0.0, 1.0);

  tau_k ~ normal(0.0, 1.0);
  tau_delta ~ normal(0.0, 0.5);

  alpha_k ~ normal(7.0, 2.0);
  alpha_delta ~ normal(0.0, 1.0);

  invphi ~ normal(0.5, 1.0);
}

// Generated Quantities
generated quantities {
  vector[n_s] p_delta;
  vector[N] log_lik;  // Log-likelihood for each observation
  vector[N] pred;     // Predictions (if needed)
  real k_eff;
  real be;
  real mu;

  for (i in 1:n_s) {
    if (delta[i] >= 0.0) {
      p_delta[i] = 1.0;
    } else {
      p_delta[i] = 0.0;
    }
  }

  // Log-likelihood computation
  for (i in 1:N) {
    if (t_scale[i] >= t_0) {
      k_eff = k[state[i]] + delta[state[i]];
    } else {
      k_eff = k[state[i]];
    }
    mu = y_0[state[i]] * exp(k_eff * t_scale[i]);
    be = invphi / mu;
    log_lik[i] = gamma_lpdf(y[i] | invphi, be);
  }

  // Predictions (if pred_ind == 1)
  if (pred_ind == 1) {
    for (i in 1:N) {
      if (t_scale[i] >= t_0) {
        k_eff = k[state[i]] + delta[state[i]];
      } else {
        k_eff = k[state[i]];
      }
      mu = y_0[state[i]] * exp(k_eff * t_scale[i]);
      pred[i] = mu;
    }
  }
}
