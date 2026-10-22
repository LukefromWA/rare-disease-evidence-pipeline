// M2a: global coefficients, independent evidence channels (guide A0).
functions {
  // Gaussian-copula log density, vectorized, per-row rho
  vector gauss_copula_logd(vector zu, vector zv, vector rho) {
    vector[rows(zu)] d = 1 - square(rho);
    return -0.5 * log(d) - (square(rho) .* (square(zu) + square(zv)) - 2 * rho .* zu .* zv) ./ (2 * d);
  }
}
data {
  // ---- training cases: candidates contiguous by case; y_pos = LOCAL index of the single causal candidate ----
  int<lower=1> N; int<lower=1> C; int<lower=1> K; int<lower=1> G;
  matrix[N, K] X;
  array[N] int<lower=1, upper=G> group_id;
  vector<lower=1e-6, upper=1 - 1e-6>[N] u_hpo;
  vector<lower=1e-6, upper=1 - 1e-6>[N] u_network;
  array[N] int<lower=0, upper=1> copula_available;
  array[C] int<lower=1, upper=N> case_start;
  array[C] int<lower=1, upper=N> case_end;
  array[C] int<lower=1> y_pos;
  // ---- priors (proper; all scales are data so prior sensitivity needs no recompile) ----
  vector<lower=0>[K] beta_prior_scale; vector<lower=0>[K] tau_beta_prior_scale;
  real<lower=0> mu_rho_prior_scale; real<lower=0> tau_rho_prior_scale;
  real<lower=0> beta_copula_prior_scale; real<lower=0> tau_copula_prior_scale;
  // ---- held-out cases scored in generated quantities (sizes 0 when fitting for strict BMA) ----
  int<lower=0> N_new; int<lower=0> C_new;
  matrix[N_new, K] X_new;
  array[N_new] int<lower=1, upper=G + 1> group_new;            // G + 1 = unseen group (effects drawn from the population)
  vector<lower=1e-6, upper=1 - 1e-6>[N_new] u_hpo_new;
  vector<lower=1e-6, upper=1 - 1e-6>[N_new] u_network_new;
  array[N_new] int<lower=0, upper=1> copula_available_new;
  array[C_new] int<lower=1, upper=N_new> case_start_new;
  array[C_new] int<lower=1, upper=N_new> case_end_new;
}
transformed data {
  vector[N] z_h = inv_Phi(u_hpo); vector[N] z_n = inv_Phi(u_network); vector[N] cop = to_vector(copula_available);
  vector[N_new] z_h_new = inv_Phi(u_hpo_new); vector[N_new] z_n_new = inv_Phi(u_network_new); vector[N_new] cop_new = to_vector(copula_available_new);
  for (c in 1:C) {
    if (case_end[c] < case_start[c]) reject("empty case ", c);
    if (y_pos[c] > case_end[c] - case_start[c] + 1) reject("y_pos out of range in case ", c);
  }
}
parameters { vector[K] beta; }
model {
  vector[N] eta = X * beta;
  target += normal_lpdf(beta | 0, beta_prior_scale);
  for (c in 1:C)
    target += categorical_logit_lpmf(y_pos[c] | segment(eta, case_start[c], case_end[c] - case_start[c] + 1));
}
generated quantities {
  vector[N_new] p_new;                                          // the ONLY saved quantity: held-out candidate probabilities
  if (N_new > 0) {
    vector[N_new] eta_new = X_new * beta;
    for (c in 1:C_new) {
      int s = case_start_new[c]; int n = case_end_new[c] - s + 1;
      p_new[s:(s + n - 1)] = softmax(segment(eta_new, s, n));
    }
  }
}
