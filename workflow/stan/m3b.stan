// M3b: hierarchical coefficients + group-level copula dependence and utility (guide A3).
functions {
  // Gaussian-copula log density, vectorized, per-row rho
  vector gauss_copula_logd(vector zu, vector zv, vector rho) {
    vector[rows(zu)] d = 1 - square(rho);
    return -0.5 * log(d) - (square(rho) .* (square(zu) + square(zv)) - 2 * rho .* zu .* zv) ./ (2 * d);
  }
  matrix make_beta_g(vector beta, vector tau, matrix z) {      // G x K partially pooled coefficients
    return rep_matrix(beta', rows(z)) + rep_matrix(tau', rows(z)) .* z;
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
parameters {
  vector[K] beta; matrix[G, K] z_beta; vector<lower=0>[K] tau_beta;
  real mu_rho_raw; real<lower=0> tau_rho; vector[G] z_rho;
  real beta_cop; real<lower=0> tau_cop; vector[G] z_cop;
}
model {
  vector[G] rho_g = tanh(mu_rho_raw + tau_rho * z_rho);
  vector[G] bcop_g = beta_cop + tau_cop * z_cop;
  vector[N] eta = rows_dot_product(X, make_beta_g(beta, tau_beta, z_beta)[group_id])
                  + bcop_g[group_id] .* (cop .* gauss_copula_logd(z_h, z_n, rho_g[group_id]));
  target += normal_lpdf(beta | 0, beta_prior_scale);
  target += std_normal_lpdf(to_vector(z_beta));
  target += normal_lpdf(tau_beta | 0, tau_beta_prior_scale) + K * log(2);
  target += normal_lpdf(mu_rho_raw | 0, mu_rho_prior_scale);
  target += normal_lpdf(tau_rho | 0, tau_rho_prior_scale) + log(2);
  target += std_normal_lpdf(z_rho);
  target += normal_lpdf(beta_cop | 0, beta_copula_prior_scale);
  target += normal_lpdf(tau_cop | 0, tau_copula_prior_scale) + log(2);
  target += std_normal_lpdf(z_cop);
  for (c in 1:C)
    target += categorical_logit_lpmf(y_pos[c] | segment(eta, case_start[c], case_end[c] - case_start[c] + 1));
}
generated quantities {
  vector[N_new] p_new;                                          // the ONLY saved quantity: held-out candidate probabilities
  if (N_new > 0) {
    matrix[G + 1, K] bg = append_row(make_beta_g(beta, tau_beta, z_beta), (beta + tau_beta .* to_vector(normal_rng(rep_vector(0, K), 1)))');
    vector[G + 1] rho_all = append_row(tanh(mu_rho_raw + tau_rho * z_rho), [tanh(mu_rho_raw + tau_rho * normal_rng(0, 1))]');
    vector[G + 1] bcop = append_row(beta_cop + tau_cop * z_cop, [beta_cop + tau_cop * normal_rng(0, 1)]');
    vector[N_new] eta_new = rows_dot_product(X_new, bg[group_new]) + bcop[group_new] .* (cop_new .* gauss_copula_logd(z_h_new, z_n_new, rho_all[group_new]));
    for (c in 1:C_new) {
      int s = case_start_new[c]; int n = case_end_new[c] - s + 1;
      p_new[s:(s + n - 1)] = softmax(segment(eta_new, s, n));
    }
  }
}
