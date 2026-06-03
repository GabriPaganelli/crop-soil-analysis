// =============================================================================
// m_sp_rirs_B.stan  —  Modello B (response-specific): M-SP-RIRS con struttura
//                       adattata per risposta
//
// Confronto rispetto a M-SP-RIRS (script 13_confronto_modelli.R):
//   - SOC: struttura RI+RS bivariate (tau_alpha_SOC, tau_beta_SOC, rho_SOC)
//   - N:   solo random intercept (tau_beta_N = 0, prior molto stretto)
//          motivazione: rho_N ≈ 0 e tau_beta_N piccolo in M-SP-RIRS
//   - P:   struttura RI+RS bivariate (tau_alpha_P, tau_beta_P, rho_P)
//
// Domanda: la variabilità di slope per N ha rilevanza predittiva?
//
// X_W (K_W=5): logBottom[1], Texture1[2], Texture2[3], BulkDensity[4], PH[5]
// X_B (K_B=4): OnFarm, Irrigate, Fertilised, N_Natural
// =============================================================================

data {

  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K_W;
  int<lower=1> K_B;

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W] X_W;
  matrix[J, K_B] X_B;

}

parameters {

  // ── SOC: bivariate RI+RS ──────────────────────────────────────────────────────
  matrix[2, J]            z_SOC;
  cholesky_factor_corr[2] L_SOC;
  vector<lower=0>[2]      tau_SOC;

  // ── N: solo random intercept (tau_beta_N → 0) ────────────────────────────────
  vector[J]     z_N_int;   // solo intercetta
  real<lower=0> tau_N_int;

  // ── P: bivariate RI+RS ──────────────────────────────────────────────────────
  matrix[2, J]            z_P;
  cholesky_factor_corr[2] L_P;
  vector<lower=0>[2]      tau_P;

  // ── Effetti fissi ──────────────────────────────────────────────────────────────
  real          alpha_SOC;
  vector[K_W]   gamma_SOC;
  vector[K_B]   beta_SOC;
  real<lower=0> sigma_SOC;

  real          alpha_N;
  vector[K_W]   gamma_N;
  vector[K_B]   beta_N;
  real<lower=0> sigma_N;

  real          alpha_P;
  vector[K_W]   gamma_P;
  vector[K_B]   beta_P;
  real<lower=0> sigma_P;

}

transformed parameters {

  matrix[2, J] u_SOC = diag_pre_multiply(tau_SOC, L_SOC) * z_SOC;
  vector[J]    u_N   = tau_N_int * z_N_int;   // solo RI per N
  matrix[2, J] u_P   = diag_pre_multiply(tau_P,   L_P)   * z_P;

}

model {

  to_vector(z_SOC) ~ std_normal();
  L_SOC ~ lkj_corr_cholesky(2);
  tau_SOC ~ normal(0, 0.5);

  z_N_int ~ std_normal();
  tau_N_int ~ normal(0, 0.5);

  to_vector(z_P) ~ std_normal();
  L_P ~ lkj_corr_cholesky(2);
  tau_P ~ normal(0, 0.5);

  alpha_SOC ~ normal(0, 2);
  gamma_SOC ~ normal(0, 1);
  beta_SOC  ~ normal(0, 1);
  sigma_SOC ~ normal(0, 0.5);

  alpha_N ~ normal(0, 2);
  gamma_N ~ normal(0, 1);
  beta_N  ~ normal(0, 1);
  sigma_N ~ normal(0, 0.5);

  alpha_P ~ normal(0, 2);
  gamma_P ~ normal(0, 1);
  beta_P  ~ normal(0, 1);
  sigma_P ~ normal(0, 0.5);

  {
    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
    vector[N] gP   = X_W * gamma_P;
    vector[J] bSOC = X_B * beta_SOC;
    vector[J] bN   = X_B * beta_N;
    vector[J] bP   = X_B * beta_P;

    for (n in 1:N) {
      int  j      = field_id[n];
      real logB_n = X_W[n, 1];

      real mu_SOC = alpha_SOC + u_SOC[1,j] + u_SOC[2,j]*logB_n + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + u_N[j]                         + gN[n]   + bN[j];
      real mu_P   = alpha_P   + u_P[1,j]   + u_P[2,j]  *logB_n + gP[n]   + bP[j];

      logSOC[n] ~ normal(mu_SOC, sigma_SOC);
      logN[n]   ~ normal(mu_N,   sigma_N);
      logP[n]   ~ normal(mu_P,   sigma_P);
    }
  }

}

generated quantities {

  corr_matrix[2] Omega_SOC = multiply_lower_tri_self_transpose(L_SOC);
  corr_matrix[2] Omega_P   = multiply_lower_tri_self_transpose(L_P);

  real rho_SOC = Omega_SOC[1, 2];
  real rho_P   = Omega_P[1, 2];

  real tau_alpha_SOC = tau_SOC[1]; real tau_beta_SOC = tau_SOC[2];
  real tau_alpha_N   = tau_N_int;  // no tau_beta_N (=0 by construction)
  real tau_alpha_P   = tau_P[1];   real tau_beta_P   = tau_P[2];

  vector[N] log_lik;
  {
    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
    vector[N] gP   = X_W * gamma_P;
    vector[J] bSOC = X_B * beta_SOC;
    vector[J] bN   = X_B * beta_N;
    vector[J] bP   = X_B * beta_P;

    for (n in 1:N) {
      int  j      = field_id[n];
      real logB_n = X_W[n, 1];
      real mu_SOC = alpha_SOC + u_SOC[1,j] + u_SOC[2,j]*logB_n + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + u_N[j]                         + gN[n]   + bN[j];
      real mu_P   = alpha_P   + u_P[1,j]   + u_P[2,j]  *logB_n + gP[n]   + bP[j];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC, sigma_SOC)
                 + normal_lpdf(logN[n]   | mu_N,   sigma_N)
                 + normal_lpdf(logP[n]   | mu_P,   sigma_P);
    }
  }

}
