// =============================================================================
// m_sp_rirs_mv.stan  —  M-SP-RIRS + Multivariate residuals
//
// Estensione di m_sp_rirs.stan:
//   I residui (eps_SOC, eps_N, eps_P) alla stessa osservazione (i,j) sono
//   modellati congiuntamente tramite MVNormal con matrice di scala Sigma.
//
//   Sigma = diag(sigma) * Omega_res * diag(sigma)
//   Omega_res ~ LKJ(2)  (prior debolmente informativo)
//
// Media (M-SP-RIRS):
//   mu_r[i,j] = alpha_r
//             + u_r[1, j]                      (random intercept per campo)
//             + u_r[2, j] * logBottom_i        (random slope per campo)
//             + gamma_r * X_W_i
//             + beta_r  * X_B_j
//
//   [u_r[1,j], u_r[2,j]] ~ MVN(0, Sigma_r)  con NCP Cholesky per risposta r.
//
// PARAMETRI AGGIUNTIVI rispetto a M-SP-RIRS:
//   sigma[3]           : vettore [sigma_SOC, sigma_N, sigma_P]
//   L_Omega_res[3,3]   : Cholesky di Omega_res (correlazioni residue)
//
// GENERATED QUANTITIES:
//   rho_SOC_N, rho_SOC_P, rho_N_P  : correlazioni residue cross-risposta
//   rho_SOC, rho_N, rho_P          : correlazioni RI–RS per risposta
//
// X_W (K_W=5): logBottom[1], Texture1[2], Texture2[3], BulkDensity[4], PH[5]
// TOTALE PARAMETRI: 6J + 45 scalari  (J=40: 285)
// =============================================================================

data {

  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K_W;   // 5: logBottom, Texture1, Texture2, BulkDensity, PH
  int<lower=1> K_B;   // 4: OnFarm, Irrigate, Fertilised, N_Natural

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W] X_W;   // colonna 1 = logBottom standardizzato
  matrix[J, K_B] X_B;

}

parameters {

  // ── Effetti casuali bivariati (NCP): [intercept, slope] per campo ─────────────
  matrix[2, J]            z_SOC;
  cholesky_factor_corr[2] L_SOC;
  vector<lower=0>[2]      tau_SOC;

  matrix[2, J]            z_N;
  cholesky_factor_corr[2] L_N;
  vector<lower=0>[2]      tau_N;

  matrix[2, J]            z_P;
  cholesky_factor_corr[2] L_P;
  vector<lower=0>[2]      tau_P;

  // ── Effetti fissi ──────────────────────────────────────────────────────────────
  real          alpha_SOC;
  vector[K_W]   gamma_SOC;
  vector[K_B]   beta_SOC;

  real          alpha_N;
  vector[K_W]   gamma_N;
  vector[K_B]   beta_N;

  real          alpha_P;
  vector[K_W]   gamma_P;
  vector[K_B]   beta_P;

  // ── Struttura multivariata dei residui ────────────────────────────────────────
  vector<lower=0>[3]      sigma;       // [sigma_SOC, sigma_N, sigma_P]
  cholesky_factor_corr[3] L_Omega_res; // Cholesky di Omega residua 3×3

}

transformed parameters {

  matrix[2, J] u_SOC = diag_pre_multiply(tau_SOC, L_SOC) * z_SOC;
  matrix[2, J] u_N   = diag_pre_multiply(tau_N,   L_N)   * z_N;
  matrix[2, J] u_P   = diag_pre_multiply(tau_P,   L_P)   * z_P;

}

model {

  // ── Prior: struttura bivariata ─────────────────────────────────────────────────
  to_vector(z_SOC) ~ std_normal();
  L_SOC ~ lkj_corr_cholesky(2);
  tau_SOC ~ normal(0, 0.5);

  to_vector(z_N) ~ std_normal();
  L_N ~ lkj_corr_cholesky(2);
  tau_N ~ normal(0, 0.5);

  to_vector(z_P) ~ std_normal();
  L_P ~ lkj_corr_cholesky(2);
  tau_P ~ normal(0, 0.5);

  // ── Prior: effetti fissi ──────────────────────────────────────────────────────
  alpha_SOC ~ normal(0, 2);
  gamma_SOC ~ normal(0, 1);
  beta_SOC  ~ normal(0, 1);

  alpha_N ~ normal(0, 2);
  gamma_N ~ normal(0, 1);
  beta_N  ~ normal(0, 1);

  alpha_P ~ normal(0, 2);
  gamma_P ~ normal(0, 1);
  beta_P  ~ normal(0, 1);

  // ── Prior: struttura residua multivariata ─────────────────────────────────────
  sigma       ~ normal(0, 0.5);
  L_Omega_res ~ lkj_corr_cholesky(2);

  // ── Likelihood ───────────────────────────────────────────────────────────────
  {
    matrix[3, 3] L_Sigma = diag_pre_multiply(sigma, L_Omega_res);

    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
    vector[N] gP   = X_W * gamma_P;

    vector[J] bSOC = X_B * beta_SOC;
    vector[J] bN   = X_B * beta_N;
    vector[J] bP   = X_B * beta_P;

    for (n in 1:N) {
      int  j      = field_id[n];
      real logB_n = X_W[n, 1];   // logBottom standardizzato

      vector[3] mu_n;
      mu_n[1] = alpha_SOC + u_SOC[1,j] + u_SOC[2,j]*logB_n + gSOC[n] + bSOC[j];
      mu_n[2] = alpha_N   + u_N[1,j]   + u_N[2,j]  *logB_n + gN[n]   + bN[j];
      mu_n[3] = alpha_P   + u_P[1,j]   + u_P[2,j]  *logB_n + gP[n]   + bP[j];

      vector[3] y_n;
      y_n[1] = logSOC[n];
      y_n[2] = logN[n];
      y_n[3] = logP[n];

      y_n ~ multi_normal_cholesky(mu_n, L_Sigma);
    }
  }

}

generated quantities {

  // ── Correlazioni RI–RS per risposta ───────────────────────────────────────────
  corr_matrix[2] Omega_SOC = multiply_lower_tri_self_transpose(L_SOC);
  corr_matrix[2] Omega_N   = multiply_lower_tri_self_transpose(L_N);
  corr_matrix[2] Omega_P   = multiply_lower_tri_self_transpose(L_P);

  real rho_SOC = Omega_SOC[1, 2];
  real rho_N   = Omega_N[1, 2];
  real rho_P   = Omega_P[1, 2];

  real tau_alpha_SOC = tau_SOC[1]; real tau_beta_SOC = tau_SOC[2];
  real tau_alpha_N   = tau_N[1];   real tau_beta_N   = tau_N[2];
  real tau_alpha_P   = tau_P[1];   real tau_beta_P   = tau_P[2];

  // ── Correlazioni residue cross-risposta ───────────────────────────────────────
  corr_matrix[3] Omega_res = multiply_lower_tri_self_transpose(L_Omega_res);
  real rho_SOC_N = Omega_res[1, 2];
  real rho_SOC_P = Omega_res[1, 3];
  real rho_N_P   = Omega_res[2, 3];

  real sigma_SOC = sigma[1];
  real sigma_N   = sigma[2];
  real sigma_P   = sigma[3];

  // ── Log-likelihood per LOO-CV ─────────────────────────────────────────────────
  vector[N] log_lik;
  {
    matrix[3, 3] L_Sigma = diag_pre_multiply(sigma, L_Omega_res);

    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
    vector[N] gP   = X_W * gamma_P;
    vector[J] bSOC = X_B * beta_SOC;
    vector[J] bN   = X_B * beta_N;
    vector[J] bP   = X_B * beta_P;

    for (n in 1:N) {
      int  j      = field_id[n];
      real logB_n = X_W[n, 1];

      vector[3] mu_n;
      mu_n[1] = alpha_SOC + u_SOC[1,j] + u_SOC[2,j]*logB_n + gSOC[n] + bSOC[j];
      mu_n[2] = alpha_N   + u_N[1,j]   + u_N[2,j]  *logB_n + gN[n]   + bN[j];
      mu_n[3] = alpha_P   + u_P[1,j]   + u_P[2,j]  *logB_n + gP[n]   + bP[j];

      vector[3] y_n;
      y_n[1] = logSOC[n];
      y_n[2] = logN[n];
      y_n[3] = logP[n];

      log_lik[n] = multi_normal_cholesky_lpdf(y_n | mu_n, L_Sigma);
    }
  }

}
