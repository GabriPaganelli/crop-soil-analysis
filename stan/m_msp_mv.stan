// =============================================================================
// m_msp_mv.stan  —  M-SP + Multivariate residuals (M-SP-MV)
//
// Estensione di m4rr_v2_ri_slope_mu.stan:
//   I residui (eps_SOC, eps_N, eps_P) alla stessa osservazione (i,j) sono
//   modellati congiuntamente tramite MVNormal con matrice di scala Sigma.
//
//   Sigma = diag(sigma) * Omega * diag(sigma)
//   Omega ~ LKJ(2)  (prior debolmente informativo: favorisce indipendenza)
//
// PARAMETRI AGGIUNTIVI rispetto a M-SP:
//   sigma[3]      : vettore [sigma_SOC, sigma_N, sigma_P]  (sostituisce 3 scalari)
//   L_Omega[3,3]  : fattore di Cholesky di Omega
//
// GENERATED QUANTITIES AGGIUNTIVI:
//   rho_SOC_N, rho_SOC_P, rho_N_P  : correlazioni residue cross-risposta
//
// INTERPRETAZIONE:
//   rho_SOC_N > 0 → a stessa profondità e stesso campo, se SOC > atteso
//                   allora N tende a essere > atteso (co-variazione residua)
//
// TOTALE PARAMETRI: 3J + 36 + 3 correlazioni = 159 (con J=40)
// =============================================================================

data {

  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K_W;   // 4: logBottom, Texture1, Texture2, BulkDensity
  int<lower=1> K_B;   // 4: OnFarm, Irrigate, Fertilised, N_Natural

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W] X_W;   // colonna 1 = logBottom standardizzato
  matrix[J, K_B] X_B;

}

parameters {

  // ── SOC ──────────────────────────────────────────────────────────────────────
  real          alpha_SOC;
  vector[J]     z_nu_SOC;
  real<lower=0> psi_SOC;
  real          eta_SOC;
  vector[K_W]   gamma_SOC;
  vector[K_B]   beta_SOC;

  // ── N ────────────────────────────────────────────────────────────────────────
  real          alpha_N;
  vector[J]     z_nu_N;
  real<lower=0> psi_N;
  real          eta_N;
  vector[K_W]   gamma_N;
  vector[K_B]   beta_N;

  // ── P ────────────────────────────────────────────────────────────────────────
  real          alpha_P;
  vector[J]     z_nu_P;
  real<lower=0> psi_P;
  real          eta_P;
  vector[K_W]   gamma_P;
  vector[K_B]   beta_P;

  // ── Struttura multivariata dei residui ────────────────────────────────────────
  vector<lower=0>[3]      sigma;    // [sigma_SOC, sigma_N, sigma_P]
  cholesky_factor_corr[3] L_Omega;  // Cholesky di Omega (matrice di correlazione 3x3)

}

model {

  // ── Prior: SOC ───────────────────────────────────────────────────────────────
  alpha_SOC ~ normal(0, 2);
  z_nu_SOC  ~ std_normal();
  psi_SOC   ~ normal(0, 0.5);
  eta_SOC   ~ normal(0, 0.3);
  gamma_SOC ~ normal(0, 1);
  beta_SOC  ~ normal(0, 1);

  // ── Prior: N ─────────────────────────────────────────────────────────────────
  alpha_N ~ normal(0, 2);
  z_nu_N  ~ std_normal();
  psi_N   ~ normal(0, 0.5);
  eta_N   ~ normal(0, 0.3);
  gamma_N ~ normal(0, 1);
  beta_N  ~ normal(0, 1);

  // ── Prior: P ─────────────────────────────────────────────────────────────────
  alpha_P ~ normal(0, 2);
  z_nu_P  ~ std_normal();
  psi_P   ~ normal(0, 0.5);
  eta_P   ~ normal(0, 0.3);
  gamma_P ~ normal(0, 1);
  beta_P  ~ normal(0, 1);

  // ── Prior: struttura multivariata ─────────────────────────────────────────────
  sigma   ~ normal(0, 0.5);        // half-normal: sigma > 0 per constraint
  L_Omega ~ lkj_corr_cholesky(2);  // LKJ(2): debolmente informativo, favorisce Omega ≈ I

  // ── Likelihood ───────────────────────────────────────────────────────────────
  {
    matrix[3, 3] L_Sigma = diag_pre_multiply(sigma, L_Omega);

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
      mu_n[1] = alpha_SOC + z_nu_SOC[j] * (psi_SOC + eta_SOC * logB_n) + gSOC[n] + bSOC[j];
      mu_n[2] = alpha_N   + z_nu_N[j]   * (psi_N   + eta_N   * logB_n) + gN[n]   + bN[j];
      mu_n[3] = alpha_P   + z_nu_P[j]   * (psi_P   + eta_P   * logB_n) + gP[n]   + bP[j];

      vector[3] y_n;
      y_n[1] = logSOC[n];
      y_n[2] = logN[n];
      y_n[3] = logP[n];

      y_n ~ multi_normal_cholesky(mu_n, L_Sigma);
    }
  }

}

generated quantities {

  // ── Matrice di correlazione dei residui ───────────────────────────────────────
  corr_matrix[3] Omega = multiply_lower_tri_self_transpose(L_Omega);
  real rho_SOC_N = Omega[1, 2];   // correlazione residua SOC–N
  real rho_SOC_P = Omega[1, 3];   // correlazione residua SOC–P
  real rho_N_P   = Omega[2, 3];   // correlazione residua N–P

  // ── Deviazioni standard (per compatibilità con M-SP) ──────────────────────────
  real sigma_SOC = sigma[1];
  real sigma_N   = sigma[2];
  real sigma_P   = sigma[3];

  // ── b_r = eta_r / psi_r ───────────────────────────────────────────────────────
  real b_SOC = (psi_SOC > 0.02) ? eta_SOC / psi_SOC : 0.0;
  real b_N   = (psi_N   > 0.02) ? eta_N   / psi_N   : 0.0;
  real b_P   = (psi_P   > 0.02) ? eta_P   / psi_P   : 0.0;

  // ── Log-likelihood per LOO-CV (likelihood congiunta MVN per osservazione) ─────
  vector[N] log_lik;
  {
    matrix[3, 3] L_Sigma = diag_pre_multiply(sigma, L_Omega);

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
      mu_n[1] = alpha_SOC + z_nu_SOC[j] * (psi_SOC + eta_SOC * logB_n) + gSOC[n] + bSOC[j];
      mu_n[2] = alpha_N   + z_nu_N[j]   * (psi_N   + eta_N   * logB_n) + gN[n]   + bN[j];
      mu_n[3] = alpha_P   + z_nu_P[j]   * (psi_P   + eta_P   * logB_n) + gP[n]   + bP[j];

      vector[3] y_n;
      y_n[1] = logSOC[n];
      y_n[2] = logN[n];
      y_n[3] = logP[n];

      log_lik[n] = multi_normal_cholesky_lpdf(y_n | mu_n, L_Sigma);
    }
  }

}
