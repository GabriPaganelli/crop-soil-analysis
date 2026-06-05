// =============================================================================
// m_sp_rirs_mvre_full.stan  —  MVRE + Multivariate residuals (MVRE-FULL)
//
// Estensione di m_sp_rirs_mvre.stan:
//   Aggiunge correlazioni tra residui alla stessa osservazione (i,j)
//   tramite MVNormal. Cattura due livelli di correlazione cross-risposta:
//
//   1. Livello di campo (da MVRE):
//      V[6,J]: correlazioni RI-RS e cross-risposta tra intercette/slope
//
//   2. Livello di osservazione (NUOVO):
//      Omega_res [3x3]: correlazioni residue SOC-N, SOC-P, N-P
//      Interpretazione: dato che il campo j ha un certo profilo V[.,j],
//      i tre residui alla stessa profondità i sono ancora correlati?
//
// Media (identica a MVRE):
//   mu_r[i,j] = alpha_r + V[2r-1,j] + V[2r,j]*logBottom_i
//             + gamma_r * X_W_i + beta_r * X_B_j
//
// Likelihood (MVN, vs 3 normali indipendenti in MVRE):
//   [logSOC, logN, logP]_i ~ MVN(mu_i, Sigma_res)
//   Sigma_res = diag(sigma_3) * Omega_res * diag(sigma_3)
//   Omega_res ~ LKJ(2)
//
// PARAMETRI AGGIUNTIVI rispetto a MVRE:
//   sigma[3]         : [sigma_SOC, sigma_N, sigma_P]  (come in MVRE)
//   L_Omega_res[3,3] : Cholesky di Omega_res  (+3 correlazioni libere)
//
// GENERATED QUANTITIES nuovi:
//   rho_res_SOC_N, rho_res_SOC_P, rho_res_N_P : correlazioni residue
//
// TOTALE PARAMETRI: 6J + 60 scalari  (J=40: 300)
//   Base MVRE: 6J + 57; +3 da L_Omega_res
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

  // ── Effetti casuali 6D (identico a MVRE) ──────────────────────────────────
  matrix[6, J]            z_v;
  cholesky_factor_corr[6] L_6;
  vector<lower=0>[6]      tau_6;

  // ── Effetti fissi ──────────────────────────────────────────────────────────
  real          alpha_SOC;
  vector[K_W]   gamma_SOC;
  vector[K_B]   beta_SOC;

  real          alpha_N;
  vector[K_W]   gamma_N;
  vector[K_B]   beta_N;

  real          alpha_P;
  vector[K_W]   gamma_P;
  vector[K_B]   beta_P;

  // ── Struttura MVN residuale (NUOVO rispetto a MVRE) ────────────────────────
  vector<lower=0>[3]      sigma;        // [sigma_SOC, sigma_N, sigma_P]
  cholesky_factor_corr[3] L_Omega_res;  // Cholesky di Omega_res 3×3

}

transformed parameters {

  // V[row, j]: riga 1=SOC_int, 2=SOC_slope, 3=N_int, 4=N_slope, 5=P_int, 6=P_slope
  matrix[6, J] V = diag_pre_multiply(tau_6, L_6) * z_v;

}

model {

  // ── Prior: effetti casuali 6D ──────────────────────────────────────────────
  to_vector(z_v) ~ std_normal();
  L_6   ~ lkj_corr_cholesky(2);
  tau_6 ~ normal(0, 0.5);

  // ── Prior: effetti fissi ──────────────────────────────────────────────────
  alpha_SOC ~ normal(0, 2);
  gamma_SOC ~ normal(0, 1);
  beta_SOC  ~ normal(0, 1);

  alpha_N ~ normal(0, 2);
  gamma_N ~ normal(0, 1);
  beta_N  ~ normal(0, 1);

  alpha_P ~ normal(0, 2);
  gamma_P ~ normal(0, 1);
  beta_P  ~ normal(0, 1);

  // ── Prior: struttura residua MVN ──────────────────────────────────────────
  sigma       ~ normal(0, 0.5);
  L_Omega_res ~ lkj_corr_cholesky(2);

  // ── Likelihood: MVN per osservazione ─────────────────────────────────────
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
      mu_n[1] = alpha_SOC + V[1,j] + V[2,j]*logB_n + gSOC[n] + bSOC[j];
      mu_n[2] = alpha_N   + V[3,j] + V[4,j]*logB_n + gN[n]   + bN[j];
      mu_n[3] = alpha_P   + V[5,j] + V[6,j]*logB_n + gP[n]   + bP[j];

      vector[3] y_n = [logSOC[n], logN[n], logP[n]]';

      y_n ~ multi_normal_cholesky(mu_n, L_Sigma);
    }
  }

}

generated quantities {

  // ── Matrice di correlazione 6×6 degli RE (come in MVRE) ──────────────────
  corr_matrix[6] Omega_6 = multiply_lower_tri_self_transpose(L_6);

  // ── Correlazioni RI-RS per risposta ──────────────────────────────────────
  real rho_SOC = Omega_6[1, 2];
  real rho_N   = Omega_6[3, 4];
  real rho_P   = Omega_6[5, 6];

  real tau_alpha_SOC = tau_6[1]; real tau_beta_SOC = tau_6[2];
  real tau_alpha_N   = tau_6[3]; real tau_beta_N   = tau_6[4];
  real tau_alpha_P   = tau_6[5]; real tau_beta_P   = tau_6[6];

  // ── Correlazioni cross-risposta degli RE ──────────────────────────────────
  real rho_int_SOC_N   = Omega_6[1, 3];
  real rho_int_SOC_P   = Omega_6[1, 5];
  real rho_int_N_P     = Omega_6[3, 5];
  real rho_slope_SOC_N = Omega_6[2, 4];
  real rho_slope_SOC_P = Omega_6[2, 6];
  real rho_slope_N_P   = Omega_6[4, 6];

  // ── Correlazioni residue cross-risposta (NUOVO) ───────────────────────────
  corr_matrix[3] Omega_res = multiply_lower_tri_self_transpose(L_Omega_res);
  real rho_res_SOC_N = Omega_res[1, 2];
  real rho_res_SOC_P = Omega_res[1, 3];
  real rho_res_N_P   = Omega_res[2, 3];

  real sigma_SOC = sigma[1];
  real sigma_N   = sigma[2];
  real sigma_P   = sigma[3];

  // ── Log-likelihood per LOO-CV ─────────────────────────────────────────────
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
      mu_n[1] = alpha_SOC + V[1,j] + V[2,j]*logB_n + gSOC[n] + bSOC[j];
      mu_n[2] = alpha_N   + V[3,j] + V[4,j]*logB_n + gN[n]   + bN[j];
      mu_n[3] = alpha_P   + V[5,j] + V[6,j]*logB_n + gP[n]   + bP[j];

      vector[3] y_n = [logSOC[n], logN[n], logP[n]]';

      log_lik[n] = multi_normal_cholesky_lpdf(y_n | mu_n, L_Sigma);
    }
  }

}
