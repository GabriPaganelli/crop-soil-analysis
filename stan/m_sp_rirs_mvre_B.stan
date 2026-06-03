// =============================================================================
// m_sp_rirs_mvre_B.stan  —  M-SP-RIRS-MVRE Modello B (risposta-specifico)
//
// Semplificazione di m_sp_rirs_mvre.stan:
//   - N ha solo random intercept (no random slope)
//   - V[5, J]: [SOC_int, SOC_slope, N_int, P_int, P_slope]
//   - mu_N = alpha_N + V[3,j] + gamma_N*X_W + beta_N*X_B  (no termine slope)
//   - SOC e P hanno la struttura piena RI+RS
//
// Domanda: la variabilità di slope per N è predittivamente rilevante?
// Se LOO(B) ≈ LOO(MVRE) → la slope di N è superflua.
//
// Struttura Omega_5 (5×5):
//   righe/colonne: [SOC_int, SOC_slope, N_int, P_int, P_slope]
//   rho_SOC = Omega_5[1,2], rho_P = Omega_5[4,5]
//   rho_int_SOC_N = Omega_5[1,3], rho_int_SOC_P = Omega_5[1,4], rho_int_N_P = Omega_5[3,4]
//
// TOTALE PARAMETRI: 5J + 48 scalari  (J=40: 248)
//   tau_5[5]: 5; L_5: 10; alpha_r: 3; gamma_r[5]: 15; beta_r[4]: 12; sigma_r: 3
//   + z_v[5,J]: 5J
// =============================================================================

data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K_W;   // 5
  int<lower=1> K_B;   // 4

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W] X_W;
  matrix[J, K_B] X_B;
}

parameters {
  // ── Effetti casuali 5D: [SOC_int, SOC_slope, N_int, P_int, P_slope] ────────
  matrix[5, J]            z_v;
  cholesky_factor_corr[5] L_5;
  vector<lower=0>[5]      tau_5;

  // ── Effetti fissi ─────────────────────────────────────────────────────────
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
  // V[5, J]: righe = [SOC_int, SOC_slope, N_int, P_int, P_slope]
  matrix[5, J] V = diag_pre_multiply(tau_5, L_5) * z_v;
}

model {
  to_vector(z_v) ~ std_normal();
  L_5   ~ lkj_corr_cholesky(2);
  tau_5 ~ normal(0, 0.5);

  alpha_SOC ~ normal(0, 2);  gamma_SOC ~ normal(0, 1);  beta_SOC ~ normal(0, 1);  sigma_SOC ~ normal(0, 0.5);
  alpha_N   ~ normal(0, 2);  gamma_N   ~ normal(0, 1);  beta_N   ~ normal(0, 1);  sigma_N   ~ normal(0, 0.5);
  alpha_P   ~ normal(0, 2);  gamma_P   ~ normal(0, 1);  beta_P   ~ normal(0, 1);  sigma_P   ~ normal(0, 0.5);

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

      // SOC: V[1]=int, V[2]=slope
      logSOC[n] ~ normal(alpha_SOC + V[1,j] + V[2,j]*logB_n + gSOC[n] + bSOC[j], sigma_SOC);
      // N: V[3]=int only (no slope)
      logN[n]   ~ normal(alpha_N   + V[3,j]             + gN[n]   + bN[j],   sigma_N);
      // P: V[4]=int, V[5]=slope
      logP[n]   ~ normal(alpha_P   + V[4,j] + V[5,j]*logB_n + gP[n]   + bP[j],   sigma_P);
    }
  }
}

generated quantities {
  corr_matrix[5] Omega_5 = multiply_lower_tri_self_transpose(L_5);

  // Correlazioni RI-RS (per risposta; N non ha slope)
  real rho_SOC = Omega_5[1, 2];   // SOC int-slope
  real rho_P   = Omega_5[4, 5];   // P int-slope

  // Scales
  real tau_alpha_SOC = tau_5[1]; real tau_beta_SOC = tau_5[2];
  real tau_alpha_N   = tau_5[3];
  real tau_alpha_P   = tau_5[4]; real tau_beta_P   = tau_5[5];

  // Cross-response intercept correlations
  real rho_int_SOC_N = Omega_5[1, 3];
  real rho_int_SOC_P = Omega_5[1, 4];
  real rho_int_N_P   = Omega_5[3, 4];

  // Cross-response slope correlations (only SOC and P have slopes)
  real rho_slope_SOC_P = Omega_5[2, 5];

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
      log_lik[n] = normal_lpdf(logSOC[n] | alpha_SOC + V[1,j] + V[2,j]*logB_n + gSOC[n] + bSOC[j], sigma_SOC)
                 + normal_lpdf(logN[n]   | alpha_N   + V[3,j]             + gN[n]   + bN[j],   sigma_N)
                 + normal_lpdf(logP[n]   | alpha_P   + V[4,j] + V[5,j]*logB_n + gP[n]   + bP[j],   sigma_P);
    }
  }
}
