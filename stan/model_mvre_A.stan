// =============================================================================
// m_sp_rirs_mvre_A.stan  —  M-SP-RIRS-MVRE Modello A (ridotto)
//
// Semplificazione di m_sp_rirs_mvre.stan:
//   - Senza effetti between-field (no beta_r, no X_B)
//   - K_W = 4: logBottom[1], Texture2[2], BulkDensity[3], PH[4]
//             (drop Texture1/ILR1)
//   - Struttura 6D random effects invariata
//
// Domanda: management (beta) e ILR1 (Texture1) contribuiscono predittivamente?
// Se LOO(A) ≈ LOO(MVRE) → sono superflui.
//
// TOTALE PARAMETRI: 6J + 45 scalari  (J=40: 285)
//   tau_6[6]: 6; L_6: 15; alpha_r: 3; gamma_r[4]: 12; sigma_r: 3
//   + z_v[6,J]: 6J
// =============================================================================

data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K_W;   // 4: logBottom, Texture2, BulkDensity, PH

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W] X_W;   // colonna 1 = logBottom standardizzato
}

parameters {
  // ── Effetti casuali 6D ─────────────────────────────────────────────────────
  matrix[6, J]            z_v;
  cholesky_factor_corr[6] L_6;
  vector<lower=0>[6]      tau_6;

  // ── Effetti fissi (no beta_r) ─────────────────────────────────────────────
  real          alpha_SOC;
  vector[K_W]   gamma_SOC;
  real<lower=0> sigma_SOC;

  real          alpha_N;
  vector[K_W]   gamma_N;
  real<lower=0> sigma_N;

  real          alpha_P;
  vector[K_W]   gamma_P;
  real<lower=0> sigma_P;
}

transformed parameters {
  matrix[6, J] V = diag_pre_multiply(tau_6, L_6) * z_v;
}

model {
  to_vector(z_v) ~ std_normal();
  L_6   ~ lkj_corr_cholesky(2);
  tau_6 ~ normal(0, 0.5);

  alpha_SOC ~ normal(0, 2);  gamma_SOC ~ normal(0, 1);  sigma_SOC ~ normal(0, 0.5);
  alpha_N   ~ normal(0, 2);  gamma_N   ~ normal(0, 1);  sigma_N   ~ normal(0, 0.5);
  alpha_P   ~ normal(0, 2);  gamma_P   ~ normal(0, 1);  sigma_P   ~ normal(0, 0.5);

  {
    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
    vector[N] gP   = X_W * gamma_P;

    for (n in 1:N) {
      int  j      = field_id[n];
      real logB_n = X_W[n, 1];

      logSOC[n] ~ normal(alpha_SOC + V[1,j] + V[2,j]*logB_n + gSOC[n], sigma_SOC);
      logN[n]   ~ normal(alpha_N   + V[3,j] + V[4,j]*logB_n + gN[n],   sigma_N);
      logP[n]   ~ normal(alpha_P   + V[5,j] + V[6,j]*logB_n + gP[n],   sigma_P);
    }
  }
}

generated quantities {
  corr_matrix[6] Omega_6 = multiply_lower_tri_self_transpose(L_6);

  real rho_SOC = Omega_6[1, 2];
  real rho_N   = Omega_6[3, 4];
  real rho_P   = Omega_6[5, 6];

  real tau_alpha_SOC = tau_6[1]; real tau_beta_SOC = tau_6[2];
  real tau_alpha_N   = tau_6[3]; real tau_beta_N   = tau_6[4];
  real tau_alpha_P   = tau_6[5]; real tau_beta_P   = tau_6[6];

  real rho_int_SOC_N   = Omega_6[1, 3];
  real rho_int_SOC_P   = Omega_6[1, 5];
  real rho_int_N_P     = Omega_6[3, 5];
  real rho_slope_SOC_N = Omega_6[2, 4];
  real rho_slope_SOC_P = Omega_6[2, 6];
  real rho_slope_N_P   = Omega_6[4, 6];

  vector[N] log_lik;
  {
    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
    vector[N] gP   = X_W * gamma_P;
    for (n in 1:N) {
      int  j      = field_id[n];
      real logB_n = X_W[n, 1];
      log_lik[n] = normal_lpdf(logSOC[n] | alpha_SOC + V[1,j] + V[2,j]*logB_n + gSOC[n], sigma_SOC)
                 + normal_lpdf(logN[n]   | alpha_N   + V[3,j] + V[4,j]*logB_n + gN[n],   sigma_N)
                 + normal_lpdf(logP[n]   | alpha_P   + V[5,j] + V[6,j]*logB_n + gP[n],   sigma_P);
    }
  }
}
