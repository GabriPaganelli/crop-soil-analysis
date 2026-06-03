// =============================================================================
// m_sp_rirs_mvre.stan  —  M-SP-RIRS + Multivariate Random Effects (MVRE)
//
// Estensione di m_sp_rirs.stan: modella le correlazioni cross-risposta
// tra gli effetti di campo, non solo la correlazione RI–RS per risposta.
//
// STRUTTURA RANDOM EFFECTS:
//   Per ogni campo j, il vettore 6D degli effetti casuali è:
//   v_j = [u_int_SOC, u_slope_SOC, u_int_N, u_slope_N, u_int_P, u_slope_P]
//   v_j ~ MVN(0, Sigma_6)
//   Sigma_6 = diag(tau_6) * Omega_6 * diag(tau_6)
//   Omega_6 ~ LKJ(2)  →  15 correlazioni libere
//
//   Parametrizzazione NCP:
//   z_v[6, J] ~ N(0, 1)
//   V = diag_pre_multiply(tau_6, L_6) * z_v    (6×J)
//
// MEDIA (identica a M-SP-RIRS):
//   mu_r[i,j] = alpha_r + V[2r-1, j] + V[2r, j]*logBottom_i
//             + gamma_r * X_W_i + beta_r * X_B_j
//
// GENERATED QUANTITIES:
//   rho_r           = Omega_6[2r-1, 2r]   (correlazione RI–RS per risposta r)
//   rho_int_SOC_N   = Omega_6[1, 3]        (correlazione intercette SOC–N)
//   rho_int_SOC_P   = Omega_6[1, 5]
//   rho_int_N_P     = Omega_6[3, 5]
//   rho_slope_SOC_N = Omega_6[2, 4]        (correlazione slope SOC–N)
//   rho_slope_SOC_P = Omega_6[2, 6]
//   rho_slope_N_P   = Omega_6[4, 6]
//
// X_W (K_W=5): logBottom[1], Texture1[2], Texture2[3], BulkDensity[4], PH[5]
// Ordine righe V:  [SOC_int, SOC_slope, N_int, N_slope, P_int, P_slope]
// TOTALE PARAMETRI: 6J + 57 scalari  (J=40: 297)
//   tau_6[6]: 6; L_6 (Cholesky 6x6): 15; alpha_r: 3; gamma_r[5]: 15; beta_r[4]: 12; sigma_r: 3
//   + z_v[6,J]: 6J
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

  // ── Effetti casuali 6D: [SOC_int, SOC_slope, N_int, N_slope, P_int, P_slope] ─
  matrix[6, J]            z_v;
  cholesky_factor_corr[6] L_6;
  vector<lower=0>[6]      tau_6;

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

  // V[row, j]: riga 1=SOC_int, 2=SOC_slope, 3=N_int, 4=N_slope, 5=P_int, 6=P_slope
  matrix[6, J] V = diag_pre_multiply(tau_6, L_6) * z_v;

}

model {

  // ── Prior: effetti casuali ────────────────────────────────────────────────────
  to_vector(z_v) ~ std_normal();
  L_6   ~ lkj_corr_cholesky(2);
  tau_6 ~ normal(0, 0.5);

  // ── Prior: effetti fissi ──────────────────────────────────────────────────────
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

  // ── Likelihood ───────────────────────────────────────────────────────────────
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

      real mu_SOC = alpha_SOC + V[1,j] + V[2,j]*logB_n + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + V[3,j] + V[4,j]*logB_n + gN[n]   + bN[j];
      real mu_P   = alpha_P   + V[5,j] + V[6,j]*logB_n + gP[n]   + bP[j];

      logSOC[n] ~ normal(mu_SOC, sigma_SOC);
      logN[n]   ~ normal(mu_N,   sigma_N);
      logP[n]   ~ normal(mu_P,   sigma_P);
    }
  }

}

generated quantities {

  // ── Matrice di correlazione completa 6×6 ─────────────────────────────────────
  corr_matrix[6] Omega_6 = multiply_lower_tri_self_transpose(L_6);

  // ── Correlazioni RI–RS per risposta (diagonale di blocchi 2×2) ────────────────
  real rho_SOC = Omega_6[1, 2];   // SOC: int–slope
  real rho_N   = Omega_6[3, 4];   // N:   int–slope
  real rho_P   = Omega_6[5, 6];   // P:   int–slope

  real tau_alpha_SOC = tau_6[1]; real tau_beta_SOC = tau_6[2];
  real tau_alpha_N   = tau_6[3]; real tau_beta_N   = tau_6[4];
  real tau_alpha_P   = tau_6[5]; real tau_beta_P   = tau_6[6];

  // ── Correlazioni cross-risposta degli intercetti ───────────────────────────────
  real rho_int_SOC_N = Omega_6[1, 3];
  real rho_int_SOC_P = Omega_6[1, 5];
  real rho_int_N_P   = Omega_6[3, 5];

  // ── Correlazioni cross-risposta delle slope ───────────────────────────────────
  real rho_slope_SOC_N = Omega_6[2, 4];
  real rho_slope_SOC_P = Omega_6[2, 6];
  real rho_slope_N_P   = Omega_6[4, 6];

  // ── Log-likelihood per LOO-CV ─────────────────────────────────────────────────
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
      real mu_SOC = alpha_SOC + V[1,j] + V[2,j]*logB_n + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + V[3,j] + V[4,j]*logB_n + gN[n]   + bN[j];
      real mu_P   = alpha_P   + V[5,j] + V[6,j]*logB_n + gP[n]   + bP[j];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC, sigma_SOC)
                 + normal_lpdf(logN[n]   | mu_N,   sigma_N)
                 + normal_lpdf(logP[n]   | mu_P,   sigma_P);
    }
  }

}
