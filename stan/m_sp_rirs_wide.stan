// =============================================================================
// m_sp_rirs_wide.stan  —  M-SP-RIRS con prior larga: tau_beta ~ N(0,1) invece di N(0,0.5)
//
// Per ogni risposta r ∈ {SOC, N, P} e campo j:
//   [u_int_r[j], u_slope_r[j]] ~ MVN(0, Sigma_r)
//   Sigma_r = diag(tau_r) * Omega_r * diag(tau_r),  Omega_r ~ LKJ(2)
//   NCP: u_r = diag_pre_multiply(tau_r, L_r) * z_r,  z_r[2,J] ~ N(0,I)
//
// Struttura della media:
//   mu_r[i,j] = alpha_r
//             + u_r[1, j]                            ← random intercept
//             + u_r[2, j] * logBottom_i              ← random slope su logBottom
//             + gamma_r * X_W_i                      ← effetti fissi within
//             + beta_r  * X_B_j                      ← effetti fissi between
//
// logBottom entra due volte:
//   - gamma_r[1] * logBottom = pendenza media fissa tra campi
//   - u_r[2, j] * logBottom  = deviazione specifica del campo dalla pendenza media
//
// Relazione con M-SP:
//   M-SP è il caso limite |rho_r| → 1 (intercetta e slope perfettamente correlati),
//   con tau_alpha_r = psi_r e tau_beta_r = |eta_r|.
//   M-RI è il caso limite tau_beta_r → 0 (solo intercetta casuale).
//
// Generated quantities: rho_r = Omega_r[1,2] per r ∈ {SOC, N, P}
// Parametri totali: 6J + 39 scalari = 240 + 39 = 279 (con J=40)
//   z_r[2,J] × 3 = 6J;  tau_r[2] × 3 = 6;  L_r (1 param) × 3 = 3;
//   alpha, gamma[4], beta[4], sigma × 3 = 30
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

  // ── Effetti casuali bivariati (NCP): [intercept, slope] per campo ─────────────
  matrix[2, J]            z_SOC;
  cholesky_factor_corr[2] L_SOC;
  vector<lower=0>[2]      tau_SOC;   // [tau_alpha_SOC, tau_beta_SOC]

  matrix[2, J]            z_N;
  cholesky_factor_corr[2] L_N;
  vector<lower=0>[2]      tau_N;

  matrix[2, J]            z_P;
  cholesky_factor_corr[2] L_P;
  vector<lower=0>[2]      tau_P;

  // ── SOC ──────────────────────────────────────────────────────────────────────
  real          alpha_SOC;
  vector[K_W]   gamma_SOC;
  vector[K_B]   beta_SOC;
  real<lower=0> sigma_SOC;

  // ── N ────────────────────────────────────────────────────────────────────────
  real          alpha_N;
  vector[K_W]   gamma_N;
  vector[K_B]   beta_N;
  real<lower=0> sigma_N;

  // ── P ────────────────────────────────────────────────────────────────────────
  real          alpha_P;
  vector[K_W]   gamma_P;
  vector[K_B]   beta_P;
  real<lower=0> sigma_P;

}

transformed parameters {

  // Trasformazione NCP: u_r[:, j] ~ MVN(0, Sigma_r) con Sigma_r = diag(tau)*Omega*diag(tau)
  matrix[2, J] u_SOC = diag_pre_multiply(tau_SOC, L_SOC) * z_SOC;
  matrix[2, J] u_N   = diag_pre_multiply(tau_N,   L_N)   * z_N;
  matrix[2, J] u_P   = diag_pre_multiply(tau_P,   L_P)   * z_P;
  // u_r[1, j] = random intercept del campo j per risposta r
  // u_r[2, j] = random slope su logBottom del campo j per risposta r

}

model {

  // ── Prior: struttura bivariata SOC ───────────────────────────────────────────
  to_vector(z_SOC) ~ std_normal();
  L_SOC ~ lkj_corr_cholesky(2);
  tau_SOC ~ normal(0, 1);

  // ── Prior: struttura bivariata N ─────────────────────────────────────────────
  to_vector(z_N) ~ std_normal();
  L_N ~ lkj_corr_cholesky(2);
  tau_N ~ normal(0, 1);

  // ── Prior: struttura bivariata P ─────────────────────────────────────────────
  to_vector(z_P) ~ std_normal();
  L_P ~ lkj_corr_cholesky(2);
  tau_P ~ normal(0, 1);

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
      real logB_n = X_W[n, 1];   // logBottom standardizzato

      real mu_SOC = alpha_SOC + u_SOC[1, j] + u_SOC[2, j] * logB_n + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + u_N[1, j]   + u_N[2, j]   * logB_n + gN[n]   + bN[j];
      real mu_P   = alpha_P   + u_P[1, j]   + u_P[2, j]   * logB_n + gP[n]   + bP[j];

      logSOC[n] ~ normal(mu_SOC, sigma_SOC);
      logN[n]   ~ normal(mu_N,   sigma_N);
      logP[n]   ~ normal(mu_P,   sigma_P);
    }
  }

}

generated quantities {

  // Matrici di correlazione tra RI e RS per ciascuna risposta
  corr_matrix[2] Omega_SOC = multiply_lower_tri_self_transpose(L_SOC);
  corr_matrix[2] Omega_N   = multiply_lower_tri_self_transpose(L_N);
  corr_matrix[2] Omega_P   = multiply_lower_tri_self_transpose(L_P);

  real rho_SOC = Omega_SOC[1, 2];   // correlazione RI–RS per SOC
  real rho_N   = Omega_N[1, 2];     // correlazione RI–RS per N
  real rho_P   = Omega_P[1, 2];     // correlazione RI–RS per P

  // Deviazioni standard dei RE
  real tau_alpha_SOC = tau_SOC[1];   real tau_beta_SOC = tau_SOC[2];
  real tau_alpha_N   = tau_N[1];     real tau_beta_N   = tau_N[2];
  real tau_alpha_P   = tau_P[1];     real tau_beta_P   = tau_P[2];

  // Log-likelihood per LOO-CV
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
      real mu_SOC = alpha_SOC + u_SOC[1, j] + u_SOC[2, j] * logB_n + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + u_N[1, j]   + u_N[2, j]   * logB_n + gN[n]   + bN[j];
      real mu_P   = alpha_P   + u_P[1, j]   + u_P[2, j]   * logB_n + gP[n]   + bP[j];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC, sigma_SOC)
                 + normal_lpdf(logN[n]   | mu_N,   sigma_N)
                 + normal_lpdf(logP[n]   | mu_P,   sigma_P);
    }
  }

}
