// =============================================================================
// m_sp_mv_re.stan  —  M-SP con effetti casuali multivariati (Omega_RE 3×3)
//
// Estensione di m4rr_v2_ri_slope_mu.stan in cui i tre effetti casuali di campo
// (per SOC, N, P) sono modellati congiuntamente come MVN con correlazione
// Omega_RE ~ LKJ(2).
//
// Parametrizzazione NCP:
//   z_raw[3, J] ~ N(0, I)                          (raw standard normals)
//   z_nu[:, j]  = L_RE * z_raw[:, j]               (Cholesky NCP)
//   Omega_RE    = L_RE * L_RE'  ~ LKJ(2)
//
// Struttura:
//   mu_r[i,j] = alpha_r
//             + z_nu[r, j] * (psi_r + eta_r * logBottom_i)
//             + gamma_r * X_W_i
//             + beta_r  * X_B_j
//
// Differenza cruciale da M-SP-MV:
//   M-SP-MV:    correlazione residua epsilon (livello osservazione, dentro campo)
//   M-SP-MVRE:  correlazione tra RE di campo (livello campo, tra campi)
//
// Generated quantities: rho_SOC_N_RE, rho_SOC_P_RE, rho_N_P_RE
// Parametri totali: 3J + 36 scalari + 3 correlazioni = 159 (con J=40)
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

  // Effetti casuali multivariati (NCP)
  matrix[3, J]            z_raw;   // [SOC; N; P] × J: raw standard normals
  cholesky_factor_corr[3] L_RE;    // Cholesky di Omega_RE

  // ── SOC ──────────────────────────────────────────────────────────────────────
  real          alpha_SOC;
  real<lower=0> psi_SOC;
  real          eta_SOC;
  vector[K_W]   gamma_SOC;
  vector[K_B]   beta_SOC;
  real<lower=0> sigma_SOC;

  // ── N ────────────────────────────────────────────────────────────────────────
  real          alpha_N;
  real<lower=0> psi_N;
  real          eta_N;
  vector[K_W]   gamma_N;
  vector[K_B]   beta_N;
  real<lower=0> sigma_N;

  // ── P ────────────────────────────────────────────────────────────────────────
  real          alpha_P;
  real<lower=0> psi_P;
  real          eta_P;
  vector[K_W]   gamma_P;
  vector[K_B]   beta_P;
  real<lower=0> sigma_P;

}

model {

  // ── Prior: struttura multivariata ─────────────────────────────────────────────
  to_vector(z_raw) ~ std_normal();
  L_RE ~ lkj_corr_cholesky(2);

  // ── Prior: SOC ───────────────────────────────────────────────────────────────
  alpha_SOC ~ normal(0, 2);
  psi_SOC   ~ normal(0, 0.5);
  eta_SOC   ~ normal(0, 0.3);
  gamma_SOC ~ normal(0, 1);
  beta_SOC  ~ normal(0, 1);
  sigma_SOC ~ normal(0, 0.5);

  // ── Prior: N ─────────────────────────────────────────────────────────────────
  alpha_N ~ normal(0, 2);
  psi_N   ~ normal(0, 0.5);
  eta_N   ~ normal(0, 0.3);
  gamma_N ~ normal(0, 1);
  beta_N  ~ normal(0, 1);
  sigma_N ~ normal(0, 0.5);

  // ── Prior: P ─────────────────────────────────────────────────────────────────
  alpha_P ~ normal(0, 2);
  psi_P   ~ normal(0, 0.5);
  eta_P   ~ normal(0, 0.3);
  gamma_P ~ normal(0, 1);
  beta_P  ~ normal(0, 1);
  sigma_P ~ normal(0, 0.5);

  // ── Likelihood ───────────────────────────────────────────────────────────────
  {
    matrix[3, J] z_nu = L_RE * z_raw;   // trasformazione NCP: z_nu[:, j] ~ MVN(0, Omega_RE)

    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
    vector[N] gP   = X_W * gamma_P;

    vector[J] bSOC = X_B * beta_SOC;
    vector[J] bN   = X_B * beta_N;
    vector[J] bP   = X_B * beta_P;

    for (n in 1:N) {
      int  j      = field_id[n];
      real logB_n = X_W[n, 1];

      real mu_SOC = alpha_SOC + z_nu[1, j] * (psi_SOC + eta_SOC * logB_n) + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + z_nu[2, j] * (psi_N   + eta_N   * logB_n) + gN[n]   + bN[j];
      real mu_P   = alpha_P   + z_nu[3, j] * (psi_P   + eta_P   * logB_n) + gP[n]   + bP[j];

      logSOC[n] ~ normal(mu_SOC, sigma_SOC);
      logN[n]   ~ normal(mu_N,   sigma_N);
      logP[n]   ~ normal(mu_P,   sigma_P);
    }
  }

}

generated quantities {

  // Matrice di correlazione tra gli effetti casuali di campo
  corr_matrix[3] Omega_RE   = multiply_lower_tri_self_transpose(L_RE);
  real rho_SOC_N_RE = Omega_RE[1, 2];
  real rho_SOC_P_RE = Omega_RE[1, 3];
  real rho_N_P_RE   = Omega_RE[2, 3];

  real b_SOC = (psi_SOC > 0.02) ? eta_SOC / psi_SOC : 0.0;
  real b_N   = (psi_N   > 0.02) ? eta_N   / psi_N   : 0.0;
  real b_P   = (psi_P   > 0.02) ? eta_P   / psi_P   : 0.0;

  vector[N] log_lik;
  {
    matrix[3, J] z_nu = L_RE * z_raw;

    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
    vector[N] gP   = X_W * gamma_P;
    vector[J] bSOC = X_B * beta_SOC;
    vector[J] bN   = X_B * beta_N;
    vector[J] bP   = X_B * beta_P;

    for (n in 1:N) {
      int  j      = field_id[n];
      real logB_n = X_W[n, 1];
      real mu_SOC = alpha_SOC + z_nu[1, j] * (psi_SOC + eta_SOC * logB_n) + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + z_nu[2, j] * (psi_N   + eta_N   * logB_n) + gN[n]   + bN[j];
      real mu_P   = alpha_P   + z_nu[3, j] * (psi_P   + eta_P   * logB_n) + gP[n]   + bP[j];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC, sigma_SOC)
                 + normal_lpdf(logN[n]   | mu_N,   sigma_N)
                 + normal_lpdf(logP[n]   | mu_P,   sigma_P);
    }
  }

}
