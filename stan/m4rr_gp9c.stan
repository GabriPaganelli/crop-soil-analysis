// =============================================================================
// m4rr_gp9c.stan  (sparsificazione di gp9b)
//
// Estensione di m4rr_gp9.stan: azzeramento dei gamma_N non significativi
// basato sugli intervalli di credibilità di gp9b.
//
// MOTIVAZIONE (da gp9b, elpd_loo = -514.8 vs gp8 = -509.0):
//   gamma_N[1] (logBottom → N):   -0.261, CI90 [-0.316, -0.207] — SIGN.
//   gamma_N[2] (Texture1 → N):    -0.044, CI90 [-0.098,  0.010] — include 0 → RIMOSSO
//   gamma_N[3] (Texture2 → N):    -0.222, CI90 [-0.332, -0.114] — SIGN.
//   gamma_N[4] (BulkDensity → N): -0.057, CI90 [-0.177,  0.064] — include 0 → RIMOSSO
//
// MODIFICA RISPETTO A gp9b:
//   - gamma_N è ora vector[K_W_N] con K_W_N = 2 (solo logBottom e Texture2)
//   - X_W_N = matrice [N, 2] con colonne 1 e 3 di X_W (passata da R)
//   - gamma_SOC rimane vector[K_W] = 4 elementi, invariato
//   - Nessun parametro "fissato": i termini nulli non entrano nel blocco parameters
//
// PARAMETRI vs gp9b: -2 (gamma_N passa da 4 a 2 elementi liberi)
// PARAMETRI vs gp8:  +2 netto (4 gamma_SOC + 2 gamma_N_free vs 4 gamma_org)
// =============================================================================

functions {

  matrix matern32_cov(matrix dist_mat, real sigma, real rho) {
    int J = rows(dist_mat);
    matrix[J, J] K;
    real s2 = square(sigma);
    real sqrt3_over_rho = sqrt(3.0) / rho;
    for (i in 1:J) {
      K[i, i] = s2;
      for (j_idx in (i + 1):J) {
        real r = sqrt3_over_rho * dist_mat[i, j_idx];
        real k = s2 * (1.0 + r) * exp(-r);
        K[i, j_idx] = k;
        K[j_idx, i] = k;
      }
    }
    return K;
  }

}

data {

  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K_W;         // 4: logBottom, Texture1, Texture2, BulkDensity
  int<lower=1> K_W_N;       // 2: logBottom, Texture2 (solo effetti sign. su N)
  int<lower=1> K_B;         // 4: OnFarm, Irrigate, Fertilised, N_Natural

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W]   X_W;    // design matrix completa per logSOC (4 colonne)
  matrix[N, K_W_N] X_W_N;  // design matrix ridotta per logN (2 colonne: [1] e [3])
  matrix[J, K_B]   X_B;

  matrix[J, J] dist_mat;

  real b_slope;             // da brms LMM (script 10, invariato)

}

parameters {

  // ── Coefficienti within SOC: tutti e 4 i predittori ──────────────────────
  vector[K_W]   gamma_SOC;    // prior normal(0, 1)

  // ── Coefficienti within N: solo logBottom e Texture2 ─────────────────────
  // Texture1 e BulkDensity fissati a 0 (CI90 che include 0 in gp9b):
  //   non appaiono nel blocco parameters → vera riduzione parametrica.
  vector[K_W_N] gamma_N;      // prior normal(0, 0.3) — 2 elementi

  // ── Coefficienti between --------------------------------------------------
  vector[K_B] gamma_P_B;
  vector[K_B] gamma_N_B;

  // ── Varianze/SD within ----------------------------------------------------
  real<lower=0> sigma_W_SOC;
  real<lower=0> sigma_W_N;
  real<lower=0> psi_W_P;

  // ── Fattore between organico (NCP) ----------------------------------------
  vector[J] z_eta_org_B;

  // ── Intercette between (NCP) ----------------------------------------------
  vector[J] alpha_SOC_raw;
  vector[J] alpha_N_raw;

  // ── Fattore between P (NCP) -----------------------------------------------
  vector[J] z_P_B;

  // ── Varianze between ------------------------------------------------------
  real<lower=0> sigma_B_SOC;
  real<lower=0> sigma_B_N;
  real<lower=0> psi_B_org;
  real<lower=0> sigma_P_between;

  // ── GP organico -----------------------------------------------------------
  real<lower=0> sigma_GP_org;
  real<lower=0> rho_org;

}

transformed parameters {

  // ── Kernel GP organico ----------------------------------------------------
  matrix[J, J] L_org;
  {
    matrix[J, J] K_org = matern32_cov(dist_mat, sigma_GP_org, rho_org);
    for (j in 1:J) K_org[j, j] += psi_B_org + 1e-8;
    L_org = cholesky_decompose(K_org);
  }

  // ── Fattori between -------------------------------------------------------
  vector[J] eta_org_B    = L_org * z_eta_org_B;
  vector[J] eta_P_B      = X_B * gamma_P_B + sigma_P_between * z_P_B;

  // ── Intercette between ----------------------------------------------------
  vector[J] alpha_SOC_res = sigma_B_SOC * alpha_SOC_raw;
  vector[J] alpha_N       = X_B * gamma_N_B + sigma_B_N * alpha_N_raw;

}

model {

  // ── Prior: effetti within SOC (tutti e 4) ---------------------------------
  gamma_SOC ~ normal(0, 1);

  // ── Prior: effetti within N (solo logBottom e Texture2) ───────────────────
  gamma_N ~ normal(0, 0.3);

  // ── Prior: effetti between ------------------------------------------------
  gamma_P_B ~ normal(0, 1);
  gamma_N_B ~ normal(0, 1);

  // ── Prior: varianze within ------------------------------------------------
  sigma_W_SOC ~ normal(0, 0.5);
  sigma_W_N   ~ normal(0, 0.5);
  psi_W_P     ~ normal(0, 0.5);

  // ── Prior: NCP -----------------------------------------------------------
  z_eta_org_B   ~ std_normal();
  z_P_B         ~ std_normal();
  alpha_SOC_raw ~ std_normal();
  alpha_N_raw   ~ std_normal();

  // ── Prior: varianze between -----------------------------------------------
  sigma_B_SOC     ~ normal(0, 0.5);
  sigma_B_N       ~ normal(0, 0.5);
  psi_B_org       ~ normal(0, 1.0);
  sigma_P_between ~ normal(0, 1.0);

  // ── Prior: GP Matern 3/2 -------------------------------------------------
  sigma_GP_org ~ normal(0, 1.0);
  rho_org      ~ inv_gamma(3.0, 10.0);

  // ── Likelihood ------------------------------------------------------------
  {
    vector[N] gSOC = X_W   * gamma_SOC;   // 4 predittori
    vector[N] gN   = X_W_N * gamma_N;     // 2 predittori (logBottom + Texture2)
    for (n in 1:N) {
      int j = field_id[n];
      real mu_SOC_n = eta_org_B[j] * (1.0 + b_slope * X_W[n, 1])
                    + alpha_SOC_res[j] + gSOC[n];
      real mu_N_n   = alpha_N[j] + gN[n];
      logSOC[n] ~ normal(mu_SOC_n, sigma_W_SOC);
      logN[n]   ~ normal(mu_N_n,   sigma_W_N);
      logP[n]   ~ normal(eta_P_B[j], sqrt(psi_W_P));
    }
  }

}

generated quantities {

  vector[J] eta_org_B_out     = eta_org_B;
  vector[J] eta_P_B_out       = eta_P_B;
  vector[J] alpha_N_out       = alpha_N;
  vector[J] alpha_SOC_res_out = alpha_SOC_res;

  vector[N] log_lik;
  {
    vector[N] gSOC = X_W   * gamma_SOC;
    vector[N] gN   = X_W_N * gamma_N;
    for (n in 1:N) {
      int j = field_id[n];
      real mu_SOC_n = eta_org_B[j] * (1.0 + b_slope * X_W[n, 1])
                    + alpha_SOC_res[j] + gSOC[n];
      real mu_N_n   = alpha_N[j] + gN[n];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC_n, sigma_W_SOC)
                 + normal_lpdf(logN[n]   | mu_N_n,   sigma_W_N)
                 + normal_lpdf(logP[n]   | eta_P_B[j], sqrt(psi_W_P));
    }
  }

}
