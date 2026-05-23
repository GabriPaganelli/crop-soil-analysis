// =============================================================================
// m4rr_gp10.stan  (sparsificazione di gp9c — between-field + SOC residuo)
//
// MODIFICHE RISPETTO A gp9c:
//
// 1. Between-field sparsificato:
//    In gp9c (e tutti i modelli precedenti) gamma_N_B e gamma_P_B avevano 4
//    elementi (OnFarm, Irrigate, Fertilised, N_Natural). Da gp9c:
//      gamma_N_B[2] (Irrigate):    -0.218, CI90 [-0.709, 0.261] — include 0
//      gamma_N_B[3] (Fertilised):  -0.319, CI90 [-0.887, 0.242] — include 0
//      gamma_P_B[2] (Irrigate):     0.275, CI90 [-0.306, 0.846] — include 0
//      gamma_P_B[3] (Fertilised):  -0.464, CI90 [-1.10,  0.173] — include 0
//    Rimossi: gamma_N_B e gamma_P_B diventano vector[K_B] con K_B = 2.
//    X_B ridotta a [J x 2]: solo OnFarm e N_Natural.
//    (N_Natural domina: coeff ~-2.3 to -2.6; OnFarm bordeline su N ma plausibile)
//
// 2. Rimozione di alpha_SOC_res:
//    In gp9c: sigma_B_SOC = 0.097 (mediana 0.084) — minuscola.
//    alpha_SOC_res[j] = sigma_B_SOC * alpha_SOC_raw[j] contribuisce pochissimo
//    alla media di logSOC, ma aggiunge J + 1 = 55 parametri al modello.
//    Il GP organico (eta_org_B) cattura già la variabilità between-SOC.
//    Rimossi: sigma_B_SOC, alpha_SOC_raw[J].
//    mu_SOC_n = eta_org_B[j] * (1 + b_slope * X_W[n,1]) + gSOC[n]
//
// STRUTTURA WITHIN-FIELD: identica a gp9c (invariata)
//   gamma_SOC[K_W=4]: prior normal(0, 1)
//   gamma_N[K_W_N=2]: prior normal(0, 0.3) — solo logBottom e Texture2
//   X_W [N x 4], X_W_N [N x 2]
//   rho_W_SN = 0 (justified by gp9a and GAM)
//
// PARAMETRI vs gp9c: -4 (gamma_N_B[2,3] + gamma_P_B[2,3]) -J -1 (alpha_SOC)
//   = -5 - J parametri (con J=54: circa 59 parametri in meno)
// PARAMETRI vs gp8: stesso ordine di grandezza (netto circa +2)
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
  int<lower=1> K_W_N;       // 2: logBottom, Texture2
  int<lower=1> K_B;         // 2: OnFarm, N_Natural (ridotto da 4)

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W]   X_W;    // design matrix within: logSOC (4 col.)
  matrix[N, K_W_N] X_W_N;  // design matrix within: logN (2 col.: logBottom, Texture2)
  matrix[J, K_B]   X_B;    // design matrix between: OnFarm, N_Natural (2 col.)

  matrix[J, J] dist_mat;

  real b_slope;             // da brms LMM (script 10, invariato)

}

parameters {

  // ── Coefficienti within SOC ───────────────────────────────────────────────
  vector[K_W]   gamma_SOC;    // prior normal(0, 1)

  // ── Coefficienti within N (logBottom + Texture2 only) ────────────────────
  vector[K_W_N] gamma_N;      // prior normal(0, 0.3)

  // ── Coefficienti between: OnFarm + N_Natural only ─────────────────────────
  vector[K_B] gamma_P_B;      // prior normal(0, 1)
  vector[K_B] gamma_N_B;      // prior normal(0, 1)

  // ── Varianze/SD within ────────────────────────────────────────────────────
  real<lower=0> sigma_W_SOC;
  real<lower=0> sigma_W_N;
  real<lower=0> psi_W_P;

  // ── Fattore between organico (NCP) ────────────────────────────────────────
  // Nota: alpha_SOC_res rimosso — eta_org_B cattura tutta la variabilità SOC between.
  vector[J] z_eta_org_B;

  // ── Intercetta between N (NCP) ────────────────────────────────────────────
  vector[J] alpha_N_raw;

  // ── Fattore between P (NCP) ───────────────────────────────────────────────
  vector[J] z_P_B;

  // ── Varianze between ──────────────────────────────────────────────────────
  real<lower=0> sigma_B_N;
  real<lower=0> psi_B_org;
  real<lower=0> sigma_P_between;

  // ── GP organico ───────────────────────────────────────────────────────────
  real<lower=0> sigma_GP_org;
  real<lower=0> rho_org;

}

transformed parameters {

  // ── Kernel GP organico ────────────────────────────────────────────────────
  matrix[J, J] L_org;
  {
    matrix[J, J] K_org = matern32_cov(dist_mat, sigma_GP_org, rho_org);
    for (j in 1:J) K_org[j, j] += psi_B_org + 1e-8;
    L_org = cholesky_decompose(K_org);
  }

  // ── Fattori between ───────────────────────────────────────────────────────
  vector[J] eta_org_B = L_org * z_eta_org_B;
  vector[J] eta_P_B   = X_B * gamma_P_B + sigma_P_between * z_P_B;

  // ── Intercetta between N ──────────────────────────────────────────────────
  vector[J] alpha_N = X_B * gamma_N_B + sigma_B_N * alpha_N_raw;

}

model {

  // ── Prior: effetti within ─────────────────────────────────────────────────
  gamma_SOC ~ normal(0, 1);
  gamma_N   ~ normal(0, 0.3);

  // ── Prior: effetti between (OnFarm + N_Natural) ───────────────────────────
  gamma_P_B ~ normal(0, 1);
  gamma_N_B ~ normal(0, 1);

  // ── Prior: varianze within ────────────────────────────────────────────────
  sigma_W_SOC ~ normal(0, 0.5);
  sigma_W_N   ~ normal(0, 0.5);
  psi_W_P     ~ normal(0, 0.5);

  // ── Prior: NCP ────────────────────────────────────────────────────────────
  z_eta_org_B ~ std_normal();
  z_P_B       ~ std_normal();
  alpha_N_raw ~ std_normal();

  // ── Prior: varianze between ───────────────────────────────────────────────
  sigma_B_N       ~ normal(0, 0.5);
  psi_B_org       ~ normal(0, 1.0);
  sigma_P_between ~ normal(0, 1.0);

  // ── Prior: GP Matern 3/2 ─────────────────────────────────────────────────
  sigma_GP_org ~ normal(0, 1.0);
  rho_org      ~ inv_gamma(3.0, 10.0);

  // ── Likelihood ────────────────────────────────────────────────────────────
  // mu_SOC_n: alpha_SOC_res rimosso — solo GP + pendenza depth
  {
    vector[N] gSOC = X_W   * gamma_SOC;
    vector[N] gN   = X_W_N * gamma_N;
    for (n in 1:N) {
      int j = field_id[n];
      real mu_SOC_n = eta_org_B[j] * (1.0 + b_slope * X_W[n, 1]) + gSOC[n];
      real mu_N_n   = alpha_N[j] + gN[n];
      logSOC[n] ~ normal(mu_SOC_n, sigma_W_SOC);
      logN[n]   ~ normal(mu_N_n,   sigma_W_N);
      logP[n]   ~ normal(eta_P_B[j], sqrt(psi_W_P));
    }
  }

}

generated quantities {

  vector[J] eta_org_B_out = eta_org_B;
  vector[J] eta_P_B_out   = eta_P_B;
  vector[J] alpha_N_out   = alpha_N;

  vector[N] log_lik;
  {
    vector[N] gSOC = X_W   * gamma_SOC;
    vector[N] gN   = X_W_N * gamma_N;
    for (n in 1:N) {
      int j = field_id[n];
      real mu_SOC_n = eta_org_B[j] * (1.0 + b_slope * X_W[n, 1]) + gSOC[n];
      real mu_N_n   = alpha_N[j] + gN[n];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC_n, sigma_W_SOC)
                 + normal_lpdf(logN[n]   | mu_N_n,   sigma_W_N)
                 + normal_lpdf(logP[n]   | eta_P_B[j], sqrt(psi_W_P));
    }
  }

}
