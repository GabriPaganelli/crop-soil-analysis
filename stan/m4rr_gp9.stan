// =============================================================================
// m4rr_gp9.stan  (versione finale — "gp9b" rispetto al primo tentativo)
//
// Estensione di m4rr_gp8.stan: rimozione del fattore latente within-field
// con regolarizzazione corretta.
//
// STORIA:
// Prima versione (gp9 iniziale): gamma_N liberi con prior normal(0,1).
//   Risultato: elpd_loo = -515.7, peggiore di gp8 (-509.0).
//   I 4 parametri gamma_N extra non erano sufficientemente regolarizzati.
//   rho_W_SN = 0.041 [-0.085, 0.165]: correlazione within praticamente zero.
//
// CORREZIONI RISPETTO ALLA PRIMA VERSIONE:
//   1. rho_W_SN = 0 (fissato): giustificato da rho≈0.04 in gp9a e dal GAM
//      (script 8: correlazione residua SOC-N = 0.078 dopo conditioning).
//      Rimozione di L_Omega_W (Cholesky) e passaggio a normali indipendenti.
//   2. Prior su gamma_N: normal(0, 1) → normal(0, 0.3).
//      Prior più stretta agisce come regolarizzazione: shrinka gamma_N verso
//      zero quando i dati sono deboli (come per Texture1 e BulkDensity su N).
//      Giustificazione: su scala log standardizzata, effetti > 0.3 log/SD
//      per N within-field sono biologicamente improbabili senza segnale forte.
//
// STRUTTURA FINALE:
//   Within-field: logSOC e logN INDIPENDENTI (rho = 0 giustificato)
//     logSOC[n] ~ normal(mu_SOC[n], sigma_W_SOC)
//     logN[n]   ~ normal(mu_N[n],   sigma_W_N)
//     logP[n]   ~ normal(eta_P_B[j], sqrt(psi_W_P))
//
//   gamma_SOC[K_W]: effetti within su logSOC, prior normal(0, 1)
//   gamma_N[K_W]:   effetti within su logN,   prior normal(0, 0.3) — regolarizzato
//
//   Between-field: invariata da gp8
//     mu_SOC = eta_org_B[j]*(1+b_slope*X_W[n,1]) + alpha_SOC_res[j] + gamma_SOC'*X_W
//     mu_N   = alpha_N[j] + gamma_N'*X_W
//
// Parametri vs gp8: +4 (gamma_SOC) +4 (gamma_N) -4 (gamma_org) -1 (psi_W_org) = +3
//   (rho_W_SN rimosso rispetto a gp9a → netto +3 vs gp8, +0 vs gp9a con rho)
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
  int<lower=1> K_B;         // 4: OnFarm, Irrigate, Fertilised, N_Natural

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W] X_W;
  matrix[J, K_B] X_B;

  matrix[J, J] dist_mat;

  real b_slope;             // da brms LMM (script 10, invariato)

}

parameters {

  // ── Coefficienti within: SOC e N indipendenti -----------------------------
  vector[K_W] gamma_SOC;    // prior normal(0, 1)
  vector[K_W] gamma_N;      // prior normal(0, 0.3) — regolarizzato

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

  // ── Prior: effetti within -------------------------------------------------
  gamma_SOC ~ normal(0, 1);
  // Prior stretto su gamma_N: regolarizza verso 0 quando il segnale e' debole.
  // Equivalente funzionale della regolarizzazione che forniva lambda_N in gp8,
  // senza imporre un rapporto fisso tra SOC e N.
  gamma_N   ~ normal(0, 0.3);

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

  // ── Likelihood: logSOC e logN indipendenti within-field -------------------
  // rho_W_SN = 0 giustificato: gp9a stima 0.041 [-0.085, 0.165]; GAM = 0.078.
  {
    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
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
    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
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
