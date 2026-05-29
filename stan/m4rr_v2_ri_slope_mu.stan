// =============================================================================
// m4rr_v2_ri_slope_mu.stan  —  RI + slope proporzionale + intercetta globale
//
// STRUTTURA:
//   Come m4rr_v2_ri_slope.stan, con aggiunta di intercette globali alpha_r.
//
//   Per ogni risposta r ∈ {SOC, N, P} e osservazione i nel campo j:
//
//   y_r[i,j] ~ N(mu_r[i,j], sigma_r²)
//
//   mu_r[i,j] = alpha_r                                         ← intercetta globale (NUOVA)
//             + z_nu_r[j] · (psi_r + eta_r · logBottom_i)      ← RI + slope prop.
//             + gamma_r · X_W_i                                 ← effetti within fissi
//             + beta_r  · X_B_j                                 ← effetti between fissi
//
// INTERPRETAZIONE:
//   alpha_r = media globale di log-r (assorbe ~-2 per logN, logP)
//   psi_r   = SD dell'effetto casuale alla profondità media (logBottom = 0)
//   eta_r   = come quella SD cambia con la profondità
//   b_r     = eta_r / psi_r  → coefficiente di proporzionalità slope-intercetta
//
// MOTIVAZIONE:
//   Identica a m4rr_v2_no_gp_mu.stan: il random intercept z_nu[j] ha media
//   zero per costruzione; senza alpha_r i coefficienti management esplodono.
//
// DIFFERENZA DA v2_ri_slope:
//   + alpha_r (3 parametri) → totale 3J + 36 scalari = 156 (con J=40)
//
// PARAMETRI: 3J + 36 scalari = 120 + 36 = 156 (con J=40)
//   NCP: z_nu_r[J] × 3 = 3J
//   Scalari: alpha_r, psi_r, eta_r, gamma_r[4], beta_r[4], sigma_r → 12 × 3 = 36
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
  real          alpha_SOC;      // intercetta globale (NUOVA)
  vector[J]     z_nu_SOC;
  real<lower=0> psi_SOC;
  real          eta_SOC;
  vector[K_W]   gamma_SOC;
  vector[K_B]   beta_SOC;
  real<lower=0> sigma_SOC;

  // ── N ────────────────────────────────────────────────────────────────────────
  real          alpha_N;
  vector[J]     z_nu_N;
  real<lower=0> psi_N;
  real          eta_N;
  vector[K_W]   gamma_N;
  vector[K_B]   beta_N;
  real<lower=0> sigma_N;

  // ── P ────────────────────────────────────────────────────────────────────────
  real          alpha_P;
  vector[J]     z_nu_P;
  real<lower=0> psi_P;
  real          eta_P;
  vector[K_W]   gamma_P;
  vector[K_B]   beta_P;
  real<lower=0> sigma_P;

}

model {

  // ── Prior: SOC ───────────────────────────────────────────────────────────────
  alpha_SOC ~ normal(0, 2);
  z_nu_SOC  ~ std_normal();
  psi_SOC   ~ normal(0, 0.5);
  eta_SOC   ~ normal(0, 0.3);
  gamma_SOC ~ normal(0, 1);
  beta_SOC  ~ normal(0, 1);
  sigma_SOC ~ normal(0, 0.5);

  // ── Prior: N ─────────────────────────────────────────────────────────────────
  alpha_N ~ normal(0, 2);
  z_nu_N  ~ std_normal();
  psi_N   ~ normal(0, 0.5);
  eta_N   ~ normal(0, 0.3);
  gamma_N ~ normal(0, 1);
  beta_N  ~ normal(0, 1);
  sigma_N ~ normal(0, 0.5);

  // ── Prior: P ─────────────────────────────────────────────────────────────────
  alpha_P ~ normal(0, 2);
  z_nu_P  ~ std_normal();
  psi_P   ~ normal(0, 0.5);
  eta_P   ~ normal(0, 0.3);
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
      int j = field_id[n];
      real logB_n = X_W[n, 1];   // logBottom standardizzato

      real mu_SOC = alpha_SOC
                  + z_nu_SOC[j] * (psi_SOC + eta_SOC * logB_n)
                  + gSOC[n] + bSOC[j];

      real mu_N   = alpha_N
                  + z_nu_N[j]   * (psi_N   + eta_N   * logB_n)
                  + gN[n]   + bN[j];

      real mu_P   = alpha_P
                  + z_nu_P[j]   * (psi_P   + eta_P   * logB_n)
                  + gP[n]   + bP[j];

      logSOC[n] ~ normal(mu_SOC, sigma_SOC);
      logN[n]   ~ normal(mu_N,   sigma_N);
      logP[n]   ~ normal(mu_P,   sigma_P);
    }
  }

}

generated quantities {

  // b_r = eta_r / psi_r: coefficiente di proporzionalità slope-intercetta
  real b_SOC = (psi_SOC > 0.02) ? eta_SOC / psi_SOC : 0.0;
  real b_N   = (psi_N   > 0.02) ? eta_N   / psi_N   : 0.0;
  real b_P   = (psi_P   > 0.02) ? eta_P   / psi_P   : 0.0;

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
      int j = field_id[n];
      real logB_n = X_W[n, 1];
      real mu_SOC = alpha_SOC + z_nu_SOC[j] * (psi_SOC + eta_SOC * logB_n) + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + z_nu_N[j]   * (psi_N   + eta_N   * logB_n) + gN[n]   + bN[j];
      real mu_P   = alpha_P   + z_nu_P[j]   * (psi_P   + eta_P   * logB_n) + gP[n]   + bP[j];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC, sigma_SOC)
                 + normal_lpdf(logN[n]   | mu_N,   sigma_N)
                 + normal_lpdf(logP[n]   | mu_P,   sigma_P);
    }
  }

}
