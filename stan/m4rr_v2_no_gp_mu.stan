// =============================================================================
// m4rr_v2_no_gp_mu.stan  —  v2_no_gp con intercetta globale (serie v2)
//
// STRUTTURA:
//   Come m4rr_v2_no_gp.stan, con aggiunta di intercette globali alpha_r
//   che assorbono la media globale di ogni log-risposta.
//
//   Per ogni risposta r ∈ {SOC, N, P} e osservazione i nel campo j:
//
//   y_r[i,j] ~ N(mu_r[i,j], sigma_r²)
//
//   mu_r[i,j] = alpha_r                  ← intercetta globale (NUOVA)
//             + psi_r · z_nu_r[j]        ← random intercept i.i.d. (deviazione dal globale)
//             + gamma_r · X_W_i          ← effetti within-field fissi
//             + beta_r  · X_B_j          ← effetti between-field (management)
//
// MOTIVAZIONE:
//   In v2_no_gp il random intercept z_nu_r[j] ha prior N(0,1) con media zero,
//   e tutte le covariate sono standardizzate o binarie centrate vicino a zero.
//   Senza alpha_r, la media globale delle log-risposte (~-2 per logN, logP)
//   non ha parametro dove appoggiarsi → beta_r esplode (beta_N[N_Natural]=-2.3).
//   In gp_v2_full il GP degenere (rho >> area) assolveva accidentalmente
//   questa funzione. Qui la modelliamo esplicitamente.
//
// DIFFERENZA DA v2_no_gp:
//   + alpha_r (3 parametri, prior normal(0,2)) → totale 3J + 33 scalari = 153
//
// PARAMETRI: 3J + 33 scalari = 120 + 33 = 153 (con J=40)
//   NCP: z_nu_r[J] × 3 = 3J
//   Scalari: alpha_r, psi_r, gamma_r[4], beta_r[4], sigma_r → 11 × 3 = 33
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

  matrix[N, K_W] X_W;
  matrix[J, K_B] X_B;

}

parameters {

  // ── SOC ──────────────────────────────────────────────────────────────────────
  real          alpha_SOC;      // intercetta globale (NUOVA)
  vector[J]     z_nu_SOC;      // NCP: random intercept i.i.d.
  real<lower=0> psi_SOC;       // ampiezza random intercept
  vector[K_W]   gamma_SOC;     // effetti within-field
  vector[K_B]   beta_SOC;      // effetti between-field (management)
  real<lower=0> sigma_SOC;     // SD residua

  // ── N ────────────────────────────────────────────────────────────────────────
  real          alpha_N;
  vector[J]     z_nu_N;
  real<lower=0> psi_N;
  vector[K_W]   gamma_N;
  vector[K_B]   beta_N;
  real<lower=0> sigma_N;

  // ── P ────────────────────────────────────────────────────────────────────────
  real          alpha_P;
  vector[J]     z_nu_P;
  real<lower=0> psi_P;
  vector[K_W]   gamma_P;
  vector[K_B]   beta_P;
  real<lower=0> sigma_P;

}

model {

  // ── Prior: SOC ───────────────────────────────────────────────────────────────
  alpha_SOC ~ normal(0, 2);
  z_nu_SOC  ~ std_normal();
  psi_SOC   ~ normal(0, 0.5);
  gamma_SOC ~ normal(0, 1);
  beta_SOC  ~ normal(0, 1);
  sigma_SOC ~ normal(0, 0.5);

  // ── Prior: N ─────────────────────────────────────────────────────────────────
  alpha_N ~ normal(0, 2);
  z_nu_N  ~ std_normal();
  psi_N   ~ normal(0, 0.5);
  gamma_N ~ normal(0, 1);
  beta_N  ~ normal(0, 1);
  sigma_N ~ normal(0, 0.5);

  // ── Prior: P ─────────────────────────────────────────────────────────────────
  alpha_P ~ normal(0, 2);
  z_nu_P  ~ std_normal();
  psi_P   ~ normal(0, 0.5);
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

      real mu_SOC = alpha_SOC + psi_SOC * z_nu_SOC[j] + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + psi_N   * z_nu_N[j]   + gN[n]   + bN[j];
      real mu_P   = alpha_P   + psi_P   * z_nu_P[j]   + gP[n]   + bP[j];

      logSOC[n] ~ normal(mu_SOC, sigma_SOC);
      logN[n]   ~ normal(mu_N,   sigma_N);
      logP[n]   ~ normal(mu_P,   sigma_P);
    }
  }

}

generated quantities {

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
      real mu_SOC = alpha_SOC + psi_SOC * z_nu_SOC[j] + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + psi_N   * z_nu_N[j]   + gN[n]   + bN[j];
      real mu_P   = alpha_P   + psi_P   * z_nu_P[j]   + gP[n]   + bP[j];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC, sigma_SOC)
                 + normal_lpdf(logN[n]   | mu_N,   sigma_N)
                 + normal_lpdf(logP[n]   | mu_P,   sigma_P);
    }
  }

}
