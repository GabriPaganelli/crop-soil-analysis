// =============================================================================
// m_sp_shared_re.stan  —  M-SP con effetto casuale condiviso (unico z_nu[j])
//
// Variante di m4rr_v2_ri_slope_mu.stan in cui le tre risposte (SOC, N, P)
// condividono un unico effetto casuale di campo z_nu[j] ~ N(0, 1).
//
// Struttura:
//   mu_r[i,j] = alpha_r
//             + z_nu[j] * (psi_r + eta_r * logBottom_i)
//             + gamma_r * X_W_i
//             + beta_r  * X_B_j
//
// Interpretazione:
//   z_nu[j] = fattore latente unico di "qualità del suolo" del campo j.
//   I ranking dei campi sono identici per SOC, N, P (correlazione RE = 1).
//   psi_r e eta_r rimangono risposta-specifici: la scala e il profilo
//   verticale restano liberi; solo il fattore di campo è condiviso.
//
// Più parsimonioso di M-SP: J invece di 3J effetti casuali.
// Parametri totali: J + 36 scalari = 40 + 36 = 76 (con J=40)
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

  // Effetto casuale CONDIVISO tra SOC, N, P
  vector[J]     z_nu;

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

  z_nu ~ std_normal();

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
    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
    vector[N] gP   = X_W * gamma_P;

    vector[J] bSOC = X_B * beta_SOC;
    vector[J] bN   = X_B * beta_N;
    vector[J] bP   = X_B * beta_P;

    for (n in 1:N) {
      int  j      = field_id[n];
      real logB_n = X_W[n, 1];

      real mu_SOC = alpha_SOC + z_nu[j] * (psi_SOC + eta_SOC * logB_n) + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + z_nu[j] * (psi_N   + eta_N   * logB_n) + gN[n]   + bN[j];
      real mu_P   = alpha_P   + z_nu[j] * (psi_P   + eta_P   * logB_n) + gP[n]   + bP[j];

      logSOC[n] ~ normal(mu_SOC, sigma_SOC);
      logN[n]   ~ normal(mu_N,   sigma_N);
      logP[n]   ~ normal(mu_P,   sigma_P);
    }
  }

}

generated quantities {

  real b_SOC = (psi_SOC > 0.02) ? eta_SOC / psi_SOC : 0.0;
  real b_N   = (psi_N   > 0.02) ? eta_N   / psi_N   : 0.0;
  real b_P   = (psi_P   > 0.02) ? eta_P   / psi_P   : 0.0;

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
      real mu_SOC = alpha_SOC + z_nu[j] * (psi_SOC + eta_SOC * logB_n) + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + z_nu[j] * (psi_N   + eta_N   * logB_n) + gN[n]   + bN[j];
      real mu_P   = alpha_P   + z_nu[j] * (psi_P   + eta_P   * logB_n) + gP[n]   + bP[j];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC, sigma_SOC)
                 + normal_lpdf(logN[n]   | mu_N,   sigma_N)
                 + normal_lpdf(logP[n]   | mu_P,   sigma_P);
    }
  }

}
