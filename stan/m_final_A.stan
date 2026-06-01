// =============================================================================
// m_final_A.stan  —  Opzione A: M-SP ridotto (no management, no Texture1)
//
// Identico a m4rr_v2_ri_slope_mu.stan con due sole differenze:
//   1. X_B e beta_r rimossi (nessun effetto between-field di management)
//   2. K_W = 3: logBottom (col 1), Texture2 (col 2), BulkDensity (col 3)
//      — Texture1 rimossa sulla base della selezione projpred
//
// Struttura per r ∈ {SOC, N, P}:
//   mu_r[i,j] = alpha_r
//             + z_nu_r[j] · (psi_r + eta_r · X_W[i,1])   ← slope proporzionale
//             + gamma_r · X_W[i,]                          ← effetti fissi within
//
// Stessi prior di m4rr_v2_ri_slope_mu.stan.
// =============================================================================

data {

  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K_W;   // 3: logBottom, Texture2, BulkDensity

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W] X_W;   // col 1 = logBottom standardizzato

}

parameters {

  // ── SOC ──────────────────────────────────────────────────────────────────────
  real          alpha_SOC;
  vector[J]     z_nu_SOC;
  real<lower=0> psi_SOC;
  real          eta_SOC;
  vector[K_W]   gamma_SOC;
  real<lower=0> sigma_SOC;

  // ── N ────────────────────────────────────────────────────────────────────────
  real          alpha_N;
  vector[J]     z_nu_N;
  real<lower=0> psi_N;
  real          eta_N;
  vector[K_W]   gamma_N;
  real<lower=0> sigma_N;

  // ── P ────────────────────────────────────────────────────────────────────────
  real          alpha_P;
  vector[J]     z_nu_P;
  real<lower=0> psi_P;
  real          eta_P;
  vector[K_W]   gamma_P;
  real<lower=0> sigma_P;

}

model {

  // ── Prior: SOC ───────────────────────────────────────────────────────────────
  alpha_SOC ~ normal(0, 2);
  z_nu_SOC  ~ std_normal();
  psi_SOC   ~ normal(0, 0.5);
  eta_SOC   ~ normal(0, 0.3);
  gamma_SOC ~ normal(0, 1);
  sigma_SOC ~ normal(0, 0.5);

  // ── Prior: N ─────────────────────────────────────────────────────────────────
  alpha_N ~ normal(0, 2);
  z_nu_N  ~ std_normal();
  psi_N   ~ normal(0, 0.5);
  eta_N   ~ normal(0, 0.3);
  gamma_N ~ normal(0, 1);
  sigma_N ~ normal(0, 0.5);

  // ── Prior: P ─────────────────────────────────────────────────────────────────
  alpha_P ~ normal(0, 2);
  z_nu_P  ~ std_normal();
  psi_P   ~ normal(0, 0.5);
  eta_P   ~ normal(0, 0.3);
  gamma_P ~ normal(0, 1);
  sigma_P ~ normal(0, 0.5);

  // ── Likelihood ───────────────────────────────────────────────────────────────
  {
    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
    vector[N] gP   = X_W * gamma_P;

    for (n in 1:N) {
      int  j      = field_id[n];
      real logB_n = X_W[n, 1];

      real mu_SOC = alpha_SOC + z_nu_SOC[j] * (psi_SOC + eta_SOC * logB_n) + gSOC[n];
      real mu_N   = alpha_N   + z_nu_N[j]   * (psi_N   + eta_N   * logB_n) + gN[n];
      real mu_P   = alpha_P   + z_nu_P[j]   * (psi_P   + eta_P   * logB_n) + gP[n];

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

    for (n in 1:N) {
      int  j      = field_id[n];
      real logB_n = X_W[n, 1];
      real mu_SOC = alpha_SOC + z_nu_SOC[j] * (psi_SOC + eta_SOC * logB_n) + gSOC[n];
      real mu_N   = alpha_N   + z_nu_N[j]   * (psi_N   + eta_N   * logB_n) + gN[n];
      real mu_P   = alpha_P   + z_nu_P[j]   * (psi_P   + eta_P   * logB_n) + gP[n];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC, sigma_SOC)
                 + normal_lpdf(logN[n]   | mu_N,   sigma_N)
                 + normal_lpdf(logP[n]   | mu_P,   sigma_P);
    }
  }

}
