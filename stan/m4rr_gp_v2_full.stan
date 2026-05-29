// =============================================================================
// m4rr_gp_v2_full.stan  —  Modello simmetrico a tre risposte (serie v2)
//
// STRUTTURA CONCETTUALE:
//   Per ogni risposta r ∈ {SOC, N, P} e osservazione i nel campo j:
//
//   y_r[i,j] ~ N(mu_r[i,j], sigma_r²)
//
//   mu_r[i,j] = f_sp_r[j] · (sigma_GP_r + eta_r · logB_i)   ← GP × profondità
//             + psi_r · z_nu_r[j]                            ← random intercept i.i.d.
//             + gamma_r · X_W_i                              ← within-field fissi
//             + beta_r · X_B_j                               ← between-field management
//
//   f_sp_r ~ GP(0, C_Matérn(rho_r))   [correlazione spaziale, ampiezza = 1]
//
// REPARAMETRIZZAZIONE CHIAVE:
//   Invece di stimare b_r (bilineare con sigma_GP_r), si stimano:
//     sigma_GP_r  = ampiezza GP alla profondità media (logB = 0)
//     eta_r       = b_r * sigma_GP_r  [variazione dell'ampiezza per SD di profondità]
//   da cui: b_r = eta_r / sigma_GP_r  [quantità derivata, in generated quantities]
//
//   I due parametri sono quasi ortogonali nel posteriore:
//     sigma_GP_r: identificato dalla varianza tra campi alla profondità media
//     eta_r:      identificato da come quella varianza cambia con la profondità
//
// NOVITÀ RISPETTO A gp9c:
//   1. GP per N e P (prima assenti)
//   2. eta_r stimato per tutte le risposte (b_SOC non più fisso da LMM)
//   3. beta_SOC: X_B → SOC (prima assente)
//   4. gamma_P: within-field per P (prima assente)
//   5. gamma_N completo (4 pred., non sparsificato a priori)
//   6. psi_r: nugget i.i.d. separato dal GP per ogni risposta
//
// DEGENERAZIONI ATTESE (da preview):
//   - sigma_GP_N → 0: N non ha struttura spaziale → z_sp_N, eta_N non identificati
//   - sigma_GP_P → 0: P probabilmente analogo a N
//   - eta_N, eta_P → 0: profilo verticale di N, P non dipende dal campo
//   - beta_SOC → 0: management non predice SOC beyond GP
//   - gamma_P → 0: P non varia with-field con le covariate fisiche
//   Se sigma_GP_r → 0: ESS basso per z_sp_r, eta_r → li si rimuove nel modello sparse.
//
// PARAMETRI: 6J + 39 scalari = 6×54 + 39 = 363 (con J=54)
//   Random: z_sp_r[J] + z_nu_r[J] per r ∈ {SOC,N,P} → 6J
//   Scalari per risposta: sigma_GP, eta, psi, rho, gamma[4], beta[4], sigma → 13 × 3 = 39
//
// CONFRONTO: gp8 / gp9c / gp10 / gp_v2_full → poi gp_v2_sparse
// =============================================================================

functions {

  // Matrice di correlazione Matérn 3/2 (ampiezza = 1, no nugget)
  matrix matern32_corr(matrix dist_mat, real rho) {
    int J = rows(dist_mat);
    matrix[J, J] C;
    real sqrt3_over_rho = sqrt(3.0) / rho;
    for (i in 1:J) {
      C[i, i] = 1.0;
      for (j in (i + 1):J) {
        real r = sqrt3_over_rho * dist_mat[i, j];
        real c = (1.0 + r) * exp(-r);
        C[i, j] = c;
        C[j, i] = c;
      }
    }
    return C;
  }

}

data {

  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K_W;   // 4: logBottom, Texture1, Texture2, BulkDensity
  int<lower=1> K_B;   // 4: OnFarm, Irrigate, Fertilised, N_Natural

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W] X_W;   // colonna 1 = logBottom (standardizzato), usato per eta
  matrix[J, K_B] X_B;

  matrix[J, J] dist_mat;

}

parameters {

  // ── SOC ──────────────────────────────────────────────────────────────────────
  vector[J] z_sp_SOC;           // NCP spaziale:  f_sp_SOC = L_corr_SOC * z_sp_SOC
  vector[J] z_nu_SOC;           // NCP i.i.d.:    contributo = psi_SOC * z_nu_SOC
  real<lower=0> sigma_GP_SOC;   // ampiezza GP a profondità media
  real         eta_SOC;         // = b_SOC * sigma_GP_SOC (può essere neg. in principio)
  real<lower=0> psi_SOC;        // ampiezza random intercept i.i.d.
  real<lower=0> rho_SOC;        // range correlazione spaziale (km)
  vector[K_W]  gamma_SOC;       // effetti within-field
  vector[K_B]  beta_SOC;        // effetti between-field (management → SOC)
  real<lower=0> sigma_SOC;      // SD residua within

  // ── N ────────────────────────────────────────────────────────────────────────
  vector[J] z_sp_N;
  vector[J] z_nu_N;
  real<lower=0> sigma_GP_N;
  real         eta_N;
  real<lower=0> psi_N;
  real<lower=0> rho_N;
  vector[K_W]  gamma_N;
  vector[K_B]  beta_N;
  real<lower=0> sigma_N;

  // ── P ────────────────────────────────────────────────────────────────────────
  vector[J] z_sp_P;
  vector[J] z_nu_P;
  real<lower=0> sigma_GP_P;
  real         eta_P;
  real<lower=0> psi_P;
  real<lower=0> rho_P;
  vector[K_W]  gamma_P;
  vector[K_B]  beta_P;
  real<lower=0> sigma_P;

}

transformed parameters {

  // ── Fattori spaziali standardizzati (ampiezza = 1) ────────────────────────────
  // f_sp_r[j] ~ GP(0, C_Matérn(rho_r))
  // Il contributo al predittore lineare è: f_sp_r[j] * (sigma_GP_r + eta_r * logB_i)

  vector[J] f_sp_SOC;
  {
    matrix[J, J] C = matern32_corr(dist_mat, rho_SOC);
    for (j in 1:J) C[j, j] += 1e-8;
    f_sp_SOC = cholesky_decompose(C) * z_sp_SOC;
  }

  vector[J] f_sp_N;
  {
    matrix[J, J] C = matern32_corr(dist_mat, rho_N);
    for (j in 1:J) C[j, j] += 1e-8;
    f_sp_N = cholesky_decompose(C) * z_sp_N;
  }

  vector[J] f_sp_P;
  {
    matrix[J, J] C = matern32_corr(dist_mat, rho_P);
    for (j in 1:J) C[j, j] += 1e-8;
    f_sp_P = cholesky_decompose(C) * z_sp_P;
  }

}

model {

  // ── Prior: SOC ───────────────────────────────────────────────────────────────
  z_sp_SOC      ~ std_normal();
  z_nu_SOC      ~ std_normal();
  sigma_GP_SOC  ~ normal(0, 0.5);
  eta_SOC       ~ normal(0, 0.3);    // atteso > 0: campo organico → decadimento più ripido
  psi_SOC       ~ normal(0, 0.5);
  rho_SOC       ~ inv_gamma(3.0, 10.0);
  gamma_SOC     ~ normal(0, 1);
  beta_SOC      ~ normal(0, 1);
  sigma_SOC     ~ normal(0, 0.5);

  // ── Prior: N ─────────────────────────────────────────────────────────────────
  z_sp_N        ~ std_normal();
  z_nu_N        ~ std_normal();
  sigma_GP_N    ~ normal(0, 0.5);
  eta_N         ~ normal(0, 0.3);    // atteso piccolo: struttura spaziale N debole
  psi_N         ~ normal(0, 0.5);
  rho_N         ~ inv_gamma(3.0, 10.0);
  gamma_N       ~ normal(0, 1);
  beta_N        ~ normal(0, 1);
  sigma_N       ~ normal(0, 0.5);

  // ── Prior: P ─────────────────────────────────────────────────────────────────
  z_sp_P        ~ std_normal();
  z_nu_P        ~ std_normal();
  sigma_GP_P    ~ normal(0, 0.5);
  eta_P         ~ normal(0, 0.3);    // ignoto: mai stimato
  psi_P         ~ normal(0, 0.5);
  rho_P         ~ inv_gamma(3.0, 10.0);
  gamma_P       ~ normal(0, 1);
  beta_P        ~ normal(0, 1);
  sigma_P       ~ normal(0, 0.5);

  // ── Likelihood ───────────────────────────────────────────────────────────────
  {
    // Effetti within fissi pre-computati
    vector[N] gSOC = X_W * gamma_SOC;
    vector[N] gN   = X_W * gamma_N;
    vector[N] gP   = X_W * gamma_P;

    // Effetti between fissi pre-computati
    vector[J] bSOC = X_B * beta_SOC;
    vector[J] bN   = X_B * beta_N;
    vector[J] bP   = X_B * beta_P;

    for (n in 1:N) {
      int j = field_id[n];
      real logB_n = X_W[n, 1];   // logBottom standardizzato (colonna 1 di X_W)

      // mu_r = f_sp[j] * (sigma_GP + eta * logB)   ← GP modulato da profondità
      //      + psi * z_nu[j]                        ← random intercept i.i.d.
      //      + gamma * X_W                          ← effetti within fissi
      //      + beta * X_B                           ← effetti between fissi
      real mu_SOC = f_sp_SOC[j] * (sigma_GP_SOC + eta_SOC * logB_n)
                  + psi_SOC * z_nu_SOC[j]
                  + gSOC[n] + bSOC[j];

      real mu_N   = f_sp_N[j] * (sigma_GP_N + eta_N * logB_n)
                  + psi_N * z_nu_N[j]
                  + gN[n] + bN[j];

      real mu_P   = f_sp_P[j] * (sigma_GP_P + eta_P * logB_n)
                  + psi_P * z_nu_P[j]
                  + gP[n] + bP[j];

      logSOC[n] ~ normal(mu_SOC, sigma_SOC);
      logN[n]   ~ normal(mu_N,   sigma_N);
      logP[n]   ~ normal(mu_P,   sigma_P);
    }
  }

}

generated quantities {

  // ── b_r derivato (con guardia numerica per sigma_GP → 0) ─────────────────────
  // b_r = eta_r / sigma_GP_r:  rapporto pendenza/intercetta del GP con la profondità
  real b_SOC = (sigma_GP_SOC > 0.02) ? eta_SOC / sigma_GP_SOC : 0.0;
  real b_N   = (sigma_GP_N   > 0.02) ? eta_N   / sigma_GP_N   : 0.0;
  real b_P   = (sigma_GP_P   > 0.02) ? eta_P   / sigma_GP_P   : 0.0;

  // ── Factor scores in unità originali ─────────────────────────────────────────
  vector[J] u_SOC_out = sigma_GP_SOC * f_sp_SOC;   // random intercept SOC at logB=0
  vector[J] u_N_out   = sigma_GP_N   * f_sp_N;
  vector[J] u_P_out   = sigma_GP_P   * f_sp_P;

  // ── Log-likelihood per LOO-CV ────────────────────────────────────────────────
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
      real mu_SOC = f_sp_SOC[j] * (sigma_GP_SOC + eta_SOC * logB_n)
                  + psi_SOC * z_nu_SOC[j] + gSOC[n] + bSOC[j];
      real mu_N   = f_sp_N[j]   * (sigma_GP_N   + eta_N   * logB_n)
                  + psi_N   * z_nu_N[j]   + gN[n]   + bN[j];
      real mu_P   = f_sp_P[j]   * (sigma_GP_P   + eta_P   * logB_n)
                  + psi_P   * z_nu_P[j]   + gP[n]   + bP[j];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC, sigma_SOC)
                 + normal_lpdf(logN[n]   | mu_N,   sigma_N)
                 + normal_lpdf(logP[n]   | mu_P,   sigma_P);
    }
  }

}
