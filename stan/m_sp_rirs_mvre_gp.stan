// =============================================================================
// m_sp_rirs_mvre_gp.stan  —  M-SP-RIRS + GP spaziale sugli intercetti
//
// Variante robustezza di m_sp_rirs_mvre.stan: sostituisce il prior iid sugli
// intercetti di campo con un Gaussian Process esponenziale separato per
// risposta. Le slope rimangono iid (nessuna struttura MVRE sulle slope,
// dato che rho_slope erano vicine a zero nel modello MVRE).
//
// STRUTTURA INTERCETTI:
//   u_int_r[j] ~ MVN(0, K_r)
//   K_r[i,j]   = alpha_gp_r^2 * exp(-D_ij / ell_r)  +  (i==j) * nu_r^2
//
//   K_r è il kernel esponenziale (Matérn 1/2):
//     - alpha_gp_r: SD spaziale (ampiezza)
//     - ell_r     : range in km (distanza a cui la correlazione scende a 1/e ≈ 0.37)
//     - nu_r      : nugget (variazione non-spaziale, non-strutturata)
//
// STRUTTURA SLOPE:
//   u_slope_r[j] ~ N(0, tau_slope_r^2)   (iid per campo)
//
// MEDIA (identica a MVRE):
//   mu_r[i,j] = alpha_r + u_int_r[j] + u_slope_r[j]*logBottom_i
//             + gamma_r * X_W_i + beta_r * X_B_j
//
// INPUT AGGIUNTIVO RISPETTO A MVRE:
//   Dmat[J,J]: matrice distanze euclidee in km tra campi (simmetrica, diag=0)
//
// NOTA: questo modello NON stima rho_int_SOC_N — la correlazione cross-risposta
// degli intercetti è ora assorbita dalla struttura spaziale condivisa.
// Un'eventuale correlazione residua dopo il GP può essere calcolata ex post
// dai draws di u_int_r.
//
// PARAMETRI TOTALI: 9 (GP iper) + 3J (z_int) + 3J (z_slope) + 3 + 33 = 285 (J=40)
//   GP iper: alpha_gp_r[3], ell_r[3], nu_r[3]
//   NCP:     z_int_r[J×3], z_slope_r[J×3]
//   Slope:   tau_slope_r[3]
//   Fixed:   alpha_r[3], gamma_r[5×3], beta_r[4×3], sigma_r[3]
// =============================================================================

data {

  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K_W;   // 5: logBottom, Texture1, Texture2, BulkDensity, PH
  int<lower=1> K_B;   // 4: OnFarm, Irrigate, Fertilised, N_Natural

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W] X_W;   // colonna 1 = logBottom standardizzato
  matrix[J, K_B] X_B;

  matrix[J, J] Dmat;    // distanze euclidee tra campi in km (diag = 0)

}

parameters {

  // ── GP sugli intercetti: iperparametri ────────────────────────────────────
  real<lower=0> alpha_gp_SOC;   // SD spaziale intercetti SOC
  real<lower=0> ell_SOC;        // range GP SOC (km)
  real<lower=0> nu_SOC;         // nugget SD SOC

  real<lower=0> alpha_gp_N;
  real<lower=0> ell_N;
  real<lower=0> nu_N;

  real<lower=0> alpha_gp_P;
  real<lower=0> ell_P;
  real<lower=0> nu_P;

  // ── GP sugli intercetti: realizzazioni NCP ────────────────────────────────
  vector[J] z_int_SOC;
  vector[J] z_int_N;
  vector[J] z_int_P;

  // ── Slope iid ─────────────────────────────────────────────────────────────
  real<lower=0> tau_slope_SOC;
  real<lower=0> tau_slope_N;
  real<lower=0> tau_slope_P;

  vector[J] z_slope_SOC;
  vector[J] z_slope_N;
  vector[J] z_slope_P;

  // ── Effetti fissi ──────────────────────────────────────────────────────────
  real          alpha_SOC;
  vector[K_W]   gamma_SOC;
  vector[K_B]   beta_SOC;
  real<lower=0> sigma_SOC;

  real          alpha_N;
  vector[K_W]   gamma_N;
  vector[K_B]   beta_N;
  real<lower=0> sigma_N;

  real          alpha_P;
  vector[K_W]   gamma_P;
  vector[K_B]   beta_P;
  real<lower=0> sigma_P;

}

transformed parameters {

  // ── Kernel GP esponenziale con nugget ─────────────────────────────────────
  // K[i,j] = alpha_gp^2 * exp(-D_ij / ell) + (i==j) * nu^2 + jitter
  matrix[J,J] K_SOC = add_diag(
    square(alpha_gp_SOC) * exp(-Dmat / ell_SOC),
    square(nu_SOC) + 1e-8
  );
  matrix[J,J] K_N = add_diag(
    square(alpha_gp_N) * exp(-Dmat / ell_N),
    square(nu_N) + 1e-8
  );
  matrix[J,J] K_P = add_diag(
    square(alpha_gp_P) * exp(-Dmat / ell_P),
    square(nu_P) + 1e-8
  );

  // ── Cholesky dei kernel ───────────────────────────────────────────────────
  matrix[J,J] L_K_SOC = cholesky_decompose(K_SOC);
  matrix[J,J] L_K_N   = cholesky_decompose(K_N);
  matrix[J,J] L_K_P   = cholesky_decompose(K_P);

  // ── Effetti casuali ───────────────────────────────────────────────────────
  vector[J] u_int_SOC   = L_K_SOC * z_int_SOC;   // GP intercetti SOC
  vector[J] u_int_N     = L_K_N   * z_int_N;
  vector[J] u_int_P     = L_K_P   * z_int_P;

  vector[J] u_slope_SOC = tau_slope_SOC * z_slope_SOC;   // slope iid
  vector[J] u_slope_N   = tau_slope_N   * z_slope_N;
  vector[J] u_slope_P   = tau_slope_P   * z_slope_P;

}

model {

  // ── Prior: GP iperparametri ────────────────────────────────────────────────
  // alpha_gp: SD spaziale — stesso ordine di grandezza di tau_6 in MVRE
  alpha_gp_SOC ~ normal(0, 0.5);
  alpha_gp_N   ~ normal(0, 0.5);
  alpha_gp_P   ~ normal(0, 0.5);

  // ell: range in km — half-normal(0, 40) consente range da sub-km a ~80km
  ell_SOC ~ normal(0, 40);
  ell_N   ~ normal(0, 40);
  ell_P   ~ normal(0, 40);

  // nu: nugget — variazione non-spaziale residua
  nu_SOC ~ normal(0, 0.5);
  nu_N   ~ normal(0, 0.5);
  nu_P   ~ normal(0, 0.5);

  // ── Prior: realizzazioni GP (NCP) ─────────────────────────────────────────
  z_int_SOC ~ std_normal();
  z_int_N   ~ std_normal();
  z_int_P   ~ std_normal();

  // ── Prior: slope ──────────────────────────────────────────────────────────
  z_slope_SOC ~ std_normal();
  z_slope_N   ~ std_normal();
  z_slope_P   ~ std_normal();

  tau_slope_SOC ~ normal(0, 0.5);
  tau_slope_N   ~ normal(0, 0.5);
  tau_slope_P   ~ normal(0, 0.5);

  // ── Prior: effetti fissi ──────────────────────────────────────────────────
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

  // ── Likelihood ────────────────────────────────────────────────────────────
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

      real mu_SOC = alpha_SOC + u_int_SOC[j] + u_slope_SOC[j]*logB_n + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + u_int_N[j]   + u_slope_N[j]*logB_n   + gN[n]   + bN[j];
      real mu_P   = alpha_P   + u_int_P[j]   + u_slope_P[j]*logB_n   + gP[n]   + bP[j];

      logSOC[n] ~ normal(mu_SOC, sigma_SOC);
      logN[n]   ~ normal(mu_N,   sigma_N);
      logP[n]   ~ normal(mu_P,   sigma_P);
    }
  }

}

generated quantities {

  // ── Log-likelihood per LOO-CV ─────────────────────────────────────────────
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
      real mu_SOC = alpha_SOC + u_int_SOC[j] + u_slope_SOC[j]*logB_n + gSOC[n] + bSOC[j];
      real mu_N   = alpha_N   + u_int_N[j]   + u_slope_N[j]*logB_n   + gN[n]   + bN[j];
      real mu_P   = alpha_P   + u_int_P[j]   + u_slope_P[j]*logB_n   + gP[n]   + bP[j];
      log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC, sigma_SOC)
                 + normal_lpdf(logN[n]   | mu_N,   sigma_N)
                 + normal_lpdf(logP[n]   | mu_P,   sigma_P);
    }
  }

}
