// =============================================================================
// m4rr_gp8.stan
//
// Estensione di m4rr_gp7.stan con random slope proporzionale all'intercetta
// casuale per logSOC.
//
// MODIFICA RISPETTO A gp7:
// In gp7, l'effetto between organico (eta_org_B[j]) entra come intercetta pura:
//   mu_SOC_n = alpha_SOC[j] + gamma_org' x_W_n
//   alpha_SOC[j] = eta_org_B[j] + sigma_B_SOC * alpha_SOC_raw[j]
//
// In gp8, eta_org_B[j] modula anche la pendenza rispetto a logBottom (X_W[n,1]):
//   mu_SOC_n = eta_org_B[j] * (1 + b_slope * X_W[n,1])
//            + alpha_SOC_res[j]
//            + gamma_org' x_W_n
//
// dove:
//   b_slope     = coefficiente fisso (stima mediana da brms LMM, passato come data)
//   alpha_SOC_res[j] = sigma_B_SOC * alpha_SOC_raw[j]  (residuo i.i.d. di SOC)
//
// Interpretazione:
//   Campi con eta_org_B[j] alto (SOC elevato) hanno anche un profilo di
//   decadimento con la profondità più pronunciato, nella misura stimata da b_slope.
//   b_slope è fissato (non stimato) per evitare il problema bilineare
//   eta_org_B[j] * b_slope * X_W[n,1] che causerebbe non-identificabilità.
//   La stima di b_slope proviene da un LMM bayesiano separato (brms, script 10).
//
// Tutto il resto del modello e' identico a gp7.
//
// RISULTATI run gp8 (4 catene x 5000 sampling + 3000 warmup, adapt_delta=0.97):
//   b_slope = 0.467 (fissato da brms LMM, CI 90% [0.21, 0.70])
//   Divergenze: 0 | Rhat > 1.05: 0 | ESS bulk min ~3500 | E-BFMI: 0.68-0.71
//
//   LOO-CV: elpd_loo(gp8) = -509.0 vs elpd_loo(gp7) = -520.3
//   DELTA_ELPD = +11.3 (SE 5.6) — rapporto 2.0 SE: gp8 migliore di gp7.
//
//   Effetto della modifica:
//     theta_W_SOC:  0.302 -> 0.261  (meno varianza within SOC inspiegata)
//     sigma_B_SOC:  0.156 -> 0.100  (meno residuo between SOC i.i.d.)
//     psi_B_org:    0.076 -> 0.133  (nugget GP riassorbe varianza liberata)
//     sigma_GP_org: 0.234 -> 0.182  (ampiezza GP lievemente piu' piccola)
//     rho_org:      5.12  -> 5.21 km (invariato — struttura spaziale stabile)
//   Tutti gli altri parametri stabili tra gp7 e gp8.
//
//   gp8 e' il modello finale del progetto.
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

  matrix[J, J] dist_mat;    // distanze geografiche (km), usate per GP organico

  real<lower=0> lambda_N;   // loading within fisso (MLE M4rr = 0.636)

  real b_slope;             // pendenza casuale proporzionale a eta_org_B
                            // stimata da LMM brms (script 10, sezione 1)
                            // b = rho * sigma_slp / sigma_int (mediana posteriore)

}

parameters {

  // ── Coefficienti within: fattore organico ----------------------------------
  vector[K_W] gamma_org;

  // ── Coefficienti between: fattore P e fattore N ---------------------------
  vector[K_B] gamma_P_B;
  vector[K_B] gamma_N_B;    // [1]=OnFarm [2]=Irrigate [3]=Fertilised [4]=N_Natural

  // ── Varianze within --------------------------------------------------------
  real<lower=0> psi_W_org;
  real<lower=0> psi_W_P;
  real<lower=0> theta_W_SOC;
  real<lower=0> theta_W_N;

  // ── Fattore between organico (NCP) -----------------------------------------
  vector[J] z_eta_org_B;

  // ── Intercetta residua between SOC (NCP) -----------------------------------
  // alpha_SOC_res[j] = sigma_B_SOC * alpha_SOC_raw[j]
  // Cattura la variabilità between SOC non spiegata da eta_org_B.
  vector[J] alpha_SOC_raw;

  // ── Intercetta between N (NCP) ---------------------------------------------
  vector[J] alpha_N_raw;

  // ── Fattore between P (NCP) ------------------------------------------------
  vector[J] z_P_B;

  // ── Varianze between -------------------------------------------------------
  real<lower=0> sigma_B_SOC;
  real<lower=0> sigma_B_N;
  real<lower=0> psi_B_org;
  real<lower=0> sigma_P_between;

  // ── GP organico ------------------------------------------------------------
  real<lower=0> sigma_GP_org;
  real<lower=0> rho_org;

}

transformed parameters {

  // ── Kernel GP organico -----------------------------------------------------
  matrix[J, J] L_org;
  {
    matrix[J, J] K_org = matern32_cov(dist_mat, sigma_GP_org, rho_org);
    for (j in 1:J) K_org[j, j] += psi_B_org + 1e-8;
    L_org = cholesky_decompose(K_org);
  }

  // ── Fattori between --------------------------------------------------------
  vector[J] eta_org_B    = L_org * z_eta_org_B;
  vector[J] eta_P_B      = X_B * gamma_P_B + sigma_P_between * z_P_B;

  // ── Intercetta residua SOC between (i.i.d., dopo GP) ----------------------
  vector[J] alpha_SOC_res = sigma_B_SOC * alpha_SOC_raw;

  // ── Intercetta between N ---------------------------------------------------
  vector[J] alpha_N = X_B * gamma_N_B + sigma_B_N * alpha_N_raw;

  // ── Covarianza within (SOC, N) ---------------------------------------------
  matrix[2, 2] Sigma_W_SN;
  matrix[2, 2] L_W_SN;
  {
    Sigma_W_SN[1, 1] = psi_W_org + theta_W_SOC;
    Sigma_W_SN[1, 2] = lambda_N * psi_W_org;
    Sigma_W_SN[2, 1] = lambda_N * psi_W_org;
    Sigma_W_SN[2, 2] = square(lambda_N) * psi_W_org + theta_W_N;
    L_W_SN = cholesky_decompose(Sigma_W_SN);
  }

}

model {

  // ── Prior: gamma_org -------------------------------------------------------
  gamma_org ~ normal(0, 1);

  // ── Prior: gamma_P_B -------------------------------------------------------
  gamma_P_B ~ normal(0, 1);

  // ── Prior: gamma_N_B -------------------------------------------------------
  gamma_N_B ~ normal(0, 1);

  // ── Prior: varianze within -------------------------------------------------
  psi_W_org   ~ normal(0, 0.5);
  psi_W_P     ~ normal(0, 0.5);
  theta_W_SOC ~ normal(0, 0.5);
  theta_W_N   ~ normal(0, 0.5);

  // ── Prior: NCP fattori between ---------------------------------------------
  z_eta_org_B ~ std_normal();
  z_P_B       ~ std_normal();

  // ── Prior: NCP intercette between -----------------------------------------
  alpha_SOC_raw ~ std_normal();
  alpha_N_raw   ~ std_normal();

  // ── Prior: varianze between ------------------------------------------------
  sigma_B_SOC     ~ normal(0, 0.5);
  sigma_B_N       ~ normal(0, 0.5);
  psi_B_org       ~ normal(0, 1.0);
  sigma_P_between ~ normal(0, 1.0);

  // ── Prior: GP organico (Matern 3/2, invariato) ----------------------------
  sigma_GP_org ~ normal(0, 1.0);
  rho_org      ~ inv_gamma(3.0, 10.0);

  // ── Likelihood ------------------------------------------------------------
  {
    vector[N] gxW = X_W * gamma_org;
    for (n in 1:N) {
      int j = field_id[n];
      vector[2] mu_SN;
      // mu_SOC: eta_org_B modula sia l'intercetta sia la pendenza di logBottom
      // X_W[n, 1] = logBottom (standardizzato/centrato come in gp7)
      mu_SN[1] = eta_org_B[j] * (1.0 + b_slope * X_W[n, 1])
               + alpha_SOC_res[j]
               + gxW[n];
      mu_SN[2] = alpha_N[j] + lambda_N * gxW[n];
      [logSOC[n], logN[n]]' ~ multi_normal_cholesky(mu_SN, L_W_SN);
      logP[n] ~ normal(eta_P_B[j], sqrt(psi_W_P));
    }
  }

}

generated quantities {

  vector[J] eta_org_B_out  = eta_org_B;
  vector[J] eta_P_B_out    = eta_P_B;
  vector[J] alpha_N_out    = alpha_N;    // intercetta N between (trend + residuo)
  vector[J] alpha_SOC_res_out = alpha_SOC_res;  // residuo SOC between

  vector[N] log_lik;
  {
    vector[N] gxW = X_W * gamma_org;
    for (n in 1:N) {
      int j = field_id[n];
      vector[2] mu_SN;
      mu_SN[1] = eta_org_B[j] * (1.0 + b_slope * X_W[n, 1])
               + alpha_SOC_res[j]
               + gxW[n];
      mu_SN[2] = alpha_N[j] + lambda_N * gxW[n];
      log_lik[n] = multi_normal_cholesky_lpdf(
                     [logSOC[n], logN[n]]' | mu_SN, L_W_SN
                   )
                   + normal_lpdf(logP[n] | eta_P_B[j], sqrt(psi_W_P));
    }
  }

}
