// =============================================================================
// m4rr_gp7.stan
//
// Evoluzione di m4rr_gp3.stan. Rispetto a gp3:
// il trend lineare sulle coordinate (beta_N_x, beta_N_y — risultato nullo,
// CI [-0.96, +0.77] e [-0.94, +0.85] rispettivamente) viene sostituito da
// covariati di management between per logN: gamma_N_B * X_B.
//
// Motivazione:
// gp2 (GP Matern per N): rho_N = 20.9 km, CI [8.8, 71.0] — degenere (area 9.5x19.8 km).
// gp3 (trend lineare): beta_N_x e beta_N_y entrambi non significativi — nessun
//   gradiente direzionale. sigma_B_N risale a 2.10 (come nel caso i.i.d. prima del GP).
// Conclusione: la variazione between di logN non ha struttura spaziale semplice
// (ne' Matern ne' gradiente lineare). L'interpretazione piu' plausibile e' che
// la variabilita' sia guidata da fattori di management (OnFarm, Irrigate,
// Fertilised, N_Natural), gia' inclusi come predittori between per logP.
// Si aggiunge quindi gamma_N_B[K_B] con la stessa struttura di gamma_P_B.
// Rimossi: beta_N_x, beta_N_y, x_km_sc, y_km_sc.
//
// RISULTATI run gp7:
// sigma_B_N = 0.94 (gp2/GP: 0.096 [degenere]; gp3/trend: 2.10; gp7: 0.94).
// I covariati spiegano ~80% della varianza between-field di logN (1-(0.94/2.10)^2),
// quasi interamente attraverso N_Natural:
//   gamma_N_B[1] OnFarm:     +0.364, CI [-0.284, +1.00]  — attraversa zero
//   gamma_N_B[2] Irrigate:   -0.228, CI [-0.745, +0.285] — attraversa zero
//   gamma_N_B[3] Fertilised: -0.238, CI [-0.809, +0.339] — attraversa zero
//   gamma_N_B[4] N_Natural:  -2.26,  CI [-2.67, -1.83]   — chiaramente != 0
// Il sigma_B_N residuo = 0.94 e' irriducibile con i dati disponibili.
// Togliere OnFarm/Irrigate/Fertilised non e' conveniente: la prior normal(0,1)
// le shrinka gia' verso zero, e la simmetria con gamma_P_B (stessa struttura
// X_B) e' preferibile a un'asimmetria non motivata. gp7 e' il modello finale.
// GP per SOC stabile: rho_org = 3.85 km (identico a gp2, robusto tra i run).
//
// NOTA — possibile struttura a 2 cluster per logN:
// Il correlogramma di Moran's I (positivo 0-10 km, negativo 10-21 km) e'
// compatibile con due cluster spaziali: SW (J=5, N_Natural=0, dist. intra 2.3 km)
// vs resto. Questa struttura e' quasi completamente catturata da N_Natural come
// covariata (i 5 field SW coincidono con N_Natural=0). Se dopo questo run
// sigma_B_N rimane grande, si valutera' una struttura di covarianza a blocchi.
//
// Struttura del modello:
//
//   WITHIN (N=220 osservazioni, indice n, campo j[n]):
//     (logSOC_n, logN_n) ~ MVN2(mu_SN_n, Sigma_W_SN)
//       mu_SOC = alpha_SOC[j] + gamma_org' x_W_n
//       mu_N   = alpha_N[j]   + lambda_N * gamma_org' x_W_n
//       Sigma_W_SN = [[psi_W_org + theta_W_SOC,  lambda_N * psi_W_org        ],
//                     [lambda_N * psi_W_org,      lambda_N^2*psi_W_org+theta_W_N]]
//     logP_n ~ N(eta_P_B[j], sqrt(psi_W_P))
//
//   BETWEEN (J=40 field, indice j):
//     eta_org_B ~ MVN(0, sigma_GP_org^2 * K_matern32(D, rho_org) + psi_B_org*I)
//     alpha_SOC[j] = eta_org_B[j] + sigma_B_SOC * raw_SOC[j]
//     alpha_N[j]   = X_B[j] * gamma_N_B + sigma_B_N * raw_N[j]   [management + i.i.d.]
//     eta_P_B[j]   = X_B[j] * gamma_P_B + sigma_P_between * z_P_B[j]
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

  // ── Intercette between SOC e N (NCP) --------------------------------------
  vector[J] alpha_SOC_raw;
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
  vector[J] eta_org_B = L_org * z_eta_org_B;
  vector[J] eta_P_B   = X_B * gamma_P_B + sigma_P_between * z_P_B;

  // ── Intercette between -----------------------------------------------------
  vector[J] alpha_SOC = eta_org_B + sigma_B_SOC * alpha_SOC_raw;
  vector[J] alpha_N   = X_B * gamma_N_B + sigma_B_N * alpha_N_raw;

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
  // gamma_P_B[4] (N_Natural): confounding geografico, non interpretabile
  // causalmente (5 field SW, dist. intra 2.3 km, differenza grezza logP = 0.11).
  gamma_P_B ~ normal(0, 1);

  // ── Prior: gamma_N_B -------------------------------------------------------
  // Stessa prior di gamma_P_B: weakly informative, sd=1 su scala log.
  // N_Natural (gamma_N_B[4]): stesso avvertimento di gamma_P_B[4].
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
      mu_SN[1] = alpha_SOC[j] + gxW[n];
      mu_SN[2] = alpha_N[j]   + lambda_N * gxW[n];
      [logSOC[n], logN[n]]' ~ multi_normal_cholesky(mu_SN, L_W_SN);
      logP[n] ~ normal(eta_P_B[j], sqrt(psi_W_P));
    }
  }

}

generated quantities {

  vector[J] eta_org_B_out = eta_org_B;
  vector[J] eta_P_B_out   = eta_P_B;
  vector[J] alpha_N_out   = alpha_N;   // intercetta N between (trend + residuo)

  vector[N] log_lik;
  {
    vector[N] gxW = X_W * gamma_org;
    for (n in 1:N) {
      int j = field_id[n];
      vector[2] mu_SN;
      mu_SN[1] = alpha_SOC[j] + gxW[n];
      mu_SN[2] = alpha_N[j]   + lambda_N * gxW[n];
      log_lik[n] = multi_normal_cholesky_lpdf(
                     [logSOC[n], logN[n]]' | mu_SN, L_W_SN
                   )
                   + normal_lpdf(logP[n] | eta_P_B[j], sqrt(psi_W_P));
    }
  }

}
