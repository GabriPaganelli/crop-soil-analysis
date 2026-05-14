// =============================================================================
// m4rr_gp3.stan
//
// Evoluzione di m4rr_gp2.stan. Unica modifica rispetto a gp2:
// il GP Matern 3/2 per eta_N_B viene sostituito da un trend lineare spaziale
// sulle coordinate dei field (x_km_sc, y_km_sc).
//
// Motivazione:
// Il run di gp2 ha restituito rho_N (mediana) = 20.9 km, CI [8.8, 71.0 km],
// su un'area di 9.5 km (E-W) x 19.8 km (N-S), con max inter-field = 20.3 km.
// Stesso problema di degenerazione gia' visto per rho_P in gp1 (run1: rho_P=18 km).
// Il correlogramma di Moran's I su logN mostra autocorrelazione positiva a
// distanze 0-10 km e negativa a 10-21 km: firma di un gradiente spaziale globale,
// non di clustering locale. Il kernel Matern 3/2 produce solo correlazioni
// non-negative e non puo' rappresentare questo pattern.
// Un trend lineare su coordinate standardizzate (2 parametri: beta_N_x, beta_N_y)
// cattura il gradiente con parsimonia: parametri 221 -> ~181 (-43 GP, +2 trend).
// Il GP Matern 3/2 per eta_org_B (SOC) rimane invariato (rho_org = 3.85 km,
// stabile e ben identificato, area coperta 9.5x19.8 km).
//
// RISULTATO run gp3 (motivazione abbandono):
// beta_N_x = -0.106, CI [-0.963, +0.769] — attraversa zero.
// beta_N_y = -0.040, CI [-0.939, +0.853] — attraversa zero.
// sigma_B_N risale a 2.10 (identico al caso i.i.d. prima del GP in gp2).
// Il trend lineare non spiega nulla della variazione between di logN.
// Questo e' coerente con la struttura spaziale osservata nell'area:
// la variazione between di logN e' organizzata per cluster di tipo di suolo
// (terreni di tipo A in una zona, terreni di tipo B in un'altra), non come
// gradiente continuo. Ne' Matern (correlazione locale) ne' un piano inclinato
// (gradiente direzionale) possono catturare cluster discreti.
// Modello successivo (gp7): covariati di management between per logN.
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
//     alpha_N[j]   = beta_N_x * x_km_sc[j] + beta_N_y * y_km_sc[j]
//                    + sigma_B_N * raw_N[j]                [trend lineare + i.i.d.]
//     eta_P_B[j]   ~ N(X_B[j] * gamma_P_B, sigma_P_between)   [i.i.d.]
//
// Rispetto a gp2, rimossi: z_eta_N_B[J], sigma_GP_N, rho_N, psi_B_N, L_N, K_N.
// Aggiunti: beta_N_x, beta_N_y (2 parametri, prior normal(0,1)).
// x_km_sc e y_km_sc passati come dati (coordinate standardizzate, SD=1).
// =============================================================================

functions {

  // Kernel Matern 3/2 da matrice delle distanze gia' calcolata (in km)
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

  int<lower=1> N;           // osservazioni totali (220)
  int<lower=1> J;           // field (40)
  int<lower=1> K_W;         // predittori within (4: logBottom,Texture1,Texture2,BulkDensity)
  int<lower=1> K_B;         // predittori between (4: OnFarm,Irrigate,Fertilised,N_Natural)

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  matrix[N, K_W] X_W;       // within: senza PH
  matrix[J, K_B] X_B;       // between: OnFarm, Irrigate, Fertilised, N_Natural

  matrix[J, J] dist_mat;    // distanze geografiche (km), usate per GP organico

  // Coordinate standardizzate dei field (SD=1, calcolate in R)
  // Usate per il trend lineare spaziale di logN between
  vector[J] x_km_sc;
  vector[J] y_km_sc;

  // Loading within fisso (MLE da M4rr, lambda_N = 0.636)
  real<lower=0> lambda_N;

}

parameters {

  // ── Coefficienti within per fattore organico (prior normal(0,1)) -----------
  // [1]=logBottom [2]=Texture1 [3]=Texture2 [4]=BulkDensity
  vector[K_W] gamma_org;

  // ── Coefficienti between per fattore P -------------------------------------
  vector[K_B] gamma_P_B;

  // ── Trend lineare spaziale per logN between --------------------------------
  // Sostituisce il GP Matern 3/2 di gp2 (rho_N degenere a 20.9 km > area).
  // beta_N_x: gradiente E-W; beta_N_y: gradiente N-S.
  // Interpretazione: aumento atteso di logN per +1 SD di coordinata.
  real beta_N_x;
  real beta_N_y;

  // ── Varianze within --------------------------------------------------------
  real<lower=0> psi_W_org;
  real<lower=0> psi_W_P;
  real<lower=0> theta_W_SOC;
  real<lower=0> theta_W_N;

  // ── Fattore between organico (NCP) -----------------------------------------
  vector[J] z_eta_org_B;   // NCP: eta_org_B = L_org * z_eta_org_B

  // ── Intercette between SOC e N (NCP residui dopo GP / dopo trend) ----------
  vector[J] alpha_SOC_raw;
  vector[J] alpha_N_raw;

  // ── Fattore between P (NCP, random intercept i.i.d.) ----------------------
  vector[J] z_P_B;           // NCP: eta_P_B = X_B*gamma_P_B + sigma_P_between*z_P_B

  // ── Varianze between -------------------------------------------------------
  real<lower=0> sigma_B_SOC;
  real<lower=0> sigma_B_N;     // SD residuo i.i.d. di logN dopo trend lineare
  real<lower=0> psi_B_org;     // nugget GP organico
  real<lower=0> sigma_P_between;

  // ── GP parameters: solo organico -------------------------------------------
  real<lower=0> sigma_GP_org;
  real<lower=0> rho_org;

}

transformed parameters {

  // ── Kernel GP organico con nugget ------------------------------------------
  matrix[J, J] L_org;
  {
    matrix[J, J] K_org = matern32_cov(dist_mat, sigma_GP_org, rho_org);
    for (j in 1:J) K_org[j, j] += psi_B_org + 1e-8;
    L_org = cholesky_decompose(K_org);
  }

  // ── Fattori between --------------------------------------------------------
  vector[J] eta_org_B = L_org * z_eta_org_B;
  vector[J] eta_P_B   = X_B * gamma_P_B + sigma_P_between * z_P_B;

  // ── Intercette between (NCP) -----------------------------------------------
  // alpha_SOC[j] = eta_org_B[j] + sigma_B_SOC * raw_SOC[j]
  // alpha_N[j]   = beta_N_x * x_km_sc[j] + beta_N_y * y_km_sc[j]
  //                + sigma_B_N * raw_N[j]
  vector[J] alpha_SOC = eta_org_B + sigma_B_SOC * alpha_SOC_raw;
  vector[J] alpha_N   = beta_N_x * x_km_sc + beta_N_y * y_km_sc
                        + sigma_B_N * alpha_N_raw;

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

  // ── Prior: gamma_org (normal(0,1) — no horseshoe) -------------------------
  gamma_org ~ normal(0, 1);

  // ── Prior: gamma_P_B -------------------------------------------------------
  // NOTA: gamma_P_B[4] (N_Natural) non e' interpretabile causalmente —
  // 5 field N_Natural=0 tutti nel cluster SW (dist. intra-gruppo 2.3 km),
  // la differenza grezza di logP tra gruppi e' 0.11 log-unita', non ~2.5.
  gamma_P_B ~ normal(0, 1);

  // ── Prior: trend lineare spaziale per logN ---------------------------------
  // Coordinate standardizzate (SD=1): prior normal(0,1) lascia spazio a
  // gradienti di dimensione realistica (~ 1 log-unita' per SD di coordinata).
  beta_N_x ~ normal(0, 1);
  beta_N_y ~ normal(0, 1);

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

  // ── Prior: GP organico (Matern 3/2, invariato da gp2) ---------------------
  sigma_GP_org ~ normal(0, 1.0);
  rho_org      ~ inv_gamma(3.0, 10.0);  // media 5 km, moda 2.5 km

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

  // Trend lineare stimato per logN between (senza il residuo i.i.d.)
  vector[J] trend_N_B = beta_N_x * x_km_sc + beta_N_y * y_km_sc;

  // Log-likelihood per LOO-CV
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
