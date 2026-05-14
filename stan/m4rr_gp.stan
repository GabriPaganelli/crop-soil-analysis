// =============================================================================
// m4rr_gp.stan
//
// MIMIC multilevel 2 fattori — versione bayesiana con:
//   - horseshoe regolarizzato (Piironen & Vehtari 2017) su gamma_W_org (5 pred.)
//   - GP Matern 3/2 sui fattori between, range separati rho_org e rho_P
//   - hard zeros da M4rr: phi^W = phi^B = 0, gamma^W_P = 0, gamma^B_org = 0
//
// NOTA su lambda_N (loading logN sul fattore organico)
// -------------------------------------------------------
// lambda_N è trattato come DATO fisso (non campionato).
// Motivazione: nella versione con lambda_N ~ normal(0.63, 0.25) il termine
// lambda_N * eta_org_B[j] crea un accoppiamento bilineare in 41 dimensioni
// (1 loading × 40 field), che genera una posterior a forma di cresta iperbolica.
// NUTS non riesce a percorrerla: nel run precedente ESS~7 e Rhat~1.53 per tutti
// i parametri del fattore organico, pur con E-BFMI>0.6 e solo 4 divergenze.
// La stima frequentista da M4rr (MLR, lavaan) è lambda_N = 0.636, SE = 0.059,
// molto precisa; fissarla perde pochissima informazione rispetto al guadagno
// computazionale. Il valore è passato da R in stan_data$lambda_N.
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
//     alpha_SOC[j] | eta_org_B[j] ~ N(eta_org_B[j], sigma_B_SOC)
//     alpha_N[j]   | eta_org_B[j] ~ N(lambda_N * eta_org_B[j], sigma_B_N)
//     eta_P_B  ~ MVN(X_B * gamma_P_B, sigma_GP_P^2 * K_matern32(D, rho_P) + psi_B_P*I)
//
// Parametri liberi: ~190 (tra cui 80 random intercepts between + 80 NCP GP)
// =============================================================================

functions {

  // Kernel Matern 3/2 da matrice delle distanze gia' calcolata (in km)
  // K[i,j] = sigma^2 * (1 + sqrt(3)*d/rho) * exp(-sqrt(3)*d/rho)
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

  // Dimensioni
  int<lower=1> N;           // osservazioni totali (220)
  int<lower=1> J;           // field (40)
  int<lower=1> K_W;         // predittori within (5)
  int<lower=1> K_B;         // predittori between (4)

  // Indici
  array[N] int<lower=1, upper=J> field_id;

  // Variabili risposta (log-trasformate)
  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  // Covariate
  matrix[N, K_W] X_W;        // within: logBottom, PH, Texture1, Texture2, BulkDensity (std)
  matrix[J, K_B] X_B;        // between: OnFarm, Irrigate, Fertilised, N_Natural (0/1)

  // Distanze geografiche (km, euclidee su coordinate proiettate)
  matrix[J, J] dist_mat;

  // Horseshoe hyperparametri (calcolati in R)
  real<lower=0> tau0;         // scala del global shrinkage
  real<lower=0> slab_scale;   // larghezza dello slab (default 2)
  real<lower=0> slab_df;      // df dello slab (default 4)

  // Loading fisso (vedi nota in testa al file)
  real<lower=0> lambda_N;     // loading logN su fattore organico (MLE da M4rr)

}

transformed data {
  // Nessuna trasformazione necessaria
  int<lower=0> KW_m1 = K_W;  // usato per documentazione
}

parameters {

  // ── Horseshoe su gamma_W_org (predittori within per fattore organico) ----
  // (lambda_N rimosso dai parameters: ora e' un dato fisso, vedi blocco data)
  vector[K_W] z_gamma_org;            // non-centered: N(0,1)
  vector<lower=0>[K_W] lambda_hs;     // local scales (half-Cauchy)
  real<lower=0> tau_hs;               // global scale
  real<lower=0> c2_hs;                // slab regularization parameter

  // ── Coefficienti between per fattore P -----------------------------------
  vector[K_B] gamma_P_B;

  // ── Varianze within ------------------------------------------------------
  real<lower=0> psi_W_org;    // varianza del fattore organico within (residuo)
  real<lower=0> psi_W_P;      // varianza del fattore P within (= logP^W)
  real<lower=0> theta_W_SOC;  // varianza residua misura logSOC within
  real<lower=0> theta_W_N;    // varianza residua misura logN within

  // ── Fattore between organico (non-centrato) ------------------------------
  vector[J] z_eta_org_B;      // NCP: eta_org_B = L_org * z_eta_org_B

  // ── Intercette between per SOC e N — parametrizzazione non-centrata -----
  // NCP: alpha = media + SD * raw,  raw ~ N(0,1)
  // Elimina il funnel tra (alpha, sigma_B) che causava divergenze e E-BFMI<0.3
  vector[J] alpha_SOC_raw;    // NCP raw per logSOC between
  vector[J] alpha_N_raw;      // NCP raw per logN between

  // ── Fattore between P (non-centrato) ------------------------------------
  vector[J] z_eta_P_B;        // NCP: eta_P_B = mu_P + L_P * z_eta_P_B

  // ── Varianze between ----------------------------------------------------
  real<lower=0> sigma_B_SOC;  // SD residua between per logSOC (oltre fattore)
  real<lower=0> sigma_B_N;    // SD residua between per logN (oltre fattore)
  real<lower=0> psi_B_org;    // nugget fattore organico between
  real<lower=0> psi_B_P;      // nugget fattore P between

  // ── GP parameters -------------------------------------------------------
  real<lower=0> sigma_GP_org; // ampiezza GP organico
  real<lower=0> sigma_GP_P;   // ampiezza GP fosforo
  real<lower=0> rho_org;      // range GP organico (km)
  real<lower=0> rho_P;        // range GP fosforo (km)

}

transformed parameters {

  // ── Horseshoe: gamma_org ------------------------------------------------
  vector[K_W] gamma_org;
  {
    // lambda_tilde: local scale regolarizzato (evita slab troppo largo)
    vector[K_W] lambda2    = square(lambda_hs);
    vector[K_W] lambda_til = sqrt(
      c2_hs * lambda2 ./ (c2_hs + square(tau_hs) * lambda2)
    );
    gamma_org = z_gamma_org .* lambda_til * tau_hs;
  }

  // ── Kernel GP con nugget ------------------------------------------------
  matrix[J, J] K_org;
  matrix[J, J] K_P;
  {
    K_org = matern32_cov(dist_mat, sigma_GP_org, rho_org);
    K_P   = matern32_cov(dist_mat, sigma_GP_P,   rho_P);
    // Nugget sulla diagonale (fattore-level variance + jitter numerico)
    for (j in 1:J) {
      K_org[j, j] += psi_B_org + 1e-8;
      K_P[j, j]   += psi_B_P   + 1e-8;
    }
  }

  // ── Cholesky dei kernel -------------------------------------------------
  matrix[J, J] L_org = cholesky_decompose(K_org);
  matrix[J, J] L_P   = cholesky_decompose(K_P);

  // ── Fattori between (da NCP a scala originale) -------------------------
  vector[J] eta_org_B = L_org * z_eta_org_B;
  vector[J] eta_P_B   = X_B * gamma_P_B + L_P * z_eta_P_B;

  // ── Intercette between: da NCP a scala originale ------------------------
  // alpha_SOC[j] = eta_org_B[j] + sigma_B_SOC * raw[j]
  // alpha_N[j]   = lambda_N * eta_org_B[j] + sigma_B_N * raw[j]
  vector[J] alpha_SOC = eta_org_B + sigma_B_SOC * alpha_SOC_raw;
  vector[J] alpha_N   = lambda_N * eta_org_B + sigma_B_N * alpha_N_raw;

  // ── Matrice di covarianza within (SOC, N) — stessa per tutte le oss. ---
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

  // ── Prior: horseshoe regolarizzato su gamma_org -------------------------
  // (lambda_N non ha prior: e' un dato fisso, vedi blocco data)
  z_gamma_org ~ std_normal();
  lambda_hs   ~ cauchy(0, 1);
  tau_hs      ~ normal(0, tau0);
  c2_hs       ~ inv_gamma(slab_df / 2.0, slab_df / 2.0 * square(slab_scale));

  // ── Prior: gamma_P_B (weakly informative) --------------------------------
  gamma_P_B ~ normal(0, 1);

  // ── Prior: varianze within (half-normal) ---------------------------------
  psi_W_org   ~ normal(0, 0.5);
  psi_W_P     ~ normal(0, 0.5);
  theta_W_SOC ~ normal(0, 0.5);
  theta_W_N   ~ normal(0, 0.5);

  // ── Prior: NCP fattori between ------------------------------------------
  z_eta_org_B ~ std_normal();
  z_eta_P_B   ~ std_normal();

  // ── Prior: NCP raw per intercette between --------------------------------
  // Equivalente matematico di alpha ~ normal(mu, sigma) ma senza funnel
  alpha_SOC_raw ~ std_normal();
  alpha_N_raw   ~ std_normal();

  // ── Prior: varianze between (half-normal) --------------------------------
  sigma_B_SOC ~ normal(0, 0.5);
  sigma_B_N   ~ normal(0, 0.5);
  psi_B_org   ~ normal(0, 1.0);
  psi_B_P     ~ normal(0, 1.0);

  // ── Prior: GP amplitudes e ranges ----------------------------------------
  // Amplitudes: half-normal debolmente informativo
  sigma_GP_org ~ normal(0, 1.0);
  sigma_GP_P   ~ normal(0, 1.0);
  // Range: inv_gamma(3, 10) → media 5 km, moda 2.5 km, copre 0.5-20 km
  rho_org ~ inv_gamma(3.0, 10.0);
  rho_P   ~ inv_gamma(3.0, 10.0);

  // ── Likelihood -----------------------------------------------------------
  {
    // Prepariamo il prodotto gamma_org' * x_W per tutte le N osservazioni
    vector[N] gxW = X_W * gamma_org;

    for (n in 1:N) {
      int j = field_id[n];

      // (logSOC, logN): MVN bivariata con covarianza within
      vector[2] mu_SN;
      mu_SN[1] = alpha_SOC[j] + gxW[n];
      mu_SN[2] = alpha_N[j]   + lambda_N * gxW[n];
      [logSOC[n], logN[n]]' ~ multi_normal_cholesky(mu_SN, L_W_SN);

      // logP: indipendente da (logSOC, logN) within
      logP[n] ~ normal(eta_P_B[j], sqrt(psi_W_P));
    }
  }

}

generated quantities {

  // ── Factor scores between (per spatial plots e diagnostiche) -------------
  vector[J] eta_org_B_out = eta_org_B;
  vector[J] eta_P_B_out   = eta_P_B;

  // ── lambda_N fissato: lo riecoheggamo per tracciabilita' nei draw --------
  real lambda_N_out = lambda_N;

  // ── Log-likelihood per LOO-CV (se necessario) ---------------------------
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
