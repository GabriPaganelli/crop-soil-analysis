// =============================================================================
// m4rr_gp2.stan
//
// Evoluzione di m4rr_gp.stan con tre modifiche motivate dai risultati del
// primo run bayesiano (see scripts/5. bayesian_gp2.R per i dettagli):
//
// MODIFICA 1 — Loading between libero (lambda_N_B) [SOPRASSEDUTA DA MODIFICA 4]
// ------------------------------------------------------------------------------
// In m4rr_gp.stan il loading di logN sul fattore organico era lo stesso a
// entrambi i livelli (lambda_N fisso = 0.636, MLE within da M4rr).
// Il run precedente ha mostrato sigma_B_N = 2.03 >> sigma_B_SOC = 0.16:
// il fattore organico between spiega quasi tutta la variazione between di
// logSOC ma quasi nulla di logN. Questo indica che il vincolo di invarianza
// del loading tra livelli (measurement invariance) non regge per logN.
// Tentativo: aggiungere lambda_N_B come parametro libero per il livello
// between. Il run con questa modifica ha però restituito gli stessi
// Rhat=1.53 e ESS=7 del run con lambda_N libero (vedi MODIFICA 4).
// Nella versione attuale lambda_N_B non è campionato: alpha_N[j] è i.i.d.
//
// MODIFICA 2 — GP per P rimosso, sostituito da random intercept
// --------------------------------------------------------------
// In m4rr_gp.stan il fattore P between seguiva un GP Matern 3/2 con
// parametri sigma_GP_P e rho_P. Il run precedente ha mostrato
// rho_P (mediana) = 18 km > distanza massima tra field (20.4 km):
// il GP degenerava verso una covarianza costante (intercetta globale).
// Tenere rho_P e sigma_GP_P era quindi un costo senza beneficio.
// Soluzione: sostituire il GP con un semplice random intercept i.i.d.:
//   eta_P_B = X_B * gamma_P_B + sigma_P_between * z_P_B
// Si rimuovono sigma_GP_P, rho_P, la matrice K_P e il Cholesky L_P.
// Si aggiunge sigma_P_between ~ normal(0, 1) (half-normal).
// dist_mat rimane in stan_data (ancora usato per il GP organico).
//
// MODIFICA 3 — PH rimosso da X_W
// --------------------------------
// In m4rr_gp.stan PH era incluso in X_W (K_W=5) e il horseshoe lo
// shrinkava a circa zero (gamma_org[2] = -0.032, CI [-0.099, +0.016]),
// confermando la restrizione di M4rr. Nel modello bayesiano non c'è
// il vincolo tecnico dell'LRT che richiedeva di tenerlo nel set delle
// variabili osservate; è più pulito rimuoverlo direttamente da X_W.
// K_W diventa 4: logBottom, Texture1, Texture2, BulkDensity.
//
// MODIFICA 4 — lambda_N_B rimosso; alpha_N[j] con GP spaziale proprio
// ---------------------------------------------------------------------
// Il secondo run (m4rr_gp2 con lambda_N_B libero) ha restituito
// Rhat=1.53 e ESS=7 su tutti i parametri del fattore organico between,
// identici al run iniziale con lambda_N libero. La causa è la stessa:
// il termine lambda_N_B * eta_org_B[j] è un accoppiamento bilineare tra
// due parametri campionati su 40 field → cresta iperbolica → NUTS bloccato.
// Il problema non è specifico di lambda_N: qualsiasi loading libero in
// posizione moltiplicativa con un fattore latente between genera la stessa
// geometria patologica. lambda_N_B rimosso definitivamente.
//
// Invece di lasciare alpha_N[j] puramente i.i.d., il test di Moran's I
// sulle medie per-field ha mostrato autocorrelazione spaziale di logN
// (I = +0.652, p < 0.0001) più forte di logSOC (I = +0.406). La variazione
// between di logN ha struttura spaziale propria, indipendente dal fattore
// organico. Soluzione: GP Matern 3/2 separato per eta_N_B, con parametri
// sigma_GP_N e rho_N propri. Nessun accoppiamento bilineare: eta_N_B non
// è moltiplicato per alcun loading libero.
//   alpha_N[j] = eta_N_B[j] + sigma_B_N * raw[j]
//   eta_N_B ~ MVN(0, sigma_GP_N^2 * K_matern32(D, rho_N) + psi_B_N * I)
//
// MODIFICA 5 — Horseshoe rimosso da gamma_org; prior normal(0,1) uniformi
// -------------------------------------------------------------------------
// Il horseshoe era stato introdotto per shrinkare automaticamente i predittori
// within non informativi (K_W=5 in m4rr_gp.stan). Con K_W=4 e subject-matter
// knowledge che tutti e quattro i predittori hanno meccanismo plausibile
// (logBottom: profilo verticale; Texture: ritenzione organica; BulkDensity:
// compattazione), il horseshoe è sovrabbondante. In particolare, la collinearità
// BulkDensity-Texture2 (r=0.552) portava il horseshoe a shrinkare BulkDensity
// verso zero (-0.107, CI [-0.251, +0.005]) nonostante il frequentista restituisse
// z=-4.01. Il prior normal(0,1) allinea BulkDensity agli altri predittori senza
// penalizzazione arbitraria per la collinearità.
// Rimossi: z_gamma_org, lambda_hs, tau_hs, c2_hs, tau0, slab_scale, slab_df.
// gamma_org ora dichiarato direttamente nei parameters con prior normal(0,1).
//
// NOTA — N_Natural come predittore di P between (gamma_P_B[4])
// -------------------------------------------------------------
// gamma_P_B[4] (N_Natural) = -2.58 (CI [-3.03, -2.12]) vs frequentista -0.309.
// La discrepanza è spuria: N_Natural=0 conta solo 5 field, tutti concentrati
// nell'angolo SW dell'area (distanza media intra-gruppo 2.3 km vs 10 km
// inter-gruppo). La differenza grezza di logP tra i due gruppi è 0.11 log-unità
// (~12%), non 2.58. Il random intercept i.i.d. non assorbe la struttura spaziale
// del cluster SW, e tutta la variazione geografica "SW ha P leggermente più alto"
// confluisce nel coefficiente di N_Natural. Limitazione del disegno: con J=5
// controlli tutti in un angolo non è possibile separare effetto N_Natural da
// effetto geografico. gamma_P_B[4] non è interpretabile causalmente.
//
// Struttura del modello (versione corrente):
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
//     eta_N_B   ~ MVN(0, sigma_GP_N^2  * K_matern32(D, rho_N)   + psi_B_N*I)
//     alpha_SOC[j] = eta_org_B[j] + sigma_B_SOC * raw_SOC[j]
//     alpha_N[j]   = eta_N_B[j]   + sigma_B_N   * raw_N[j]
//     eta_P_B[j]   ~ N(X_B[j] * gamma_P_B, sigma_P_between)   [i.i.d.]
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

  matrix[N, K_W] X_W;       // within: senza PH (modifica 3)
  matrix[J, K_B] X_B;       // between: OnFarm, Irrigate, Fertilised, N_Natural

  matrix[J, J] dist_mat;    // distanze geografiche (km), usate per GP organico e GP N

  // (tau0, slab_scale, slab_df rimossi — modifica 5: horseshoe rimosso)

  // Loading within fisso (MLE da M4rr, vedi modifica 1)
  real<lower=0> lambda_N;

}

parameters {

  // ── (lambda_N_B rimosso — modifica 4: accoppiamento bilineare, vedi header) --

  // ── Coefficienti within per fattore organico (prior normal(0,1) — modifica 5) --
  // [1]=logBottom [2]=Texture1 [3]=Texture2 [4]=BulkDensity
  vector[K_W] gamma_org;

  // ── Coefficienti between per fattore P -------------------------------------
  vector[K_B] gamma_P_B;

  // ── Varianze within --------------------------------------------------------
  real<lower=0> psi_W_org;
  real<lower=0> psi_W_P;
  real<lower=0> theta_W_SOC;
  real<lower=0> theta_W_N;

  // ── Fattori between: organico e N (NCP) ------------------------------------
  vector[J] z_eta_org_B;   // NCP: eta_org_B = L_org * z_eta_org_B
  vector[J] z_eta_N_B;     // NCP: eta_N_B   = L_N   * z_eta_N_B

  // ── Intercette between SOC e N (NCP residui dopo GP) ----------------------
  vector[J] alpha_SOC_raw;
  vector[J] alpha_N_raw;

  // ── Fattore between P (NCP, random intercept i.i.d.) ----------------------
  vector[J] z_P_B;           // NCP: eta_P_B = X_B*gamma_P_B + sigma_P_between*z_P_B

  // ── Varianze between -------------------------------------------------------
  real<lower=0> sigma_B_SOC;
  real<lower=0> sigma_B_N;     // residuo i.i.d. di logN after GP (atteso ~0.1-0.3)
  real<lower=0> psi_B_org;     // nugget GP organico
  real<lower=0> psi_B_N;       // nugget GP N
  real<lower=0> sigma_P_between;

  // ── GP parameters: organico e N (separati, processi indipendenti) ----------
  real<lower=0> sigma_GP_org;
  real<lower=0> rho_org;
  real<lower=0> sigma_GP_N;
  real<lower=0> rho_N;

}

transformed parameters {

  // ── (gamma_org in parameters, prior normal(0,1)) ----------------------------

  // ── Kernel GP organico con nugget ------------------------------------------
  matrix[J, J] L_org;
  {
    matrix[J, J] K_org = matern32_cov(dist_mat, sigma_GP_org, rho_org);
    for (j in 1:J) K_org[j, j] += psi_B_org + 1e-8;
    L_org = cholesky_decompose(K_org);
  }

  // ── Kernel GP N con nugget -------------------------------------------------
  matrix[J, J] L_N;
  {
    matrix[J, J] K_N = matern32_cov(dist_mat, sigma_GP_N, rho_N);
    for (j in 1:J) K_N[j, j] += psi_B_N + 1e-8;
    L_N = cholesky_decompose(K_N);
  }

  // ── Fattori between --------------------------------------------------------
  vector[J] eta_org_B = L_org * z_eta_org_B;
  vector[J] eta_N_B   = L_N   * z_eta_N_B;
  vector[J] eta_P_B   = X_B * gamma_P_B + sigma_P_between * z_P_B;

  // ── Intercette between (NCP residui dopo GP) --------------------------------
  // alpha_SOC[j] = eta_org_B[j] + sigma_B_SOC * raw[j]
  // alpha_N[j]   = eta_N_B[j]   + sigma_B_N   * raw[j]
  vector[J] alpha_SOC = eta_org_B + sigma_B_SOC * alpha_SOC_raw;
  vector[J] alpha_N   = eta_N_B   + sigma_B_N   * alpha_N_raw;

  // ── Covarianza within (SOC, N) — usa lambda_N fisso (within) ---------------
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

  // ── (lambda_N_B rimosso — modifica 4, nessun prior necessario) ---------------

  // ── Prior: gamma_org (modifica 5: normal(0,1), horseshoe rimosso) ----------
  // Weakly informative: SD=1 su covariate standardizzate lascia ampio spazio
  // a effetti di dimensione realistica (-0.5 a -0.1 su scala log).
  gamma_org ~ normal(0, 1);

  // ── Prior: gamma_P_B -------------------------------------------------------
  // NOTA: gamma_P_B[4] (N_Natural) non è interpretabile causalmente —
  // vedere MODIFICA 5 / NOTE nell'header per la spiegazione del confounding.
  gamma_P_B ~ normal(0, 1);

  // ── Prior: varianze within -------------------------------------------------
  psi_W_org   ~ normal(0, 0.5);
  psi_W_P     ~ normal(0, 0.5);
  theta_W_SOC ~ normal(0, 0.5);
  theta_W_N   ~ normal(0, 0.5);

  // ── Prior: NCP fattori between ---------------------------------------------
  z_eta_org_B ~ std_normal();
  z_eta_N_B   ~ std_normal();
  z_P_B       ~ std_normal();

  // ── Prior: NCP intercette between (residui dopo GP) -----------------------
  alpha_SOC_raw ~ std_normal();
  alpha_N_raw   ~ std_normal();

  // ── Prior: varianze between ------------------------------------------------
  sigma_B_SOC     ~ normal(0, 0.5);
  sigma_B_N       ~ normal(0, 0.5);   // residuo i.i.d. dopo GP; prior stretto
  psi_B_org       ~ normal(0, 1.0);
  psi_B_N         ~ normal(0, 1.0);
  sigma_P_between ~ normal(0, 1.0);

  // ── Prior: GP organico e GP N (stessa struttura prior, range separati) ----
  sigma_GP_org ~ normal(0, 1.0);
  sigma_GP_N   ~ normal(0, 1.0);
  rho_org      ~ inv_gamma(3.0, 10.0);  // media 5 km, moda 2.5 km
  rho_N        ~ inv_gamma(3.0, 10.0);  // stessa prior; rho_N stimato indipendentemente

  // ── Likelihood ------------------------------------------------------------
  {
    vector[N] gxW = X_W * gamma_org;
    for (n in 1:N) {
      int j = field_id[n];
      vector[2] mu_SN;
      mu_SN[1] = alpha_SOC[j] + gxW[n];
      mu_SN[2] = alpha_N[j]   + lambda_N * gxW[n];   // within usa lambda_N fisso
      [logSOC[n], logN[n]]' ~ multi_normal_cholesky(mu_SN, L_W_SN);
      logP[n] ~ normal(eta_P_B[j], sqrt(psi_W_P));
    }
  }

}

generated quantities {

  vector[J] eta_org_B_out = eta_org_B;
  vector[J] eta_N_B_out   = eta_N_B;
  vector[J] eta_P_B_out   = eta_P_B;

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
