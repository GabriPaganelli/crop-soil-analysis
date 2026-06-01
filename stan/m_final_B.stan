// =============================================================================
// m_final_B.stan  —  Opzione B: modello finale response-specific (da projpred)
//
// Struttura fedele alla selezione projpred (alpha=0.10):
//
//   SOC: mu = alpha_SOC + z_nu_SOC[j]*(psi_SOC + eta_SOC*logBottom)
//                       + gamma_SOC[1]*logBottom + gamma_SOC[2]*Texture2
//        → slope proporzionale + {logBottom, Texture2} fissi
//
//   N:   mu = alpha_N + z_nu_N[j]*psi_N
//                     + gamma_N[1]*logBottom + gamma_N[2]*Texture2 + gamma_N[3]*BulkDensity
//        → SOLO random intercept (no eta_N) + {logBottom, Texture2, BulkDensity} fissi
//
//   P:   mu = alpha_P + gamma_P_logB*logBottom + z_nu_P[j]*(psi_P + eta_P*logBottom)
//        → slope proporzionale + logBottom fisso (necessario: E[z_nu]=0, senza
//          gamma_P_logB non ci sarebbe trend medio di profondità per P)
//
// Note:
//   - N ha psi_N ma non eta_N: il random intercept è la sola componente casuale
//   - P non ha gamma_P: la profondità entra solo tramite la struttura proporzionale
//   - Nessun effetto between-field (management rimosso per tutte le risposte)
//   - Stessi prior di m4rr_v2_ri_slope_mu.stan su tutti i parametri condivisi
//
// Parametri totali: 3J + 27 scalari = 120 + 27 = 147 (con J=40)
//   SOC: alpha, psi, eta, gamma[2], sigma, z_nu[J]  = 5 + 2 + J
//   N:   alpha, psi,      gamma[3], sigma, z_nu[J]  = 5 + 3 + J  (no eta)
//   P:   alpha, psi, eta,           sigma, z_nu[J]  = 4 + J       (no gamma)
// =============================================================================

data {

  int<lower=1> N;
  int<lower=1> J;

  array[N] int<lower=1, upper=J> field_id;

  vector[N] logSOC;
  vector[N] logN;
  vector[N] logP;

  // Covariate (tutte scalate come in script 03)
  vector[N] logBottom_s;    // logBottom standardizzato — usato in SOC, N, P
  vector[N] Texture2_s;     // Texture2 (ILR2) standardizzata — usato in SOC, N
  vector[N] BulkDensity_s;  // BulkDensity standardizzata — usato solo in N

}

parameters {

  // ── SOC ──────────────────────────────────────────────────────────────────────
  real          alpha_SOC;
  vector[J]     z_nu_SOC;
  real<lower=0> psi_SOC;
  real          eta_SOC;
  vector[2]     gamma_SOC;   // [1]=logBottom, [2]=Texture2
  real<lower=0> sigma_SOC;

  // ── N (solo random intercept) ─────────────────────────────────────────────
  real          alpha_N;
  vector[J]     z_nu_N;
  real<lower=0> psi_N;       // SD random intercept per N (no eta_N)
  vector[3]     gamma_N;     // [1]=logBottom, [2]=Texture2, [3]=BulkDensity
  real<lower=0> sigma_N;

  // ── P (slope proporzionale + logBottom fisso) ────────────────────────────
  // logBottom fisso necessario: projpred lo seleziona come termine 2 per P.
  // Senza di esso, il modello non avrebbe nessun trend medio di profondità
  // (E[z_nu_P[j]] = 0 per costruzione NCP).
  real          alpha_P;
  vector[J]     z_nu_P;
  real<lower=0> psi_P;
  real          eta_P;
  real          gamma_P_logB;   // effetto fisso di logBottom per P
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
  psi_N   ~ normal(0, 0.5);   // stesso prior di psi in M-SP
  gamma_N ~ normal(0, 1);
  sigma_N ~ normal(0, 0.5);

  // ── Prior: P ─────────────────────────────────────────────────────────────────
  alpha_P      ~ normal(0, 2);
  z_nu_P       ~ std_normal();
  psi_P        ~ normal(0, 0.5);
  eta_P        ~ normal(0, 0.3);
  gamma_P_logB ~ normal(0, 1);
  sigma_P      ~ normal(0, 0.5);

  // ── Likelihood ───────────────────────────────────────────────────────────────
  for (n in 1:N) {
    int j = field_id[n];

    real mu_SOC = alpha_SOC
                + z_nu_SOC[j] * (psi_SOC + eta_SOC * logBottom_s[n])
                + gamma_SOC[1] * logBottom_s[n]
                + gamma_SOC[2] * Texture2_s[n];

    real mu_N = alpha_N
              + z_nu_N[j] * psi_N
              + gamma_N[1] * logBottom_s[n]
              + gamma_N[2] * Texture2_s[n]
              + gamma_N[3] * BulkDensity_s[n];

    real mu_P = alpha_P
              + gamma_P_logB * logBottom_s[n]
              + z_nu_P[j] * (psi_P + eta_P * logBottom_s[n]);

    logSOC[n] ~ normal(mu_SOC, sigma_SOC);
    logN[n]   ~ normal(mu_N,   sigma_N);
    logP[n]   ~ normal(mu_P,   sigma_P);
  }

}

generated quantities {

  // b_r = eta_r / psi_r (solo dove eta esiste)
  real b_SOC = (psi_SOC > 0.02) ? eta_SOC / psi_SOC : 0.0;
  real b_P   = (psi_P   > 0.02) ? eta_P   / psi_P   : 0.0;
  // b_N non esiste in questo modello (eta_N assente)

  vector[N] log_lik;
  for (n in 1:N) {
    int j = field_id[n];

    real mu_SOC = alpha_SOC
                + z_nu_SOC[j] * (psi_SOC + eta_SOC * logBottom_s[n])
                + gamma_SOC[1] * logBottom_s[n]
                + gamma_SOC[2] * Texture2_s[n];

    real mu_N = alpha_N
              + z_nu_N[j] * psi_N
              + gamma_N[1] * logBottom_s[n]
              + gamma_N[2] * Texture2_s[n]
              + gamma_N[3] * BulkDensity_s[n];

    real mu_P = alpha_P
              + gamma_P_logB * logBottom_s[n]
              + z_nu_P[j] * (psi_P + eta_P * logBottom_s[n]);

    log_lik[n] = normal_lpdf(logSOC[n] | mu_SOC, sigma_SOC)
               + normal_lpdf(logN[n]   | mu_N,   sigma_N)
               + normal_lpdf(logP[n]   | mu_P,   sigma_P);
  }

}
