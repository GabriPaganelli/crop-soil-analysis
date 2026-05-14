# ── 0. SETUP ────────────────────────────────────────────────────────────────────

library(tidyverse)
library(lavaan)
library(galamm)
library(here)

dati <- readRDS(here("data", "dati.rds"))


# ── 1. PREPARAZIONE DEI DATI ─────────────────────────────────────────────────────

# Converti dummy da factor a intero 0/1
dati <- dati |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x))))

# Log-trasformazione delle risposte e di Bottom
dati <- dati |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom)
  )

# Standardizzazione covariate continue (media 0, sd 1)
# Le dummy rimangono 0/1; i gamma saranno confrontabili tra loro
dati <- dati |>
  mutate(across(c(logBottom, PH, Texture1, Texture2, BulkDensity),
                ~ c(scale(.x))))

# Dataset completo e senza outlier di PercTotPhos
dati_full   <- dati
dati_no_out <- dati |> filter(PercTotPhos <= 0.5)

cat("Outlier PercTotPhos rimossi:", nrow(dati_full) - nrow(dati_no_out), "\n")

# Shorthand per le covariate (usato in tutti i modelli lavaan)
covar_all <- "logBottom + PH + Texture1 + Texture2 + BulkDensity +
              OnFarm + Irrigate + Fertilised + N_Natural"
covar_W   <- "logBottom + PH + Texture1 + Texture2 + BulkDensity"
covar_B   <- "OnFarm + Irrigate + Fertilised + N_Natural"


# ── HELPER: estrazione indici di fit (lavaan) ────────────────────────────────────

estrai_fit_lavaan <- function(fit, nome) {
  fm <- fitMeasures(fit, c("chisq", "df", "pvalue",
                           "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
                           "cfi", "srmr", "aic", "bic"))
  tibble(
    modello  = nome,
    chisq    = round(fm["chisq"], 2),
    df       = fm["df"],
    p_chisq  = round(fm["pvalue"], 3),
    rmsea    = round(fm["rmsea"], 3),
    rmsea_lo = round(fm["rmsea.ci.lower"], 3),
    rmsea_hi = round(fm["rmsea.ci.upper"], 3),
    cfi      = round(fm["cfi"], 3),
    srmr     = round(fm["srmr"], 3),
    aic      = round(fm["aic"], 1),
    bic      = round(fm["bic"], 1)
  )
}

estrai_fit_galamm <- function(fit, nome) {
  ll  <- logLik(fit)
  k   <- attr(ll, "df")
  n   <- nobs(fit)
  aic <- AIC(fit)
  bic <- BIC(fit)
  tibble(
    modello = nome,
    loglik  = round(as.numeric(ll), 2),
    n_par   = k,
    aic     = round(aic, 1),
    bic     = round(bic, 1)
  )
}


# ── MODELLO 1: MIMIC flat, 1 fattore ────────────────────────────────────────────
#
# Specifica statistica:
#   Misura:     y_ji = lambda_j * eta_i + eps_ji,  eps_ji ~ N(0, theta_j)
#               j in {SOC, N, P};  lambda_SOC = 1 (vincolo di scala)
#   Struttura:  eta_i = gamma' x_i + zeta_i,  zeta_i ~ N(0, psi)
#               x_i = (logBottom, PH, Texture1, Texture2, BulkDensity,
#                      OnFarm, Irrigate, Fertilised, N_Natural)
#   Stimatore:  MLR (robusto a non-normalità e dati cluster)
#   Assunzione: tutte le n osservazioni sono indipendenti (no clustering).
#               Questo modello ignora la struttura gerarchica — serve da prototipo.

mod1 <- paste0('
  fertilita =~ logSOC + logN + logP
  fertilita ~ ', covar_all)

fit1 <- sem(mod1, data = dati_full, estimator = "MLR")
summary(fit1, fit.measures = TRUE, standardized = TRUE)


# ── MODELLO 2: MIMIC flat, 2 fattori ────────────────────────────────────────────
#
# Specifica statistica:
#   Misura:
#     fertil_organica =~ logSOC + logN
#       logSOC = 1 * eta1 + eps1   [lambda_SOC = 1 fissato]
#       logN   = lambda_N * eta1 + eps2
#     fertil_P: singolo indicatore con misura perfetta
#       logP   = 1 * eta2 + 0      [lambda_P = 1, theta_P = 0 fissati => eta2 = logP]
#   Struttura:
#     eta1 = gamma1' x + zeta1
#     eta2 = gamma2' x + zeta2
#     Cov(zeta1, zeta2) = phi12  (libera: cattura correlazione residua fra i fattori)
#   Stimatore: MLR
#   Note: eta2 = logP per costruzione (fattore = indicatore). I due modelli
#         non sono strettamente nested => confronto via AIC/BIC, non Delta-chi2.

mod2 <- paste0('
  fertil_organica =~ logSOC + logN
  fertil_P        =~ 1*logP
  logP            ~~ 0*logP

  fertil_organica ~ ', covar_all, '
  fertil_P        ~ ', covar_all, '

  fertil_organica ~~ fertil_P')

fit2 <- sem(mod2, data = dati_full, estimator = "MLR")
summary(fit2, fit.measures = TRUE, standardized = TRUE)


# ── CONFRONTO FLAT: 1 vs 2 fattori ──────────────────────────────────────────────

tab_fit_flat <- bind_rows(
  estrai_fit_lavaan(fit1, "M1: flat 1F"),
  estrai_fit_lavaan(fit2, "M2: flat 2F")
)
print(tab_fit_flat)

# Delta-chi2 solo se i modelli sono nested (qui NON lo sono => solo AIC/BIC)
cat("\nDelta AIC (M1 - M2):", tab_fit_flat$aic[1] - tab_fit_flat$aic[2], "\n")
cat("Delta BIC (M1 - M2):", tab_fit_flat$bic[1] - tab_fit_flat$bic[2], "\n")


# ── MODIFICATION INDICES (su Modello 1 flat) ────────────────────────────────────
#
# Verifica due percorsi biologicamente motivati non inclusi nel modello:
#   1. logBottom → logP diretto: il fosforo e' poco mobile e si concentra in
#      superficie per ragioni fisiche, non solo attraverso la fertilita' generale.
#   2. logSOC ~~ logN correlazione residua: il rapporto C:N e' quasi fisso (~12);
#      puo' restare covarianza residua specifica alla coppia C-N anche dopo eta.
#
# Soglia: MI > 3.84 (chi2(1) al 5%). Si liberano solo percorsi con giustificazione
# biologica esplicita — non applicare meccanicamente tutti gli MI.

mi1 <- modindices(fit1, sort. = TRUE, maximum.number = 15)
print(mi1)

# Estratti specifici di interesse
mi1 |> filter(
  (lhs == "logBottom" & rhs == "logP")   |
  (lhs == "logP"      & rhs == "logBottom") |
  (lhs == "logSOC"    & rhs == "logN")   |
  (lhs == "logN"      & rhs == "logSOC")
)

# Modello 1 arricchito con gli effetti diretti se MI > 3.84
mod1_mi <- paste0('
  fertilita =~ logSOC + logN + logP
  fertilita ~ ', covar_all, '

  # effetti diretti motivati biologicamente
  logP   ~ logBottom
  logSOC ~~ logN')

fit1_mi <- sem(mod1_mi, data = dati_full, estimator = "MLR")
summary(fit1_mi, fit.measures = TRUE, standardized = TRUE)

# Confronto con e senza MI
tab_fit_mi <- bind_rows(
  estrai_fit_lavaan(fit1,    "M1: flat 1F"),
  estrai_fit_lavaan(fit1_mi, "M1mi: flat 1F + MI")
)
print(tab_fit_mi)


# ── MODELLO 3: MIMIC multilevel, 1 fattore (lavaan) ─────────────────────────────
#
# Specifica statistica:
#   Decomposizione a due livelli: Y_ij = Y^W_ij + Y^B_j
#   (Y^W = deviazione within-Field; Y^B = media between-Field)
#
#   Livello 1 (within-Field, varia con la profondita'):
#     eta^W_ij = gamma^W' x^W_ij + zeta^W_ij
#     logSOC^W  = 1 * eta^W + eps^W_1   [loading fissato per scala]
#     logN^W    = a * eta^W + eps^W_2   [a stimato, vincolato = livello 2]
#     logP^W    = b * eta^W + eps^W_3   [b stimato, vincolato = livello 2]
#     x^W = (logBottom, PH, Texture1, Texture2, BulkDensity)
#
#   Livello 2 (between-Field, varia tra siti):
#     eta^B_j   = gamma^B' x^B_j + zeta^B_j
#     logSOC^B  = 1 * eta^B
#     logN^B    = a * eta^B             [loading = stesso del livello 1]
#     logP^B    = b * eta^B             [loading = stesso del livello 1]
#     x^B = (OnFarm, Irrigate, Fertilised, N_Natural)
#     Var(zeta^B) = sigma^2_u: varianza del random intercept di Field non
#                  spiegata dalle dummy di gestione.
#
#   Vincolo di invarianza dei loadings tra livelli (cross-level invariance):
#   necessario per identificazione con soli 40 Field (6 momenti between vs
#   10 parametri liberi senza vincolo).
#   Stimatore: MLR

mod3 <- '
  level: 1
    eta_W =~ logSOC + a*logN + b*logP
    eta_W ~ logBottom + PH + Texture1 + Texture2 + BulkDensity

  level: 2
    eta_B =~ logSOC + a*logN + b*logP
    eta_B ~ OnFarm + Irrigate + Fertilised + N_Natural
'

fit3 <- sem(mod3, data = dati_full, cluster = "Field", estimator = "MLR")
summary(fit3, fit.measures = TRUE, standardized = TRUE)


# ── MODELLO 4: MIMIC multilevel, 2 fattori (lavaan) ─────────────────────────────
#
# Specifica statistica:
#   Stessa decomposizione a due livelli del Modello 3.
#
#   Livello 1 (within-Field):
#     eta1^W_ij = gamma1^W' x^W_ij + zeta1^W_ij  [fertil_organica within]
#     eta2^W_ij ≡ logP^W_ij                       [fertil_P within = logP^W]
#     logSOC^W  = 1 * eta1^W + eps^W_1
#     logN^W    = a * eta1^W + eps^W_2
#     logP^W    = 1 * eta2^W + 0
#     eta1^W ~~ eta2^W: covarianza residua within tra i due fattori
#
#   Livello 2 (between-Field):
#     eta1^B_j = gamma1^B' x^B_j + zeta1^B_j
#     eta2^B_j ≡ logP^B_j
#     Struttura dei loadings identica al livello 1 (invarianza).
#     Covarianza eta1^B ~~ eta2^B: quanto fertil_organica e P covariano tra Field.
#
#   Interpretazione chiave: se Cov(eta1^B, eta2^B) e' alta per Natural e bassa
#   per Farm fertilizzate, stai separando due regimi del ciclo del fosforo.

mod4 <- '
  level: 1
    fertil_org_W =~ logSOC + a*logN
    fertil_P_W   =~ 1*logP
    logP         ~~ 0*logP
    fertil_org_W ~~ fertil_P_W
    fertil_org_W ~ logBottom + PH + Texture1 + Texture2 + BulkDensity
    fertil_P_W   ~ logBottom + PH + Texture1 + Texture2 + BulkDensity

  level: 2
    fertil_org_B =~ logSOC + a*logN
    fertil_P_B   =~ 1*logP
    logP         ~~ 0*logP
    fertil_org_B ~~ fertil_P_B
    fertil_org_B ~ OnFarm + Irrigate + Fertilised + N_Natural
    fertil_P_B   ~ OnFarm + Irrigate + Fertilised + N_Natural
'

fit4 <- sem(mod4, data = dati_full, cluster = "Field", estimator = "MLR")
summary(fit4, fit.measures = TRUE, standardized = TRUE)


# ── CONFRONTO MULTILEVEL: 1 vs 2 fattori ────────────────────────────────────────

tab_fit_ml <- bind_rows(
  estrai_fit_lavaan(fit3, "M3: ML 1F"),
  estrai_fit_lavaan(fit4, "M4: ML 2F")
)
print(tab_fit_ml)

cat("\nDelta AIC (M3 - M4):", tab_fit_ml$aic[1] - tab_fit_ml$aic[2], "\n")
cat("Delta BIC (M3 - M4):", tab_fit_ml$bic[1] - tab_fit_ml$bic[2], "\n")


# ── MODELLO 4r: MIMIC multilevel 2F — versione ristretta ────────────────────────
#
# Specifica statistica:
#   Identico a M4 tranne per le seguenti restrizioni (motivate empiricamente):
#
#   1. phi^W = 0: covarianza within tra fertil_org_W e fertil_P_W fissata a zero.
#      Motivazione: in M7 la correlazione within SOC-N era 0.068 e SOC-P -0.002;
#      i due fattori sono indipendenti dentro il profilo verticale.
#
#   2. phi^B = 0: covarianza between tra fertil_org_B e fertil_P_B fissata a zero.
#      Motivazione: z = 0.77 in M4, non significativo.
#
#   3. gamma^W(P) = 0: nessun predittore within per fertil_P_W.
#      Motivazione: tutti i coefficienti non significativi in M4 (max |z| = 1.52);
#      il fosforo non mostra profilo verticale sistematico dopo il random Field.
#
#   4. gamma^B(org) = 0: nessun predittore between per fertil_org_B.
#      Motivazione: tutti i coefficienti non significativi in M4 (max |z| = 0.78);
#      la gestione attuale non spiega la variazione between del fattore organico.
#      Interpretazione biologica: la fertilita' organica tra Field e' determinata
#      da storia d'uso e pedologia, non dallo status gestionale corrente.
#
#   Parametri rimossi: 11 (2 covarianze + 5 gamma^W_P + 4 gamma^B_org)
#   Parametri liberi: 32 - 11 = 21
#   df residui: df_M4 + 11 = 21
#   M4r e' nested in M4 -> confronto via LRT Satorra-Bentler valido.

mod4r <- '
  level: 1
    fertil_org_W =~ logSOC + a*logN
    fertil_P_W   =~ 1*logP
    logP         ~~ 0*logP
    fertil_org_W ~~ 0*fertil_P_W
    fertil_org_W ~ logBottom + PH + Texture1 + Texture2 + BulkDensity

  level: 2
    fertil_org_B =~ logSOC + a*logN
    fertil_P_B   =~ 1*logP
    logP         ~~ 0*logP
    fertil_org_B ~~ 0*fertil_P_B
    fertil_P_B   ~ OnFarm + Irrigate + Fertilised + N_Natural
'

fit4r <- sem(mod4r, data = dati_full, cluster = "Field", estimator = "MLR")
summary(fit4r, fit.measures = TRUE, standardized = TRUE)


# ── CONFRONTO M4 vs M4r ─────────────────────────────────────────────────────────

tab_m4_m4r <- bind_rows(
  estrai_fit_lavaan(fit4,  "M4:  ML 2F completo"),
  estrai_fit_lavaan(fit4r, "M4r: ML 2F ristretto")
)
cat("\n── Confronto M4 vs M4r ──────────────────────────────────────────────\n")
print(tab_m4_m4r)
cat("\nDelta AIC (M4r - M4):", round(tab_m4_m4r$aic[2] - tab_m4_m4r$aic[1], 1), "\n")
cat("Delta BIC (M4r - M4):", round(tab_m4_m4r$bic[2] - tab_m4_m4r$bic[1], 1), "\n")

# LRT Satorra-Bentler: M4r nested in M4 (11 restrizioni)
cat("\n── LRT Satorra-Bentler (M4r vs M4) ─────────────────────────────────\n")
print(lavTestLRT(fit4r, fit4))


# ── MODELLO 4rr: MIMIC multilevel 2F — versione ulteriormente ristretta ──────────
#
# Rispetto a M4r, si aggiunge:
#   5. gamma^W_org(PH) = 0: z = -0.83 in M4r, std.all = -0.04.
#      PH non ha effetto sul fattore organico within; fissato esplicitamente (0*PH)
#      per mantenere PH nell'insieme delle variabili osservate e rendere il
#      confronto LRT valido (stesso saturated model di M4r).
#
#   Candidate non rimosse:
#   - gamma^W_org(Texture1): z = -1.53 in M4r (p=0.127): borderline; LRT
#     M4rr_b vs M4rr_a ha p=0.154 e AIC peggiora di +1.7 -> si mantiene.
#
#   Parametri liberi: 21 - 1 = 20  |  df residui: 22
#   M4rr nested in M4r -> LRT Satorra-Bentler vs M4r.

mod4rr <- '
  level: 1
    fertil_org_W =~ logSOC + a*logN
    fertil_P_W   =~ 1*logP
    logP         ~~ 0*logP
    fertil_org_W ~~ 0*fertil_P_W
    fertil_org_W ~ logBottom + 0*PH + Texture1 + Texture2 + BulkDensity

  level: 2
    fertil_org_B =~ logSOC + a*logN
    fertil_P_B   =~ 1*logP
    logP         ~~ 0*logP
    fertil_org_B ~~ 0*fertil_P_B
    fertil_P_B   ~ OnFarm + Irrigate + Fertilised + N_Natural
'

fit4rr <- sem(mod4rr, data = dati_full, cluster = "Field", estimator = "MLR")
summary(fit4rr, fit.measures = TRUE, standardized = TRUE)
saveRDS(fit4rr, here("stan", "fit4rr_freq.rds"))


# ── CONFRONTO M4r vs M4rr ────────────────────────────────────────────────────────

tab_m4r_m4rr <- bind_rows(
  estrai_fit_lavaan(fit4r,  "M4r:  ML 2F ristretto (21 par)"),
  estrai_fit_lavaan(fit4rr, "M4rr: ML 2F min        (20 par)")
)
cat("\n── Confronto M4r vs M4rr ─────────────────────────────────────────────\n")
print(tab_m4r_m4rr)
cat("\nDelta AIC (M4rr - M4r):", round(tab_m4r_m4rr$aic[2] - tab_m4r_m4rr$aic[1], 1), "\n")
cat("Delta BIC (M4rr - M4r):", round(tab_m4r_m4rr$bic[2] - tab_m4r_m4rr$bic[1], 1), "\n")

# LRT Satorra-Bentler: M4rr nested in M4r (1 restrizione: PH=0)
cat("\n── LRT Satorra-Bentler (M4rr vs M4r) ───────────────────────────────\n")
print(lavTestLRT(fit4rr, fit4r))


# ── MODELLO 5: galamm, 1 fattore ────────────────────────────────────────────────
#
# Specifica statistica:
#   Riformulazione come LMM multivariato con struttura fattoriale (Roy & Lin 2002).
#   Formato long: una riga per (osservazione, indicatore).
#
#   Modello marginale per il Field j (integrato su u_j):
#     y_j ~ N(Lambda * gamma' X_j,
#             sigma^2_u * Lambda*Lambda' ⊗ 11' +
#             sigma^2_zeta * Lambda*Lambda' ⊗ I +
#             Theta ⊗ I)
#
#   Lambda = [1, lambda_N, lambda_P]' (lambda_SOC = 1 fissato)
#   u_j ~ N(0, sigma^2_u): random intercept di Field sul fattore latente
#   Theta = diag(theta_SOC, theta_N, theta_P): varianze specifiche degli indicatori
#
#   Nota: in galamm i predittori entrano con effetto DIRETTO su ogni indicatore
#   (non vincolati a essere proporzionali ai loadings come nel MIMIC lavaan).
#   Questo e' un modello piu' generale; se gli effetti stimati risultano
#   approssimativamente proporzionali ai loadings, il vincolo MIMIC e' supportato.

# Formato long
vars_id   <- c("Field", "OnFarm", "Irrigate", "Fertilised", "N_Natural",
               "logBottom", "PH", "Texture1", "Texture2", "BulkDensity")

dati_long <- dati_full |>
  select(all_of(c(vars_id, "logSOC", "logN", "logP"))) |>
  pivot_longer(
    cols      = c(logSOC, logN, logP),
    names_to  = "indicatore",
    values_to = "valore"
  ) |>
  mutate(indicatore = factor(indicatore, levels = c("logSOC", "logN", "logP")))

# Matrice dei loadings: 3 indicatori × 1 fattore
# logSOC: fissato a 1 (referenza), logN e logP: stimati (NA)
lambda_1f <- matrix(c(1, NA, NA), ncol = 1)

fit5 <- galamm(
  formula  = valore ~ 0 + indicatore +
               indicatore:(logBottom + PH + Texture1 + Texture2 +
                           BulkDensity + OnFarm + Irrigate +
                           Fertilised + N_Natural) +
               (0 + loading | Field),
  data     = dati_long,
  factor   = "loading",
  load_var = "indicatore",
  lambda   = lambda_1f
)

summary(fit5)


# ── MODELLO 6: galamm, 2 fattori ────────────────────────────────────────────────
#
# Specifica statistica:
#   Come Modello 5 ma con due fattori latenti:
#     Fattore 1 (fertil_organica): carica su logSOC (lambda=1) e logN (lambda stimato)
#     Fattore 2 (fertil_P):        carica su logP (lambda=1, singolo indicatore)
#   Matrice dei loadings Lambda (3 indicatori × 2 fattori):
#     [ 1   0  ]   <- logSOC: solo su fattore 1
#     [ NA  0  ]   <- logN:   solo su fattore 1
#     [ 0   1  ]   <- logP:   solo su fattore 2
#   u_j = [u1_j, u2_j]': random intercepts per Field su ciascun fattore,
#         con struttura di covarianza libera (cattura correlazione tra fattori).

lambda_2f <- matrix(
  c(1,  0,
    NA, 0,
    0,  1),
  nrow = 3, ncol = 2, byrow = TRUE
)

fit6 <- galamm(
  formula  = valore ~ 0 + indicatore +
               indicatore:(logBottom + PH + Texture1 + Texture2 +
                           BulkDensity + OnFarm + Irrigate +
                           Fertilised + N_Natural) +
               (0 + loading1 | Field) + (0 + loading2 | Field),
  data     = dati_long,
  factor   = c("loading1", "loading2"),
  load_var = "indicatore",
  lambda   = lambda_2f
)

summary(fit6)


# ── TABELLA RIASSUNTIVA FIT ──────────────────────────────────────────────────────

tab_lavaan <- bind_rows(
  estrai_fit_lavaan(fit1,    "M1:  flat 1F"),
  estrai_fit_lavaan(fit1_mi, "M1mi:flat 1F + MI"),
  estrai_fit_lavaan(fit2,    "M2:  flat 2F"),
  estrai_fit_lavaan(fit3,    "M3:  ML 1F"),
  estrai_fit_lavaan(fit4,    "M4:  ML 2F"),
  estrai_fit_lavaan(fit4r,   "M4r: ML 2F ristretto"),
  estrai_fit_lavaan(fit7,    "M7:  ML 3F (multivariata)")
)

tab_galamm <- bind_rows(
  estrai_fit_galamm(fit5, "M5: galamm 1F"),
  estrai_fit_galamm(fit6, "M6: galamm 2F")
)

cat("\n── Fit lavaan ──────────────────────────────────────────────────────\n")
print(tab_lavaan)

cat("\n── Fit galamm (no RMSEA/CFI: non sono modelli SEM) ────────────────\n")
print(tab_galamm)

# Confronto AIC/BIC tra lavaan e galamm (stessa unita' di misura se stesso dataset)
cat("\nAIC confrontabile tra lavaan e galamm (stesso n, stesso likelihood):\n")
bind_rows(
  tab_lavaan |> select(modello, aic, bic),
  tab_galamm |> select(modello, aic, bic)
) |> arrange(aic) |> print()


# ══════════════════════════════════════════════════════════════════════════════
# MODELLO 7: BASELINE MULTIVARIATA (ML, 3 fattori = nessuna struttura latente)
# ══════════════════════════════════════════════════════════════════════════════
#
# Specifica statistica:
#   Ogni indicatore ha il proprio fattore con misura perfetta (lambda=1, theta=0):
#     eta_SOC =~ 1*logSOC;  logSOC ~~ 0*logSOC  =>  eta_SOC = logSOC esattamente
#     eta_N   =~ 1*logN;    logN   ~~ 0*logN
#     eta_P   =~ 1*logP;    logP   ~~ 0*logP
#   Covarianza libera tra tutti e tre i fattori a entrambi i livelli.
#   Struttura equivalente a un LMM multivariato senza compressione latente:
#     ogni risposta e' modellata indipendentemente, con matrice di covarianza libera.
#
# M4 e' una versione vincolata di M7: impone (a) covarianza residua SOC-N = 0
# tramite il fattore condiviso, e (b) proporzionalita' dei coefficienti strutturali
# su SOC e N (effetto_N = lambda_N * effetto_SOC per ogni covariata).
#
# Tre confronti:
#   1. AIC/BIC: M4 vs M7 — la struttura fattoriale paga in parsimonia?
#   2. Correlazione SOC-N nei fattori di M7: se grande, il vincolo di M4 e' rigido
#   3. Proporzionalita' dei gamma in M7: se gamma_N/gamma_SOC ~ lambda_N (M4),
#      il fattore condiviso e' una buona sintesi dei dati

mod7 <- '
  level: 1
    eta_SOC_W =~ 1*logSOC
    eta_N_W   =~ 1*logN
    eta_P_W   =~ 1*logP
    logSOC    ~~ 0*logSOC
    logN      ~~ 0*logN
    logP      ~~ 0*logP
    eta_SOC_W ~~ eta_N_W
    eta_SOC_W ~~ eta_P_W
    eta_N_W   ~~ eta_P_W
    eta_SOC_W ~ logBottom + PH + Texture1 + Texture2 + BulkDensity
    eta_N_W   ~ logBottom + PH + Texture1 + Texture2 + BulkDensity
    eta_P_W   ~ logBottom + PH + Texture1 + Texture2 + BulkDensity

  level: 2
    eta_SOC_B =~ 1*logSOC
    eta_N_B   =~ 1*logN
    eta_P_B   =~ 1*logP
    logSOC    ~~ 0*logSOC
    logN      ~~ 0*logN
    logP      ~~ 0*logP
    eta_SOC_B ~~ eta_N_B
    eta_SOC_B ~~ eta_P_B
    eta_N_B   ~~ eta_P_B
    eta_SOC_B ~ OnFarm + Irrigate + Fertilised + N_Natural
    eta_N_B   ~ OnFarm + Irrigate + Fertilised + N_Natural
    eta_P_B   ~ OnFarm + Irrigate + Fertilised + N_Natural
'

fit7 <- sem(mod7, data = dati_full, cluster = "Field", estimator = "MLR")
summary(fit7, fit.measures = TRUE, standardized = TRUE)


# ── 1. AIC/BIC: M4 vs M7 ────────────────────────────────────────────────────

tab_m4_m7 <- bind_rows(
  estrai_fit_lavaan(fit4, "M4: ML 2F (fattoriale)"),
  estrai_fit_lavaan(fit7, "M7: ML 3F (multivariata)")
)
cat("\n── Confronto M4 vs M7 ──────────────────────────────────────────────\n")
print(tab_m4_m7)
cat("\nDelta AIC (M4 - M7):", round(tab_m4_m7$aic[1] - tab_m4_m7$aic[2], 1), "\n")
cat("Delta BIC (M4 - M7):", round(tab_m4_m7$bic[1] - tab_m4_m7$bic[2], 1), "\n")

# Test chi-quadrato (Satorra-Bentler scaled): M4 e M7 non sono nested in modo
# standard (parameterizzazioni diverse), ma lavaan prova il confronto se i df
# sono compatibili. Se restituisce errore, si usa solo AIC/BIC.
tryCatch(
  { cat("\n── lavTestLRT (Satorra-Bentler) ──\n"); print(lavTestLRT(fit4, fit7)) },
  error = function(e) cat("LRT non applicabile: modelli non nested nella stessa parametrizzazione.\n")
)


# ── 2. Correlazione residua SOC-N in M7 ─────────────────────────────────────
# In M7 la covarianza tra eta_SOC e eta_N e' libera: e' la correlazione tra logSOC
# e logN che RIMANE dopo aver rimosso i predittori e la struttura gerarchica.
# Se e' alta => M4 (che la forza a zero tramite il fattore) e' troppo rigido.

vc_m7 <- lavInspect(fit7, "est")

# Covarianza within (livello 1)
psi_W <- vc_m7$within$psi
cat("\n── Matrice di covarianza fattori within (M7, livello 1) ──\n")
print(round(psi_W, 4))

# Correlazione within
sd_W    <- sqrt(diag(psi_W))
cor_W   <- psi_W / outer(sd_W, sd_W)
cat("\n── Correlazione tra fattori within (M7) ──\n")
print(round(cor_W, 3))

# Covarianza between (livello 2)
psi_B <- vc_m7$Field$psi
cat("\n── Matrice di covarianza fattori between (M7, livello 2) ──\n")
print(round(psi_B, 4))

sd_B  <- sqrt(diag(psi_B))
cor_B <- psi_B / outer(sd_B, sd_B)
cat("\n── Correlazione tra fattori between (M7) ──\n")
print(round(cor_B, 3))


# ── 3. Proporzionalita' dei gamma strutturali M7 vs M4 ──────────────────────
# In M4: effetto di ogni covariata su logN = lambda_N * effetto su logSOC.
# Verifica: il rapporto gamma_N / gamma_SOC in M7 e' circa costante (= lambda_N)?
# Se si => il fattore condiviso e' giustificato.
# Se no => le due risposte reagiscono diversamente alle covariate => 3F preferibile.

lambda_N_m4 <- coef(fit4)["a"]  # loading stimato in M4

pe_m7 <- parameterEstimates(fit7, standardized = FALSE) |>
  filter(op == "~", level == 1) |>
  select(lhs, rhs, est, se) |>
  filter(lhs %in% c("eta_SOC_W", "eta_N_W")) |>
  tidyr::pivot_wider(names_from = lhs,
                     values_from = c(est, se),
                     names_sep   = "_") |>
  mutate(
    rapporto   = round(est_eta_N_W / est_eta_SOC_W, 3),
    lambda_N   = round(lambda_N_m4, 3),
    delta      = round(rapporto - lambda_N_m4, 3)
  )

cat("\n── Proporzionalita' gamma_N / gamma_SOC vs lambda_N di M4 (within) ──\n")
cat("lambda_N stimato in M4:", round(lambda_N_m4, 3), "\n\n")
print(as.data.frame(pe_m7 |> select(rhs, est_eta_SOC_W, est_eta_N_W,
                                     rapporto, lambda_N, delta)))


# ── SENSITIVITY ANALYSIS: con e senza outlier PercTotPhos ────────────────────────
#
# Il modello finale (ML 2 fattori, fit4) viene re-stimato sul dataset senza outlier.
# Si confrontano: loadings, gamma delle covariate, covarianza tra fattori.
# Se i coefficienti cambiano sostanzialmente, l'outlier e' influente e va segnalato.

dati_long_no_out <- dati_no_out |>
  select(all_of(c(vars_id, "logSOC", "logN", "logP"))) |>
  pivot_longer(
    cols      = c(logSOC, logN, logP),
    names_to  = "indicatore",
    values_to = "valore"
  ) |>
  mutate(indicatore = factor(indicatore, levels = c("logSOC", "logN", "logP")))

# Modello 3 (lavaan ML 1F) senza outlier
fit3_no_out <- sem(mod3, data = dati_no_out, cluster = "Field", estimator = "MLR")

# Modello 4 (lavaan ML 2F) senza outlier
fit4_no_out <- sem(mod4, data = dati_no_out, cluster = "Field", estimator = "MLR")

# Modello 5 (galamm 1F) senza outlier
fit5_no_out <- galamm(
  formula  = valore ~ 0 + indicatore +
               indicatore:(logBottom + PH + Texture1 + Texture2 +
                           BulkDensity + OnFarm + Irrigate +
                           Fertilised + N_Natural) +
               (0 + loading | Field),
  data     = dati_long_no_out,
  factor   = "loading",
  load_var = "indicatore",
  lambda   = lambda_1f
)

# Confronto coefficienti lavaan ML 2F: con vs senza outlier
coef_full   <- parameterEstimates(fit4,        standardized = TRUE)
coef_no_out <- parameterEstimates(fit4_no_out, standardized = TRUE)

sensitivity_tab <- coef_full |>
  select(lhs, op, rhs, level, est_full = std.all) |>
  left_join(
    coef_no_out |> select(lhs, op, rhs, level, est_no_out = std.all),
    by = c("lhs", "op", "rhs", "level")
  ) |>
  mutate(delta = round(est_no_out - est_full, 3)) |>
  filter(abs(delta) > 0.05)  # mostra solo differenze > 0.05 in scala standardizzata

cat("\n── Sensitivity: parametri che cambiano > 0.05 std rimuovendo outlier ──\n")
print(sensitivity_tab)


# ══════════════════════════════════════════════════════════════════════════════
# SEZIONE B — FACTOR SCORES (BLUP del Modello 4)
# ══════════════════════════════════════════════════════════════════════════════
#
# Il Modello 4 (ML 2F lavaan) decompone ogni osservazione in:
#   Y_ij = Y^W_ij + Y^B_j
# I "factor scores" sono le stime a posteriori (Empirical Bayes / BLUP) dei
# fattori latenti condizionatamente agli indicatori osservati.
#
# Level 2 (between-Field): eta^B_j — un valore per Field, riflette la fertilità
#   media del sito dopo aver rimosso l'effetto delle covariate between (gestione).
#
# Level 1 (within-Field): eta^W_ij — un valore per osservazione (Field × depth),
#   riflette la deviazione dalla media di Field lungo il profilo verticale.
#
# Nota: lavPredict() usa le equazioni di predizione MAP/EB, non campiona dal
#   posterior — l'incertezza non è automaticamente propagata. Per l'incertezza
#   completa servirebbe un approccio bayesiano (blavaan o brms).

# Coordinate geografiche e Landuse dei Field
coords_field <- crop_raw |>
  group_by(Field) |>
  summarise(Lat  = mean(Lat,  na.rm = TRUE),
            Long = mean(Long, na.rm = TRUE),
            .groups = "drop") |>
  left_join(crop_raw |> distinct(Field, Landuse), by = "Field")


# ── B1. MAPPA GEOGRAFICA DEGLI SCORE BETWEEN-FIELD ───────────────────────────
#
# Estrae lo score di fertilità organica (eta1^B) e di fosforo (eta2^B) per ogni
# Field e li visualizza sulla mappa lat/long.
# Pattern spaziale residuo → variabili geografiche mancanti nel modello (es.
# temperatura, precipitazione, storia d'uso). Vedi SEZIONE C per la correzione.

# lavPredict(level=2) restituisce UNA RIGA PER CLUSTER (40 righe, non 220).
# I cluster sono nell'ordine in cui appaiono per la prima volta in dati_full.
scores_B_raw <- lavPredict(fit4, level = 2)   # 40 × 2
field_ids    <- unique(dati_full$Field)        # stesso ordine visto da lavaan

scores_B <- data.frame(
  Field     = field_ids,
  eta_org_B = scores_B_raw[, "fertil_org_B"],
  eta_P_B   = scores_B_raw[, "fertil_P_B"]
) |>
  left_join(coords_field, by = "Field") |>
  mutate(
    z_org = (eta_org_B - mean(eta_org_B)) / sd(eta_org_B),
    z_P   = (eta_P_B   - mean(eta_P_B))   / sd(eta_P_B)
  )

p_map_org <- ggplot(scores_B, aes(x = Long, y = Lat,
                                   colour = eta_org_B, size = abs(z_org))) +
  geom_point(alpha = 0.85) +
  scale_colour_gradient2(low = "#d73027", mid = "gray90", high = "#1a9850",
                         midpoint = 0, name = "Score\norganico") +
  scale_size_continuous(range = c(2, 7), guide = "none") +
  labs(title    = "Score between-Field: fertilita' organica (eta1^B)",
       subtitle = "Verde = sopra media; Rosso = sotto media",
       x = "Longitudine", y = "Latitudine") +
  theme_minimal()

p_map_P <- ggplot(scores_B, aes(x = Long, y = Lat,
                                 colour = eta_P_B, size = abs(z_P))) +
  geom_point(alpha = 0.85) +
  scale_colour_gradient2(low = "#d73027", mid = "gray90", high = "#1a9850",
                         midpoint = 0, name = "Score\nfosforo") +
  scale_size_continuous(range = c(2, 7), guide = "none") +
  labs(title    = "Score between-Field: fosforo (eta2^B)",
       subtitle = "Verde = sopra media; Rosso = sotto media",
       x = "Longitudine", y = "Latitudine") +
  theme_minimal()

p_map_org + p_map_P


# ── B2. PROFILI VERTICALI DEGLI SCORE WITHIN-FIELD ───────────────────────────
#
# Per ogni Field, eta1^W varia con la profondità.
# Profili decrescenti → fertilità organica più alta in superficie (atteso).
# Profili piatti o positivi → anomalia (Field con accumulo in profondità).

scores_W_raw <- lavPredict(fit4, level = 1)  # matrice n_obs × 2
scores_W <- data.frame(
  Field     = dati_full$Field,
  Bottom    = dati_full$Bottom,
  eta_org_W = scores_W_raw[, "fertil_org_W"],
  eta_P_W   = scores_W_raw[, "fertil_P_W"]
) |>
  left_join(crop_raw |> distinct(Field, Landuse), by = "Field")

# Tutti i Field, colorati per Landuse
ggplot(scores_W, aes(x = Bottom, y = eta_org_W,
                     group = Field, colour = Landuse)) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  geom_point(size = 1.2, alpha = 0.8) +
  scale_x_continuous(breaks = c(20, 30, 40, 50, 60, 80)) +
  labs(title    = "Profili verticali: score within-Field (eta1^W, fertilita' organica)",
       subtitle = "Modello 4 (ML 2F) — ogni linea e' un Field",
       x = "Profondita' (cm)", y = "Score eta^W") +
  theme_minimal()

# Facet per singolo Field (campione di 16 per leggibilita')
set.seed(42)
campi_sample <- sample(unique(scores_W$Field), size = 16)

scores_W |>
  filter(Field %in% campi_sample) |>
  ggplot(aes(x = Bottom, y = eta_org_W, colour = Landuse, group = Field)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ Field, ncol = 4, scales = "free_y") +
  labs(title = "Profili within-Field per campione di 16 Field",
       x = "Profondita' (cm)", y = "Score eta^W (fertilita' org.)") +
  theme_minimal(base_size = 8)


# ── B3. OUTLIER DETECTION TRAMITE SCORE ESTREMI ──────────────────────────────
#
# Field con |z| > 2 su almeno un fattore sono siti "anomali":
#   eta^B > 0: fertilità residua superiore al previsto dalla gestione/tessitura
#   eta^B < 0: fertilità inferiore al previsto
# Questi Field meritano ispezione: errore di misura? gestione non catturata?

outlier_B <- scores_B |>
  filter(abs(z_org) > 2 | abs(z_P) > 2) |>
  arrange(desc(z_org^2 + z_P^2)) |>
  select(Field, Landuse, Lat, Long, eta_org_B, eta_P_B, z_org, z_P)

cat("\n── Field con score between estremi (|z| > 2 su almeno un fattore) ──\n")
print(outlier_B)

# Scatter biplot degli score between
ggplot(scores_B, aes(x = eta_org_B, y = eta_P_B)) +
  geom_point(aes(colour = abs(z_org) > 2 | abs(z_P) > 2,
                 shape  = Landuse), size = 2.5) +
  ggrepel::geom_text_repel(
    data   = filter(scores_B, abs(z_org) > 2 | abs(z_P) > 2),
    aes(label = Field), size = 3, max.overlaps = 20
  ) +
  scale_colour_manual(values = c("FALSE" = "steelblue", "TRUE" = "#d73027"),
                      name   = "|z| > 2") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray60") +
  labs(title    = "BLUP between-Field: eta^B organico vs fosforo",
       subtitle = "Rosso = outlier su almeno un fattore",
       x = "eta^B (fertilita' organica)", y = "eta^B (fosforo)") +
  theme_minimal()


# ── B4. REGRESSIONE DEGLI SCORE SU VARIABILI OMESSE ──────────────────────────
#
# Se i residui dei BLUPs sono correlati con Lat/Long, il modello non ha catturato
# tutta la variazione geografica. Questo è il segnale diagnostico chiave per
# giustificare l'aggiunta di un GP sui random effects (SEZIONE C).
#
# Un R² significativo di eta^B ~ Lat + Long significa che c'è variazione
# spaziale non spiegata dalle covariate di gestione e tessitura.

lm_org_geo <- lm(eta_org_B ~ Lat + Long, data = scores_B)
lm_P_geo   <- lm(eta_P_B   ~ Lat + Long, data = scores_B)

cat("\n── Regressione eta^B organico ~ Lat + Long ──\n")
print(summary(lm_org_geo))
cat("\n── Regressione eta^B fosforo  ~ Lat + Long ──\n")
print(summary(lm_P_geo))

# Visualizzazione: eta^B vs Lat, colorato per Landuse
scores_B |>
  pivot_longer(cols      = c(eta_org_B, eta_P_B),
               names_to  = "fattore",
               values_to = "score") |>
  ggplot(aes(x = Lat, y = score, colour = Landuse)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              colour = "gray30", linewidth = 0.8) +
  facet_wrap(~ fattore, scales = "free_y",
             labeller = labeller(
               fattore = c(eta_org_B = "Fertilita' organica (eta1^B)",
                           eta_P_B   = "Fosforo (eta2^B)"))) +
  labs(title    = "Score between-Field vs latitudine",
       subtitle = "Pendenza significativa => struttura spaziale residua non catturata",
       x = "Latitudine", y = "Score BLUP") +
  theme_minimal()

# Moran's I (test formale di autocorrelazione spaziale) se ape è disponibile
if (requireNamespace("ape", quietly = TRUE)) {
  library(ape)
  # Matrice di distanza geografica tra Field
  coord_mat  <- as.matrix(scores_B[, c("Long", "Lat")])
  dist_geo   <- as.matrix(dist(coord_mat))
  # Matrice di pesi spaziali: inverso della distanza (escluso diagonale)
  w_mat            <- 1 / dist_geo
  diag(w_mat)      <- 0
  w_mat            <- w_mat / rowSums(w_mat)   # normalizzazione row-standard

  mi_org <- Moran.I(scores_B$eta_org_B, w_mat)
  mi_P   <- Moran.I(scores_B$eta_P_B,   w_mat)

  cat("\n── Moran's I (autocorrelazione spaziale dei BLUP) ──\n")
  cat("eta^B organico: I =", round(mi_org$observed, 3),
      "  p =", round(mi_org$p.value, 4), "\n")
  cat("eta^B fosforo:  I =", round(mi_P$observed,   3),
      "  p =", round(mi_P$p.value,   4), "\n")
} else {
  cat("\nInstalla il pacchetto 'ape' per il test di Moran's I.\n")
}
