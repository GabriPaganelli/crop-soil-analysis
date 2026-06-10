# =============================================================================
# 06_lmm_esplorativi.R  —  LMM multivariati frequentisti (esplorativi)
#
# Tre varianti del modello LMM multivariato, usate come motivazione per il
# passaggio all'approccio bayesiano gerarchico (script 07-09).
#
# ── Variante A (Sezione 1) ── modello base con covariabili binarie (0/1)
#    Risposta log-scale; predittori: N_Natural, OnFarm, Irrigate, Fertilised,
#    logBottom_c, PH, BulkDensity, Texture1, Texture2
#
# ── Variante B (Sezione 2) ── usa Landuse (7 livelli) come between-field
#    Al posto dei 4 binari, una singola variabile categorica.
#    Classe di riferimento: Landuse = 4.
#
# ── Variante C (Sezione 3) ── risposta su scala naturale × BulkDensity
#    y_r = Perc_r × BulkDensity  (grammi per unità di volume)
#    Testa se la scelta della scala (log vs naturale) cambia le conclusioni.
#
# PERCHÉ ABBANDONATI A FAVORE DI STAN:
#   1. MCMCglmm non modella esplicitamente η_r (correlazione slope–intercetta
#      per campo): la matrice G 6×6 per int+slope è al limite dell'identificabilità
#      con J=40 campi e 21 parametri liberi.
#   2. Nessun parametro interpretabile analogo a η_r (o ρ_r in M-SP-RIRS)
#      che quantifichi direttamente la dipendenza slope–livello-campo.
#   3. Difficoltà numeriche di convergenza con sommer per la struttura G 6×6.
#   → M-SP-RIRS (script 09) risolve questi problemi con parameterizzazione NCP
#     e prior LKJ su Omega_r.
#
# Dipende da: data/crop_analytic.rds, data/crop_full.rds, scripts/00_utilities.R
# =============================================================================

library(MCMCglmm)
library(sommer)
library(lme4)
library(lmerTest)
library(tidyverse)
library(patchwork)
library(here)
source(here("scripts", "00_utilities.R"))

# MCMCglmm carica MASS, che maschera dplyr::select, dplyr::filter, dplyr::rename.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename


# ──────────────────────────────────────────────────────────────────────────────
# SEZIONE 1: Variante base — binari di gestione
# ──────────────────────────────────────────────────────────────────────────────

dati <- readRDS(here("data", "crop_analytic.rds"))
dati$logBottom_c <- log(dati$Bottom) - mean(log(dati$Bottom))

dati <- dati |>
  mutate(
    logSOC = log(PercSOC),
    logN   = log(PercTotNitro),
    logP   = log(PercTotPhos),
    Field  = factor(Field)
  )

target_vars <- c("logSOC", "logN", "logP")
n_risposte  <- length(target_vars)
n_fields    <- length(unique(dati$Field))
n_obs       <- nrow(dati)

# Feasibility check
params_G_int   <- n_risposte * (n_risposte + 1) / 2
params_G_slope <- (2 * n_risposte) * (2 * n_risposte + 1) / 2

cat("═══ FEASIBILITY CHECK ═══════════════════════════════════════════\n")
cat(sprintf("  G 3×3 (intercept only): %d param,  %.1f gruppi/param → FATTIBILE\n",
            params_G_int, n_fields / params_G_int))
cat(sprintf("  G 6×6 (int + slope):   %d param,  %.1f gruppi/param → AL LIMITE\n",
            params_G_slope, n_fields / params_G_slope))

# Correlazione BLUP empirica
formula_fx <- "~ N_Natural + OnFarm + Irrigate + Fertilised +
                 logBottom_c + PH + BulkDensity + Texture1 + Texture2"

dati_sc <- dati |> mutate(
  across(c(logBottom_c, BulkDensity, PH, Texture1, Texture2), scale)
)

m_uni <- setNames(
  lapply(target_vars, function(y) {
    lmer(as.formula(paste(y, formula_fx, "+ (1 | Field)")),
         data = dati_sc, REML = TRUE)
  }),
  target_vars
)

blup_mat <- sapply(m_uni, function(m) ranef(m)$Field[[1]])
colnames(blup_mat) <- target_vars

cat("\n═══ CORRELAZIONE BLUP TRA-FIELD (variante base) ════════════════\n")
print(round(cor(blup_mat), 3))

# MCMCglmm — random intercept
dati_long <- dati |>
  select(Field, N_Natural, OnFarm, Irrigate, Fertilised,
         logBottom_c, PH, BulkDensity, Texture1, Texture2,
         logSOC, logN, logP) |>
  pivot_longer(cols = all_of(target_vars),
               names_to  = "componente",
               values_to = "y") |>
  mutate(componente = factor(componente, levels = target_vars)) |>
  as.data.frame()

prior_int <- list(
  G = list(G1 = list(V = diag(n_risposte) * 0.1, nu = n_risposte + 2)),
  R = list(       V = diag(n_risposte) * 0.1, nu = n_risposte + 2)
)

formula_mcmc <- y ~ componente - 1 +
  componente:N_Natural + componente:OnFarm + componente:Irrigate + componente:Fertilised +
  componente:logBottom_c +
  componente:PH + componente:BulkDensity + componente:Texture1 + componente:Texture2

cat("\nAvvio MCMCglmm — variante base (G 3×3)...\n")
set.seed(42)
m_A <- MCMCglmm(
  formula_mcmc,
  random  = ~ us(componente):Field,
  rcov    = ~ us(componente):units,
  family  = "gaussian",
  prior   = prior_int,
  data    = dati_long,
  nitt    = 100000, burnin = 20000, thin = 80,
  verbose = FALSE
)
cat("Completato.\n"); summary(m_A)

# Estrazione correlazioni G e R
vcv_G <- m_A$VCV[, grep("Field", colnames(m_A$VCV))]
vcv_R <- m_A$VCV[, grep("units", colnames(m_A$VCV))]
G_mat <- matrix(colMeans(vcv_G), n_risposte, n_risposte,
                dimnames = list(c("SOC","N","P"), c("SOC","N","P")))
R_mat <- matrix(colMeans(vcv_R), n_risposte, n_risposte,
                dimnames = list(c("SOC","N","P"), c("SOC","N","P")))

G_cor <- cov2cor(G_mat)
R_cor <- cov2cor(R_mat)

cat("\n── G (tra-Field):\n"); print(round(G_cor, 3))
cat("\n── R (residua):\n");   print(round(R_cor, 3))

# DIC confronto int vs int+slope
prior_slope <- list(
  G = list(G1 = list(V = diag(2 * n_risposte) * 0.1, nu = 2 * n_risposte + 5)),
  R = list(       V = diag(n_risposte) * 0.1, nu = n_risposte + 2)
)
set.seed(43)
m_A_slope <- MCMCglmm(
  formula_mcmc,
  random  = ~ us(componente + componente:logBottom_c):Field,
  rcov    = ~ us(componente):units,
  family  = "gaussian",
  prior   = prior_slope,
  data    = dati_long,
  nitt    = 150000, burnin = 30000, thin = 120,
  verbose = FALSE
)
cat(sprintf("\nDIC int:          %.2f\n", m_A$DIC))
cat(sprintf("DIC int + slope:  %.2f  (Delta = %.2f)\n",
            m_A_slope$DIC, m_A$DIC - m_A_slope$DIC))
cat("Delta DIC > 10 → variazione slope rilevante; < 0 → int sufficiente\n")

# sommer (verifica frequentista)
dati_s <- dati |>
  mutate(across(c(N_Natural, OnFarm, Irrigate, Fertilised),
                ~ as.numeric(as.character(.x))))
m_sommer <- tryCatch(
  mmer(cbind(logSOC, logN, logP) ~
         N_Natural + OnFarm + Irrigate + Fertilised +
         logBottom_c + PH + BulkDensity + Texture1 + Texture2,
       random  = ~ vsr(Field, Gtc = unsm(3)),
       rcov    = ~ vsr(units, Gtc = unsm(3)),
       data    = as.data.frame(dati_s), verbose = FALSE),
  error = function(e) { cat("sommer:", e$message, "\n"); NULL }
)
if (!is.null(m_sommer)) {
  cat("\n── sommer G:\n"); print(round(cov2cor(m_sommer$sigma$`u:Field`), 3))
}


# ──────────────────────────────────────────────────────────────────────────────
# SEZIONE 2: Variante Landuse — variabile categorica between-field
# ──────────────────────────────────────────────────────────────────────────────

cat("\n\n═══ SEZIONE 2: Variante Landuse (riferimento = 4) ══════════════\n")

crop <- readRDS(here("data", "crop_full.rds"))
dati2 <- dati |>
  mutate(Landuse = relevel(factor(crop$Landuse), ref = "4"))

dati_long2 <- dati2 |>
  select(Field, Landuse, logBottom_c, PH, BulkDensity, Texture1, Texture2,
         logSOC, logN, logP) |>
  pivot_longer(cols = all_of(target_vars),
               names_to  = "componente",
               values_to = "y") |>
  mutate(componente = factor(componente, levels = target_vars)) |>
  as.data.frame()

formula_B <- y ~ componente - 1 +
  componente:Landuse +
  componente:logBottom_c +
  componente:PH + componente:BulkDensity + componente:Texture1 + componente:Texture2

set.seed(44)
m_B <- MCMCglmm(
  formula_B,
  random  = ~ us(componente):Field,
  rcov    = ~ us(componente):units,
  family  = "gaussian",
  prior   = prior_int,
  data    = dati_long2,
  nitt    = 100000, burnin = 20000, thin = 80,
  verbose = FALSE
)
cat("Completato — variante Landuse.\n"); summary(m_B)

landuse_cols <- grep("componente.*:Landuse", colnames(m_B$Sol), value = TRUE)
landuse_levels <- unique(sub("componente.*:Landuse", "", landuse_cols))
cat("Livelli Landuse attivi:", paste(landuse_levels, collapse = ", "),
    "  (Landuse 4 = baseline)\n")


# ──────────────────────────────────────────────────────────────────────────────
# SEZIONE 3: Variante scala naturale × BulkDensity
# ──────────────────────────────────────────────────────────────────────────────

cat("\n\n═══ SEZIONE 3: Variante scala naturale × BulkDensity ════════════\n")
cat("  y_r = Perc_r * BulkDensity  (g/cm³ × %)  [~grammi per volume di suolo]\n\n")

dati3 <- dati |>
  mutate(
    ySOC = PercSOC * BulkDensity,
    yN   = PercTotNitro * BulkDensity,
    yP   = PercTotPhos * BulkDensity
  )

dati_long3 <- dati3 |>
  select(Field, N_Natural, OnFarm, Irrigate, Fertilised,
         logBottom_c, PH, BulkDensity, Texture1, Texture2,
         ySOC, yN, yP) |>
  pivot_longer(cols = c(ySOC, yN, yP),
               names_to  = "componente",
               values_to = "y") |>
  mutate(componente = factor(componente, levels = c("ySOC","yN","yP"))) |>
  as.data.frame()

formula_C <- y ~ componente - 1 +
  componente:N_Natural + componente:OnFarm + componente:Irrigate + componente:Fertilised +
  componente:logBottom_c +
  componente:PH + componente:Texture1 + componente:Texture2

set.seed(45)
m_C <- MCMCglmm(
  formula_C,
  random  = ~ us(componente):Field,
  rcov    = ~ us(componente):units,
  family  = "gaussian",
  prior   = prior_int,
  data    = dati_long3,
  nitt    = 100000, burnin = 20000, thin = 80,
  verbose = FALSE
)
cat("Completato — variante scala naturale.\n"); summary(m_C)

cat("\n── Conclusione (scala naturale vs log):\n")
cat("  Le conclusioni qualitative sulle correlazioni G e R devono essere simili.\n")
cat("  La scala log è preferita: distribuzioni più simmetriche, residui più normali.\n")


# ──────────────────────────────────────────────────────────────────────────────
# MOTIVAZIONE PER IL PASSAGGIO AL MODELLO BAYESIANO STAN
# ──────────────────────────────────────────────────────────────────────────────

cat("\n\n═══ MOTIVAZIONE PER M-SP-RIRS (stan/model_rirs.stan) ════════════\n\n")
cat("  1. PARAMETRO ETA/RHO NON IDENTIFICABILE VIA MCMCglmm:\n")
cat("     η_r (o ρ_r in RIRS) quantifica la correlazione slope-intercetta per campo.\n")
cat("     In MCMCglmm, la matrice G 6×6 ha 21 param con J=40 (1.9 gruppi/param).\n")
cat("     Il prior LKJ(2) di Stan identifica meglio ρ_r con la NCP.\n\n")
cat("  2. DIFFICOLTÀ DI CONVERGENZA:\n")
cat("     sommer con G 6×6 spesso non converge o restituisce stime degenerate.\n")
cat("     M-SP originale (Script 08) già mostra η_SOC=0.21 [0.12,0.30] — robusto.\n\n")
cat("  3. INFERENZA DIRETTA SU SLOPE:\n")
cat("     M-SP-RIRS stima τ_β_r (SD slope tra-campo) e ρ_r (correlazione RI-RS)\n")
cat("     come parametri diretti della posteriore, con CI credibili interpretabili.\n\n")
cat("  → Continuare con script 07 (M-RI) e 09 (M-SP-RIRS).\n")

cat("\n── Fine script 06 ──────────────────────────────────────────────\n")
