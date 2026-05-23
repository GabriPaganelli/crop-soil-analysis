# =============================================================================
# 8. gam_comparison.R
#
# Due analisi indipendenti come confronto esterno ai modelli MIMIC (script 3-7):
#
# SEZIONE A — GAM univariati (mgcv)
#   Per ciascuna risposta (logSOC, logN, logP): GAMM con smooth di profondità
#   e random intercept per Field. Domanda: la forma di s(Bottom) è compatibile
#   con la parametrizzazione log(Bottom) usata in gp7, o il decadimento ha
#   una forma diversa?
#
# SEZIONE B — GAM multivariato (brms, opzionale)
#   Tre risposte modellate insieme con correlazione residua libera.
#   Permette di confrontare gli effetti fissi stimati con quelli di gp7.
#   Richiede Stan/brms — commentato di default, decommenta per eseguire.
#
# Note:
#   - Nessuna componente spaziale (collinearità Lat/Long con dummy di gestione)
#   - Texture1/Texture2 entrano come lineari (ILR già calcolato in dati.rds)
#   - PH incluso per N e P (diversamente da M4rr dove era fisso a 0 per SOC)
# =============================================================================


# ── 0. SETUP ─────────────────────────────────────────────────────────────────

library(tidyverse)
library(mgcv)
library(here)

dati_raw <- readRDS(here("data", "dati.rds"))

dati <- dati_raw |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(
    logSOC = log(PercSOC),
    logN   = log(PercTotNitro),
    logP   = log(PercTotPhos),
    Field  = factor(Field)
  ) |>
  mutate(across(c(PH, Texture1, Texture2, BulkDensity), ~ c(scale(.x))))

cat("N =", nrow(dati), "| J =", nlevels(dati$Field), "\n\n")


# ── HELPER: indici di fit GAMM ────────────────────────────────────────────────

estrai_fit_gamm <- function(fit, nome) {
  s <- summary(fit)
  tibble(
    modello  = nome,
    dev_expl = round(s$dev.expl * 100, 1),
    r2_adj   = round(s$r.sq, 3),
    aic      = round(AIC(fit), 1),
    gcv      = round(fit$gcv.ubre, 3)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SEZIONE A — GAM UNIVARIATI
# ══════════════════════════════════════════════════════════════════════════════
#
# Per ogni risposta, tre modelli:
#   GA0: effetti lineari, random intercept Field (equivalente a lme4)
#   GA1: s(Bottom) smooth + random intercept Field  ← principale
#   GA2: s(logBottom) smooth + random intercept Field (confronto parametrizzazione)
#
# s(Bottom, bs="cr", k=5): cubic regression spline, k=5 nodi
#   k=5 è il massimo sensato con 6 profondità per Field.
#   La penalizzazione regolerà i gradi di libertà effettivi (tipicamente 2-3).
#
# s(Field, bs="re"): random intercept per Field (equivalente a (1|Field) in lme4)

risposte <- list(
  logSOC = "logSOC",
  logN   = "logN",
  logP   = "logP"
)

cov_lin <- "Texture1 + Texture2 + BulkDensity + PH +
             OnFarm + Irrigate + Fertilised + N_Natural"

fit_gam <- list()
tab_fit  <- list()

for (resp in names(risposte)) {

  cat("\n", strrep("═", 55), "\n")
  cat(" SEZIONE A —", resp, "\n")
  cat(strrep("═", 55), "\n")

  # GA0: lineare su logBottom (baseline compatibile con M4rr/gp7)
  f_ga0 <- as.formula(paste(
    resp, "~ log(Bottom) +", cov_lin, "+ s(Field, bs='re')"
  ))

  # GA1: smooth su Bottom grezzo
  f_ga1 <- as.formula(paste(
    resp, "~ s(Bottom, bs='cr', k=5) +", cov_lin, "+ s(Field, bs='re')"
  ))

  # GA2: smooth su log(Bottom) — più flessibile di GA0 ma stessa scala
  f_ga2 <- as.formula(paste(
    resp, "~ s(log(Bottom), bs='cr', k=5) +", cov_lin, "+ s(Field, bs='re')"
  ))

  ga0 <- gam(f_ga0, data = dati, method = "REML")
  ga1 <- gam(f_ga1, data = dati, method = "REML")
  ga2 <- gam(f_ga2, data = dati, method = "REML")

  fit_gam[[resp]] <- list(ga0 = ga0, ga1 = ga1, ga2 = ga2)

  cat("\n── Fit ───────────────────────────────────\n")
  tab_r <- bind_rows(
    estrai_fit_gamm(ga0, paste0(resp, " GA0: log(Bottom) lineare")),
    estrai_fit_gamm(ga1, paste0(resp, " GA1: s(Bottom)")),
    estrai_fit_gamm(ga2, paste0(resp, " GA2: s(log(Bottom))"))
  )
  print(as.data.frame(tab_r))
  tab_fit[[resp]] <- tab_r

  # ── Gradi di libertà effettivi dello smooth ─────────────────────────────
  cat("\n── EDF smooth (GA1 e GA2) ────────────────\n")
  cat(sprintf("  GA1 s(Bottom):      EDF = %.2f\n",
              summary(ga1)$s.table["s(Bottom)", "edf"]))
  cat(sprintf("  GA2 s(log(Bottom)): EDF = %.2f\n",
              summary(ga2)$s.table["s(log(Bottom))", "edf"]))

  # ── Effetti fissi (GA1) ─────────────────────────────────────────────────
  cat("\n── Effetti fissi GA1 ────────────────────\n")
  ptab <- summary(ga1)$p.table
  print(round(ptab, 3))
}


# ── A1. PLOT SMOOTH DI PROFONDITÀ ─────────────────────────────────────────────

# Predizioni del solo smooth s(Bottom) da GA1 e GA2, con CI,
# sovrapposte alla regressione lineare di log(Bottom) da GA0.

for (resp in names(risposte)) {

  ga0 <- fit_gam[[resp]]$ga0
  ga1 <- fit_gam[[resp]]$ga1
  ga2 <- fit_gam[[resp]]$ga2

  depth_seq <- seq(18, 82, length.out = 100)

  # Predizioni dei soli termini lisci (esclude random Field)
  # Usiamo predict con type="terms" e excl per escludere il random effect

  new_base <- data.frame(
    Bottom    = depth_seq,
    Texture1  = 0, Texture2 = 0, BulkDensity = 0, PH = 0,
    OnFarm = 0, Irrigate = 0, Fertilised = 0, N_Natural = 0,
    Field  = factor(levels(dati$Field)[1], levels = levels(dati$Field))
  )

  # GA0: lineare log(Bottom)
  p0 <- predict(ga0, newdata = new_base, type = "link",
                exclude = "s(Field)", se.fit = TRUE)

  # GA1: s(Bottom)
  p1 <- predict(ga1, newdata = new_base, type = "link",
                exclude = "s(Field)", se.fit = TRUE)

  # GA2: s(log(Bottom))
  p2 <- predict(ga2, newdata = new_base, type = "link",
                exclude = "s(Field)", se.fit = TRUE)

  df_pred <- bind_rows(
    tibble(Bottom = depth_seq, fit = p0$fit, se = p0$se.fit, modello = "GA0: log(Bottom) lineare"),
    tibble(Bottom = depth_seq, fit = p1$fit, se = p1$se.fit, modello = "GA1: s(Bottom)"),
    tibble(Bottom = depth_seq, fit = p2$fit, se = p2$se.fit, modello = "GA2: s(log(Bottom))")
  ) |>
    mutate(
      lo = fit - 1.96 * se,
      hi = fit + 1.96 * se,
      # centra sullo stesso valore medio per comparabilità visiva
      intercept_adj = fit[Bottom == depth_seq[1]][1]
    )

  # Centratura per confronto visivo
  df_pred <- df_pred |>
    group_by(modello) |>
    mutate(fit = fit - mean(fit)) |>
    ungroup()

  p <- ggplot(df_pred, aes(x = Bottom, y = fit,
                            colour = modello, fill = modello)) +
    geom_ribbon(aes(ymin = fit - 1.96 * se,
                    ymax = fit + 1.96 * se), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
    scale_x_continuous(breaks = c(20, 30, 40, 50, 60, 80)) +
    scale_colour_manual(values = c("gray30", "steelblue", "firebrick")) +
    scale_fill_manual(values   = c("gray30", "steelblue", "firebrick")) +
    labs(
      title    = paste("Profilo di profondità —", resp),
      subtitle = "Centrato sulla media; CI 95%. Confronto forme funzionali",
      x = "Profondità Bottom (cm)", y = "Effetto su scala log (centrato)",
      colour = NULL, fill = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(p)
}


# ── A2. CONFRONTO EFFETTI FISSI: GAM vs MIMIC ────────────────────────────────
#
# Carica fit4rr (lavaan M4rr, modello frequentista di riferimento) se disponibile
# e affianca i coefficienti delle covariate comuni.

fit4rr_path <- here("stan", "fit4rr_freq.rds")

if (file.exists(fit4rr_path) && requireNamespace("lavaan", quietly = TRUE)) {
  library(lavaan)
  fit4rr <- readRDS(fit4rr_path)

  freq_coef <- parameterEstimates(fit4rr, standardized = TRUE) |>
    filter(op == "~") |>
    select(lhs, rhs, est_mimic = std.all, z_mimic = z) |>
    filter(rhs %in% c("logBottom", "Texture1", "Texture2", "BulkDensity",
                       "OnFarm", "Irrigate", "Fertilised", "N_Natural"))

  cat("\n── A2. Effetti fissi standardizzati: GAM (GA1) vs MIMIC (M4rr) ──\n")

  for (resp in names(risposte)) {
    ga1_coef <- summary(fit_gam[[resp]]$ga1)$p.table |>
      as.data.frame() |>
      rownames_to_column("rhs") |>
      filter(rhs != "(Intercept)") |>
      select(rhs, est_gam = Estimate, se_gam = `Std. Error`)

    cat("\n", resp, ":\n")
    print(as.data.frame(ga1_coef))
  }

  cat("\nMIMIC M4rr (effetti sul fattore organico within):\n")
  print(as.data.frame(freq_coef |> filter(lhs == "fertil_org_W")))

} else {
  cat("\nfit4rr_freq.rds non trovato: esegui 3. mimic.R per il confronto.\n")
}


# ── A3. RESIDUI E CORRELAZIONE TRA RISPOSTE ──────────────────────────────────
#
# Correlazione residua tra logSOC, logN, logP dopo aver rimosso l'effetto
# delle covariate e del random intercept Field.
# Se alta => la struttura fattoriale di M4rr è giustificata.
# Se bassa per SOC-P => i due fattori separati sono giustificati.

res_df <- dati |>
  mutate(
    res_SOC = residuals(fit_gam$logSOC$ga1),
    res_N   = residuals(fit_gam$logN$ga1),
    res_P   = residuals(fit_gam$logP$ga1)
  ) |>
  select(Field, Bottom, res_SOC, res_N, res_P)

cat("\n── A3. Correlazione residua tra risposte (dopo GAM) ──\n")
cor_res <- cor(res_df[, c("res_SOC", "res_N", "res_P")], use = "complete.obs")
print(round(cor_res, 3))

cat("\nInterpretazione:\n")
cat("  Corr(SOC, N): se alta => fattore condiviso giustificato\n")
cat("  Corr(SOC, P): se bassa => fattore separato per P giustificato\n")


# ── A4. SD RANDOM INTERCEPT PER FIELD E ICC ───────────────────────────────────

cat("\n── A4. SD random intercept per Field (GA1) ──\n")
# Varianza random intercept: sigma^2_Field = sigma^2_residua / lambda_Field
# dove lambda_Field è il parametro di smoothing associato a s(Field, bs='re')
for (resp in names(risposte)) {
  ga1        <- fit_gam[[resp]]$ga1
  sigma2_e   <- ga1$sig2
  # sp: vettore dei parametri di smoothing nell'ordine dei termini s()
  # Il termine s(Field, bs='re') è l'ultimo smooth → ultimo elemento di sp
  lambda_re  <- ga1$sp[length(ga1$sp)]
  sigma2_re  <- sigma2_e / lambda_re
  icc        <- sigma2_re / (sigma2_re + sigma2_e)
  cat(sprintf("  %-8s: SD(Field) = %.3f | SD(residuo) = %.3f | ICC = %.3f\n",
              resp, sqrt(sigma2_re), sqrt(sigma2_e), icc))
}


# ══════════════════════════════════════════════════════════════════════════════
# SEZIONE B — GAM MULTIVARIATO (brms)
# ══════════════════════════════════════════════════════════════════════════════
#
# Tre risposte modellate congiuntamente con correlazione residua libera.
# Smooth s(Bottom, k=5) per ciascuna risposta.
# Random intercept per Field su ciascuna risposta.
# Prior debolmente informativi (student_t(3, 0, 2.5)).
#
# Decommenta per eseguire (richiede ~15-40 min su 4 catene).

if (FALSE) {

  library(brms)

  bf_SOC <- bf(logSOC ~ s(Bottom, k = 5) + Texture1 + Texture2 +
                  BulkDensity + PH + OnFarm + Irrigate + Fertilised +
                  N_Natural + (1 | Field))
  bf_N   <- bf(logN   ~ s(Bottom, k = 5) + Texture1 + Texture2 +
                  BulkDensity + PH + OnFarm + Irrigate + Fertilised +
                  N_Natural + (1 | Field))
  bf_P   <- bf(logP   ~ s(Bottom, k = 5) + Texture1 + Texture2 +
                  BulkDensity + PH + OnFarm + Irrigate + Fertilised +
                  N_Natural + (1 | Field))

  fit_mv_gam <- brm(
    mvbf(bf_SOC, bf_N, bf_P) + set_rescor(TRUE),
    data    = dati,
    prior   = c(
      prior(student_t(3, 0, 2.5), class = b),
      prior(student_t(3, 0, 2.5), class = Intercept),
      prior(student_t(3, 0, 2.5), class = sd),
      prior(student_t(3, 0, 2.5), class = sigma),
      prior(lkj(2), class = rescor)       # prior sulla matrice di correlazione residua
    ),
    chains  = 4, cores = 4,
    iter    = 3000, warmup = 1000,
    seed    = 2024,
    file    = here("stan", "fit_mv_gam")
  )

  summary(fit_mv_gam)

  # Confronto correlazione residua: GAM multivariato vs M7 (lavaan)
  posterior_summary(fit_mv_gam, variable = "^rescor", regex = TRUE)

  # Confronto effetti fissi: GAM multivariato vs gp7 (bayesiano)
  # I coefficienti di OnFarm, Irrigate, ecc. dovrebbero essere comparabili
  # a gamma_org e gamma_P_B di gp7 (tenendo conto dei loadings)
  posterior_summary(fit_mv_gam,
                    variable = c("b_logSOC_OnFarm", "b_logN_OnFarm", "b_logP_OnFarm"))
}

cat("\n── Fine script 8 ────────────────────────────────────────────────────\n")
