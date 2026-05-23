# =============================================================================
# 9. random_slopes.R
#
# Esplorazione random slopes per profondità su logSOC, logN, logP.
# Due domande:
#   A. Il tasso di decadimento verticale varia tra Field? (varianza della slope)
#   B. La correlazione intercetta-slope è negativa? (siti ricchi decadono più rapido)
#
# Due parametrizzazioni della profondità a confronto:
#   M1: slope su logBottom (scala log-log, power law)  → come in gp7
#   M2: slope su scale(Bottom) (scala log-lin, esponenziale) → alternativa
#
# Modelli frequentisti (lme4, ML). Risultati chiave:
#   - logSOC: ΔAIC M1 = -12.3 (p=0.0003), ma singolarità (Corr=1.0) con 5-6 obs/Field.
#     Il segnale è reale; la singolarità è un limite del ML con pochi livelli.
#     Estensione naturale: random slope su logBottom in Stan (potenziale gp8).
#   - logN: ΔAIC M1 = -3.3 (p=0.025), nessuna singolarità (Corr=-0.25). Segnale marginale.
#   - logP: ΔAIC M1 = -0.8 (p=0.09), singolarità. Nessun segnale.
#   - M2 sempre peggiore di M1: parametrizzazione log-lin (esponenziale) esclusa.
#   - Correlazione tra slopes dei tre fattori per Field: SOC-N=0.24, SOC-P=0.07, N-P=-0.21.
#     Non c'è un unico "tasso di decadimento del sito" condiviso tra le risposte.
# =============================================================================

library(tidyverse)
library(lme4)
library(here)

# ── 0. DATI ──────────────────────────────────────────────────────────────────

dati_raw <- readRDS(here("data", "dati.rds"))

dati <- dati_raw |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom),
    Bottom_sc = c(scale(Bottom))      # centrato E scalato (era solo /100: non centrato)
  ) |>
  mutate(across(c(logBottom, PH, Texture1, Texture2, BulkDensity), ~ c(scale(.x))))

cat("N =", nrow(dati), "| J =", length(unique(dati$Field)),
    "| profondità:", paste(sort(unique(dati$Bottom)), collapse = ", "), "cm\n\n")


# ── 1. SETUP FORMULE ─────────────────────────────────────────────────────────

# Covariate fisse: stesse di M4rr + PH incluso per N e P
# (in M4rr PH era fissato a 0 per l'organico, ma qui esploriamo le risposte separatamente)
fixed_cov <- "Texture1 + Texture2 + BulkDensity + PH +
              OnFarm + Irrigate + Fertilised + N_Natural"

risposte <- c("logSOC", "logN", "logP")

risultati <- list()


# ── 2. LOOP SU RISPOSTE ───────────────────────────────────────────────────────

for (resp in risposte) {

  cat("\n", strrep("═", 60), "\n")
  cat(" Risposta:", resp, "\n")
  cat(strrep("═", 60), "\n")

  # M0: random intercept solo
  f0 <- as.formula(paste(resp, "~ logBottom +", fixed_cov, "+ (1 | Field)"))

  # M1: random intercept + slope su logBottom (power law)
  f1 <- as.formula(paste(resp, "~ logBottom +", fixed_cov, "+ (1 + logBottom | Field)"))

  # M2: random intercept + slope su Bottom_sc (centrato e scalato)
  # Usa Bottom_sc (scale(Bottom)) invece di logBottom: testa se la relazione
  # è log-lin (esponenziale) anziché log-log (power law)
  f2 <- as.formula(paste(resp, "~ Bottom_sc +", fixed_cov, "+ (1 + Bottom_sc | Field)"))

  m0 <- lmer(f0, data = dati, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa"))
  m1 <- lmer(f1, data = dati, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa"))
  m2 <- lmer(f2, data = dati, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa"))

  # ── AIC ──────────────────────────────────────────────────────────────────
  cat("\n── AIC ──────────────────────────────────────\n")
  cat(sprintf("  M0 (RI solo):          %8.1f\n", AIC(m0)))
  cat(sprintf("  M1 (+ slope logBottom): %8.1f  (Δ = %+.1f vs M0)\n",
              AIC(m1), AIC(m1) - AIC(m0)))
  cat(sprintf("  M2 (+ slope Bottom_sc): %8.1f  (Δ = %+.1f vs M0)\n",
              AIC(m2), AIC(m2) - AIC(m0)))

  # ── Varianze random effects M1 ─────────────────────────────────────────
  vc1 <- VarCorr(m1)
  re1 <- as.data.frame(vc1)
  var_int1   <- re1[re1$var1 == "(Intercept)" & is.na(re1$var2), "vcov"]
  var_slp1   <- re1[re1$var1 == "logBottom"   & is.na(re1$var2), "vcov"]
  cov_is1    <- re1[!is.na(re1$var2), "vcov"]
  cor_is1    <- cov_is1 / sqrt(var_int1 * var_slp1)
  var_eps1   <- sigma(m1)^2

  cat("\n── Varianze M1 (slope = logBottom) ─────────\n")
  cat(sprintf("  SD intercetta per Field: %.3f\n", sqrt(var_int1)))
  cat(sprintf("  SD slope logBottom:      %.3f\n", sqrt(var_slp1)))
  cat(sprintf("  Corr intercetta-slope:   %.3f\n", cor_is1))
  cat(sprintf("  SD residua:              %.3f\n", sqrt(var_eps1)))

  # ── Varianze random effects M2 ─────────────────────────────────────────
  vc2 <- VarCorr(m2)
  re2 <- as.data.frame(vc2)
  var_int2   <- re2[re2$var1 == "(Intercept)" & is.na(re2$var2), "vcov"]
  var_slp2   <- re2[re2$var1 == "Bottom_sc"  & is.na(re2$var2), "vcov"]
  cov_is2    <- re2[!is.na(re2$var2), "vcov"]
  cor_is2    <- cov_is2 / sqrt(var_int2 * var_slp2)

  cat("\n── Varianze M2 (slope = scale(Bottom)) ─────\n")
  cat(sprintf("  SD intercetta per Field: %.3f\n", sqrt(var_int2)))
  cat(sprintf("  SD slope Bottom_sc:      %.3f\n", sqrt(var_slp2)))
  cat(sprintf("  Corr intercetta-slope:   %.3f\n", cor_is2))

  # ── LRT: M0 vs M1 ────────────────────────────────────────────────────
  cat("\n── LRT M0 vs M1 ────────────────────────────\n")
  print(anova(m0, m1))

  # ── Effetti fissi M1 ──────────────────────────────────────────────────
  cat("\n── Effetti fissi M1 ────────────────────────\n")
  fe <- fixef(m1)
  se <- sqrt(diag(vcov(m1)))
  tab_fe <- data.frame(
    stima = round(fe, 3),
    se    = round(se, 3),
    z     = round(fe / se, 2)
  )
  print(tab_fe)

  risultati[[resp]] <- list(m0 = m0, m1 = m1, m2 = m2,
                             aic = c(M0 = AIC(m0), M1 = AIC(m1), M2 = AIC(m2)),
                             cor_is_logBottom = cor_is1,
                             cor_is_Bottom    = cor_is2,
                             sd_slope_logBottom = sqrt(var_slp1),
                             sd_slope_Bottom    = sqrt(var_slp2))
}


# ── 3. TABELLA RIASSUNTIVA ───────────────────────────────────────────────────

cat("\n\n", strrep("═", 60), "\n")
cat(" TABELLA RIASSUNTIVA\n")
cat(strrep("═", 60), "\n\n")

tab <- map_dfr(risposte, function(r) {
  res <- risultati[[r]]
  tibble(
    risposta          = r,
    AIC_M0            = round(res$aic["M0"], 1),
    AIC_M1            = round(res$aic["M1"], 1),
    delta_M1_M0       = round(res$aic["M1"] - res$aic["M0"], 1),
    AIC_M2            = round(res$aic["M2"], 1),
    delta_M2_M0       = round(res$aic["M2"] - res$aic["M0"], 1),
    SD_slope_logBottom = round(res$sd_slope_logBottom, 3),
    Corr_int_slp_log  = round(res$cor_is_logBottom, 3),
    SD_slope_Bottom   = round(res$sd_slope_Bottom, 3),
    Corr_int_slp_lin  = round(res$cor_is_Bottom, 3)
  )
})
print(as.data.frame(tab))


# ── 4. GRAFICI PROFILI PER FIELD ─────────────────────────────────────────────

for (resp in risposte) {

  m1  <- risultati[[resp]]$m1
  re  <- ranef(m1)$Field |>
    rownames_to_column("Field") |>
    rename(u_int = `(Intercept)`, u_slp = logBottom)

  # Scatter intercetta vs slope (segno della correlazione)
  p_scatter <- ggplot(re, aes(x = u_int, y = u_slp)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, colour = "firebrick", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "gray60") +
    labs(
      title    = paste("Random effects M1 —", resp),
      subtitle = paste("Corr intercetta-slope:",
                       round(risultati[[resp]]$cor_is_logBottom, 3)),
      x = "u_intercetta (deviazione da media Field)",
      y = "u_slope (deviazione tasso decadimento)"
    ) +
    theme_minimal()
  print(p_scatter)

  # Profili predetti per Field (componente random: intercetta + slope)
  depth_seq <- seq(min(dati$logBottom), max(dati$logBottom), length.out = 50)
  beta_int  <- fixef(m1)["(Intercept)"]
  beta_slp  <- fixef(m1)["logBottom"]

  profili <- re |>
    mutate(Field = as.integer(Field)) |>
    rowwise() |>
    mutate(
      data = list(
        tibble(
          logBottom_std = depth_seq,
          pred = (beta_int + u_int) + (beta_slp + u_slp) * depth_seq
        )
      )
    ) |>
    unnest(data)

  p_profili <- ggplot(profili, aes(x = logBottom_std, y = pred, group = Field)) +
    geom_line(alpha = 0.5, linewidth = 0.5, colour = "steelblue") +
    geom_abline(intercept = beta_int, slope = beta_slp,
                colour = "firebrick", linewidth = 1.2, linetype = "dashed") +
    labs(
      title    = paste("Profili per Field (random intercept + slope) —", resp),
      subtitle = "Linea rossa = effetto medio fisso",
      x = "logBottom (standardizzato)",
      y = paste("Predetto:", resp)
    ) +
    theme_minimal()
  print(p_profili)
}


# ── 5. CORRELAZIONE TRA SLOPE DEI TRE FATTORI ───────────────────────────────
# Se la slope di SOC e la slope di N per Field sono correlate, la variazione
# del tasso di decadimento è una caratteristica del sito (non della risposta).

slopes_df <- map_dfr(risposte, function(r) {
  ranef(risultati[[r]]$m1)$Field |>
    rownames_to_column("Field") |>
    select(Field, slope = logBottom) |>
    mutate(risposta = r)
}) |>
  pivot_wider(names_from = risposta, values_from = slope)

cat("\n── Correlazione tra slopes per Field (M1) ──\n")
print(round(cor(slopes_df[, risposte], use = "complete.obs"), 3))

cat("\n── Fine esplorazione random slopes ──────────────────────────────────\n")
