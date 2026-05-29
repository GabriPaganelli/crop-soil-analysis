# =============================================================================
# 07_eda_spatial_gp.R  —  Struttura spaziale e motivazione per abbandonare il GP
#
# DOMANDE:
#   A. C'è autocorrelazione spaziale nei residui between-field? (variogramma,
#      Moran's I) → SOC: sì, debole; N e P: no → GP degenere per N e P.
#   B. La relazione intercetta-slope è visibile nei LMM? → conferma eta_SOC forte.
#   C. Correlazioni cross-risposta a livello di campo.
#   D. Range GP stimato (M-GPS) vs range empirico dal variogramma SOC.
#
# Dipende da: data/dati.rds, data/crop.rds
# =============================================================================

suppressMessages({
  library(tidyverse)
  library(here)
  library(lme4)
})

# ── 0. DATI ──────────────────────────────────────────────────────────────────

dati <- readRDS(here("data", "dati.rds")) |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom)
  ) |>
  mutate(across(c(logBottom, Texture1, Texture2, BulkDensity),
                ~ c(scale(.x)))) |>
  mutate(Field = factor(Field))

# Coordinate da crop.rds (ha Lat/Long)
crop_raw <- readRDS(here("data", "crop.rds"))
coords <- crop_raw |>
  group_by(Field) |>
  summarise(Lat = mean(Lat, na.rm = TRUE), Long = mean(Long, na.rm = TRUE)) |>
  mutate(Field = factor(Field, levels = levels(dati$Field)))

ref_lat  <- mean(coords$Lat)
ref_long <- mean(coords$Long)
coords <- coords |>
  mutate(
    x_km = (Long - ref_long) * 111.32 * cos(ref_lat * pi / 180),
    y_km = (Lat  - ref_lat)  * 110.54
  ) |>
  arrange(Field)
rm(crop_raw)

cat("════════════════════════════════════════════════════════════════\n")
cat(" OVERVIEW\n")
cat("════════════════════════════════════════════════════════════════\n")
cat(sprintf("N = %d | J = %d campi | profondità: %s cm\n",
            nrow(dati), length(unique(dati$Field)),
            paste(sort(unique(dati$Bottom)), collapse = ", ")))
cat(sprintf("Area: %.1f × %.1f km (x: [%.1f, %.1f]  y: [%.1f, %.1f])\n",
            diff(range(coords$x_km)), diff(range(coords$y_km)),
            min(coords$x_km), max(coords$x_km),
            min(coords$y_km), max(coords$y_km)))
cat(sprintf("Distanza mediana tra campi: %.1f km | max: %.1f km\n",
            median(dist(coords[, c("x_km", "y_km")])),
            max(dist(coords[, c("x_km", "y_km")]))))


# ═══════════════════════════════════════════════════════════════════════════════
# PARTE A — STRUTTURA SPAZIALE DEGLI EFFETTI BETWEEN-FIELD
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n════════════════════════════════════════════════════════════════\n")
cat(" PARTE A — STRUTTURA SPAZIALE\n")
cat("════════════════════════════════════════════════════════════════\n")

risposte <- c("logSOC", "logN", "logP")

# LMM con random intercept → BLUP per campo
lmm_ri <- list()
blups_df <- coords |> select(Field, x_km, y_km)

for (resp in risposte) {
  f <- as.formula(paste(resp,
    "~ logBottom + Texture1 + Texture2 + BulkDensity +
       OnFarm + Irrigate + Fertilised + N_Natural + (1 | Field)"))
  m <- lmer(f, data = dati, REML = TRUE,
            control = lmerControl(optimizer = "bobyqa"))
  lmm_ri[[resp]] <- m

  re <- ranef(m)$Field |>
    rownames_to_column("Field") |>
    rename(!!paste0("blup_", resp) := `(Intercept)`) |>
    mutate(Field = factor(Field, levels = levels(dati$Field)))

  blups_df <- blups_df |>
    left_join(re |> select(Field, starts_with("blup")), by = "Field")
}

# Varianza between-field
cat("\n── Varianza between-field (random intercept dopo covariate) ──\n")
for (resp in risposte) {
  vc  <- as.data.frame(VarCorr(lmm_ri[[resp]]))
  v_b <- vc[vc$grp == "Field", "vcov"]
  v_w <- sigma(lmm_ri[[resp]])^2
  cat(sprintf("  %-8s: SD_between = %.3f | SD_within = %.3f | ICC = %.3f\n",
              resp, sqrt(v_b), sqrt(v_w), v_b / (v_b + v_w)))
}

# Matrice delle distanze
D <- as.matrix(dist(blups_df[, c("x_km", "y_km")]))

# ── A1. Correlazione empirica per bin di distanza ──────────────────────────
cat("\n── Correlazione tra BLUP di campi a distanza h ──\n")
cat("  (se positiva e decrescente → struttura spaziale)\n\n")

breaks <- c(0, 1, 2, 4, 7, 11, 16, 22, Inf)
labs   <- c("<1", "1-2", "2-4", "4-7", "7-11", "11-16", "16-22", ">22")

cat(sprintf("  %-8s | %s\n", "Risposta",
            paste(sprintf("%7s", labs), collapse = " ")))
cat(sprintf("  %-8s | %s\n", "", paste(sprintf("%7s", "km"), collapse = " ")))
cat(strrep("-", 74), "\n")

for (resp in risposte) {
  z <- scale(blups_df[[paste0("blup_", resp)]])[, 1]
  cors <- sapply(seq_along(labs), function(k) {
    idx <- which(D > breaks[k] & D <= breaks[k+1], arr.ind = TRUE)
    idx <- idx[idx[,1] < idx[,2], , drop = FALSE]
    if (nrow(idx) < 4) return(NA_real_)
    cor(z[idx[,1]], z[idx[,2]])
  })
  n_pairs <- sapply(seq_along(labs), function(k) {
    idx <- which(D > breaks[k] & D <= breaks[k+1], arr.ind = TRUE)
    sum(idx[,1] < idx[,2])
  })
  cat(sprintf("  %-8s | %s\n", resp,
              paste(sprintf("%+7.3f", cors), collapse = " ")))
  cat(sprintf("  %-8s | %s\n", paste0("(n=)"),
              paste(sprintf("%7d", n_pairs), collapse = " ")))
  cat("\n")
}

# ── A2. Variogramma empirico ───────────────────────────────────────────────
cat("── Variogramma empirico (semivarianza) ──\n")
cat("  Sill = varianza totale. Se γ̂(h) cresce con h → struttura spaziale.\n\n")

vario_res <- list()
for (resp in risposte) {
  z <- blups_df[[paste0("blup_", resp)]]
  sill_tot <- var(z)

  n_bins <- 10
  max_d  <- 22
  breaks_v <- seq(0, max_d, length.out = n_bins + 1)
  mids   <- (breaks_v[-1] + breaks_v[-(n_bins+1)]) / 2

  gamma <- sapply(seq_len(n_bins), function(k) {
    idx <- which(D > breaks_v[k] & D <= breaks_v[k+1], arr.ind = TRUE)
    idx <- idx[idx[,1] < idx[,2], , drop = FALSE]
    if (nrow(idx) < 4) return(NA_real_)
    mean(0.5 * (z[idx[,1]] - z[idx[,2]])^2)
  })

  vario_res[[resp]] <- tibble(dist = mids, gamma = gamma)
  cat(sprintf("  %s  (sill totale = %.4f):\n", resp, sill_tot))
  cat(sprintf("    dist(km): %s\n",
              paste(sprintf("%6.1f", mids), collapse = " ")))
  cat(sprintf("    γ̂(h):    %s\n",
              paste(sprintf("%6.4f", gamma), collapse = " ")))
  # Frazione del sill raggiunta a diverse distanze
  frac <- gamma / sill_tot
  cat(sprintf("    γ/sill:   %s\n\n",
              paste(sprintf("%6.2f", frac), collapse = " ")))
}

# ── A3. Moran's I ──────────────────────────────────────────────────────────
cat("── Moran's I (pesi = 1/distanza per d > 0.3 km) ──\n")
W <- ifelse(D > 0.3, 1/D, 0)
diag(W) <- 0
W <- W / rowSums(W)

for (resp in risposte) {
  z <- blups_df[[paste0("blup_", resp)]]
  z <- z - mean(z)
  I <- (t(z) %*% W %*% z / t(z) %*% z) |> as.numeric()
  cat(sprintf("  %-8s: Moran's I = %+.3f  (>0 = autocorr. positiva)\n", resp, I))
}


# ═══════════════════════════════════════════════════════════════════════════════
# PARTE B — RELAZIONE INTERCETTA–SLOPE (logBottom)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n════════════════════════════════════════════════════════════════\n")
cat(" PARTE B — INTERCETTA–SLOPE per logBottom\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("  Domanda: i campi con intercetta alta decadono più in fretta?\n")
cat("  Se cor(intercetta, slope) ≈ -1 → struttura u_j*(1+b*logBottom) OK.\n\n")

slope_fits <- list()
for (resp in risposte) {
  f <- as.formula(paste(resp,
    "~ logBottom + Texture1 + Texture2 + BulkDensity +
       OnFarm + Irrigate + Fertilised + N_Natural + (1 + logBottom | Field)"))
  suppressMessages(
    m <- tryCatch(
      lmer(f, data = dati, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 3e5))),
      error = function(e) NULL
    )
  )
  slope_fits[[resp]] <- m
  if (is.null(m)) { cat(sprintf("  %s: modello fallito\n", resp)); next }

  vc  <- as.data.frame(VarCorr(m))
  v_i    <- vc[vc$grp == "Field" & vc$var1 == "(Intercept)" & is.na(vc$var2), "vcov"]
  v_s    <- vc[vc$grp == "Field" & vc$var1 == "logBottom"   & is.na(vc$var2), "vcov"]
  cov_is <- vc[vc$grp == "Field" & !is.na(vc$var2), "vcov"]
  if (length(v_i) == 0) v_i <- 0
  if (length(v_s) == 0) v_s <- 0
  if (length(cov_is) == 0) cov_is <- 0
  v_i <- v_i[1]; v_s <- v_s[1]; cov_is <- cov_is[1]
  cor_is <- if (v_i > 0 && v_s > 0) cov_is / sqrt(v_i * v_s) else NA_real_
  sing   <- isSingular(m)

  re    <- ranef(m)$Field
  b_ols <- coef(lm(re[, "logBottom"] ~ re[, "(Intercept)"]))[2]
  r2    <- cor(re[, "(Intercept)"], re[, "logBottom"])^2

  delta_aic <- AIC(lmm_ri[[resp]]) - AIC(m)

  cat(sprintf("  %s:\n", resp))
  cat(sprintf("    SD intercetta per campo: %.3f\n", sqrt(v_i)))
  cat(sprintf("    SD slope logBottom:      %.3f\n", sqrt(v_s)))
  cat(sprintf("    Corr(intercetta, slope): %+.3f  %s\n",
              cor_is, if (sing) "[SINGOLARE]" else ""))
  cat(sprintf("    b_OLS (slope ~ intercetta): %.3f  (R² = %.2f)\n", b_ols, r2))
  cat(sprintf("    ΔAIC vs random-intercept:   %+.1f\n\n", delta_aic))
}


# ── B2. Scatter BLUP intercetta vs slope ──────────────────────────────────
cat("── BLUP intercetta vs slope per campo (correlazione empirica BLUPs) ──\n")
for (resp in risposte) {
  m <- slope_fits[[resp]]
  if (is.null(m)) next
  re <- ranef(m)$Field
  r  <- cor(re[,1], re[,2])
  cat(sprintf("  %-8s: cor(BLUP_int, BLUP_slope) = %+.3f\n", resp, r))
}


# ── B3. Varianza spiegata dall'interazione intercetta-slope ───────────────
cat("\n── Quanto varia il decadimento tra campi? ──\n")
cat("  (varianza della slope casuale / varianza totale within-field)\n")
for (resp in risposte) {
  m <- slope_fits[[resp]]
  if (is.null(m)) next
  vc  <- as.data.frame(VarCorr(m))
  v_s <- vc[vc$var1 == "logBottom" & is.na(vc$var2), "vcov"]
  v_w <- sigma(m)^2
  # Varianza aggiuntiva dovuta alla slope casuale (su scala logBottom)
  # var(logBottom) * v_s
  v_lb <- var(dati$logBottom)
  extra_var <- v_s * v_lb
  cat(sprintf("  %-8s: SD slope = %.3f | SD residua = %.3f | slope spiega %.1f%% var. within\n",
              resp, sqrt(v_s), sqrt(v_w), 100 * extra_var / (extra_var + v_w)))
}


# ═══════════════════════════════════════════════════════════════════════════════
# PARTE C — CORRELAZIONI CROSS-RISPOSTA A LIVELLO DI CAMPO
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n════════════════════════════════════════════════════════════════\n")
cat(" PARTE C — CORRELAZIONI TRA RISPOSTE A LIVELLO DI CAMPO\n")
cat("════════════════════════════════════════════════════════════════\n")

blups_wide <- blups_df |>
  select(starts_with("blup_")) |>
  rename_with(~ str_remove(.x, "blup_"))

cat("\n  Correlazione tra effetti casuali (BLUP_intercetta) per risposta:\n")
print(round(cor(blups_wide, use = "complete.obs"), 3))

cat("\n  Interpretazione: se SOC-N alta → campi ricchi in SOC tendono\n")
cat("  ad avere anche più N → possibile fattore comune.\n")


# ═══════════════════════════════════════════════════════════════════════════════
# PARTE D — DOMANDA GP: MATÉRN STIMATO VS CORRELAZIONE FISSA
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n════════════════════════════════════════════════════════════════\n")
cat(" PARTE D — GP NECESSARIO? CONFRONTO RANGE STIMATI\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("\n  Dal modello gp_v2_full (già girato):\n")
cat("    SOC: rho ≈ 3.8 km, sigma_GP ≈ 0.26 → range credibile, GP utile\n")
cat("    N:   rho ≈ 20  km, sigma_GP ≈ 1.19 → rho >> area, GP degenere\n")
cat("    P:   rho ≈ 20  km, sigma_GP ≈ 1.14 → rho >> area, GP degenere\n")
cat("\n  Implicazione per SOC: il GP stima rho; potremmo fissare rho\n")
cat("  a un valore ragionevole (~3-5 km) e stimare solo sigma.\n")
cat("  Questo riduce 2 parametri a 1, ma perde flessibilità.\n")
cat("  Con J=40 campi e area limitata, la stima di rho è possibile ma incerta.\n\n")

# Range stimato empiricamente dal variogramma di SOC
vario_soc <- vario_res[["logSOC"]]
sill_soc  <- var(blups_df$blup_logSOC)
# range ≈ prima distanza in cui γ > 0.63 * sill
r_approx <- vario_soc |>
  filter(!is.na(gamma), gamma > 0.63 * sill_soc) |>
  slice(1) |> pull(dist)

cat(sprintf("  Range empirico SOC dal variogramma: %s km\n",
            if (length(r_approx) == 0) ">22" else sprintf("%.1f", r_approx)))
cat(sprintf("  Range stimato dal GP (gp_v2_full mediana): 3.8 km\n"))
cat(sprintf("  → Coerenti? Questo conferma/smentisce l'utilità del GP.\n"))

cat("\n── Fine esplorazione ────────────────────────────────────────────\n")
