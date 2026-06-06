# =============================================================================
# 20_spatial_confounding.R  —  Diagnostica spatial confounding
#
# DOMANDA: le covariate del modello (X_B e X_W) sono spazialmente clusterizzate?
# Se sì, gli effetti stimati beta_r e gamma_r potrebbero essere distorti perché
# la gestione agricola e le proprietà del suolo confondono con la posizione.
#
# ANALISI:
#   1. Moran's I su ciascuna variabile X_B (gestione) e X_W (suolo) a livello campo
#   2. Correlazione di Pearson tra covariata e coordinate geografiche (Lat/Long)
#   3. Moran's I sugli intercetti posteriori MVRE (dal fit finale)
#   4. Scatter: intercetti posteriori MVRE vs coordinate → visualizzazione confounding
#   5. Tabella riassuntiva
#
# Dipende da: data/dati.rds, data/crop.rds, stan/fit_msp_rirs_mvre.rds
# =============================================================================

library(tidyverse)
library(here)
library(spdep)

if (requireNamespace("posterior", quietly = TRUE)) library(posterior)


# ── 0. DATI E COORDINATE ──────────────────────────────────────────────────────

dati <- readRDS(here("data", "dati.rds")) |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(logBottom = log(Bottom), Field = factor(Field))

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

cat(sprintf("Area: %.1f × %.1f km | J = %d campi\n",
            diff(range(coords$x_km)), diff(range(coords$y_km)), nrow(coords)))


# ── 1. DATI A LIVELLO CAMPO ───────────────────────────────────────────────────

# Per X_B: un valore per campo (binarie → media = frequenza)
# Per X_W: media per campo (suolo, tessitura, ecc.)
field_df <- dati |>
  mutate(across(c(logBottom, Texture1, Texture2, BulkDensity, PH),
                ~ c(scale(.x)))) |>
  group_by(Field) |>
  summarise(
    OnFarm     = mean(OnFarm),
    Irrigate   = mean(Irrigate),
    Fertilised = mean(Fertilised),
    N_Natural  = mean(N_Natural),
    Texture1   = mean(Texture1),
    Texture2   = mean(Texture2),
    BulkDensity = mean(BulkDensity),
    PH         = mean(PH),
    .groups    = "drop"
  ) |>
  left_join(coords, by = "Field") |>
  arrange(Field)

cat(sprintf("Campi nel dataset field_df: %d\n", nrow(field_df)))


# ── 2. MATRICE PESI SPAZIALI ──────────────────────────────────────────────────

coords_mat <- as.matrix(field_df[, c("x_km", "y_km")])
D <- as.matrix(dist(coords_mat))

# Pesi 1/distanza con soglia 0.3 km (identico a script 05)
W <- ifelse(D > 0.3, 1 / D, 0)
diag(W) <- 0
W_std <- W / rowSums(W)   # row-standardized

lw <- mat2listw(W_std, style = "W")


# ── 3. MORAN'S I PER CIASCUNA COVARIATA ──────────────────────────────────────

cat("\n═══ MORAN'S I SULLE COVARIATE DI CAMPO ══════════════════════════\n")
cat(sprintf("%-14s | %8s | %8s | %8s\n", "Variabile", "Moran I", "Z-score", "p-value"))
cat(strrep("-", 46), "\n")

covariates <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural",
                "Texture1", "Texture2", "BulkDensity", "PH")

moran_tab <- map_dfr(covariates, function(v) {
  x  <- field_df[[v]]
  mt <- moran.test(x, lw, alternative = "greater")
  z  <- (mt$estimate["Moran I statistic"] - mt$estimate["Expectation"]) /
        sqrt(mt$estimate["Variance"])
  tibble(
    variabile = v,
    moran_I   = mt$estimate["Moran I statistic"],
    z_score   = z,
    p_value   = mt$p.value,
    tipo      = if_else(v %in% c("OnFarm","Irrigate","Fertilised","N_Natural"),
                        "X_B (gestione)", "X_W (suolo)")
  )
})

for (i in seq_len(nrow(moran_tab))) {
  r <- moran_tab[i, ]
  cat(sprintf("%-14s | %+8.4f | %+8.2f | %.4f  %s\n",
              r$variabile, r$moran_I, r$z_score, r$p_value,
              if_else(r$p_value < 0.05, "***", if_else(r$p_value < 0.10, "*", ""))))
}

cat("\nNote: *** p<0.05, * p<0.10  (test unilaterale alternativa 'greater')\n")


# ── 4. CORRELAZIONE CON COORDINATE ───────────────────────────────────────────

cat("\n═══ CORRELAZIONE CON COORDINATE GEOGRAFICHE ══════════════════════\n")
cat(sprintf("%-14s | %8s %8s | %8s %8s\n",
            "Variabile", "cor(Lat)", "cor(Long)", "p(Lat)", "p(Long)"))
cat(strrep("-", 56), "\n")

cor_tab <- map_dfr(covariates, function(v) {
  x  <- field_df[[v]]
  ct_lat  <- cor.test(x, field_df$Lat,  method = "spearman", exact = FALSE)
  ct_long <- cor.test(x, field_df$Long, method = "spearman", exact = FALSE)
  tibble(
    variabile = v,
    cor_lat   = ct_lat$estimate,
    p_lat     = ct_lat$p.value,
    cor_long  = ct_long$estimate,
    p_long    = ct_long$p.value
  )
})

for (i in seq_len(nrow(cor_tab))) {
  r <- cor_tab[i, ]
  cat(sprintf("%-14s | %+8.3f %8.3f | %8.4f %8.4f\n",
              r$variabile, r$cor_lat, r$cor_long, r$p_lat, r$p_long))
}

cat("\n(Spearman; p-value non corretti per test multipli)\n")


# ── 5. INTERCETTI POSTERIORI MVRE ─────────────────────────────────────────────

fit_path <- here("stan", "fit_msp_rirs_mvre.rds")
if (!file.exists(fit_path)) {
  cat("\n[SKIP] fit_msp_rirs_mvre.rds non trovato — salto analisi intercetti.\n")
} else {
  cat("\n═══ MORAN SU INTERCETTI POSTERIORI MVRE ══════════════════════════\n")

  fit_mvre <- readRDS(fit_path)

  field_levels <- sort(unique(as.integer(as.character(dati$Field))))
  J            <- length(field_levels)

  draws <- fit_mvre$draws(
    variables = c(paste0("V[1,", 1:J, "]"),
                  paste0("V[3,", 1:J, "]"),
                  paste0("V[5,", 1:J, "]")),
    format = "matrix"
  )

  re_geo <- coords |>
    mutate(
      u_int_SOC = colMeans(draws[, paste0("V[1,", 1:J, "]")]),
      u_int_N   = colMeans(draws[, paste0("V[3,", 1:J, "]")]),
      u_int_P   = colMeans(draws[, paste0("V[5,", 1:J, "]")])
    )

  cat(sprintf("%-12s | %8s | %8s | %8s\n",
              "Intercetta", "Moran I", "Z-score", "p-value"))
  cat(strrep("-", 46), "\n")

  for (v in c("u_int_SOC", "u_int_N", "u_int_P")) {
    mt <- moran.test(re_geo[[v]], lw)
    z  <- (mt$estimate["Moran I statistic"] - mt$estimate["Expectation"]) /
          sqrt(mt$estimate["Variance"])
    cat(sprintf("%-12s | %+8.4f | %+8.2f | %.6f\n",
                v, mt$estimate["Moran I statistic"], z, mt$p.value))
  }

  # Correlazione intercetti posteriori con coordinate
  cat("\n── Correlazione intercetti posteriori con Lat/Long ─────────────\n")
  cat(sprintf("%-12s | %8s %8s\n", "", "cor(Lat)", "cor(Long)"))
  for (v in c("u_int_SOC", "u_int_N", "u_int_P")) {
    cat(sprintf("%-12s | %+8.3f %+8.3f\n", v,
                cor(re_geo[[v]], re_geo$Lat),
                cor(re_geo[[v]], re_geo$Long)))
  }

  # Correlazione intercetti posteriori con gestione (spatial confounding key check)
  cat("\n── Correlazione intercetti posteriori con X_B ────────────────────\n")
  cat(sprintf("%-12s |", ""))
  for (xb in c("OnFarm","Irrigate","Fertilised","N_Natural"))
    cat(sprintf(" %11s", xb))
  cat("\n")
  cat(strrep("-", 60), "\n")

  re_geo2 <- re_geo |> left_join(field_df |> select(Field, all_of(covariates)), by = "Field")

  for (v in c("u_int_SOC", "u_int_N", "u_int_P")) {
    cat(sprintf("%-12s |", v))
    for (xb in c("OnFarm","Irrigate","Fertilised","N_Natural"))
      cat(sprintf(" %+11.3f", cor(re_geo2[[v]], re_geo2[[xb]])))
    cat("\n")
  }

  cat("\n(Se |cor| alto → il campo casuale è correlato con la gestione → confounding)\n")

  rm(fit_mvre, draws); gc()
}


# ── 6. RIEPILOGO ──────────────────────────────────────────────────────────────

cat("\n═══ RIEPILOGO SPATIAL CONFOUNDING ════════════════════════════════\n")
cat("X_B (gestione):\n")
moran_tab |> filter(tipo == "X_B (gestione)") |>
  mutate(flag = if_else(p_value < 0.05, "AUTOCORR. SIGNIFICATIVA",
                        if_else(p_value < 0.10, "BORDERLINE", "ok"))) |>
  select(variabile, moran_I, z_score, p_value, flag) |>
  print(n = Inf)

cat("\nX_W (suolo, medie per campo):\n")
moran_tab |> filter(tipo == "X_W (suolo)") |>
  mutate(flag = if_else(p_value < 0.05, "AUTOCORR. SIGNIFICATIVA",
                        if_else(p_value < 0.10, "BORDERLINE", "ok"))) |>
  select(variabile, moran_I, z_score, p_value, flag) |>
  print(n = Inf)

cat("\n── Fine script 20 ───────────────────────────────────────────────\n")
