# =============================================================================
# 22_robustezza_latlong.R  —  Lat/Long come effetti fissi in X_B
#
# Baseline spaziale semplice: aggiunge le coordinate geografiche standardizzate
# (Lat_std, Long_std) come covariate between-field, usando lo stesso Stan file
# di MVRE (model_mvre.stan) con K_B = 6 invece di 4.
#
# Cattura la tendenza spaziale lineare (trend di primo ordine), NON la struttura
# spaziale locale. Atteso: effetti di gestione beta_r cambiano (se confounding),
# LOO quasi invariato o peggiore (aggiunge parametri senza struttura propria).
#
# DOMANDE:
#   A. Gli effetti di gestione beta_r si spostano dopo controllo Lat/Long?
#   B. Lat e Long sono predittori utili (CI lontani da zero)?
#   C. Il LOO migliora o peggiora rispetto a MVRE (base)?
#
# Stan file riusato: stan/model_mvre.stan  (K_B è variabile data → K_B=6 ok)
# Fit file: stan/fit_msp_rirs_latlong.rds
# =============================================================================


# ── 0. SETUP ──────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
library(cmdstanr)
source(here("scripts", "00_utilities.R"))
setup_rtools()

if (requireNamespace("posterior",  quietly = TRUE)) library(posterior)
if (requireNamespace("loo",        quietly = TRUE)) library(loo)


# ── 1. DATI ────────────────────────────────────────────────────────────────────

dati <- carica_dati()

field_levels <- sort(unique(as.integer(as.character(dati$Field))))
J <- length(field_levels)
N <- nrow(dati)

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))


# ── 2. COORDINATE ─────────────────────────────────────────────────────────────

crop_raw <- readRDS(here("data", "crop.rds"))
coords <- crop_raw |>
  group_by(Field) |>
  summarise(Lat = mean(Lat, na.rm = TRUE), Long = mean(Long, na.rm = TRUE)) |>
  mutate(Field = factor(Field, levels = levels(dati$Field))) |>
  arrange(Field)
rm(crop_raw); gc()

# Standardizzazione su scala di dati (non x_km/y_km ma gradi — viene scalato uguale)
coords <- coords |>
  mutate(
    Lat_std  = c(scale(Lat)),
    Long_std = c(scale(Long))
  )

cat(sprintf("Lat range: [%.4f, %.4f] → std [%.2f, %.2f]\n",
            min(coords$Lat), max(coords$Lat),
            min(coords$Lat_std), max(coords$Lat_std)))
cat(sprintf("Long range: [%.4f, %.4f] → std [%.2f, %.2f]\n",
            min(coords$Long), max(coords$Long),
            min(coords$Long_std), max(coords$Long_std)))


# ── 3. STAN DATA CON LAT/LONG IN X_B ─────────────────────────────────────────

X_W_cols    <- c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")
X_B_base    <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")   # 4 originali
X_B_latlong <- c(X_B_base, "Lat_std", "Long_std")                  # 6 con coordinate

field_xb <- dati_int |>
  distinct(field_int, across(all_of(X_B_base))) |>
  arrange(field_int) |>
  left_join(
    coords |> mutate(field_int = as.integer(factor(Field, levels = levels(Field)))) |>
      select(field_int, Lat_std, Long_std),
    by = "field_int"
  ) |>
  select(all_of(X_B_latlong))

stan_data <- list(
  N        = N,
  J        = J,
  K_W      = length(X_W_cols),
  K_B      = length(X_B_latlong),   # 6
  field_id = dati_int$field_int,
  logSOC   = dati_int$logSOC,
  logN     = dati_int$logN,
  logP     = dati_int$logP,
  X_W      = as.matrix(dati_int[, X_W_cols]),
  X_B      = as.matrix(field_xb)
)

cat(sprintf("N = %d | J = %d | K_W = %d | K_B = %d (4 gestione + 2 coordinate)\n",
            N, J, length(X_W_cols), length(X_B_latlong)))


# ── 4. FIT ─────────────────────────────────────────────────────────────────────

fit_path  <- here("stan", "fit_msp_rirs_latlong.rds")
stan_file <- here("stan", "model_mvre.stan")   # stesso Stan file, K_B diverso

if (!file.exists(fit_path)) {
  mod <- cmdstan_model(stan_file, compile = TRUE)
  cat("Avvio MCMC LatLong (4 catene × 5000 sampling + 3000 warmup)...\n\n")
  fit <- mod$sample(
    data            = stan_data,
    seed            = 2024,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 3000,
    iter_sampling   = 5000,
    adapt_delta     = 0.97,
    max_treedepth   = 11,
    refresh         = 200,
    output_dir      = here("stan")
  )
  fit$save_object(fit_path)
  invisible(file.remove(c(
    list.files(here("stan"), pattern = "\\.csv$",  full.names = TRUE),
    list.files(here("stan"), pattern = "\\.json$", full.names = TRUE)
  )))
} else {
  cat("Carico fit LatLong da:", fit_path, "\n")
  fit <- readRDS(fit_path)
}


# ── 5. DIAGNOSTICA ─────────────────────────────────────────────────────────────

cat("\n═══ DIAGNOSTICA MVRE + Lat/Long ════════════════════════════════════\n")

smry <- fit$summary()
np   <- fit$sampler_diagnostics()

n_divs <- sum(np[, , "divergent__"])
max_td  <- max(np[, , "treedepth__"])
ebfmi   <- apply(np[, , "energy__"], 2, function(e) mean(diff(e)^2) / var(e))

cat("Divergenze:", n_divs, "| Max treedepth:", max_td, "\n")
cat("E-BFMI per catena:", round(ebfmi, 3), "\n")

bad_rhat <- smry |> filter(rhat > 1.05) |> arrange(desc(rhat))
bad_ess  <- smry |> filter(ess_bulk < 400) |> arrange(ess_bulk)
cat("Rhat > 1.05:", nrow(bad_rhat))
if (nrow(bad_rhat) > 0) { cat("\n"); print(head(bad_rhat, 10)) } else cat(" — OK\n")
cat("ESS bulk < 400:", nrow(bad_ess))
if (nrow(bad_ess) > 0) { cat("\n"); print(head(bad_ess, 10)) } else cat(" — OK\n")


# ── 6. EFFETTI DI GESTIONE e COORDINATE ──────────────────────────────────────

cat("\n═══ EFFETTI BETWEEN-FIELD beta_r ════════════════════════════════════\n")
cat("(X_B ordine: OnFarm, Irrigate, Fertilised, N_Natural, Lat_std, Long_std)\n\n")

xb_labels <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural", "Lat_std", "Long_std")

cat(sprintf("%-14s | %-12s | %8s | %8s %8s\n",
            "Risposta", "Variabile", "Mediana", "CI5", "CI95"))
cat(strrep("-", 58), "\n")

for (r in c("SOC", "N", "P")) {
  for (k in seq_along(xb_labels)) {
    nm  <- sprintf("beta_%s[%d]", r, k)
    row <- smry |> filter(variable == nm)
    if (nrow(row) == 0) next
    flag <- if_else(sign(row$q5) == sign(row$q95), "(*)", "   ")
    cat(sprintf("log%-8s | %-12s | %+8.3f | %+8.3f %+8.3f  %s\n",
                r, xb_labels[k], row$median, row$q5, row$q95, flag))
  }
  cat("\n")
}

cat("(*) CI5/CI95 stesso segno → effetto non trascurabile\n")


# ── 7. CONFRONTO BETA MVRE vs LATLONG ─────────────────────────────────────────

cat("\n═══ CONFRONTO beta MVRE vs MVRE+LatLong ════════════════════════════\n")
cat("(Spostamento = spatial confounding sugli effetti di gestione)\n\n")

mvre_path <- here("stan", "fit_msp_rirs_mvre.rds")
if (file.exists(mvre_path)) {
  fit_mvre <- readRDS(mvre_path)
  smry_mvre <- fit_mvre$summary(
    variables = c(paste0("beta_SOC[", 1:4, "]"),
                  paste0("beta_N[",   1:4, "]"),
                  paste0("beta_P[",   1:4, "]"))
  )
  rm(fit_mvre); gc()

  cat(sprintf("%-10s | %-12s | %10s | %10s | %8s\n",
              "Risposta", "Variabile", "MVRE med", "LL med", "Δ"))
  cat(strrep("-", 58), "\n")

  for (r in c("SOC", "N", "P")) {
    for (k in 1:4) {
      nm_mvre <- sprintf("beta_%s[%d]", r, k)
      nm_ll   <- sprintf("beta_%s[%d]", r, k)
      row_m <- smry_mvre  |> filter(variable == nm_mvre)
      row_l <- smry        |> filter(variable == nm_ll)
      if (nrow(row_m) == 0 || nrow(row_l) == 0) next
      delta <- row_l$median - row_m$median
      cat(sprintf("log%-7s | %-12s | %+10.3f | %+10.3f | %+8.3f\n",
                  r, xb_labels[k], row_m$median, row_l$median, delta))
    }
  }
  cat("\n(Δ grande → spatial confounding; Δ ≈ 0 → gestione indipendente dalla posizione)\n")
}


# ── 8. LOO-CV ──────────────────────────────────────────────────────────────────

cat("\n═══ LOO-CV: MVRE vs MVRE+LatLong ═══════════════════════════════════\n")

cat("LOO MVRE+LatLong...\n")
loo_ll <- loo(fit$draws("log_lik", format = "matrix"), cores = 4)
cat("\nMVRE+LatLong:\n"); print(loo_ll)

if (file.exists(mvre_path)) {
  cat("\nLOO M-SP-RIRS-MVRE (base)...\n")
  fit_mvre <- readRDS(mvre_path)
  loo_mvre <- loo(fit_mvre$draws("log_lik", format = "matrix"), cores = 4)
  rm(fit_mvre); gc()

  cmp <- loo_compare(list(`MVRE+LatLong` = loo_ll, `MVRE` = loo_mvre))
  cat("\n── Confronto LOO ─────────────────────────────────────────────────\n")
  print(cmp)

  delta <- loo_ll$estimates["elpd_loo","Estimate"] -
           loo_mvre$estimates["elpd_loo","Estimate"]
  cat(sprintf("\nΔELPD (LatLong − MVRE) = %+.1f\n", delta))
  cat("  > 0 → coordinate migliorano il fit\n")
  cat("  < 0 → coordinate solo rumore aggiuntivo (atteso)\n")
}

cat("\n── Fine script 22 ───────────────────────────────────────────────\n")
