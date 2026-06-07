# =============================================================================
# 10_model_final_mvre.R  —  M-SP-RIRS-MVRE: MODELLO FINALE (RE 6D cross-risposta)
#
# Modello finale dell'analisi. Estende M-SP-RIRS modellando tutti gli effetti
# casuali di campo (intercette e slope per SOC, N, P) come MVN 6D congiunta.
# Omega_6 ~ LKJ(2) → 15 correlazioni libere.
#
# Risultato chiave:
#   rho_int_SOC_N ≈ +0.386 [0.038, 0.667]: i campi con alto SOC tendono ad
#   avere alto N (legame biologico identificabile a livello di campo).
#   rho_SOC = +0.676: M-SP-RIRS confermato per SOC (slope proporzionale a intercetta).
#
# Struttura V[6,J]: [SOC_int, SOC_slope, N_int, N_slope, P_int, P_slope]
# mu_r[i,j] = alpha_r + V[2r-1,j] + V[2r,j]*logBottom_i + gamma_r*X_W + beta_r*X_B
#
# Stan: stan/model_mvre.stan
# Fit:  stan/fit_msp_rirs_mvre.rds
# TOTALE PARAMETRI: 6J + 57 = 297 (J=40)
# =============================================================================


# ── 0. SETUP ──────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
library(cmdstanr)
source(here("scripts", "00_utilities.R"))
setup_rtools()

if (requireNamespace("posterior",  quietly = TRUE)) library(posterior)
if (requireNamespace("bayesplot",  quietly = TRUE)) library(bayesplot)
if (requireNamespace("loo",        quietly = TRUE)) library(loo)


# ── 1. DATI ───────────────────────────────────────────────────────────────────

dati <- carica_dati()

field_levels <- sort(unique(as.integer(as.character(dati$Field))))
J <- length(field_levels)
N <- nrow(dati)

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))

X_W_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")
X_B_cols <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")

stan_data <- list(
  N        = N,
  J        = J,
  K_W      = length(X_W_cols),   # 5
  K_B      = length(X_B_cols),
  field_id = dati_int$field_int,
  logSOC   = dati_int$logSOC,
  logN     = dati_int$logN,
  logP     = dati_int$logP,
  X_W      = as.matrix(dati_int[, X_W_cols]),
  X_B      = dati_int |> distinct(field_int, across(all_of(X_B_cols))) |>
               arrange(field_int) |> select(all_of(X_B_cols)) |> as.matrix()
)

cat(sprintf("N = %d | J = %d | K_W = %d\n", N, J, length(X_W_cols)))


# ── 2. FIT ────────────────────────────────────────────────────────────────────

fit_path  <- here("stan", "fit_msp_rirs_mvre.rds")
stan_file <- here("stan", "model_mvre.stan")

if (!file.exists(fit_path)) {
  mod <- cmdstan_model(stan_file, compile = TRUE)
  cat("Avvio MCMC M-SP-MVRE (4 catene × 5000 sampling + 3000 warmup)...\n\n")
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
  cat("Carico fit da:", fit_path, "\n")
  fit <- readRDS(fit_path)
}


# ── 3. DIAGNOSTICA ────────────────────────────────────────────────────────────

cat("\n═══ DIAGNOSTICA M-SP-MVRE ════════════════════════════════════════\n")

smry <- fit$summary()
np   <- fit$sampler_diagnostics()

n_divs <- sum(np[, , "divergent__"])
max_td  <- max(np[, , "treedepth__"])
ebfmi   <- apply(np[, , "energy__"], 2, function(e) mean(diff(e)^2) / var(e))

bad_rhat <- smry |> filter(rhat > 1.05) |> arrange(desc(rhat))
bad_ess  <- smry |> filter(ess_bulk < 400) |> arrange(ess_bulk)

cat("Divergenze:", n_divs, "| Max treedepth:", max_td, "\n")
cat("E-BFMI per catena:", round(ebfmi, 3), "\n")
cat("\nRhat > 1.05:", nrow(bad_rhat))
if (nrow(bad_rhat) > 0) { cat("\n"); print(head(bad_rhat, 10)) } else cat(" — OK\n")
cat("\nESS bulk < 400:", nrow(bad_ess))
if (nrow(bad_ess)  > 0) { cat("\n"); print(head(bad_ess,  10)) } else cat(" — OK\n")


# ── 4. PARAMETRI CHIAVE ───────────────────────────────────────────────────────

params_chiave <- c(
  "alpha_SOC", "tau_alpha_SOC", "tau_beta_SOC", "rho_SOC", "sigma_SOC",
  "alpha_N",   "tau_alpha_N",   "tau_beta_N",   "rho_N",   "sigma_N",
  "alpha_P",   "tau_alpha_P",   "tau_beta_P",   "rho_P",   "sigma_P",
  # correlazioni cross-risposta degli intercetti
  "rho_int_SOC_N", "rho_int_SOC_P", "rho_int_N_P",
  # correlazioni cross-risposta delle slope
  "rho_slope_SOC_N", "rho_slope_SOC_P", "rho_slope_N_P"
)

smry_key <- fit$summary(params_chiave) |>
  select(variable, median, sd, q5, q95, rhat, ess_bulk) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\n═══ STRUTTURA BIVARIATA RI+RS ════════════════════════════════════\n")
cat(sprintf("%-8s | %10s %8s | %10s %8s | %8s %8s | %8s\n",
            "Risposta", "tau_alpha", "sd", "tau_beta", "sd", "rho", "sd", "sigma"))
cat(strrep("-", 80), "\n")
for (r in c("SOC", "N", "P")) {
  ta  <- smry_key |> filter(variable == paste0("tau_alpha_", r))
  tb  <- smry_key |> filter(variable == paste0("tau_beta_",  r))
  rho <- smry_key |> filter(variable == paste0("rho_",       r))
  sg  <- smry_key |> filter(variable == paste0("sigma_",     r))
  cat(sprintf("%-8s | %10.3f %8.3f | %10.3f %8.3f | %8.3f %8.3f | %8.3f\n",
              paste0("log", r), ta$median, ta$sd, tb$median, tb$sd,
              rho$median, rho$sd, sg$median))
}

cat("\n═══ CORRELAZIONI CROSS-RISPOSTA DEGLI INTERCETTI ════════════════\n")
cat(sprintf("%-20s | %8s %8s | %8s %8s\n", "Correlazione", "Mediana","SD","CI5","CI95"))
cat(strrep("-", 60), "\n")
for (v in c("rho_int_SOC_N","rho_int_SOC_P","rho_int_N_P")) {
  r <- smry_key |> filter(variable == v)
  cat(sprintf("%-20s | %8.3f %8.3f | %8.3f %8.3f\n", v, r$median,r$sd,r$q5,r$q95))
}

cat("\n═══ CORRELAZIONI CROSS-RISPOSTA DELLE SLOPE ═════════════════════\n")
cat(sprintf("%-22s | %8s %8s | %8s %8s\n", "Correlazione", "Mediana","SD","CI5","CI95"))
cat(strrep("-", 62), "\n")
for (v in c("rho_slope_SOC_N","rho_slope_SOC_P","rho_slope_N_P")) {
  r <- smry_key |> filter(variable == v)
  cat(sprintf("%-22s | %8.3f %8.3f | %8.3f %8.3f\n", v, r$median,r$sd,r$q5,r$q95))
}


# ── 5. LOO-CV COMPLETO: MVRE vs M-SP-RIRS vs M-SP vs M-RI ──────────────────

cat("\n═══ LOO-CV COMPLETO: M-SP-RIRS-MVRE vs tutti i modelli ════════════\n")

cat("LOO M-SP-RIRS-MVRE (modello finale)...\n")
loo_mvre <- loo(fit$draws("log_lik", format = "matrix"), cores = 4)
cat("\nM-SP-RIRS-MVRE:\n"); print(loo_mvre)

loo_list <- list(`M-SP-RIRS-MVRE` = loo_mvre)

model_map <- c(
  `M-SP-RIRS` = "fit_msp_rirs",
  `M-SP`      = "fit_msp",
  `M-RI`      = "fit_mri"
)
for (lbl in names(model_map)) {
  p <- here("stan", paste0(model_map[[lbl]], ".rds"))
  if (file.exists(p)) {
    cat(sprintf("LOO %s...\n", lbl))
    f <- readRDS(p)
    loo_list[[lbl]] <- loo(f$draws("log_lik", format = "matrix"), cores = 4)
    rm(f); gc()
  } else {
    cat(sprintf("  SKIP %s: file non trovato.\n", lbl))
  }
}

if (length(loo_list) > 1) {
  cat("\n═══ CONFRONTO LOO (migliore in cima) ══════════════════════════════\n")
  cmp <- loo_compare(loo_list)
  print(cmp)
  cat("\n")
  for (lbl in names(model_map)) {
    if (!is.null(loo_list[[lbl]])) {
      delta <- loo_mvre$estimates["elpd_loo","Estimate"] -
               loo_list[[lbl]]$estimates["elpd_loo","Estimate"]
      cat(sprintf("  ΔELPD MVRE vs %s = %+.1f\n", lbl, delta))
    }
  }
  cat("\n  (> 0 → MVRE migliore; < 0 → modello più semplice preferibile)\n")
}

cat("\n── Fine script 10 ───────────────────────────────────────────────\n")
