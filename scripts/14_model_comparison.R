# =============================================================================
# 14_model_comparison.R  —  Confronto modelli A / B / M-SP-RIRS-MVRE (LOO + parametri)
#
# Testa due varianti semplificate di M-SP-RIRS-MVRE (modello finale):
#
#   Modello A (model_mvre_A.stan):
#     MVRE senza effetti between-field (beta) e senza Texture1 (ILR1).
#     X_W = {logBottom[1], Texture2[2], BulkDensity[3], PH[4]}  (K_W = 4)
#     Domanda: management e ILR1 contribuiscono predittivamente?
#
#   Modello B (model_mvre_B.stan):
#     Response-specific: SOC/P con RI+RS; N con solo RI (struttura 5D).
#     X_W come MVRE (K_W = 5), X_B presente.
#     Domanda: la variabilità di slope per N è predittivamente rilevante?
#
# Output:
#   stan/fit_mvre_A.rds              → fit Stan modello A
#   stan/fit_mvre_B.rds              → fit Stan modello B
#   output/cache/loo_mvre_AB.rds     → lista LOO: M-SP-RIRS-MVRE, A, B
#   output/tables/tab_14_params.csv  → tabella confronto parametri
#   output/tables/tab_14_loo.csv     → tabella confronto LOO
#
# MCMC: 4 catene × 3000wu + 5000samp, adapt_delta=0.97, seed=2024
# Dipende da: stan/fit_msp_rirs_mvre.rds, data/crop_analytic.rds
# =============================================================================


# ── 0. SETUP ──────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
library(cmdstanr)
library(posterior)
library(loo)
source(here("scripts", "00_utilities.R"))
setup_rtools()

dir.create(here("output", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "cache"),  recursive = TRUE, showWarnings = FALSE)
tab_dir   <- here("output", "tables")
cache_dir <- here("output", "cache")

cat("CmdStan version:", cmdstan_version(), "\n")

smry_draws <- function(fit, vars) {
  fit$summary(variables = vars) |>
    select(variable, median, q5, q95) |>
    rename(q05 = q5) |>
    mutate(across(c(median, q05, q95), ~round(.x, 3)))
}

fmt_row <- function(smry_df, var_name) {
  row <- smry_df |> filter(variable == var_name)
  if (nrow(row) == 0) return("—")
  sprintf("%.3f [%.3f, %.3f]", row$median[1], row$q05[1], row$q95[1])
}


# ── 1. DATI ───────────────────────────────────────────────────────────────────

cat("\nPreparazione dati...\n")
dati <- readRDS(here("data", "crop_analytic.rds")) |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom)
  ) |>
  mutate(across(c(logBottom, Texture1, Texture2, BulkDensity, PH),
                ~ c(scale(.x)))) |>
  mutate(Field = factor(Field))

field_levels <- sort(unique(as.integer(as.character(dati$Field))))
J <- length(field_levels)
N <- nrow(dati)

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))
cat(sprintf("N = %d | J = %d\n", N, J))

X_B_cols <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")
X_B <- dati_int |>
  distinct(field_int, across(all_of(X_B_cols))) |>
  arrange(field_int) |>
  select(all_of(X_B_cols)) |>
  as.matrix()

# Stan data per modello A (senza Texture1, senza X_B)
# X_W_A: logBottom[1], Texture2[2], BulkDensity[3], PH[4]
X_W_A_cols <- c("logBottom", "Texture2", "BulkDensity", "PH")
stan_data_A <- list(
  N        = N,  J = J,
  K_W      = length(X_W_A_cols),
  field_id = dati_int$field_int,
  logSOC   = dati_int$logSOC,
  logN     = dati_int$logN,
  logP     = dati_int$logP,
  X_W      = as.matrix(dati_int[, X_W_A_cols])
)

# Stan data per modello B (X_W completo, X_B presente)
X_W_B_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")
stan_data_B <- list(
  N        = N,  J = J,
  K_W      = length(X_W_B_cols),
  K_B      = length(X_B_cols),
  field_id = dati_int$field_int,
  logSOC   = dati_int$logSOC,
  logN     = dati_int$logN,
  logP     = dati_int$logP,
  X_W      = as.matrix(dati_int[, X_W_B_cols]),
  X_B      = X_B
)
gc()


# ── 2. MODELLO A ──────────────────────────────────────────────────────────────

cat("\n═══ MODELLO A ═══════════════════════════════════════════════════\n")
fit_A_path <- here("stan", "fit_mvre_A.rds")

if (!file.exists(fit_A_path)) {
  cat("Compilazione model_mvre_A.stan...\n")
  mod_A <- cmdstan_model(here("stan", "model_mvre_A.stan"), compile = TRUE)
  cat("Compilazione OK.\n\n")
  fit_A <- mod_A$sample(
    data            = stan_data_A,
    seed            = 2024,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 3000,
    iter_sampling   = 5000,
    adapt_delta     = 0.97,
    max_treedepth   = 11,
    refresh         = 300,
    output_dir      = here("stan")
  )
  fit_A$save_object(fit_A_path)
  cat("\nFit A salvato in:", fit_A_path, "\n")
  invisible(file.remove(list.files(here("stan"), pattern = "\\.csv$", full.names = TRUE)))
} else {
  cat("Carico fit A da:", fit_A_path, "\n")
  fit_A <- readRDS(fit_A_path)
}

np_A     <- fit_A$sampler_diagnostics()
n_divs_A <- sum(np_A[, , "divergent__"])
max_td_A <- max(np_A[, , "treedepth__"])
smry_A_full  <- fit_A$summary()
bad_rhat_A   <- sum(smry_A_full$rhat > 1.05, na.rm = TRUE)
bad_ess_A    <- sum(smry_A_full$ess_bulk < 400, na.rm = TRUE)
cat(sprintf("  Divergenze: %d | Max treedepth: %d | Rhat>1.05: %d | ESS<400: %d\n",
            n_divs_A, max_td_A, bad_rhat_A, bad_ess_A))

smry_A <- smry_draws(fit_A, c(
  "alpha_SOC", "tau_alpha_SOC", "tau_beta_SOC", "rho_SOC", "sigma_SOC",
  "alpha_N",   "tau_alpha_N",   "tau_beta_N",   "rho_N",   "sigma_N",
  "alpha_P",   "tau_alpha_P",   "tau_beta_P",   "rho_P",   "sigma_P",
  paste0("gamma_SOC[", 1:4, "]"),
  paste0("gamma_N[",   1:4, "]"),
  paste0("gamma_P[",   1:4, "]")
))
print(smry_A, n = 35)
rm(smry_A_full); gc()


# ── 3. MODELLO B ──────────────────────────────────────────────────────────────

cat("\n═══ MODELLO B ═══════════════════════════════════════════════════\n")
fit_B_path <- here("stan", "fit_mvre_B.rds")

if (!file.exists(fit_B_path)) {
  cat("Compilazione model_mvre_B.stan...\n")
  mod_B <- cmdstan_model(here("stan", "model_mvre_B.stan"), compile = TRUE)
  cat("Compilazione OK.\n\n")
  fit_B <- mod_B$sample(
    data            = stan_data_B,
    seed            = 2024,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 3000,
    iter_sampling   = 5000,
    adapt_delta     = 0.97,
    max_treedepth   = 11,
    refresh         = 300,
    output_dir      = here("stan")
  )
  fit_B$save_object(fit_B_path)
  cat("\nFit B salvato in:", fit_B_path, "\n")
  invisible(file.remove(list.files(here("stan"), pattern = "\\.csv$", full.names = TRUE)))
} else {
  cat("Carico fit B da:", fit_B_path, "\n")
  fit_B <- readRDS(fit_B_path)
}

np_B     <- fit_B$sampler_diagnostics()
n_divs_B <- sum(np_B[, , "divergent__"])
max_td_B <- max(np_B[, , "treedepth__"])
smry_B_full  <- fit_B$summary()
bad_rhat_B   <- sum(smry_B_full$rhat > 1.05, na.rm = TRUE)
bad_ess_B    <- sum(smry_B_full$ess_bulk < 400, na.rm = TRUE)
cat(sprintf("  Divergenze: %d | Max treedepth: %d | Rhat>1.05: %d | ESS<400: %d\n",
            n_divs_B, max_td_B, bad_rhat_B, bad_ess_B))

smry_B <- smry_draws(fit_B, c(
  "alpha_SOC", "tau_alpha_SOC", "tau_beta_SOC", "rho_SOC", "sigma_SOC",
  "alpha_N",   "tau_alpha_N",                              "sigma_N",
  "alpha_P",   "tau_alpha_P",   "tau_beta_P",   "rho_P",   "sigma_P"
))
print(smry_B, n = 20)
rm(smry_B_full); gc()


# ── 4. M-SP-RIRS: LOO + SUMMARY ──────────────────────────────────────────────

cat("\n═══ M-SP-RIRS-MVRE (reference) ══════════════════════════════════\n")
loo_cache_path <- file.path(cache_dir, "loo_mvre_AB.rds")

cat("Carico M-SP-RIRS-MVRE...\n")
fit_rirs <- readRDS(here("stan", "fit_msp_rirs_mvre.rds"))

smry_rirs <- smry_draws(fit_rirs, c(
  "alpha_SOC", "tau_alpha_SOC", "tau_beta_SOC", "rho_SOC", "sigma_SOC",
  "alpha_N",   "tau_alpha_N",   "tau_beta_N",   "rho_N",   "sigma_N",
  "alpha_P",   "tau_alpha_P",   "tau_beta_P",   "rho_P",   "sigma_P",
  paste0("gamma_SOC[", 1:5, "]"),
  paste0("gamma_N[",   1:5, "]"),
  paste0("gamma_P[",   1:5, "]")
))

if (!file.exists(loo_cache_path)) {
  cat("Calcolo LOO M-SP-RIRS...\n")
  ll_rirs  <- fit_rirs$draws("log_lik", format = "matrix")
  loo_rirs <- loo(ll_rirs, cores = 4)
  rm(ll_rirs); gc()
}
rm(fit_rirs); gc()


# ── 5. LOO ────────────────────────────────────────────────────────────────────

cat("\n═══ LOO-CV ══════════════════════════════════════════════════════\n")

if (!file.exists(loo_cache_path)) {
  cat("Calcolo LOO modello A...\n")
  ll_A   <- fit_A$draws("log_lik", format = "matrix")
  loo_A  <- loo(ll_A, cores = 4)
  rm(ll_A); gc()

  cat("Calcolo LOO modello B...\n")
  ll_B   <- fit_B$draws("log_lik", format = "matrix")
  loo_B  <- loo(ll_B, cores = 4)
  rm(ll_B); gc()

  loo_list <- list(`M-SP-RIRS-MVRE` = loo_rirs, `A (ridotto)` = loo_A, `B (resp-spec)` = loo_B)
  saveRDS(loo_list, loo_cache_path)
  cat("Cache LOO salvata:", loo_cache_path, "\n")
  rm(loo_rirs, loo_A, loo_B); gc()
} else {
  cat("Carico LOO da cache:", loo_cache_path, "\n")
}

loo_list <- readRDS(loo_cache_path)

for (nm in names(loo_list)) { cat(sprintf("\n%s:\n", nm)); print(loo_list[[nm]]) }

cat("\n--- Confronto LOO (migliore in cima) ---\n")
loo_comp <- loo_compare(loo_list)
print(loo_comp)

loo_tab <- as.data.frame(loo_comp) |>
  tibble::rownames_to_column("modello") |>
  mutate(across(where(is.numeric), ~round(.x, 2)))
write.csv(loo_tab, file.path(tab_dir, "tab_14_loo.csv"), row.names = FALSE)
cat("  Salvato: output/tables/tab_14_loo.csv\n")


# ── 6. TABELLA CONFRONTO PARAMETRI ────────────────────────────────────────────

cat("\n═══ TABELLA CONFRONTO PARAMETRI ═════════════════════════════════\n")
cat("Legenda:\n  M-SP-RIRS-MVRE gamma[1..5]: logBottom, Texture1, Texture2, BulkDensity, PH\n")
cat("  A gamma[1..4]:              logBottom, Texture2, BulkDensity, PH\n\n")

tab_rows <- bind_rows(lapply(c("SOC", "N", "P"), function(r) {
  rows_struct <- tribble(
    ~Parametro, ~`M-SP-RIRS-MVRE`, ~`A (ridotto)`, ~`B (resp-spec)`,
    paste0("alpha_", r),
      fmt_row(smry_rirs, paste0("alpha_", r)),
      fmt_row(smry_A,    paste0("alpha_", r)),
      fmt_row(smry_B,    paste0("alpha_", r)),
    paste0("tau_alpha_", r),
      fmt_row(smry_rirs, paste0("tau_alpha_", r)),
      fmt_row(smry_A,    paste0("tau_alpha_", r)),
      fmt_row(smry_B,    paste0("tau_alpha_", r)),
    paste0("tau_beta_", r),
      fmt_row(smry_rirs, paste0("tau_beta_", r)),
      fmt_row(smry_A,    paste0("tau_beta_", r)),
      if (r == "N") "— (solo RI in B)" else fmt_row(smry_B, paste0("tau_beta_", r)),
    paste0("rho_", r),
      fmt_row(smry_rirs, paste0("rho_", r)),
      fmt_row(smry_A,    paste0("rho_", r)),
      if (r == "N") "— (solo RI in B)" else fmt_row(smry_B, paste0("rho_", r)),
    paste0("sigma_", r),
      fmt_row(smry_rirs, paste0("sigma_", r)),
      fmt_row(smry_A,    paste0("sigma_", r)),
      fmt_row(smry_B,    paste0("sigma_", r))
  )
  rows_gamma <- tribble(
    ~Parametro, ~`M-SP-RIRS-MVRE`, ~`A (ridotto)`, ~`B (resp-spec)`,
    paste0("γ_", r, " logBottom"),
      fmt_row(smry_rirs, paste0("gamma_", r, "[1]")),
      fmt_row(smry_A,    paste0("gamma_", r, "[1]")),
      fmt_row(smry_B,    paste0("gamma_", r, "[1]")),
    paste0("γ_", r, " Texture2"),
      fmt_row(smry_rirs, paste0("gamma_", r, "[3]")),
      fmt_row(smry_A,    paste0("gamma_", r, "[2]")),
      fmt_row(smry_B,    paste0("gamma_", r, "[3]")),
    paste0("γ_", r, " BulkDensity"),
      fmt_row(smry_rirs, paste0("gamma_", r, "[4]")),
      fmt_row(smry_A,    paste0("gamma_", r, "[3]")),
      fmt_row(smry_B,    paste0("gamma_", r, "[4]")),
    paste0("γ_", r, " PH"),
      fmt_row(smry_rirs, paste0("gamma_", r, "[5]")),
      fmt_row(smry_A,    paste0("gamma_", r, "[4]")),
      fmt_row(smry_B,    paste0("gamma_", r, "[5]"))
  )
  bind_rows(rows_struct, rows_gamma) |> mutate(Risposta = r)
}))

for (r in c("SOC", "N", "P")) {
  cat(sprintf("\n--- %s ---\n", r))
  print(tab_rows |> filter(Risposta == r) |> select(-Risposta) |> as.data.frame(),
        row.names = FALSE)
}

write.csv(tab_rows |> select(Risposta, everything()),
          file.path(tab_dir, "tab_14_params.csv"), row.names = FALSE)
cat("\n  Salvato: output/tables/tab_14_params.csv\n")


# ── 7. SINTESI DIAGNOSTICA ────────────────────────────────────────────────────

cat("\n═══ SINTESI DIAGNOSTICA ══════════════════════════════════════════\n")
cat(sprintf("%-14s | %10s | %12s | %9s | %8s\n",
            "Modello", "Divergenze", "MaxTreedepth", "Rhat>1.05", "ESS<400"))
cat(strrep("-", 62), "\n")
cat(sprintf("%-14s | %10d | %12d | %9d | %8d\n",
            "A (ridotto)",   n_divs_A, max_td_A, bad_rhat_A, bad_ess_A))
cat(sprintf("%-14s | %10d | %12d | %9d | %8d\n",
            "B (resp-spec)", n_divs_B, max_td_B, bad_rhat_B, bad_ess_B))

cat("\n── Fine script 14 ──────────────────────────────────────────────\n")
