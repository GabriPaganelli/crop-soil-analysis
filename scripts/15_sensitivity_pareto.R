# =============================================================================
# 15_sensitivity_pareto.R  —  Sensitivity analysis: rimozione osservazioni influenti
#
# Usa i Pareto k dal LOO di M-SP-RIRS-MVRE per identificare le osservazioni
# più influenti (k > 0.7), le rimuove dal dataset e rifita M-SP-RIRS-MVRE.
# Confronta i parametri chiave (rho_r, tau_beta_r, gamma_r, sigma_r) tra
# MVRE-full e MVRE-no-influential per valutare la robustezza.
#
# Logica:
#   - k > 0.7: PSIS-LOO inaffidabile per quell'osservazione → influente
#   - k in [0.5, 0.7]: potenzialmente problematica → riportata ma non rimossa
#   - Rimozione per singola osservazione (non per campo intero)
#
# Output:
#   stan/fit_msp_rirs_mvre_noinfl.rds      → fit M-SP-RIRS-MVRE senza obs influenti
#   output/cache/pareto_k_mvre.rds         → vettore k per M-SP-RIRS-MVRE
#   output/tables/tab_15_influential.csv   → lista osservazioni con k > 0.5
#   output/tables/tab_15_sensitivity.csv   → confronto parametri full vs no-infl
#
# MCMC: identici a MVRE (4 catene × 3000wu + 5000samp, adapt_delta=0.97, seed=2024)
# Dipende da: stan/fit_msp_rirs_mvre.rds, data/dati.rds
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


# ── 1. DATI ───────────────────────────────────────────────────────────────────

cat("\nCarico dati...\n")
dati <- readRDS(here("data", "dati.rds")) |>
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
J_full <- length(field_levels)
N_full <- nrow(dati)

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))
cat(sprintf("Dataset completo: N = %d | J = %d\n", N_full, J_full))


# ── 2. PARETO K DA M-SP ───────────────────────────────────────────────────────

pareto_cache <- file.path(cache_dir, "pareto_k_mvre.rds")

if (!file.exists(pareto_cache)) {
  cat("\nCalcolo Pareto k da M-SP-RIRS-MVRE (carico fit_msp_rirs_mvre.rds)...\n")
  fit_msp_full <- readRDS(here("stan", "fit_msp_rirs_mvre.rds"))
  ll_full      <- fit_msp_full$draws("log_lik", format = "matrix")
  loo_full     <- loo(ll_full, cores = 4)
  pareto_k     <- loo_full$diagnostics$pareto_k
  saveRDS(pareto_k, pareto_cache)
  cat(sprintf("  Pareto k salvato: %d valori\n", length(pareto_k)))
  rm(fit_msp_full, ll_full, loo_full); gc()
} else {
  cat("Carico Pareto k da cache:", pareto_cache, "\n")
  pareto_k <- readRDS(pareto_cache)
}

stopifnot(length(pareto_k) == N_full)


# ── 3. IDENTIFICA OSSERVAZIONI INFLUENTI ─────────────────────────────────────

cat("\n═══ PARETO K — DIAGNOSI ══════════════════════════════════════════\n")
cat(sprintf("  k < 0.5  (buono):        %d osservazioni\n",  sum(pareto_k < 0.5)))
cat(sprintf("  0.5 ≤ k < 0.7 (ok):     %d osservazioni\n",  sum(pareto_k >= 0.5 & pareto_k < 0.7)))
cat(sprintf("  0.7 ≤ k < 1.0 (alto):   %d osservazioni\n",  sum(pareto_k >= 0.7 & pareto_k < 1.0)))
cat(sprintf("  k ≥ 1.0  (molto alto):  %d osservazioni\n",  sum(pareto_k >= 1.0)))

# Tabella osservazioni con k > 0.5
infl_idx  <- which(pareto_k >= 0.5)
remove_idx <- which(pareto_k >= 0.7)  # queste vengono rimosse

cat(sprintf("\n  Osservazioni da rimuovere (k ≥ 0.7): %d\n", length(remove_idx)))

infl_df <- dati_int[infl_idx, ] |>
  mutate(
    obs_idx   = infl_idx,
    pareto_k  = pareto_k[infl_idx],
    categoria = case_when(
      pareto_k >= 1.0 ~ "molto_alto",
      pareto_k >= 0.7 ~ "alto",
      TRUE            ~ "moderato"
    )
  ) |>
  select(obs_idx, pareto_k, categoria, Field, Bottom, logSOC, logN, logP,
         logBottom, Texture2, BulkDensity) |>
  arrange(desc(pareto_k))

cat("\nOsservazioni con k > 0.5:\n")
print(infl_df |> select(obs_idx, pareto_k, categoria, Field, Bottom) |>
      mutate(pareto_k = round(pareto_k, 3)) |> as.data.frame())

write.csv(infl_df |> mutate(across(where(is.numeric), ~round(.x, 3))),
          file.path(tab_dir, "tab_15_influential.csv"), row.names = FALSE)
cat("  Salvato: output/tables/tab_15_influential.csv\n")

if (length(remove_idx) == 0) {
  cat("\nNessuna osservazione con k ≥ 0.7. Sensitivity analysis non necessaria.\n")
  cat("Le conclusioni di M-SP-RIRS-MVRE sono robuste per tutti i punti del dataset.\n")
  cat("\n── Fine script 15 (nessuna rimozione) ──\n")
} else {


# ── 4. DATASET RIDOTTO ────────────────────────────────────────────────────────

cat(sprintf("\nRimuovo %d osservazioni (k ≥ 0.7) su %d totali (%.1f%%).\n",
            length(remove_idx), N_full, 100 * length(remove_idx) / N_full))

dati_red <- dati_int[-remove_idx, ]

# Ricalcola field_int sul dataset ridotto (alcuni campi potrebbero perdere obs)
# I campi con 0 osservazioni rimaste vengono eliminati automaticamente
field_levels_red <- sort(unique(as.integer(as.character(dati_red$Field))))
J_red <- length(field_levels_red)
N_red <- nrow(dati_red)

dati_red <- dati_red |>
  mutate(field_int_red = as.integer(factor(as.integer(as.character(Field)),
                                           levels = field_levels_red)))

cat(sprintf("Dataset ridotto: N = %d | J = %d\n", N_red, J_red))
if (J_red < J_full) {
  removed_fields <- setdiff(field_levels, field_levels_red)
  cat(sprintf("  Campi eliminati completamente: %s\n",
              paste(removed_fields, collapse = ", ")))
}

# Dati Stan per dataset ridotto (stessa struttura di M-SP-RIRS, K_W=5)
X_W_red <- as.matrix(dati_red[, c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")])
X_B_red <- dati_red |>
  distinct(field_int_red, OnFarm, Irrigate, Fertilised, N_Natural) |>
  arrange(field_int_red) |>
  select(OnFarm, Irrigate, Fertilised, N_Natural) |>
  as.matrix()

stan_data_red <- list(
  N        = N_red,
  J        = J_red,
  K_W      = 5L,   # logBottom, Texture1, Texture2, BulkDensity, PH
  K_B      = 4L,
  field_id = dati_red$field_int_red,
  logSOC   = dati_red$logSOC,
  logN     = dati_red$logN,
  logP     = dati_red$logP,
  X_W      = X_W_red,
  X_B      = X_B_red
)
rm(X_W_red, X_B_red); gc()


# ── 5. FIT M-SP SU DATASET RIDOTTO ────────────────────────────────────────────

cat("\n═══ FIT M-SP SENZA INFLUENTI ════════════════════════════════════\n")
fit_noinfl_path <- here("stan", "fit_msp_rirs_mvre_noinfl.rds")

if (!file.exists(fit_noinfl_path)) {
  cat("Compilazione (riusa model_mvre.stan)...\n")
  mod_msp <- cmdstan_model(here("stan", "model_mvre.stan"), compile = TRUE)
  cat("OK.\n\n")

  cat("Avvio MCMC (4 × 3000wu + 5000samp)...\n\n")
  fit_noinfl <- mod_msp$sample(
    data            = stan_data_red,
    seed            = 2024,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 3000,
    iter_sampling   = 5000,
    adapt_delta     = 0.97,
    max_treedepth   = 11,
    refresh         = 500,
    output_dir      = here("stan")
  )
  fit_noinfl$save_object(fit_noinfl_path)
  cat("Fit salvato in:", fit_noinfl_path, "\n")
  invisible(file.remove(list.files(here("stan"), pattern = "\\.csv$", full.names = TRUE)))
} else {
  cat("Carico da cache:", fit_noinfl_path, "\n")
  fit_noinfl <- readRDS(fit_noinfl_path)
}


# ── 6. DIAGNOSTICA ────────────────────────────────────────────────────────────

cat("\n--- Diagnostica M-SP-no-influential ---\n")
np_ni    <- fit_noinfl$sampler_diagnostics()
n_divs   <- sum(np_ni[, , "divergent__"])
max_td   <- max(np_ni[, , "treedepth__"])
smry_ni_full <- fit_noinfl$summary()
bad_rhat <- sum(smry_ni_full$rhat > 1.05, na.rm = TRUE)
bad_ess  <- sum(smry_ni_full$ess_bulk < 400, na.rm = TRUE)
cat(sprintf("  Divergenze: %d | Max treedepth: %d | Rhat>1.05: %d | ESS<400: %d\n",
            n_divs, max_td, bad_rhat, bad_ess))
rm(smry_ni_full); gc()


# ── 7. CONFRONTO PARAMETRI: FULL vs NO-INFLUENTIAL ───────────────────────────

cat("\n═══ CONFRONTO PARAMETRI ══════════════════════════════════════════\n")

vars_key <- c(
  "alpha_SOC", "tau_alpha_SOC", "tau_beta_SOC", "rho_SOC", "sigma_SOC",
  "alpha_N",   "tau_alpha_N",   "tau_beta_N",   "rho_N",   "sigma_N",
  "alpha_P",   "tau_alpha_P",   "tau_beta_P",   "rho_P",   "sigma_P",
  # gamma within-field (logBottom=1, Texture1=2, Texture2=3, BulkDensity=4, PH=5)
  "gamma_SOC[1]", "gamma_SOC[3]", "gamma_SOC[4]", "gamma_SOC[5]",
  "gamma_N[1]",   "gamma_N[3]",   "gamma_N[4]",   "gamma_N[5]",
  "gamma_P[1]",   "gamma_P[3]",   "gamma_P[4]",   "gamma_P[5]",
  # beta between-field (OnFarm=1, Irrigate=2, Fertilised=3, N_Natural=4)
  "beta_SOC[1]", "beta_SOC[2]",
  "beta_N[1]",   "beta_N[2]",
  "beta_P[1]",   "beta_P[2]"
)

# Summary M-SP full (caricato una sola volta)
cat("Carico M-SP-RIRS-MVRE full per confronto...\n")
fit_msp_full2 <- readRDS(here("stan", "fit_msp_rirs_mvre.rds"))
smry_full <- smry_draws(fit_msp_full2, vars_key)
rm(fit_msp_full2); gc()

# Summary M-SP no-influential
smry_noinfl <- smry_draws(fit_noinfl, vars_key)

# Confronto: differenza mediane e sovrapposizione CI90%
labels_map <- c(
  "alpha_SOC"="alpha_SOC", "tau_alpha_SOC"="tau_alpha_SOC", "tau_beta_SOC"="tau_beta_SOC",
  "rho_SOC"="rho_SOC", "sigma_SOC"="sigma_SOC",
  "alpha_N"="alpha_N", "tau_alpha_N"="tau_alpha_N", "tau_beta_N"="tau_beta_N",
  "rho_N"="rho_N", "sigma_N"="sigma_N",
  "alpha_P"="alpha_P", "tau_alpha_P"="tau_alpha_P", "tau_beta_P"="tau_beta_P",
  "rho_P"="rho_P", "sigma_P"="sigma_P",
  "gamma_SOC[1]"="γ_SOC logBottom", "gamma_SOC[3]"="γ_SOC Texture2",
  "gamma_SOC[4]"="γ_SOC BulkDensity", "gamma_SOC[5]"="γ_SOC PH",
  "gamma_N[1]"="γ_N logBottom",  "gamma_N[3]"="γ_N Texture2",
  "gamma_N[4]"="γ_N BulkDensity", "gamma_N[5]"="γ_N PH",
  "gamma_P[1]"="γ_P logBottom",  "gamma_P[3]"="γ_P Texture2",
  "gamma_P[4]"="γ_P BulkDensity", "gamma_P[5]"="γ_P PH",
  "beta_SOC[1]"="β_SOC OnFarm",  "beta_SOC[2]"="β_SOC Irrigate",
  "beta_N[1]"="β_N OnFarm",     "beta_N[2]"="β_N Irrigate",
  "beta_P[1]"="β_P OnFarm",     "beta_P[2]"="β_P Irrigate"
)

compare_tab <- smry_full |>
  rename(med_full = median, q05_full = q05, q95_full = q95) |>
  inner_join(
    smry_noinfl |> rename(med_ni = median, q05_ni = q05, q95_ni = q95),
    by = "variable"
  ) |>
  mutate(
    label    = labels_map[variable],
    diff_abs = abs(med_full - med_ni),
    # Segnale di divergenza: la mediana no-infl cade fuori dal CI90% full
    fuori_ci = med_ni < q05_full | med_ni > q95_full,
    flag     = case_when(
      fuori_ci          ~ "ATTENZIONE",
      diff_abs > 0.05   ~ "nota",
      TRUE              ~ "OK"
    )
  ) |>
  select(variable, label, med_full, q05_full, q95_full,
         med_ni, q05_ni, q95_ni, diff_abs, flag) |>
  mutate(across(where(is.numeric), ~round(.x, 3)))

# Stampa parametri chiave (eta, psi, gamma principali)
cat("\n--- Parametri strutturali (alpha, tau, rho, sigma) ---\n")
print(
  compare_tab |>
    filter(grepl("^(alpha_|tau_alpha_|tau_beta_|rho_|sigma_)", variable)) |>
    select(variable, med_full, q05_full, q95_full, med_ni, q05_ni, q95_ni, flag) |>
    as.data.frame(),
  row.names = FALSE
)

cat("\n--- Effetti fissi principali (gamma) ---\n")
print(
  compare_tab |>
    filter(grepl("^gamma", variable)) |>
    select(variable, label, med_full, med_ni, diff_abs, flag) |>
    as.data.frame(),
  row.names = FALSE
)

# Eventuali divergenze
attenzione <- compare_tab |> filter(flag == "ATTENZIONE")
if (nrow(attenzione) > 0) {
  cat(sprintf("\nATTENZIONE: %d parametri cambiano significativamente:\n",
              nrow(attenzione)))
  print(attenzione |> select(variable, label, med_full, med_ni, diff_abs) |>
        as.data.frame(), row.names = FALSE)
} else {
  cat("\nNessun parametro cambia significativamente. Conclusioni robuste.\n")
}



# ── 8. SINTESI ────────────────────────────────────────────────────────────────

cat("\n═══ SINTESI SENSITIVITY ANALYSIS ═══════════════════════════════\n")
cat(sprintf("  Osservazioni rimosse: %d / %d (k ≥ 0.7)\n",
            length(remove_idx), N_full))
cat(sprintf("  Campi rimasti: %d / %d\n", J_red, J_full))
cat(sprintf("  Parametri OK (diff < 0.05):  %d\n", sum(compare_tab$flag == "OK")))
cat(sprintf("  Parametri 'nota' (diff ≥ 0.05): %d\n", sum(compare_tab$flag == "nota")))
cat(sprintf("  Parametri ATTENZIONE (fuori CI): %d\n", sum(compare_tab$flag == "ATTENZIONE")))

cat("\n  Parametri chiave (rho_r, tau_beta_r) — confronto diretto:\n")
for (r in c("SOC", "N", "P")) {
  rho_f <- compare_tab |> filter(variable == paste0("rho_",      r))
  tb_f  <- compare_tab |> filter(variable == paste0("tau_beta_", r))
  if (nrow(rho_f) > 0)
    cat(sprintf("    rho_%s:      full = %.3f [%.3f,%.3f]  |  no-infl = %.3f [%.3f,%.3f]  → %s\n",
                r, rho_f$med_full, rho_f$q05_full, rho_f$q95_full,
                rho_f$med_ni, rho_f$q05_ni, rho_f$q95_ni, rho_f$flag))
  if (nrow(tb_f) > 0)
    cat(sprintf("    tau_beta_%s: full = %.3f [%.3f,%.3f]  |  no-infl = %.3f [%.3f,%.3f]  → %s\n",
                r, tb_f$med_full, tb_f$q05_full, tb_f$q95_full,
                tb_f$med_ni, tb_f$q05_ni, tb_f$q95_ni, tb_f$flag))
}

write.csv(compare_tab, file.path(tab_dir, "tab_15_sensitivity.csv"), row.names = FALSE)
cat("  Salvato: output/tables/tab_15_sensitivity.csv\n")

cat("\n── Fine script 15 ──────────────────────────────────────────────\n")
} # end else (remove_idx > 0)
