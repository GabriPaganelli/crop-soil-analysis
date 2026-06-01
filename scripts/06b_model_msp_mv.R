# =============================================================================
# 06b_model_msp_mv.R  —  M-SP-MV: Modello M-SP con residui multivariati
#
# Estensione di 06_model_msp.R:
#   I residui (eps_SOC, eps_N, eps_P) alla stessa osservazione (campo j,
#   profondità i) sono modellati congiuntamente tramite MVNormal con
#   matrice di correlazione Omega_3x3.
#
# PARAMETRI AGGIUNTIVI:
#   rho_SOC_N, rho_SOC_P, rho_N_P  : correlazioni residue cross-risposta
#
# INTERPRETAZIONE:
#   rho_SOC_N > 0 → a stessa profondità e campo, se SOC è sopra la media
#                   predetta allora N tende ad esserlo (co-variazione intrinseca)
#
# Tutti gli altri parametri (alpha, psi, eta, gamma, beta) restano invariati
# e possono essere confrontati direttamente con M-SP.
#
# Fit salvato in: stan/fit_msp_mv.rds
# =============================================================================


# ── 0. SETUP ──────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)

if (!requireNamespace("cmdstanr", quietly = TRUE))
  install.packages("cmdstanr",
                   repos = c("https://mc-stan.org/r-packages/",
                             "https://cloud.r-project.org"))
library(cmdstanr)
source(here("scripts", "00_utilities.R"))
setup_rtools()

if (requireNamespace("posterior",  quietly = TRUE)) library(posterior)
if (requireNamespace("bayesplot",  quietly = TRUE)) library(bayesplot)


# ── 1. DATI (identici a 06_model_msp.R) ──────────────────────────────────────

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

field_levels <- sort(unique(as.integer(as.character(dati$Field))))
J <- length(field_levels)
N <- nrow(dati)

cat(sprintf("N = %d | J = %d\n", N, J))

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))

X_W_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity")
X_W      <- as.matrix(dati_int[, X_W_cols])

X_B_cols <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")
X_B <- dati_int |>
  distinct(field_int, across(all_of(X_B_cols))) |>
  arrange(field_int) |>
  select(all_of(X_B_cols)) |>
  as.matrix()

stan_data <- list(
  N        = N,
  J        = J,
  K_W      = length(X_W_cols),
  K_B      = length(X_B_cols),
  field_id = dati_int$field_int,
  logSOC   = dati_int$logSOC,
  logN     = dati_int$logN,
  logP     = dati_int$logP,
  X_W      = X_W,
  X_B      = X_B
)

cat(sprintf("Parametri totali attesi: 3x%d + 36 + 3_corr = %d\n", J, 3 * J + 39))
rm(X_W, X_B); gc()


# ── 2. COMPILAZIONE ───────────────────────────────────────────────────────────

stan_file <- here("stan", "m_msp_mv.stan")
cat("\nCompilazione m_msp_mv.stan...\n")
mod <- cmdstan_model(stan_file, compile = TRUE)
cat("OK. CmdStan version:", cmdstan_version(), "\n")


# ── 3. MCMC ───────────────────────────────────────────────────────────────────

fit_path <- here("stan", "fit_msp_mv.rds")

if (!file.exists(fit_path)) {

  cat("\nAvvio MCMC (4 catene x 5000 sampling + 3000 warmup)...\n\n")

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
  cat("Fit salvato in:", fit_path, "\n")

  csv_files  <- list.files(here("stan"), pattern = "\\.csv$",  full.names = TRUE)
  json_files <- list.files(here("stan"), pattern = "\\.json$", full.names = TRUE)
  invisible(file.remove(c(csv_files, json_files)))

} else {
  cat("\nCarico fit da:", fit_path, "\n")
  fit <- readRDS(fit_path)
}

rm(stan_data); gc()


# ── 4. DIAGNOSTICA ────────────────────────────────────────────────────────────

cat("\n=== DIAGNOSTICA ======================================================\n")

smry_full <- fit$summary()
np        <- fit$sampler_diagnostics()

n_divs <- sum(np[, , "divergent__"])
max_td  <- max(np[, , "treedepth__"])
ebfmi   <- apply(np[, , "energy__"], 2, function(e) mean(diff(e)^2) / var(e))

bad_rhat <- smry_full |> filter(rhat > 1.05) |> arrange(desc(rhat))
bad_ess  <- smry_full |> filter(ess_bulk < 400) |> arrange(ess_bulk)

cat("Divergenze:", n_divs, "| Max treedepth:", max_td, "\n")
cat("E-BFMI per catena:", round(ebfmi, 3), "\n")
cat("\nRhat > 1.05:", nrow(bad_rhat))
if (nrow(bad_rhat) > 0) { cat("\n"); print(head(bad_rhat, 10)) } else cat(" — OK\n")
cat("\nESS bulk < 400:", nrow(bad_ess))
if (nrow(bad_ess)  > 0) { cat("\n"); print(head(bad_ess,  10)) } else cat(" — OK\n")


# ── 5. PARAMETRI CHIAVE ───────────────────────────────────────────────────────

smry <- fit$summary(c(
  "alpha_SOC", "psi_SOC", "eta_SOC", "b_SOC", "sigma_SOC",
  "alpha_N",   "psi_N",   "eta_N",   "b_N",   "sigma_N",
  "alpha_P",   "psi_P",   "eta_P",   "b_P",   "sigma_P",
  "rho_SOC_N", "rho_SOC_P", "rho_N_P",
  "gamma_SOC", "gamma_N", "gamma_P",
  "beta_SOC",  "beta_N",  "beta_P"
)) |>
  select(variable, mean, median, sd, q5, q95, rhat, ess_bulk) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\n=== PARAMETRI CHIAVE =================================================\n")
print(smry, n = 70)


# ── 6. CORRELAZIONI RESIDUE (RISULTATO PRINCIPALE DI QUESTO SCRIPT) ──────────

cat("\n=== CORRELAZIONI RESIDUE CROSS-RISPOSTA ==============================\n")
cat("(Interpretazione: correlazione tra residui (eps_SOC, eps_N, eps_P) alla\n")
cat(" stessa osservazione, dopo aver rimosso effetti di profondit'a, tessitura,\n")
cat(" management e variabilita' di campo)\n\n")

for (rho_nm in c("rho_SOC_N", "rho_SOC_P", "rho_N_P")) {
  r <- smry |> filter(variable == rho_nm)
  cat(sprintf("  %-12s : mediana = %+.3f  [CI90: %+.3f, %+.3f]  (sd = %.3f)\n",
              rho_nm, r$median, r$q5, r$q95, r$sd))
}

cat("\n  rho > 0: co-variazione positiva residua (suoli piu' ricchi in un\n")
cat("           nutriente tendono ad essere piu' ricchi nell'altro)\n")
cat("  rho ~ 0: i residui dei tre nutrienti sono indipendenti dopo la\n")
cat("           correzione per profondit'a, tessitura e campo\n")


# ── 7. CONFRONTO M-SP vs M-SP-MV (eta_r) ─────────────────────────────────────

cat("\n=== CONFRONTO eta_r: M-SP (base) vs M-SP-MV ==========================\n")
cat("(se M-SP-MV cambia eta_r, la correlazione cross-risposta assorbiva\n")
cat(" variabilita' che erroneamente compariva in eta)\n\n")

# Carica M-SP per confronto
msp_path <- here("stan", "fit_msp.rds")
if (file.exists(msp_path)) {
  fit_msp  <- readRDS(msp_path)
  smry_msp <- fit_msp$summary(c("eta_SOC", "eta_N", "eta_P")) |>
    select(variable, median, q5, q95) |>
    mutate(modello = "M-SP")
  smry_mv  <- smry |> filter(variable %in% c("eta_SOC", "eta_N", "eta_P")) |>
    select(variable, median, q5, q95) |>
    mutate(modello = "M-SP-MV")
  cat(sprintf("  %-12s | %8s [%6s, %6s]\n", "Modello/Par.", "mediana", "q5", "q95"))
  cat(strrep("-", 45), "\n")
  for (p in c("eta_SOC", "eta_N", "eta_P")) {
    r1 <- smry_msp |> filter(variable == p)
    r2 <- smry_mv  |> filter(variable == p)
    cat(sprintf("  %-12s | M-SP:    %+7.3f [%+6.3f, %+6.3f]\n", p, r1$median, r1$q5, r1$q95))
    cat(sprintf("  %-12s | M-SP-MV: %+7.3f [%+6.3f, %+6.3f]\n", "", r2$median, r2$q5, r2$q95))
  }
  rm(fit_msp, smry_msp); gc()
}


# ── 8. LOO-CV: M-SP vs M-SP-MV ───────────────────────────────────────────────

if (requireNamespace("loo", quietly = TRUE)) {
  library(loo)

  cat("\n=== LOO-CV: M-SP vs M-SP-MV ==========================================\n")
  ll_mv  <- fit$draws("log_lik", format = "matrix")
  loo_mv <- loo(ll_mv, cores = 4)
  cat("\nM-SP-MV:\n"); print(loo_mv)

  msp_path <- here("stan", "fit_msp.rds")
  if (file.exists(msp_path)) {
    fit_msp  <- readRDS(msp_path)
    ll_msp   <- fit_msp$draws("log_lik", format = "matrix")
    loo_msp  <- loo(ll_msp, cores = 4)
    cat("\nConfronto LOO (M-SP-MV vs M-SP):\n")
    print(loo_compare(list(`M-SP-MV` = loo_mv, `M-SP` = loo_msp)))
    cat("\nNota: un ΔELPD positivo per M-SP-MV indica che la correlazione residua\n")
    cat("migliora la predizione. ΔELPD < 2*SE = nessuna differenza significativa.\n")
    rm(fit_msp, ll_msp, loo_msp); gc()
  }

  rm(ll_mv); gc()
}


cat("\n-- Fine script 06b_model_msp_mv.R ------------------------------------\n")
