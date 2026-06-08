# =============================================================================
# 16_robustezza_msp_mv.R  —  MVRE-FULL: MVRE + residui multivariati (robustezza)
#
# Estensione di M-SP-RIRS-MVRE (modello finale):
#   Aggiunge correlazioni tra residui alla stessa osservazione (i,j) tramite
#   MVNormal. Verifica se, dopo aver rimosso gli effetti di campo (V[6,J]),
#   i residui SOC/N/P sono ancora correlati (stoichiometria residua).
#
# PARAMETRI AGGIUNTIVI rispetto a MVRE:
#   sigma[3]        : [sigma_SOC, sigma_N, sigma_P]
#   L_Omega_res[3x3]: Cholesky di Omega_res (correlazioni residue)
#   → rho_res_SOC_N, rho_res_SOC_P, rho_res_N_P
#
# INTERPRETAZIONE:
#   rho_res ≈ 0 e CI90% include zero → nutrienti condizionalmente indipendenti
#              dato M-SP-RIRS-MVRE: confermato che MVRE è il modello finale
#   rho_res > 0.15 e CI lontano da zero → correlazione residua non catturata da MVRE
#
# RISULTATO ATTESO (basato su test preliminare, script 20):
#   rho_res_SOC_N ≈ +0.015 [-0.121, +0.151] → zero
#   ΔELPD MVRE-FULL vs MVRE = -3.1 (SE=1.1) → MVRE preferibile
#
# Stan: stan/model_mvre_full.stan
# Fit:  stan/fit_mvre_full.rds
# TOTALE PARAMETRI: 6J + 60 = 300 (J=40)
#
# Dipende da: stan/fit_msp_rirs_mvre.rds, stan/fit_mvre_full.rds, data/dati.rds
# =============================================================================


# ── 0. SETUP ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(cmdstanr)
  library(posterior)
  library(loo)
})
source(here("scripts", "00_utilities.R"))
setup_rtools()

dir.create(here("output", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "cache"),  recursive = TRUE, showWarnings = FALSE)
tab_dir   <- here("output", "tables")
cache_dir <- here("output", "cache")

cat("CmdStan version:", cmdstan_version(), "\n")


# ── 1. DATI ───────────────────────────────────────────────────────────────────

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
J <- length(field_levels)
N <- nrow(dati)
cat(sprintf("N = %d | J = %d\n", N, J))

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))

X_W_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")
X_B_cols <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")

stan_data <- list(
  N        = N,
  J        = J,
  K_W      = length(X_W_cols),
  K_B      = length(X_B_cols),
  field_id = dati_int$field_int,
  logSOC   = dati_int$logSOC,
  logN     = dati_int$logN,
  logP     = dati_int$logP,
  X_W      = as.matrix(dati_int[, X_W_cols]),
  X_B      = dati_int |> distinct(field_int, across(all_of(X_B_cols))) |>
               arrange(field_int) |> select(all_of(X_B_cols)) |> as.matrix()
)
rm(dati_int); gc()


# ── 2. FIT MVRE-FULL ─────────────────────────────────────────────────────────

fit_path <- here("stan", "fit_mvre_full.rds")

if (!file.exists(fit_path)) {
  cat("\nCompilazione model_mvre_full.stan...\n")
  mod <- cmdstan_model(here("stan", "model_mvre_full.stan"), compile = TRUE)
  cat("Avvio MCMC (4 catene × 5000 sampling + 3000 warmup)...\n")
  t0 <- proc.time()["elapsed"]
  fit <- mod$sample(
    data            = stan_data,
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
  cat(sprintf("MCMC completato in %.1f min\n", (proc.time()["elapsed"] - t0) / 60))
  fit$save_object(fit_path)
  invisible(file.remove(c(
    list.files(here("stan"), pattern = "\\.csv$",  full.names = TRUE),
    list.files(here("stan"), pattern = "\\.json$", full.names = TRUE)
  )))
} else {
  cat("\nCarico fit MVRE-FULL da:", fit_path, "\n")
  fit <- readRDS(fit_path)
}

rm(stan_data); gc()


# ── 3. DIAGNOSTICA ────────────────────────────────────────────────────────────

cat("\n═══ DIAGNOSTICA MVRE-FULL ════════════════════════════════════\n")
smry_full <- fit$summary()
np        <- fit$sampler_diagnostics()
n_divs    <- sum(np[, , "divergent__"])
max_td    <- max(np[, , "treedepth__"])
ebfmi     <- apply(np[, , "energy__"], 2, function(e) mean(diff(e)^2) / var(e))
bad_rhat  <- smry_full |> filter(rhat > 1.05)
bad_ess   <- smry_full |> filter(ess_bulk < 400)
cat(sprintf("Divergenze: %d | Max treedepth: %d\n", n_divs, max_td))
cat("E-BFMI per catena:", round(ebfmi, 3), "\n")
cat(sprintf("Rhat > 1.05: %d %s\n", nrow(bad_rhat), if (nrow(bad_rhat) == 0) "— OK" else ""))
cat(sprintf("ESS bulk < 400: %d %s\n", nrow(bad_ess),  if (nrow(bad_ess)  == 0) "— OK" else ""))
rm(smry_full); gc()


# ── 4. CORRELAZIONI RESIDUE ────────────────────────────────────────────────────

cat("\n═══ CORRELAZIONI RESIDUE CROSS-RISPOSTA ══════════════════════\n")
cat("(Interpretazione: correlazione tra i residui alla stessa osservazione\n")
cat(" dopo aver rimosso effetti di campo V[6,J], profondità, tessitura,\n")
cat(" management. Risposta alla domanda: MVRE cattura tutta la struttura?)\n\n")

smry_res <- fit$summary(c(
  "rho_res_SOC_N", "rho_res_SOC_P", "rho_res_N_P",
  "rho_int_SOC_N", "rho_int_SOC_P", "rho_int_N_P",
  "rho_SOC", "rho_N", "rho_P",
  "sigma_SOC", "sigma_N", "sigma_P",
  "tau_alpha_SOC", "tau_alpha_N", "tau_alpha_P",
  "tau_beta_SOC",  "tau_beta_N",  "tau_beta_P"
)) |>
  select(variable, median, sd, q5, q95) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat(sprintf("  %-20s | %8s  [%6s, %6s]\n", "Parametro", "mediana", "CI5", "CI95"))
cat(strrep("-", 55), "\n")
for (i in seq_len(nrow(smry_res))) {
  r <- smry_res[i, ]
  cat(sprintf("  %-20s | %+8.3f  [%+6.3f, %+6.3f]\n",
              r$variable, r$median, r$q5, r$q95))
  if (i == 3) cat(strrep("-", 55), "\n")
}

cat("\nConclusione rho_res:\n")
rho_nomi <- c("rho_res_SOC_N", "rho_res_SOC_P", "rho_res_N_P")
for (nm in rho_nomi) {
  r <- smry_res |> filter(variable == nm)
  zero_in_ci <- r$q5 < 0 & r$q95 > 0
  cat(sprintf("  %s: mediana=%+.3f CI90=[%+.3f,%+.3f] → %s\n",
              nm, r$median, r$q5, r$q95,
              if (zero_in_ci) "CI copre zero — indipendenza condizionale" else "CI lontano da zero"))
}


# ── 5. CONFRONTO PARAMETRI STRUTTURALI: MVRE-FULL vs MVRE ────────────────────

cat("\n═══ CONFRONTO rho_r e tau_beta_r: MVRE-FULL vs MVRE ══════════\n")
cat("(verifica: MVRE-FULL cambia le stime degli effetti di campo?)\n\n")

mvre_path <- here("stan", "fit_msp_rirs_mvre.rds")
if (file.exists(mvre_path)) {
  fit_mvre  <- readRDS(mvre_path)
  vars_cmp  <- c("rho_SOC","rho_N","rho_P",
                 "tau_alpha_SOC","tau_alpha_N","tau_alpha_P",
                 "tau_beta_SOC","tau_beta_N","tau_beta_P",
                 "rho_int_SOC_N","rho_int_SOC_P","rho_int_N_P")
  smry_mvre_cmp <- fit_mvre$summary(vars_cmp) |>
    select(variable, median, q5, q95) |>
    rename(med_mvre = median, q05_mvre = q5, q95_mvre = q95)
  smry_full_cmp <- fit$summary(vars_cmp) |>
    select(variable, median, q5, q95) |>
    rename(med_full = median, q05_full = q5, q95_full = q95)
  cmp_df <- inner_join(smry_mvre_cmp, smry_full_cmp, by = "variable") |>
    mutate(diff = abs(med_mvre - med_full))

  cat(sprintf("  %-20s | %-24s | %-24s | %s\n",
              "Parametro", "MVRE med [CI90]", "MVRE-FULL med [CI90]", "Diff"))
  cat(strrep("-", 80), "\n")
  for (i in seq_len(nrow(cmp_df))) {
    r <- cmp_df[i, ]
    cat(sprintf("  %-20s | %+6.3f [%+5.3f,%+5.3f] | %+6.3f [%+5.3f,%+5.3f] | %.3f\n",
                r$variable,
                r$med_mvre, r$q05_mvre, r$q95_mvre,
                r$med_full, r$q05_full, r$q95_full,
                r$diff))
  }

  cat("\n(Se diff piccole: stime stabili — MVRE non perdeva informazione)\n")
  rm(fit_mvre, smry_mvre_cmp); gc()
}


# ── 6. LOO-CV: MVRE-FULL vs MVRE ──────────────────────────────────────────────

cat("\n═══ LOO-CV: MVRE-FULL vs MVRE ════════════════════════════════\n")

loo_cache_path <- file.path(cache_dir, "loo_mvre_vs_full.rds")

if (file.exists(loo_cache_path)) {
  cat("Carico LOO da cache:", loo_cache_path, "\n")
  loo_cache <- readRDS(loo_cache_path)
  loo_full  <- loo_cache$loo_full
  loo_mvre  <- loo_cache$loo_mvre
  cmp       <- loo_cache$comparison
  delta     <- loo_cache$delta
} else {
  cat("Calcolo LOO MVRE-FULL...\n")
  ll_full  <- fit$draws("log_lik", format = "matrix")
  loo_full <- loo(ll_full, cores = 4)
  rm(ll_full); gc()

  if (file.exists(mvre_path)) {
    cat("Calcolo LOO MVRE...\n")
    fit_mvre <- readRDS(mvre_path)
    ll_mvre  <- fit_mvre$draws("log_lik", format = "matrix")
    loo_mvre <- loo(ll_mvre, cores = 4)
    rm(fit_mvre, ll_mvre); gc()
  }

  cmp   <- loo_compare(list(`MVRE-FULL` = loo_full, `MVRE` = loo_mvre))
  delta <- loo_full$estimates["elpd_loo","Estimate"] -
           loo_mvre$estimates["elpd_loo","Estimate"]

  saveRDS(list(loo_full = loo_full, loo_mvre = loo_mvre,
               comparison = cmp, delta = delta), loo_cache_path)
}

cat("\nMVRE-FULL:\n"); print(loo_full)
cat("\n── Confronto LOO ──\n")
print(cmp)
cat(sprintf("\n  ΔELPD (MVRE-FULL − MVRE) = %+.1f  (SE = %.1f)\n",
            delta, cmp[2, "se_diff"]))

cat("\n═══ CONCLUSIONE ══════════════════════════════════════════════\n")
if (delta > 4) {
  cat("  MVRE-FULL superiore: correlazioni residue predittivamente rilevanti.\n")
} else if (delta > -4) {
  cat(sprintf("  ΔELPD = %+.1f: i due modelli sono equivalenti predittivamente.\n", delta))
  cat("  Le correlazioni residue (rho_res ≈ 0) non aggiungono informazione.\n")
  cat("  MVRE è sufficiente: la struttura cross-risposta è catturata\n")
  cat("  interamente a livello di campo (rho_int_SOC_N = +0.38).\n")
} else {
  cat(sprintf("  MVRE preferibile: MVRE-FULL perde %.1f ELPD per complessità inutile.\n",
              abs(delta)))
  cat("  Confermato: MVRE è il modello finale.\n")
}
cat("══════════════════════════════════════════════════════════════\n")

# Salva tabella LOO per figura script 19
loo_tab <- as.data.frame(cmp) |>
  tibble::rownames_to_column("modello") |>
  mutate(across(where(is.numeric), ~round(.x, 2)))
write.csv(loo_tab, file.path(tab_dir, "tab_16_loo_mv.csv"), row.names = FALSE)
cat("\n  Salvato: output/tables/tab_16_loo_mv.csv\n")

cat("\n── Diagnostica residui: autocorrelazione lag-1 (MVRE) ──────────────────\n")

if (file.exists(here("stan", "fit_msp_rirs_mvre.rds"))) {
  fit_mvre_diag <- readRDS(here("stan", "fit_msp_rirs_mvre.rds"))
  dati_diag <- readRDS(here("data", "dati.rds")) |>
    mutate(lb_sc = as.numeric(scale(log(Bottom))),
           bd_sc = as.numeric(scale(BulkDensity)),
           ph_sc = as.numeric(scale(PH)))
  fields_diag <- sort(unique(dati_diag$Field))
  dati_diag$field_idx <- match(dati_diag$Field, fields_diag)

  get_mean_par <- function(par) mean(as.vector(fit_mvre_diag$draws(par)))
  alpha_SOC <- get_mean_par("alpha_SOC"); alpha_N <- get_mean_par("alpha_N"); alpha_P <- get_mean_par("alpha_P")
  gamma_SOC <- sapply(1:5, function(k) get_mean_par(sprintf("gamma_SOC[%d]",k)))
  gamma_N   <- sapply(1:5, function(k) get_mean_par(sprintf("gamma_N[%d]",k)))
  gamma_P   <- sapply(1:5, function(k) get_mean_par(sprintf("gamma_P[%d]",k)))
  beta_SOC  <- sapply(1:4, function(k) get_mean_par(sprintf("beta_SOC[%d]",k)))
  beta_N    <- sapply(1:4, function(k) get_mean_par(sprintf("beta_N[%d]",k)))
  beta_P    <- sapply(1:4, function(k) get_mean_par(sprintf("beta_P[%d]",k)))
  V_mean_diag <- matrix(NA, 6, 40)
  for(k in 1:6) for(j in 1:40)
    V_mean_diag[k,j] <- get_mean_par(sprintf("V[%d,%d]",k,j))

  field_df_diag <- dati_diag |> group_by(Field, field_idx) |>
    summarise(OnFarm=first(OnFarm), Irrigate=first(Irrigate),
              Fertilised=first(Fertilised), N_Natural=first(N_Natural), .groups="drop") |>
    arrange(field_idx)
  X_B_d <- as.matrix(field_df_diag[, c("OnFarm","Irrigate","Fertilised","N_Natural")])
  X_W_d <- cbind(dati_diag$lb_sc, dati_diag$Texture1, dati_diag$Texture2, dati_diag$bd_sc, dati_diag$ph_sc)
  j_idx_d <- dati_diag$field_idx

  mu_SOC_d <- alpha_SOC + V_mean_diag[1,j_idx_d] + V_mean_diag[2,j_idx_d]*dati_diag$lb_sc +
              as.numeric(X_W_d %*% gamma_SOC) + as.numeric(X_B_d %*% beta_SOC)[j_idx_d]
  mu_N_d   <- alpha_N   + V_mean_diag[3,j_idx_d] + V_mean_diag[4,j_idx_d]*dati_diag$lb_sc +
              as.numeric(X_W_d %*% gamma_N)   + as.numeric(X_B_d %*% beta_N)[j_idx_d]
  mu_P_d   <- alpha_P   + V_mean_diag[5,j_idx_d] + V_mean_diag[6,j_idx_d]*dati_diag$lb_sc +
              as.numeric(X_W_d %*% gamma_P)   + as.numeric(X_B_d %*% beta_P)[j_idx_d]

  df_res_diag <- bind_rows(
    data.frame(Field=dati_diag$Field, Bottom=dati_diag$Bottom, resid=log(dati_diag$PercSOC)-mu_SOC_d, response="logSOC"),
    data.frame(Field=dati_diag$Field, Bottom=dati_diag$Bottom, resid=log(dati_diag$PercTotNitro)-mu_N_d, response="logN"),
    data.frame(Field=dati_diag$Field, Bottom=dati_diag$Bottom, resid=log(dati_diag$PercTotPhos)-mu_P_d, response="logP")
  )

  cat("Lag-1 autocorrelazione residui (mediana su campi con n>=3):\n")
  for(resp in c("logSOC","logN","logP")){
    d <- df_res_diag |> filter(response == resp) |>
      group_by(Field) |> arrange(Bottom, .by_group=TRUE) |>
      filter(n() >= 3) |>
      summarise(ac = cor(resid[-n()], resid[-1]), .groups="drop")
    cat(sprintf("  %s: mediana=%.3f  [Q10=%.3f, Q90=%.3f]  n=%d\n",
                resp, median(d$ac,na.rm=T), quantile(d$ac,.10,na.rm=T),
                quantile(d$ac,.90,na.rm=T), sum(!is.na(d$ac))))
  }
  rm(fit_mvre_diag); gc()
} else {
  cat("  fit_msp_rirs_mvre.rds non trovato — skip diagnostica residui.\n")
}

cat("\n── Fine script 16 ──────────────────────────────────────────────\n")
