# =============================================================================
# 03b_model_msp_landuse.R  —  M-SP con Landuse al posto delle binarie
#
# STRUTTURA: identica a 03_model_proportional_slope.R.
#   L'unica differenza è X_B: invece delle 4 variabili binarie
#   (OnFarm, Irrigate, Fertilised, N_Natural), si usa Landuse come
#   predittore between-field con K−1 = 6 dummy variables.
#
# MOTIVAZIONE:
#   Le 4 binarie identificano Landuse quasi perfettamente (ogni livello di
#   Landuse corrisponde a una combinazione unica di binarie). Questo confronto
#   verifica se la reparametrizzazione cambia l'inferenza su eta_r e se
#   i beta_r per Landuse sono più interpretabili di quelli per le binarie.
#
# CATEGORIE LANDUSE (da crop.rds):
#   Landuse 1: OnFarm=1, Irrigate=0, Fertilised=1, N_Natural=1
#   Landuse 2: OnFarm=1, Irrigate=1, Fertilised=1, N_Natural=1
#   Landuse 3: OnFarm=0, Irrigate=0, Fertilised=0, N_Natural=1  ← RIFERIMENTO
#   Landuse 4: OnFarm=0, Irrigate=0, Fertilised=0, N_Natural=0
#   Landuse 5: OnFarm=0, Irrigate=1, Fertilised=1, N_Natural=1
#   Landuse 6: OnFarm=0, Irrigate=1, Fertilised=0, N_Natural=1
#   Landuse 7: OnFarm=0, Irrigate=0, Fertilised=1, N_Natural=1
#
#   Riferimento = Landuse 3: nessuna gestione antropica, solo azoto naturale.
#   È il più "naturale" del dataset (foresta/incolto).
#   I beta_r stimano la differenza di ogni altro landuse rispetto a Landuse 3.
#
# CONFRONTO LOO: M-SP-Landuse vs M-SP (binarie)
# Fit salvato in: stan/fit_msp_landuse.rds
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


# ── 1. DATI ───────────────────────────────────────────────────────────────────

dati <- readRDS(here("data", "dati.rds")) |>
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


# ── 1b. X_B DA LANDUSE ────────────────────────────────────────────────────────
# Landuse non è in dati.rds → recuperato da crop.rds.
# Riferimento: Landuse 3 (nessuna gestione, solo azoto naturale).
# K_B = 6 dummy (una per ogni livello ≠ 3).

crop_raw <- readRDS(here("data", "crop.rds"))

landuse_per_field <- crop_raw |>
  distinct(Field, Landuse) |>
  mutate(
    Field     = factor(as.integer(as.character(Field)), levels = field_levels),
    field_int = as.integer(Field),
    Landuse   = factor(Landuse)
  ) |>
  arrange(field_int)

REF_LEVEL          <- "3"
non_ref_levels     <- sort(setdiff(levels(landuse_per_field$Landuse), REF_LEVEL))
X_B_cols           <- paste0("Landuse", non_ref_levels)   # "Landuse1" ... "Landuse7"

X_B <- sapply(non_ref_levels, function(lv) {
  as.integer(as.character(landuse_per_field$Landuse) == lv)
})
colnames(X_B) <- X_B_cols

cat(sprintf("K_B = %d dummy Landuse (riferimento = Landuse %s)\n",
            length(X_B_cols), REF_LEVEL))
cat("Livelli non-riferimento:", paste(non_ref_levels, collapse = ", "), "\n")
cat("Distribuzione Landuse per campo:\n")
print(table(as.character(landuse_per_field$Landuse)))

rm(crop_raw)

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

n_extra <- 3 * length(X_W_cols) + 3 * length(X_B_cols) + 12  # = 3*4 + 3*6 + 12 = 42
cat(sprintf("Parametri totali attesi: 3×%d + %d = %d\n", J, n_extra, 3 * J + n_extra))
rm(X_W, X_B); gc()


# ── 2. COMPILAZIONE ───────────────────────────────────────────────────────────
# Stesso file Stan di M-SP: K_B viene passato come dato, non hardcodato.

stan_file <- here("stan", "m4rr_v2_ri_slope_mu.stan")
cat("\nCompilazione m4rr_v2_ri_slope_mu...\n")
mod <- cmdstan_model(stan_file, compile = TRUE)
cat("OK. CmdStan version:", cmdstan_version(), "\n")


# ── 3. MCMC ───────────────────────────────────────────────────────────────────

fit_path <- here("stan", "fit_msp_landuse.rds")

if (!file.exists(fit_path)) {

  cat("\nAvvio MCMC (4 catene × 5000 sampling + 2000 warmup)...\n\n")

  fit <- mod$sample(
    data            = stan_data,
    seed            = 2024,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 2000,
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

cat("\n═══ DIAGNOSTICA ══════════════════════════════════════════════════\n")

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
  "gamma_SOC", "gamma_N", "gamma_P",
  "beta_SOC",  "beta_N",  "beta_P"
)) |>
  select(variable, mean, median, sd, q5, q95, rhat, ess_bulk) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\n═══ PARAMETRI CHIAVE ═════════════════════════════════════════════\n")
print(smry, n = 70)


# ── 6. STRUTTURA RI–SLOPE + INTERCETTA GLOBALE ────────────────────────────────

cat("\n═══ INTERCETTA GLOBALE + RI–SLOPE PER RISPOSTA ══════════════════\n")
cat(sprintf("%-8s | %8s %8s | %8s %8s | %8s %8s | %8s %8s | %8s\n",
            "Risposta", "alpha", "sd", "psi", "sd", "eta", "sd", "b=eta/psi", "sd", "sigma"))
cat(strrep("-", 95), "\n")

for (r in c("SOC", "N", "P")) {
  a_r   <- smry |> filter(variable == paste0("alpha_", r))
  psi_r <- smry |> filter(variable == paste0("psi_",   r))
  eta_r <- smry |> filter(variable == paste0("eta_",   r))
  b_r   <- smry |> filter(variable == paste0("b_",     r))
  sig_r <- smry |> filter(variable == paste0("sigma_", r))
  cat(sprintf("%-8s | %8.3f %8.3f | %8.3f %8.3f | %8.3f %8.3f | %8.3f %8.3f | %8.3f\n",
              paste0("log", r),
              a_r$median,   a_r$sd,
              psi_r$median, psi_r$sd,
              eta_r$median, eta_r$sd,
              b_r$median,   b_r$sd,
              sig_r$median))
}

cat("\n  Confronto con M-SP (binarie): eta_r atteso stabile — il segnale\n")
cat("  di decadimento non dipende dalla reparametrizzazione di X_B.\n")


# ── 7. EFFETTI WITHIN E BETWEEN ───────────────────────────────────────────────

cov_W <- c("logBottom", "Texture1", "Texture2", "BulkDensity")
# Etichette leggibili per i dummy Landuse (vs riferimento Landuse 3)
cov_B_labels <- paste0("vs Landuse3: ", X_B_cols)

cat("\n═══ EFFETTI WITHIN-FIELD (gamma_r) ══════════════════════════════\n")
cat(sprintf("%-12s | %8s %8s | %8s %8s | %8s %8s\n",
            "Covariata", "SOC", "sd", "N", "sd", "P", "sd"))
cat(strrep("-", 66), "\n")
for (k in seq_along(cov_W)) {
  gs <- smry |> filter(variable == sprintf("gamma_SOC[%d]", k))
  gn <- smry |> filter(variable == sprintf("gamma_N[%d]",   k))
  gp <- smry |> filter(variable == sprintf("gamma_P[%d]",   k))
  cat(sprintf("%-12s | %8.3f %8.3f | %8.3f %8.3f | %8.3f %8.3f\n",
              cov_W[k], gs$median, gs$sd, gn$median, gn$sd, gp$median, gp$sd))
}

cat("\n═══ EFFETTI BETWEEN-FIELD (beta_r) — Landuse vs riferimento (3) ═\n")
cat(sprintf("%-20s | %8s %8s | %8s %8s | %8s %8s\n",
            "Landuse", "SOC", "sd", "N", "sd", "P", "sd"))
cat(strrep("-", 74), "\n")
for (k in seq_along(X_B_cols)) {
  bs <- smry |> filter(variable == sprintf("beta_SOC[%d]", k))
  bn <- smry |> filter(variable == sprintf("beta_N[%d]",   k))
  bp <- smry |> filter(variable == sprintf("beta_P[%d]",   k))
  cat(sprintf("%-20s | %8.3f %8.3f | %8.3f %8.3f | %8.3f %8.3f\n",
              cov_B_labels[k],
              bs$median, bs$sd, bn$median, bn$sd, bp$median, bp$sd))
}


# ── 8. LOO-CV ─────────────────────────────────────────────────────────────────

if (requireNamespace("loo", quietly = TRUE)) {
  library(loo)

  cat("\n═══ LOO-CV ══════════════════════════════════════════════════════\n")
  ll      <- fit$draws("log_lik", format = "matrix")
  loo_cur <- loo(ll, cores = 4)
  cat("\nM-SP-Landuse:\n"); print(loo_cur)

  p_msp <- here("stan", "fit_msp.rds")
  loos  <- list(`M-SP-Landuse` = loo_cur)

  if (file.exists(p_msp)) {
    cat("LOO per M-SP (binarie)...\n")
    f_msp    <- readRDS(p_msp)
    ll_msp   <- f_msp$draws("log_lik", format = "matrix")
    loos[["M-SP"]] <- loo(ll_msp, cores = 4)
    rm(f_msp, ll_msp); gc()
  }

  cat("\n═══ CONFRONTO LOO: Landuse vs binarie ═══════════════════════════\n")
  if (length(loos) > 1) print(loo_compare(loos))
  rm(ll, loos); gc()
}


# ── 9. TRACE PLOTS ────────────────────────────────────────────────────────────

if (requireNamespace("bayesplot", quietly = TRUE)) {
  params_key <- c("alpha_SOC", "psi_SOC", "eta_SOC", "b_SOC",
                  "alpha_N",   "psi_N",   "eta_N",
                  "alpha_P",   "psi_P",   "eta_P")
  draws_key <- fit$draws(variables = params_key)
  print(mcmc_trace(draws_key) + theme_minimal() +
          ggtitle("Trace plots — M-SP Landuse (eta_r: confronto con binarie)"))
  rm(draws_key); gc()
}

cat("\n── Fine script 03b ──────────────────────────────────────────────\n")
