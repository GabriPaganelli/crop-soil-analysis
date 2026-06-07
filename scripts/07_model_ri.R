# =============================================================================
# 07_model_ri.R  —  Modello M-RI (Random Intercept, confronto baseline)
#
# STRUTTURA (stan/model_ri.stan):
#   Per r ∈ {SOC, N, P}:
#     mu_r[i,j] = alpha_r              ← intercetta globale
#               + psi_r · z_nu_r[j]   ← random intercept i.i.d. (NCP)
#               + gamma_r · X_W_i     ← effetti fissi within-field
#               + beta_r  · X_B_j     ← effetti fissi between-field
#
# RUOLO NELLA PIPELINE:
#   Modello di confronto per M-SP (script 03). Cattura la variabilità
#   between-field con intercette casuali i.i.d., senza modellare la
#   dipendenza del decadimento verticale dal livello del campo.
#   alpha_r assorbe la media globale delle log-risposte (essenziale:
#   senza di essa i beta esplodono a ±2).
#
# CONFRONTO LOO: M-RI vs M-SP vs M-GPS vs M-GP
# Fit salvato in: stan/fit_mri.rds
# =============================================================================


# ── 0. SETUP ──────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
source(here("scripts", "00_utilities.R"))

if (!requireNamespace("cmdstanr", quietly = TRUE))
  install.packages("cmdstanr",
                   repos = c("https://mc-stan.org/r-packages/",
                             "https://cloud.r-project.org"))
library(cmdstanr)
setup_rtools()

if (requireNamespace("posterior",  quietly = TRUE)) library(posterior)
if (requireNamespace("bayesplot",  quietly = TRUE)) library(bayesplot)


# ── 1. DATI ───────────────────────────────────────────────────────────────────

dati <- carica_dati()

field_levels <- sort(unique(as.integer(as.character(dati$Field))))
J <- length(field_levels)
N <- nrow(dati)

cat(sprintf("N = %d | J = %d\n", N, J))

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))

X_W_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")
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
  K_W      = length(X_W_cols),   # 5: logBottom, Texture1, Texture2, BulkDensity, PH
  K_B      = length(X_B_cols),
  field_id = dati_int$field_int,
  logSOC   = dati_int$logSOC,
  logN     = dati_int$logN,
  logP     = dati_int$logP,
  X_W      = X_W,
  X_B      = X_B
)

cat(sprintf("Parametri totali attesi: 3×%d + 36 = %d\n", J, 3 * J + 36))
rm(X_W, X_B); gc()


# ── 2. COMPILAZIONE ───────────────────────────────────────────────────────────

stan_file <- here("stan", "model_ri.stan")
cat("\nCompilazione model_ri...\n")
mod <- cmdstan_model(stan_file, compile = TRUE)
cat("OK. CmdStan version:", cmdstan_version(), "\n")


# ── 3. MCMC ───────────────────────────────────────────────────────────────────

fit_path <- here("stan", "fit_mri.rds")

if (!file.exists(fit_path)) {

  cat("\nAvvio MCMC (4 catene × 5000 sampling + 3000 warmup)...\n\n")

  fit <- mod$sample(
    data            = stan_data,
    seed            = 2024,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 3000,
    iter_sampling   = 5000,
    adapt_delta     = 0.95,
    max_treedepth   = 10,
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
  "alpha_SOC", "psi_SOC", "sigma_SOC",
  "alpha_N",   "psi_N",   "sigma_N",
  "alpha_P",   "psi_P",   "sigma_P",
  "gamma_SOC", "gamma_N", "gamma_P",
  "beta_SOC",  "beta_N",  "beta_P"
)) |>
  select(variable, mean, median, sd, q5, q95, rhat, ess_bulk) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\n═══ PARAMETRI CHIAVE ═════════════════════════════════════════════\n")
print(smry, n = 50)


# ── 6. STRUTTURA GLOBALE + RANDOM PER RISPOSTA ────────────────────────────────

cat("\n═══ INTERCETTE GLOBALI E RANDOM INTERCEPT ════════════════════════\n")
cat(sprintf("%-8s | %8s %8s | %8s %8s | %8s\n",
            "Risposta", "alpha", "sd", "psi", "sd", "sigma"))
cat(strrep("-", 55), "\n")

for (r in c("SOC", "N", "P")) {
  a_r   <- smry |> filter(variable == paste0("alpha_", r))
  psi_r <- smry |> filter(variable == paste0("psi_",   r))
  sig_r <- smry |> filter(variable == paste0("sigma_", r))
  cat(sprintf("%-8s | %8.3f %8.3f | %8.3f %8.3f | %8.3f\n",
              paste0("log", r),
              a_r$median,   a_r$sd,
              psi_r$median, psi_r$sd,
              sig_r$median))
}

cat("\n  Confronto con v2_no_gp: alpha assorbe la media globale,\n")
cat("  i coefficienti management dovrebbero tornare < |0.5|.\n")


# ── 7. EFFETTI WITHIN E BETWEEN ───────────────────────────────────────────────

cov_W <- c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")
cov_B <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")

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

cat("\n═══ EFFETTI BETWEEN-FIELD (beta_r) ══════════════════════════════\n")
cat(sprintf("%-12s | %8s %8s | %8s %8s | %8s %8s\n",
            "Covariata", "SOC", "sd", "N", "sd", "P", "sd"))
cat(strrep("-", 66), "\n")
for (k in seq_along(cov_B)) {
  bs <- smry |> filter(variable == sprintf("beta_SOC[%d]", k))
  bn <- smry |> filter(variable == sprintf("beta_N[%d]",   k))
  bp <- smry |> filter(variable == sprintf("beta_P[%d]",   k))
  cat(sprintf("%-12s | %8.3f %8.3f | %8.3f %8.3f | %8.3f %8.3f\n",
              cov_B[k], bs$median, bs$sd, bn$median, bn$sd, bp$median, bp$sd))
}


# ── 8. LOO-CV ─────────────────────────────────────────────────────────────────

if (requireNamespace("loo", quietly = TRUE)) {
  library(loo)

  cat("\n═══ LOO-CV ══════════════════════════════════════════════════════\n")
  ll      <- fit$draws("log_lik", format = "matrix")
  loo_cur <- loo(ll, cores = 4)
  cat("\nv2_no_gp_mu:\n"); print(loo_cur)

  modelli <- list(
    `M-SP`  = here("stan", "fit_msp.rds"),
    `M-GPS` = here("stan", "fit_mgps.rds"),
    `M-GP`  = here("stan", "fit_mgp.rds")
  )

  loos <- list(`M-RI` = loo_cur)
  for (nm in names(modelli)) {
    p <- modelli[[nm]]
    if (file.exists(p)) {
      cat(sprintf("LOO per %s...\n", nm))
      f    <- readRDS(p)
      ll_f <- f$draws("log_lik", format = "matrix")
      loos[[nm]] <- loo(ll_f, cores = 4)
      rm(f, ll_f); gc()
    }
  }

  cat("\n═══ CONFRONTO LOO ═══════════════════════════════════════════════\n")
  if (length(loos) > 1) print(loo_compare(loos))
  rm(ll, loos); gc()
}


# ── 9. TRACE PLOTS ────────────────────────────────────────────────────────────

if (requireNamespace("bayesplot", quietly = TRUE)) {
  params_key <- c("alpha_SOC", "psi_SOC",
                  "alpha_N",   "psi_N",
                  "alpha_P",   "psi_P")
  draws_key <- fit$draws(variables = params_key)
  print(mcmc_trace(draws_key) + theme_minimal() +
          ggtitle("Trace plots — M-RI (intercette globali + random intercept)"))
  rm(draws_key); gc()
}

cat("\n── Fine script 07 ───────────────────────────────────────────────\n")
