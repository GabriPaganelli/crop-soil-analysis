# =============================================================================
# 09_model_msp_rirs.R  —  M-SP-RIRS: random intercept + random slope bivariato
#                         ← MODELLO FINALE
#
# Per ogni risposta r e campo j, intercetta e slope su logBottom sono estratti
# congiuntamente da una distribuzione bivariata:
#   [u_int_r[j], u_slope_r[j]] ~ MVN(0, Sigma_r)
#   Sigma_r = diag(tau_r) * Omega_r * diag(tau_r),  Omega_r ~ LKJ(2)
#
# Struttura:
#   mu_r[i,j] = alpha_r
#             + u_r[1, j]                        ← random intercept
#             + u_r[2, j] * logBottom_i          ← random slope (deviazione dal fisso)
#             + gamma_r * X_W_i                  ← effetti fissi within (incl. logBottom)
#             + beta_r  * X_B_j
#
# X_W (K_W=5): logBottom, Texture1(ILR1), Texture2(ILR2), BulkDensity, PH
# X_B (K_B=4): OnFarm, Irrigate, Fertilised, N_Natural
#
# rho_r = Omega_r[1,2]: correlazione tra RI e RS per risposta r.
#   rho_r ≈ +1  →  M-SP confermato (proporzionalità regge)
#   |rho_r| < 1 →  alcuni campi escono dalla proporzionalità
#   tau_beta_r ≈ 0  →  M-RI (nessuna variazione nelle pendenze)
#
# Stan: stan/m_sp_rirs.stan
# Fit:  stan/fit_msp_rirs.rds
# Parametri totali: 6J + 42 = 282  (J=40, K_W=5)
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
  K_W      = length(X_W_cols),   # 5: logBottom, Texture1, Texture2, BulkDensity, PH
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
cat(sprintf("Parametri totali attesi: 6×%d + 42 = %d\n", J, 6 * J + 42))
rm(dati, dati_int); gc()


# ── 2. COMPILAZIONE E FIT ─────────────────────────────────────────────────────

fit_path  <- here("stan", "fit_msp_rirs.rds")
stan_file <- here("stan", "m_sp_rirs.stan")

mod <- cmdstan_model(stan_file, compile = TRUE)
cat("Compilazione m_sp_rirs: OK\n")

if (!file.exists(fit_path)) {
  cat("Avvio MCMC M-SP-RIRS (4 catene × 5000 sampling + 3000 warmup)...\n\n")
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
  invisible(file.remove(c(
    list.files(here("stan"), pattern = "\\.csv$",  full.names = TRUE),
    list.files(here("stan"), pattern = "\\.json$", full.names = TRUE)
  )))
} else {
  cat("Carico fit da:", fit_path, "\n")
  fit <- readRDS(fit_path)
}

rm(stan_data); gc()


# ── 3. DIAGNOSTICA ────────────────────────────────────────────────────────────

cat("\n═══ DIAGNOSTICA M-SP-RIRS ════════════════════════════════════════\n")

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
  "alpha_P",   "tau_alpha_P",   "tau_beta_P",   "rho_P",   "sigma_P"
)

smry_key <- fit$summary(params_chiave) |>
  select(variable, median, sd, q5, q95, rhat, ess_bulk) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\n═══ STRUTTURA BIVARIATA RI+RS ════════════════════════════════════\n")
cat(sprintf("%-8s | %10s %8s | %10s %8s | %8s %8s | %8s\n",
            "Risp.", "tau_alpha", "sd", "tau_beta", "sd", "rho", "sd", "sigma"))
cat(strrep("-", 80), "\n")
for (r in c("SOC", "N", "P")) {
  ta  <- smry_key |> filter(variable == paste0("tau_alpha_", r))
  tb  <- smry_key |> filter(variable == paste0("tau_beta_",  r))
  rho <- smry_key |> filter(variable == paste0("rho_",       r))
  sg  <- smry_key |> filter(variable == paste0("sigma_",     r))
  cat(sprintf("%-8s | %10.3f %8.3f | %10.3f %8.3f | %8.3f %8.3f | %8.3f\n",
              paste0("log", r),
              ta$median, ta$sd, tb$median, tb$sd,
              rho$median, rho$sd, sg$median))
}

cat("\n  Interpretazione rho_r:\n")
cat("    rho ≈ +1 → M-SP confermato: i campi con alto baseline decadono più lentamente\n")
cat("    |rho| < 1 → alcuni campi escono dalla proporzionalità\n")
cat("    rho ≈  0 → RI e RS indipendenti\n")

cat("\n  Confronto con M-SP (tau_alpha ≈ psi, tau_beta ≈ |eta|):\n")
cat("    M-SP:  psi_SOC = 0.330, eta_SOC = +0.204  → atteso rho_SOC ≈ +1, tau_beta_SOC ≈ 0.204\n")
cat("    M-SP:  psi_N   = 0.203, eta_N   = -0.038  → atteso rho_N   ≈ 0,  tau_beta_N   ≈ 0\n")
cat("    M-SP:  psi_P   = 0.400, eta_P   = +0.097  → atteso rho_P   ≈ +1, tau_beta_P   ≈ 0.097\n")

cat("\n═══ EFFETTI FISSI WITHIN-FIELD (gamma_r) ═════════════════════════\n")
cov_W <- c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")
cat(sprintf("%-12s | %8s %8s | %8s %8s | %8s %8s\n",
            "Covariata", "SOC", "sd", "N", "sd", "P", "sd"))
cat(strrep("-", 66), "\n")
for (k in seq_along(cov_W)) {
  gs <- smry |> filter(variable == sprintf("gamma_SOC[%d]", k)) |>
        select(median, sd) |> mutate(across(everything(), ~ round(.x, 3)))
  gn <- smry |> filter(variable == sprintf("gamma_N[%d]",   k)) |>
        select(median, sd) |> mutate(across(everything(), ~ round(.x, 3)))
  gp <- smry |> filter(variable == sprintf("gamma_P[%d]",   k)) |>
        select(median, sd) |> mutate(across(everything(), ~ round(.x, 3)))
  cat(sprintf("%-12s | %8.3f %8.3f | %8.3f %8.3f | %8.3f %8.3f\n",
              cov_W[k], gs$median, gs$sd, gn$median, gn$sd, gp$median, gp$sd))
}

rm(smry); gc()


# ── 5. LOO-CV vs M-SP ─────────────────────────────────────────────────────────

cat("\n═══ LOO-CV: M-SP-RIRS vs M-SP ═══════════════════════════════════\n")

cat("LOO M-SP-RIRS...\n")
loo_rirs <- loo(fit$draws("log_lik", format = "matrix"), cores = 4)
cat("\nM-SP-RIRS:\n");  print(loo_rirs)
cat(sprintf("  Pareto k > 0.7: %d osservazioni\n", sum(loo_rirs$diagnostics$pareto_k > 0.7)))

msp_path <- here("stan", "fit_msp.rds")
if (file.exists(msp_path)) {
  cat("\nLOO M-SP...\n")
  fit_msp <- readRDS(msp_path)
  loo_msp <- loo(fit_msp$draws("log_lik", format = "matrix"), cores = 4)
  rm(fit_msp); gc()

  cat("\nM-SP:\n"); print(loo_msp)
  cat("\nConfronto LOO:\n")
  cmp <- loo_compare(list(`M-SP` = loo_msp, `M-SP-RIRS` = loo_rirs))
  print(cmp)

  delta   <- loo_rirs$estimates["elpd_loo","Estimate"] - loo_msp$estimates["elpd_loo","Estimate"]
  se_diff <- cmp[nrow(cmp), "se_diff"]
  cat(sprintf("\n  ΔELPD = %.1f  (SE = %.1f,  |z| = %.1f)\n", delta, se_diff, abs(delta/se_diff)))
} else {
  cat("\nfit_msp.rds non trovato — confronto LOO vs M-SP saltato.\n")
  cat("Esegui 08_model_msp.R e poi ricompila questo LOO.\n")
}
rm(fit); gc()

cat("\n── Fine script 09 ───────────────────────────────────────────────\n")
