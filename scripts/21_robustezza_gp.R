# =============================================================================
# 21_robustezza_gp.R  —  M-SP-RIRS-GP: GP esponenziale sugli intercetti
#
# Robustezza rispetto a M-SP-RIRS-MVRE: sostituisce il prior iid sugli
# intercetti di campo con un GP esponenziale per risposta (kernel Matérn 1/2).
# Le slope restano iid. Non esiste più la struttura 6D MVRE (nessuna Omega_6).
#
# Stan file: stan/model_mvre_gp.stan
# Fit file:  stan/fit_msp_rirs_mvre_gp.rds
#
# DOMANDE:
#   A. Gli iperparametri GP (alpha_gp, ell, nu) sono identificati?
#   B. Il LOO cambia rispetto a MVRE? (atteso: poco — il LOO obs-level non
#      premia la struttura spaziale se i campi hanno osservazioni residue)
#   C. Gli effetti di gestione beta_r cambiano? (se sì → spatial confounding)
#   D. rho_int_SOC_N "sopravvive" come correlazione residua dei GP posteriori?
#
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
if (requireNamespace("spdep",      quietly = TRUE)) library(spdep)


# ── 1. DATI ────────────────────────────────────────────────────────────────────

dati <- carica_dati()

field_levels <- sort(unique(as.integer(as.character(dati$Field))))
J <- length(field_levels)
N <- nrow(dati)

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))

X_W_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")
X_B_cols <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")


# ── 2. COORDINATE E MATRICE DISTANZE ──────────────────────────────────────────

crop_raw <- readRDS(here("data", "crop.rds"))
coords <- crop_raw |>
  group_by(Field) |>
  summarise(Lat = mean(Lat, na.rm = TRUE), Long = mean(Long, na.rm = TRUE)) |>
  mutate(Field = factor(Field, levels = levels(dati$Field))) |>
  arrange(Field)
rm(crop_raw); gc()

ref_lat  <- mean(coords$Lat)
ref_long <- mean(coords$Long)
coords <- coords |>
  mutate(
    x_km = (Long - ref_long) * 111.32 * cos(ref_lat * pi / 180),
    y_km = (Lat  - ref_lat)  * 110.54
  )

coords_mat <- as.matrix(coords[, c("x_km", "y_km")])
Dmat <- as.matrix(dist(coords_mat))   # J×J, km

cat(sprintf("Area: %.1f × %.1f km | dist mediana: %.1f km | max: %.1f km\n",
            diff(range(coords$x_km)), diff(range(coords$y_km)),
            median(Dmat[lower.tri(Dmat)]), max(Dmat)))


# ── 3. STAN DATA ───────────────────────────────────────────────────────────────

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
               arrange(field_int) |> select(all_of(X_B_cols)) |> as.matrix(),
  Dmat     = Dmat
)

cat(sprintf("N = %d | J = %d | K_W = %d | K_B = %d | Dmat dim = %dx%d\n",
            N, J, length(X_W_cols), length(X_B_cols),
            nrow(Dmat), ncol(Dmat)))


# ── 4. FIT ─────────────────────────────────────────────────────────────────────

fit_path  <- here("stan", "fit_msp_rirs_mvre_gp.rds")
stan_file <- here("stan", "model_mvre_gp.stan")

if (!file.exists(fit_path)) {
  mod <- cmdstan_model(stan_file, compile = TRUE)
  cat("Avvio MCMC GP (4 catene × 3000 sampling + 3000 warmup)...\n\n")
  fit <- mod$sample(
    data            = stan_data,
    seed            = 2024,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 3000,
    iter_sampling   = 3000,
    adapt_delta     = 0.98,
    max_treedepth   = 12,
    refresh         = 200,
    output_dir      = here("stan")
  )
  fit$save_object(fit_path)
  invisible(file.remove(c(
    list.files(here("stan"), pattern = "\\.csv$",  full.names = TRUE),
    list.files(here("stan"), pattern = "\\.json$", full.names = TRUE)
  )))
} else {
  cat("Carico fit GP da:", fit_path, "\n")
  fit <- readRDS(fit_path)
}


# ── 5. DIAGNOSTICA ─────────────────────────────────────────────────────────────

cat("\n═══ DIAGNOSTICA M-SP-RIRS-GP ══════════════════════════════════════\n")

smry <- fit$summary()
np   <- fit$sampler_diagnostics()

n_divs <- sum(np[, , "divergent__"])
max_td  <- max(np[, , "treedepth__"])
ebfmi   <- apply(np[, , "energy__"], 2, function(e) mean(diff(e)^2) / var(e))

bad_rhat <- smry |> filter(rhat > 1.05) |> arrange(desc(rhat))
bad_ess  <- smry |> filter(ess_bulk < 400) |> arrange(ess_bulk)

cat("Divergenze:", n_divs, "| Max treedepth:", max_td, "\n")
cat("E-BFMI per catena:", round(ebfmi, 3), "\n")
cat("Rhat > 1.05:", nrow(bad_rhat))
if (nrow(bad_rhat) > 0) { cat("\n"); print(head(bad_rhat, 10)) } else cat(" — OK\n")
cat("ESS bulk < 400:", nrow(bad_ess))
if (nrow(bad_ess)  > 0) { cat("\n"); print(head(bad_ess, 10)) } else cat(" — OK\n")


# ── 6. IPERPARAMETRI GP ────────────────────────────────────────────────────────

cat("\n═══ IPERPARAMETRI GP ══════════════════════════════════════════════\n")

gp_params <- c(
  "alpha_gp_SOC", "ell_SOC", "nu_SOC",
  "alpha_gp_N",   "ell_N",   "nu_N",
  "alpha_gp_P",   "ell_P",   "nu_P"
)

gp_smry <- smry |>
  filter(variable %in% gp_params) |>
  select(variable, median, sd, q5, q95) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat(sprintf("%-16s | %8s %8s | %8s %8s\n", "Parametro", "Mediana", "SD", "CI5", "CI95"))
cat(strrep("-", 58), "\n")
for (r in c("SOC","N","P")) {
  for (p in c("alpha_gp","ell","nu")) {
    nm <- paste0(p, "_", r)
    row <- gp_smry |> filter(variable == nm)
    unit <- if_else(p == "ell", "(km)", "     ")
    cat(sprintf("%-16s | %8.3f %8.3f | %8.3f %8.3f  %s\n",
                nm, row$median, row$sd, row$q5, row$q95, unit))
  }
  cat("\n")
}

cat("ell: range dove la correlazione scende a exp(-1) ≈ 0.37\n")
cat("     range effettivo (95%): ca. 3 × ell\n")


# ── 7. CONFRONTO SLOPE e SIGMA ─────────────────────────────────────────────────

cat("\n═══ SLOPE e SIGMA (confronto con MVRE) ═══════════════════════════\n")

slope_params <- c(
  "tau_slope_SOC", "sigma_SOC",
  "tau_slope_N",   "sigma_N",
  "tau_slope_P",   "sigma_P"
)

sp_smry <- smry |>
  filter(variable %in% slope_params) |>
  select(variable, median, q5, q95) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("(GP: tau_slope_r iid; MVRE: tau_beta_r dal fit precedente)\n\n")
cat(sprintf("%-16s | %8s | %8s %8s\n", "Parametro", "Mediana", "CI5", "CI95"))
cat(strrep("-", 48), "\n")
for (i in seq_len(nrow(sp_smry)))
  cat(sprintf("%-16s | %8.3f | %8.3f %8.3f\n",
              sp_smry$variable[i], sp_smry$median[i],
              sp_smry$q5[i], sp_smry$q95[i]))


# ── 8. EFFETTI DI GESTIONE beta_r ─────────────────────────────────────────────

cat("\n═══ EFFETTI GESTIONE beta_r ═══════════════════════════════════════\n")
cat("(Domanda chiave: cambiano rispetto a MVRE dopo correzione spaziale?)\n\n")

xb_labels <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")

for (r in c("SOC", "N", "P")) {
  cat(sprintf("── log%s ──\n", r))
  for (k in seq_along(xb_labels)) {
    nm  <- sprintf("beta_%s[%d]", r, k)
    row <- smry |> filter(variable == nm)
    if (nrow(row) == 0) next
    cat(sprintf("  %-14s  median = %+.3f  [%+.3f, %+.3f]\n",
                xb_labels[k], row$median, row$q5, row$q95))
  }
  cat("\n")
}


# ── 9. CORRELAZIONE RESIDUA GP POSTERIORI (ex post rho_int_SOC_N) ─────────────

cat("\n═══ CORRELAZIONE RESIDUA INTERCETTI GP (ex post) ══════════════════\n")
cat("rho calcolata sui draws posteriori di u_int_r — analogo di rho_int MVRE\n\n")

draws_int <- fit$draws(
  variables = c(paste0("u_int_SOC[", 1:J, "]"),
                paste0("u_int_N[",   1:J, "]"),
                paste0("u_int_P[",   1:J, "]")),
  format = "matrix"
)

# Per ogni draw, calcola la correlazione tra i J valori dei tre vettori
# (correlazione cross-risposta degli intercetti di campo)
set.seed(42)
samp_idx <- sample(nrow(draws_int), min(2000, nrow(draws_int)))

rho_soc_n <- rho_soc_p <- rho_n_p <- numeric(length(samp_idx))
for (s in seq_along(samp_idx)) {
  i       <- samp_idx[s]
  int_soc <- draws_int[i, paste0("u_int_SOC[", 1:J, "]")]
  int_n   <- draws_int[i, paste0("u_int_N[",   1:J, "]")]
  int_p   <- draws_int[i, paste0("u_int_P[",   1:J, "]")]
  rho_soc_n[s] <- cor(int_soc, int_n)
  rho_soc_p[s] <- cor(int_soc, int_p)
  rho_n_p[s]   <- cor(int_n,   int_p)
}

q_ci <- function(x) quantile(x, c(0.05, 0.50, 0.95))

cat(sprintf("%-22s | %8s | %8s %8s\n", "Correlazione", "Mediana", "CI5", "CI95"))
cat(strrep("-", 54), "\n")
for (nm_pair in list(
  list("rho_int_SOC_N (GP)", rho_soc_n),
  list("rho_int_SOC_P (GP)", rho_soc_p),
  list("rho_int_N_P   (GP)", rho_n_p)
)) {
  ci <- q_ci(nm_pair[[2]])
  cat(sprintf("%-22s | %+8.3f | %+8.3f %+8.3f\n",
              nm_pair[[1]], ci[2], ci[1], ci[3]))
}
cat("\n(MVRE riferimento: rho_int_SOC_N=+0.386 [0.038,0.667])\n")
cat("  Se la correlazione residua cade → era struttura spaziale condivisa.\n")
cat("  Se persiste → segnale biologico reale, non artefatto spaziale.\n")


# ── 10. MORAN SUI RESIDUI GP ──────────────────────────────────────────────────

cat("\n═══ MORAN SU INTERCETTI GP POSTERIORI ═════════════════════════════\n")
cat("(Atteso: vicino a zero dopo GP — verifica che il GP assorba la struttura)\n\n")

D_geo  <- as.matrix(dist(coords_mat))
W_geo  <- ifelse(D_geo > 0.3, 1 / D_geo, 0)
diag(W_geo) <- 0
W_geo  <- W_geo / rowSums(W_geo)
lw_geo <- mat2listw(W_geo, style = "W")

u_mean_soc <- colMeans(draws_int[, paste0("u_int_SOC[", 1:J, "]")])
u_mean_n   <- colMeans(draws_int[, paste0("u_int_N[",   1:J, "]")])
u_mean_p   <- colMeans(draws_int[, paste0("u_int_P[",   1:J, "]")])

for (lst in list(
  list("u_int_SOC (GP)", u_mean_soc),
  list("u_int_N   (GP)", u_mean_n),
  list("u_int_P   (GP)", u_mean_p)
)) {
  mt <- moran.test(lst[[2]], lw_geo)
  z  <- (mt$estimate["Moran I statistic"] - mt$estimate["Expectation"]) /
        sqrt(mt$estimate["Variance"])
  cat(sprintf("%-18s  I = %+.4f  Z = %+.2f  p = %.4f\n",
              lst[[1]], mt$estimate["Moran I statistic"], z, mt$p.value))
}

cat("\n(MVRE iid: I ≈ 0.91, Z ≈ 22 — confronta con valori sopra)\n")
rm(draws_int); gc()


# ── 11. LOO-CV: GP vs MVRE ────────────────────────────────────────────────────

cat("\n═══ LOO-CV: GP vs M-SP-RIRS-MVRE ══════════════════════════════════\n")

cat("LOO M-SP-RIRS-GP...\n")
loo_gp <- loo(fit$draws("log_lik", format = "matrix"), cores = 4)
cat("\nM-SP-RIRS-GP:\n"); print(loo_gp)

mvre_path <- here("stan", "fit_msp_rirs_mvre.rds")
if (file.exists(mvre_path)) {
  cat("\nLOO M-SP-RIRS-MVRE...\n")
  fit_mvre <- readRDS(mvre_path)
  loo_mvre <- loo(fit_mvre$draws("log_lik", format = "matrix"), cores = 4)
  rm(fit_mvre); gc()

  cat("\n── Confronto LOO ─────────────────────────────────────────────────\n")
  cmp <- loo_compare(list(`GP` = loo_gp, `MVRE` = loo_mvre))
  print(cmp)

  delta <- loo_gp$estimates["elpd_loo","Estimate"] -
           loo_mvre$estimates["elpd_loo","Estimate"]
  cat(sprintf("\nΔELPD GP vs MVRE = %+.1f\n", delta))
  cat("  > 0 → GP ha migliore LOO obs-level\n")
  cat("  < 0 → MVRE ha migliore LOO obs-level\n")
  cat("  Nota: LOO obs-level non premia molto la struttura spaziale.\n")
  cat("  Un LOO leave-one-field-out sarebbe più discriminante.\n")
}

cat("\n── Fine script 21 ───────────────────────────────────────────────\n")
