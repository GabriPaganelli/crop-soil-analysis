# =============================================================================
# 11. gp9_no_latent_within.R
#
# Estensione di gp8: rimozione del fattore latente within-field.
#
# MOTIVAZIONE (da analisi gp8):
#   psi_W_org = 0.024 spiega solo l'8% della varianza within di logSOC e logN.
#   La correlazione residua within-field SOC-N e' ≈ 0.08 (confermata anche dal
#   GAM, script 8). Il vincolo lambda_N = 0.636 (fisso da M4rr) forza una
#   proporzionalita' rigida tra gli effetti di X_W su logSOC e logN non testata.
#
# STRUTTURA DI gp9:
#   - gamma_SOC[K_W] e gamma_N[K_W] stimati indipendentemente (no factor within)
#   - Correlazione within (SOC, N) libera via Cholesky (prior LKJ(2))
#   - Struttura between-field invariata da gp8:
#       eta_org_B ~ GP Matern3/2 + nugget (SOC)
#       mu_SOC = eta_org_B[j]*(1 + b_slope*X_W[n,1]) + alpha_SOC_res[j] + gamma_SOC'*X_W
#       mu_N   = alpha_N[j] + gamma_N'*X_W
#       logP   ~ N(eta_P_B[j], sqrt(psi_W_P))
#   - b_slope caricato da fit_lmm_b_slope.rds (già stimato in script 10)
#
# DOMANDE:
#   1. Il vincolo lambda_N era utile? (DELTA_ELPD gp8 vs gp9 risponde)
#   2. La correlazione within (SOC, N) libera: quanto differisce da zero?
#      (se rho_W_SN ≈ 0.08, il fattore within era davvero inerte)
#   3. Gli effetti di X_W differiscono tra SOC e N?
#      (confronto gamma_SOC[1] vs gamma_N[1] per logBottom, ecc.)
# =============================================================================


# ── 0. SETUP ─────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)

if (!requireNamespace("cmdstanr", quietly = TRUE))
  install.packages("cmdstanr",
                   repos = c("https://mc-stan.org/r-packages/",
                             "https://cloud.r-project.org"))
library(cmdstanr)

if (.Platform$OS.type == "windows") {
  rtools_path <- Sys.getenv("RTOOLS45_HOME", unset = "C:/rtools45")
  Sys.setenv(
    PATH = paste(file.path(rtools_path, "ucrt64/bin"),
                 file.path(rtools_path, "usr/bin"),
                 Sys.getenv("PATH"), sep = ";"),
    RTOOLS44_HOME = rtools_path
  )
}

if (requireNamespace("posterior",  quietly = TRUE)) library(posterior)
if (requireNamespace("bayesplot",  quietly = TRUE)) library(bayesplot)


# ── 1. DATI ──────────────────────────────────────────────────────────────────

dati <- readRDS(here("data", "dati.rds")) |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom)
  ) |>
  mutate(across(c(logBottom, PH, Texture1, Texture2, BulkDensity),
                ~ c(scale(.x)))) |>
  mutate(Field = factor(Field))

field_levels <- sort(unique(as.integer(as.character(dati$Field))))
J <- length(field_levels)
N <- nrow(dati)

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

crop_raw     <- readRDS(here("data", "crop.rds"))
coords_field <- crop_raw |>
  group_by(Field) |>
  summarise(Lat  = mean(Lat,  na.rm = TRUE),
            Long = mean(Long, na.rm = TRUE)) |>
  arrange(Field)
rm(crop_raw); gc()

ref_lat  <- mean(coords_field$Lat,  na.rm = TRUE)
ref_long <- mean(coords_field$Long, na.rm = TRUE)

coords_km <- coords_field |>
  arrange(Field) |>
  mutate(
    x_km = (Long - ref_long) * 111.32 * cos(ref_lat * pi / 180),
    y_km = (Lat  - ref_lat)  * 110.54
  )

dist_mat <- as.matrix(dist(coords_km[, c("x_km", "y_km")]))


# ── 2. CARICA b_slope DA gp8 (nessuna re-stima brms) ────────────────────────

brms_path <- here("stan", "fit_lmm_b_slope.rds")
if (!file.exists(brms_path))
  stop("fit_lmm_b_slope.rds non trovato: esegui prima 10. gp8_random_slope.R")

fit_lmm <- readRDS(brms_path)

b_hat <- tryCatch({
  if (!requireNamespace("brms", quietly = TRUE))
    stop("brms non installato")
  library(brms)
  draws_lmm <- as_draws_df(fit_lmm)
  sigma_int  <- draws_lmm[["sd_Field__Intercept"]]
  sigma_slp  <- draws_lmm[["sd_Field__logBottom"]]
  rho        <- draws_lmm[["cor_Field__Intercept__logBottom"]]
  median(rho * sigma_slp / sigma_int)
}, error = function(e) {
  message("Impossibile ricalcolare b da brms (", conditionMessage(e),
          "): uso b_hat = 0.4674 (da run gp8)")
  0.4674
})

cat(sprintf("b_slope caricato da gp8: %.4f\n", b_hat))


# ── 3. DATI STAN (senza lambda_N) ────────────────────────────────────────────

stan_data_gp9 <- list(
  N        = N,
  J        = J,
  K_W      = length(X_W_cols),
  K_B      = length(X_B_cols),
  field_id = dati_int$field_int,
  logSOC   = dati_int$logSOC,
  logN     = dati_int$logN,
  logP     = dati_int$logP,
  X_W      = X_W,
  X_B      = X_B,
  dist_mat = dist_mat,
  b_slope  = b_hat
  # lambda_N rimosso: gamma_SOC e gamma_N sono ora indipendenti
)

rm(X_W, X_B, dist_mat); gc()


# ── 4. COMPILAZIONE ──────────────────────────────────────────────────────────

stan_file_gp9 <- here("stan", "m4rr_gp9.stan")
cat("\nCompilazione gp9...\n")
mod_gp9 <- cmdstan_model(stan_file_gp9, compile = TRUE)
cat("OK. CmdStan version:", cmdstan_version(), "\n")


# ── 5. MCMC ──────────────────────────────────────────────────────────────────

fit_path_gp9      <- here("stan", "fit_m4rr_gp9.rds")
risultati_path_gp9 <- here("stan", "risultati_gp9.RData")

# Il modello Stan è cambiato (rho fissato, prior gamma_N): cancella fit precedenti
# per forzare ri-esecuzione. Rimuovi queste righe dopo il run definitivo.
if (file.exists(fit_path_gp9))      { file.remove(fit_path_gp9);      cat("Vecchio fit_gp9 rimosso\n") }
if (file.exists(risultati_path_gp9)) { file.remove(risultati_path_gp9); cat("Vecchi risultati_gp9 rimossi\n") }

if (!file.exists(fit_path_gp9)) {

  cat("\nAvvio MCMC gp9 (4 catene × 5000 sampling + 3000 warmup)...\n")

  fit_gp9 <- mod_gp9$sample(
    data            = stan_data_gp9,
    seed            = 2024,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 3000,
    iter_sampling   = 5000,
    adapt_delta     = 0.97,
    max_treedepth   = 12,
    refresh         = 200,
    output_dir      = here("stan")
  )

  fit_gp9$save_object(fit_path_gp9)
  cat("Fit salvato in:", fit_path_gp9, "\n")

  csv_files  <- list.files(here("stan"), pattern = "\\.csv$",  full.names = TRUE)
  json_files <- list.files(here("stan"), pattern = "\\.json$", full.names = TRUE)
  invisible(file.remove(c(csv_files, json_files)))

} else {
  cat("\nCarico gp9 da:", fit_path_gp9, "\n")
  fit_gp9 <- readRDS(fit_path_gp9)
}

rm(stan_data_gp9); gc()


# ── 6. DIAGNOSTICA ───────────────────────────────────────────────────────────

if (!file.exists(risultati_path_gp9)) {

  smry_full <- fit_gp9$summary()
  np        <- fit_gp9$sampler_diagnostics()

  diag_list <- list(
    bad_rhat = smry_full |> filter(rhat > 1.05) |> arrange(desc(rhat)),
    bad_ess  = smry_full |> filter(ess_bulk < 400) |> arrange(ess_bulk),
    n_divs   = sum(np[, , "divergent__"]),
    max_td   = max(np[, , "treedepth__"]),
    ebfmi    = apply(np[, , "energy__"], 2,
                     function(e) mean(diff(e)^2) / var(e))
  )

  cat("\n── Rhat > 1.05 ──\n");    print(diag_list$bad_rhat)
  cat("\n── ESS bulk < 400 ──\n"); print(head(diag_list$bad_ess, 20))
  cat("Divergenze:", diag_list$n_divs,
      "| Max treedepth:", diag_list$max_td, "\n")
  cat("E-BFMI per catena:", round(diag_list$ebfmi, 3), "\n")

  rm(smry_full, np); gc()


  # ── 7. SUMMARY PARAMETRI ─────────────────────────────────────────────────────

  # rho_W_SN_out rimosso: correlazione within fissata a 0 in questa versione
  # (gp9a aveva stimato rho = 0.041, CI90 [-0.085, 0.165] — non diverso da 0)
  vars_key <- c(
    "gamma_SOC", "gamma_N",
    "gamma_P_B", "gamma_N_B",
    "sigma_W_SOC", "sigma_W_N", "psi_W_P",
    "sigma_B_SOC", "sigma_B_N",
    "psi_B_org", "sigma_P_between",
    "sigma_GP_org", "rho_org"
  )

  smry_gp9 <- fit_gp9$summary(vars_key) |>
    select(variable, mean, median, sd, q5, q95, rhat, ess_bulk) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  cat("\n── Parametri gp9 ──\n")
  print(smry_gp9, n = 40)


  # ── 8. CONFRONTO gamma_SOC vs gamma_N ────────────────────────────────────────
  # Risponde alla domanda: il vincolo lambda_N era giustificato?
  # Se gamma_SOC[1] ≈ 0.636 * gamma_N[1] (e idem per gli altri), il fattore era ok.
  # Se differiscono, il fattore imponeva una struttura non supportata dai dati.

  cat("\n── Confronto gamma_SOC vs gamma_N (effetti within-field) ──────────────\n")
  cat(sprintf("%-15s | %8s %8s | %8s %8s | %8s\n",
              "Predittore", "gamma_SOC", "sd", "gamma_N", "sd", "ratio"))
  cat(strrep("-", 65), "\n")

  gSOC <- smry_gp9 |> filter(startsWith(variable, "gamma_SOC"))
  gN   <- smry_gp9 |> filter(startsWith(variable, "gamma_N["))

  for (i in seq_len(nrow(gSOC))) {
    pred_lab <- X_W_cols[i]
    gs <- gSOC$mean[i]; ss <- gSOC$sd[i]
    gn <- gN$mean[i];   sn <- gN$sd[i]
    ratio <- if (abs(gn) > 0.01) gs / gn else NA_real_
    cat(sprintf("%-15s | %8.3f %8.3f | %8.3f %8.3f | %8.3f\n",
                pred_lab, gs, ss, gn, sn, ratio))
  }
  cat(sprintf("\n  lambda_N fisso in gp7/gp8 = 0.636 — confronta con le ratio sopra\n"))

  cat("\n  rho_W_SN fissato a 0 (gp9a: 0.041, CI90 [-0.085, 0.165]; GAM: 0.078)\n")


  # ── 9. CONFRONTO LOO gp7, gp8, gp9 ───────────────────────────────────────────

  if (requireNamespace("loo", quietly = TRUE)) {
    library(loo)

    ll_gp9  <- fit_gp9$draws("log_lik", format = "matrix")
    loo_gp9 <- loo(ll_gp9)
    cat("\n── LOO-CV gp9 ────────────────────────────────────────────────────\n")
    print(loo_gp9)
    rm(ll_gp9); gc()

    gp8_path <- here("stan", "fit_m4rr_gp8.rds")
    gp7_path <- here("stan", "fit_m4rr_gp7.rds")

    loo_list  <- list(gp9 = loo_gp9)

    if (file.exists(gp8_path)) {
      fit_gp8 <- readRDS(gp8_path)
      ll_gp8  <- fit_gp8$draws("log_lik", format = "matrix")
      loo_list$gp8 <- loo(ll_gp8)
      rm(fit_gp8, ll_gp8); gc()
    }
    if (file.exists(gp7_path)) {
      fit_gp7 <- readRDS(gp7_path)
      ll_gp7  <- fit_gp7$draws("log_lik", format = "matrix")
      loo_list$gp7 <- loo(ll_gp7)
      rm(fit_gp7, ll_gp7); gc()
    }

    if (length(loo_list) > 1) {
      cat("\n── Confronto LOO gp7 / gp8 / gp9 ────────────────────────────────\n")
      # unname() necessario: do.call con lista nominata passerebbe i nomi come
      # argomenti di loo_compare, ma il primo argomento si chiama "x" e non
      # "gp9"/"gp8"/"gp7" — causerebbe "argomento x assente"
      cmp <- do.call(loo_compare, unname(loo_list))
      rownames(cmp) <- names(loo_list)[order(sapply(loo_list, function(l)
        l$estimates["elpd_loo", "Estimate"]), decreasing = TRUE)]
      print(cmp)
    }
  } else {
    cat("\nInstalla 'loo' per il confronto predittivo.\n")
  }


  # ── 10. FACTOR SCORES SPAZIALI ───────────────────────────────────────────────

  org_draws <- fit_gp9$draws("eta_org_B_out", format = "matrix")
  P_draws   <- fit_gp9$draws("eta_P_B_out",   format = "matrix")
  N_draws   <- fit_gp9$draws("alpha_N_out",   format = "matrix")

  spatial_df_gp9 <- coords_km |>
    mutate(
      eta_org_mean = colMeans(org_draws),
      eta_org_sd   = apply(org_draws, 2, sd),
      alpha_N_mean = colMeans(N_draws),
      eta_P_mean   = colMeans(P_draws),
      eta_P_sd     = apply(P_draws, 2, sd)
    )

  rm(org_draws, P_draws, N_draws); gc()

  save(diag_list, smry_gp9, spatial_df_gp9, file = risultati_path_gp9)
  cat("\nRisultati gp9 salvati in:", risultati_path_gp9, "\n")

} else {
  cat("\nCarico risultati gp9 da:", risultati_path_gp9, "\n")
  load(risultati_path_gp9)
}


# ── 11. TRACE PLOTS ──────────────────────────────────────────────────────────

if (requireNamespace("bayesplot", quietly = TRUE)) {
  params_key <- c("gamma_SOC[1]", "gamma_N[1]",
                  "gamma_SOC[3]", "gamma_N[3]",
                  "sigma_W_SOC",  "sigma_W_N",
                  "sigma_GP_org", "rho_org",
                  "sigma_B_SOC",  "sigma_B_N")
  draws_key <- fit_gp9$draws(variables = params_key)
  print(mcmc_trace(draws_key) + theme_minimal())
  rm(draws_key); gc()
}


# ── 12. PLOT: gamma_SOC vs gamma_N ───────────────────────────────────────────

df_gamma <- bind_rows(
  smry_gp9 |>
    filter(startsWith(variable, "gamma_SOC")) |>
    mutate(risposta = "logSOC", predittore = X_W_cols),
  smry_gp9 |>
    filter(startsWith(variable, "gamma_N[")) |>
    mutate(risposta = "logN", predittore = X_W_cols)
)

p_gamma <- ggplot(df_gamma, aes(x = predittore, y = mean,
                                 colour = risposta, shape = risposta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
  geom_pointrange(aes(ymin = q5, ymax = q95),
                  position = position_dodge(0.4), linewidth = 0.7, size = 0.9) +
  scale_colour_manual(values = c("logSOC" = "steelblue", "logN" = "firebrick")) +
  labs(
    title    = "Effetti within-field: gamma_SOC vs gamma_N (gp9)",
    subtitle = paste0("CI 90%. Se ratio ≈ 0.636 il vincolo lambda_N era giustificato.",
                      "  b_slope = ", round(b_hat, 3)),
    x = "Predittore within-field (X_W)",
    y = "Coefficiente (scala log, X_W standardizzato)",
    colour = NULL, shape = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_gamma)


# ── 13. PLOT SPAZIALE ────────────────────────────────────────────────────────

if (exists("spatial_df_gp9") && nrow(spatial_df_gp9) > 0) {

  p_org <- ggplot(spatial_df_gp9,
                  aes(x = x_km, y = y_km, fill = eta_org_mean)) +
    geom_point(shape = 21, size = 4, stroke = 0.4) +
    scale_fill_viridis_c(option = "D", name = "Fattore\norganico") +
    labs(title    = "Fattore organico between (gp9, senza fattore latente within)",
         subtitle = sprintf("b_slope = %.3f | rho_org ≈ %.1f km", b_hat,
                            smry_gp9$median[smry_gp9$variable == "rho_org"]),
         x = "Est (km)", y = "Nord (km)") +
    theme_minimal()

  print(p_org)
}

cat("\n── Fine script 11 ───────────────────────────────────────────────────\n")
