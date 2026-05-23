# =============================================================================
# 12. gp9c_sparse_within.R
#
# Sparsificazione di gp9b: rimozione dei gamma_N non significativi.
#
# MOTIVAZIONE (da gp9b):
#   gamma_N[2] (Texture1 → N):      -0.044, CI90 [-0.098,  0.010] — include 0 → RIMOSSO
#   gamma_N[4] (BulkDensity → N):   -0.057, CI90 [-0.177,  0.064] — include 0 → RIMOSSO
#   gamma_N[1] (logBottom → N):     -0.261, CI90 [-0.316, -0.207] — significativo
#   gamma_N[3] (Texture2 → N):      -0.222, CI90 [-0.332, -0.114] — significativo
#
# STRUTTURA DI gp9c:
#   - gamma_SOC[4]: effetti within su logSOC, prior normal(0, 1)  — invariato
#   - gamma_N[2]:   effetti within su logN (solo logBottom + Texture2), prior normal(0, 0.3)
#   - X_W_N = matrice [N, 2] con colonne {logBottom, Texture2} — passata a Stan
#   - Texture1 e BulkDensity su N: fissati a 0 (rimossi come parametri)
#   - Between-field, GP, b_slope: invariati rispetto a gp9b / gp8
#
# DOMANDE RISPOSTE:
#   1. La sparsificazione migliora LOO rispetto a gp9b?
#   2. gp9c recupera almeno in parte il gap con gp8?
#   3. Quanti parametri effettivi: gp9c ha +2 parametri netti vs gp8
#      (4 gamma_SOC + 2 gamma_N liberi, vs 4 gamma_org in gp8)
#
# CONFRONTO FINALE: gp7 vs gp8 vs gp9b vs gp9c
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

# X_W completa per logSOC (4 predittori)
X_W_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity")
X_W      <- as.matrix(dati_int[, X_W_cols])

# X_W_N ridotta per logN: solo logBottom e Texture2 (significativi in gp9b)
X_W_N_cols <- c("logBottom", "Texture2")
X_W_N      <- as.matrix(dati_int[, X_W_N_cols])

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


# ── 3. DATI STAN ─────────────────────────────────────────────────────────────

stan_data_gp9c <- list(
  N        = N,
  J        = J,
  K_W      = length(X_W_cols),    # 4: tutte le covariate per logSOC
  K_W_N    = length(X_W_N_cols),  # 2: solo logBottom + Texture2 per logN
  K_B      = length(X_B_cols),
  field_id = dati_int$field_int,
  logSOC   = dati_int$logSOC,
  logN     = dati_int$logN,
  logP     = dati_int$logP,
  X_W      = X_W,    # [N × 4] per logSOC
  X_W_N    = X_W_N,  # [N × 2] per logN (Texture1 e BulkDensity → N esclusi)
  X_B      = X_B,
  dist_mat = dist_mat,
  b_slope  = b_hat
)

rm(X_W, X_W_N, X_B, dist_mat); gc()


# ── 4. COMPILAZIONE ──────────────────────────────────────────────────────────

stan_file_gp9c <- here("stan", "m4rr_gp9c.stan")
cat("\nCompilazione gp9c...\n")
mod_gp9c <- cmdstan_model(stan_file_gp9c, compile = TRUE)
cat("OK. CmdStan version:", cmdstan_version(), "\n")


# ── 5. MCMC ──────────────────────────────────────────────────────────────────

fit_path_gp9c       <- here("stan", "fit_m4rr_gp9c.rds")
risultati_path_gp9c <- here("stan", "risultati_gp9c.RData")

if (file.exists(fit_path_gp9c))       { file.remove(fit_path_gp9c);       cat("Vecchio fit_gp9c rimosso\n") }
if (file.exists(risultati_path_gp9c)) { file.remove(risultati_path_gp9c); cat("Vecchi risultati_gp9c rimossi\n") }

if (!file.exists(fit_path_gp9c)) {

  cat("\nAvvio MCMC gp9c (4 catene × 5000 sampling + 3000 warmup)...\n")

  fit_gp9c <- mod_gp9c$sample(
    data            = stan_data_gp9c,
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

  fit_gp9c$save_object(fit_path_gp9c)
  cat("Fit salvato in:", fit_path_gp9c, "\n")

  csv_files  <- list.files(here("stan"), pattern = "\\.csv$",  full.names = TRUE)
  json_files <- list.files(here("stan"), pattern = "\\.json$", full.names = TRUE)
  invisible(file.remove(c(csv_files, json_files)))

} else {
  cat("\nCarico gp9c da:", fit_path_gp9c, "\n")
  fit_gp9c <- readRDS(fit_path_gp9c)
}

rm(stan_data_gp9c); gc()


# ── 6. DIAGNOSTICA ───────────────────────────────────────────────────────────

if (!file.exists(risultati_path_gp9c)) {

  smry_full <- fit_gp9c$summary()
  np        <- fit_gp9c$sampler_diagnostics()

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


  # ── 7. SUMMARY PARAMETRI ─────────────────────────────────────────────────

  vars_key <- c(
    "gamma_SOC", "gamma_N",
    "gamma_P_B", "gamma_N_B",
    "sigma_W_SOC", "sigma_W_N", "psi_W_P",
    "sigma_B_SOC", "sigma_B_N",
    "psi_B_org", "sigma_P_between",
    "sigma_GP_org", "rho_org"
  )

  smry_gp9c <- fit_gp9c$summary(vars_key) |>
    select(variable, mean, median, sd, q5, q95, rhat, ess_bulk) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  cat("\n── Parametri gp9c ──\n")
  print(smry_gp9c, n = 40)


  # ── 8. CONFRONTO gamma_SOC vs gamma_N sparsificato ───────────────────────

  cat("\n── Effetti within-field gp9c ──────────────────────────────────────────\n")
  cat(sprintf("%-15s | %8s %8s\n", "Predittore", "gamma_SOC", "sd"))
  cat(strrep("-", 40), "\n")
  gSOC <- smry_gp9c |> filter(startsWith(variable, "gamma_SOC"))
  for (i in seq_len(nrow(gSOC))) {
    cat(sprintf("%-15s | %8.3f %8.3f\n",
                X_W_cols[i], gSOC$mean[i], gSOC$sd[i]))
  }

  cat(sprintf("\n%-15s | %8s %8s | %s\n", "Predittore", "gamma_N", "sd", "Note"))
  cat(strrep("-", 55), "\n")
  gN <- smry_gp9c |> filter(startsWith(variable, "gamma_N["))
  for (i in seq_along(X_W_N_cols)) {
    cat(sprintf("%-15s | %8.3f %8.3f | (stimato)\n",
                X_W_N_cols[i], gN$mean[i], gN$sd[i]))
  }
  cat(sprintf("%-15s | %8s %8s | FISSATO A 0 (da gp9b)\n", "Texture1",     "", ""))
  cat(sprintf("%-15s | %8s %8s | FISSATO A 0 (da gp9b)\n", "BulkDensity",  "", ""))


  # ── 9. CONFRONTO LOO gp7 / gp8 / gp9b / gp9c ────────────────────────────

  if (requireNamespace("loo", quietly = TRUE)) {
    library(loo)

    ll_gp9c  <- fit_gp9c$draws("log_lik", format = "matrix")
    loo_gp9c <- loo(ll_gp9c)
    cat("\n── LOO-CV gp9c ──────────────────────────────────────────────────\n")
    print(loo_gp9c)
    rm(ll_gp9c); gc()

    # gp9b è salvato come fit_m4rr_gp9.rds (script 11)
    gp9b_path <- here("stan", "fit_m4rr_gp9.rds")
    gp8_path  <- here("stan", "fit_m4rr_gp8.rds")
    gp7_path  <- here("stan", "fit_m4rr_gp7.rds")

    loo_list <- list(gp9c = loo_gp9c)

    if (file.exists(gp9b_path)) {
      fit_gp9b <- readRDS(gp9b_path)
      ll_gp9b  <- fit_gp9b$draws("log_lik", format = "matrix")
      loo_list$gp9b <- loo(ll_gp9b)
      rm(fit_gp9b, ll_gp9b); gc()
    }
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
      cat("\n── Confronto LOO gp7 / gp8 / gp9b / gp9c ──────────────────────\n")
      cmp <- do.call(loo_compare, unname(loo_list))
      rownames(cmp) <- names(loo_list)[order(
        sapply(loo_list, function(l) l$estimates["elpd_loo", "Estimate"]),
        decreasing = TRUE)]
      print(cmp)
    }
  } else {
    cat("\nInstalla 'loo' per il confronto predittivo.\n")
  }


  # ── 10. FACTOR SCORES SPAZIALI ───────────────────────────────────────────

  org_draws <- fit_gp9c$draws("eta_org_B_out", format = "matrix")
  P_draws   <- fit_gp9c$draws("eta_P_B_out",   format = "matrix")
  N_draws   <- fit_gp9c$draws("alpha_N_out",   format = "matrix")

  spatial_df_gp9c <- coords_km |>
    mutate(
      eta_org_mean = colMeans(org_draws),
      eta_org_sd   = apply(org_draws, 2, sd),
      alpha_N_mean = colMeans(N_draws),
      eta_P_mean   = colMeans(P_draws),
      eta_P_sd     = apply(P_draws, 2, sd)
    )

  rm(org_draws, P_draws, N_draws); gc()

  save(diag_list, smry_gp9c, spatial_df_gp9c, file = risultati_path_gp9c)
  cat("\nRisultati gp9c salvati in:", risultati_path_gp9c, "\n")

} else {
  cat("\nCarico risultati gp9c da:", risultati_path_gp9c, "\n")
  load(risultati_path_gp9c)
}


# ── 11. TRACE PLOTS ──────────────────────────────────────────────────────────

if (requireNamespace("bayesplot", quietly = TRUE)) {
  params_key <- c("gamma_SOC[1]", "gamma_N[1]",
                  "gamma_SOC[3]", "gamma_N[2]",   # gamma_N[2] = Texture2
                  "sigma_W_SOC",  "sigma_W_N",
                  "sigma_GP_org", "rho_org",
                  "sigma_B_SOC",  "sigma_B_N")
  draws_key <- fit_gp9c$draws(variables = params_key)
  print(mcmc_trace(draws_key) + theme_minimal())
  rm(draws_key); gc()
}


# ── 12. PLOT: gamma_SOC vs gamma_N sparsificato ──────────────────────────────

# Costruisce data frame per il confronto, includendo le entry fissate a 0
gSOC_df <- smry_gp9c |>
  filter(startsWith(variable, "gamma_SOC")) |>
  mutate(risposta = "logSOC", predittore = X_W_cols, stimato = TRUE)

gN_df <- smry_gp9c |>
  filter(startsWith(variable, "gamma_N[")) |>
  mutate(risposta = "logN", predittore = X_W_N_cols, stimato = TRUE)

# Aggiunge le righe "zero" per Texture1 e BulkDensity su logN
gN_zero <- tibble(
  variable = c("gamma_N_fix[2]", "gamma_N_fix[4]"),
  mean = 0, median = 0, sd = 0, q5 = 0, q95 = 0, rhat = NA, ess_bulk = NA,
  risposta = "logN",
  predittore = c("Texture1", "BulkDensity"),
  stimato = FALSE
)

df_gamma <- bind_rows(gSOC_df, gN_df, gN_zero) |>
  mutate(predittore = factor(predittore, levels = X_W_cols))

p_gamma <- ggplot(df_gamma, aes(x = predittore, y = mean,
                                 colour = risposta, shape = stimato)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
  geom_pointrange(aes(ymin = q5, ymax = q95),
                  position = position_dodge(0.4), linewidth = 0.7, size = 0.9) +
  scale_colour_manual(values = c("logSOC" = "steelblue", "logN" = "firebrick")) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "stimato", "FALSE" = "fissato = 0")) +
  labs(
    title    = "Effetti within-field: gamma_SOC vs gamma_N (gp9c sparsificato)",
    subtitle = paste0("CI 90%. Texture1 e BulkDensity su N fissati a 0.  b_slope = ",
                      round(b_hat, 3)),
    x = "Predittore within-field (X_W)",
    y = "Coefficiente (scala log, X_W standardizzato)",
    colour = NULL, shape = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_gamma)


# ── 13. PLOT SPAZIALE ────────────────────────────────────────────────────────

if (exists("spatial_df_gp9c") && nrow(spatial_df_gp9c) > 0) {

  p_org <- ggplot(spatial_df_gp9c,
                  aes(x = x_km, y = y_km, fill = eta_org_mean)) +
    geom_point(shape = 21, size = 4, stroke = 0.4) +
    scale_fill_viridis_c(option = "D", name = "Fattore\norganico") +
    labs(title    = "Fattore organico between (gp9c, gamma_N sparsificato)",
         subtitle = sprintf("b_slope = %.3f | rho_org ≈ %.1f km", b_hat,
                            smry_gp9c$median[smry_gp9c$variable == "rho_org"]),
         x = "Est (km)", y = "Nord (km)") +
    theme_minimal()

  print(p_org)
}

cat("\n── Fine script 12 (gp9c) ────────────────────────────────────────────\n")
