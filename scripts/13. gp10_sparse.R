# =============================================================================
# 13. gp10_sparse.R
#
# Sparsificazione di gp9c: rimozione effetti between-field non significativi
# e del residuo SOC between-field (alpha_SOC_res).
#
# MODIFICHE RISPETTO A gp9c:
#
# 1. Between-field sparsificato (da CI90 in gp9c):
#    gamma_N_B[2] (Irrigate → N_between):    -0.218, CI90 [-0.709, 0.261] → RIMOSSO
#    gamma_N_B[3] (Fertilised → N_between):  -0.319, CI90 [-0.887, 0.242] → RIMOSSO
#    gamma_P_B[2] (Irrigate → P_between):     0.275, CI90 [-0.306, 0.846] → RIMOSSO
#    gamma_P_B[3] (Fertilised → P_between):  -0.464, CI90 [-1.10,  0.173] → RIMOSSO
#    → X_B ridotta a [J x 2] con sole colonne {OnFarm, N_Natural}
#
# 2. Rimozione alpha_SOC_res:
#    sigma_B_SOC = 0.097 (mediana 0.084): ~0 → alpha_SOC_raw[J] + sigma_B_SOC
#    aggiungono 55 parametri con contributo trascurabile.
#    Il GP organico è sufficiente per la variabilità between-SOC.
#
# STRUTTURA WITHIN-FIELD: identica a gp9c.
#
# CONFRONTO FINALE: gp8 / gp9b / gp9c / gp10
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

# ── Design matrix within ──────────────────────────────────────────────────────

X_W_cols  <- c("logBottom", "Texture1", "Texture2", "BulkDensity")
X_W       <- as.matrix(dati_int[, X_W_cols])

X_W_N_cols <- c("logBottom", "Texture2")          # solo effetti sign. su N
X_W_N      <- as.matrix(dati_int[, X_W_N_cols])

# ── Design matrix between (RIDOTTA: solo OnFarm e N_Natural) ─────────────────
# Irrigate e Fertilised rimossi: CI90 che include 0 per sia gamma_N_B sia gamma_P_B

X_B_cols <- c("OnFarm", "N_Natural")              # 2 colonne, non più 4
X_B <- dati_int |>
  distinct(field_int, across(all_of(X_B_cols))) |>
  arrange(field_int) |>
  select(all_of(X_B_cols)) |>
  as.matrix()

# ── Coordinate e matrice distanze ─────────────────────────────────────────────

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


# ── 2. CARICA b_slope DA gp8 ─────────────────────────────────────────────────

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

stan_data_gp10 <- list(
  N        = N,
  J        = J,
  K_W      = length(X_W_cols),    # 4
  K_W_N    = length(X_W_N_cols),  # 2
  K_B      = length(X_B_cols),    # 2 (OnFarm + N_Natural)
  field_id = dati_int$field_int,
  logSOC   = dati_int$logSOC,
  logN     = dati_int$logN,
  logP     = dati_int$logP,
  X_W      = X_W,
  X_W_N    = X_W_N,
  X_B      = X_B,
  dist_mat = dist_mat,
  b_slope  = b_hat
)

rm(X_W, X_W_N, X_B, dist_mat); gc()


# ── 4. COMPILAZIONE ──────────────────────────────────────────────────────────

stan_file_gp10 <- here("stan", "m4rr_gp10.stan")
cat("\nCompilazione gp10...\n")
mod_gp10 <- cmdstan_model(stan_file_gp10, compile = TRUE)
cat("OK. CmdStan version:", cmdstan_version(), "\n")


# ── 5. MCMC ──────────────────────────────────────────────────────────────────

fit_path_gp10       <- here("stan", "fit_m4rr_gp10.rds")
risultati_path_gp10 <- here("stan", "risultati_gp10.RData")

if (file.exists(fit_path_gp10))       { file.remove(fit_path_gp10);       cat("Vecchio fit_gp10 rimosso\n") }
if (file.exists(risultati_path_gp10)) { file.remove(risultati_path_gp10); cat("Vecchi risultati_gp10 rimossi\n") }

if (!file.exists(fit_path_gp10)) {

  cat("\nAvvio MCMC gp10 (4 catene × 5000 sampling + 3000 warmup)...\n")

  fit_gp10 <- mod_gp10$sample(
    data            = stan_data_gp10,
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

  fit_gp10$save_object(fit_path_gp10)
  cat("Fit salvato in:", fit_path_gp10, "\n")

  csv_files  <- list.files(here("stan"), pattern = "\\.csv$",  full.names = TRUE)
  json_files <- list.files(here("stan"), pattern = "\\.json$", full.names = TRUE)
  invisible(file.remove(c(csv_files, json_files)))

} else {
  cat("\nCarico gp10 da:", fit_path_gp10, "\n")
  fit_gp10 <- readRDS(fit_path_gp10)
}

rm(stan_data_gp10); gc()


# ── 6. DIAGNOSTICA ───────────────────────────────────────────────────────────

if (!file.exists(risultati_path_gp10)) {

  smry_full <- fit_gp10$summary()
  np        <- fit_gp10$sampler_diagnostics()

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
    "sigma_B_N",
    "psi_B_org", "sigma_P_between",
    "sigma_GP_org", "rho_org"
  )

  smry_gp10 <- fit_gp10$summary(vars_key) |>
    select(variable, mean, median, sd, q5, q95, rhat, ess_bulk) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  cat("\n── Parametri gp10 ──\n")
  print(smry_gp10, n = 40)


  # ── 8. TABELLA EFFETTI ───────────────────────────────────────────────────

  cat("\n── Effetti within-field (gp10) ─────────────────────────────────────\n")
  cat(sprintf("%-15s | %8s %8s | %8s %8s\n",
              "Predittore", "γ_SOC", "sd", "γ_N", "sd"))
  cat(strrep("-", 58), "\n")

  gSOC <- smry_gp10 |> filter(startsWith(variable, "gamma_SOC"))
  gN   <- smry_gp10 |> filter(startsWith(variable, "gamma_N["))

  for (i in seq_along(X_W_cols)) {
    gs <- gSOC$mean[i]; ss <- gSOC$sd[i]
    ni <- which(X_W_N_cols == X_W_cols[i])
    if (length(ni) == 1) {
      gn <- gN$mean[ni]; sn <- gN$sd[ni]; note <- "(stimato)"
    } else {
      gn <- 0; sn <- NA; note <- "FISSATO A 0"
    }
    cat(sprintf("%-15s | %8.3f %8.3f | %8.3f %8s  %s\n",
                X_W_cols[i], gs, ss, gn,
                if (!is.na(sn)) sprintf("%.3f", sn) else "    —",
                note))
  }

  cat("\n── Effetti between-field (gp10, sparsificati) ───────────────────────\n")
  cat(sprintf("%-15s | %8s %8s | %8s %8s\n",
              "Predittore", "γ_P_B", "sd", "γ_N_B", "sd"))
  cat(strrep("-", 58), "\n")

  gPB <- smry_gp10 |> filter(startsWith(variable, "gamma_P_B"))
  gNB <- smry_gp10 |> filter(startsWith(variable, "gamma_N_B"))

  for (i in seq_along(X_B_cols)) {
    cat(sprintf("%-15s | %8.3f %8.3f | %8.3f %8.3f\n",
                X_B_cols[i],
                gPB$mean[i], gPB$sd[i],
                gNB$mean[i], gNB$sd[i]))
  }
  cat(sprintf("%-15s | %8s %8s | %8s %8s  FISSATO A 0\n",
              "Irrigate",   "", "", "", ""))
  cat(sprintf("%-15s | %8s %8s | %8s %8s  FISSATO A 0\n",
              "Fertilised", "", "", "", ""))


  # ── 9. LOO: gp8 / gp9b / gp9c / gp10 ───────────────────────────────────

  if (requireNamespace("loo", quietly = TRUE)) {
    library(loo)

    ll_gp10  <- fit_gp10$draws("log_lik", format = "matrix")
    loo_gp10 <- loo(ll_gp10)
    cat("\n── LOO-CV gp10 ──────────────────────────────────────────────────\n")
    print(loo_gp10)
    rm(ll_gp10); gc()

    gp9c_path <- here("stan", "fit_m4rr_gp9c.rds")
    gp9b_path <- here("stan", "fit_m4rr_gp9.rds")   # gp9b salvato come gp9
    gp8_path  <- here("stan", "fit_m4rr_gp8.rds")

    loo_list <- list(gp10 = loo_gp10)

    if (file.exists(gp9c_path)) {
      fit_gp9c <- readRDS(gp9c_path)
      ll_gp9c  <- fit_gp9c$draws("log_lik", format = "matrix")
      loo_list$gp9c <- loo(ll_gp9c)
      rm(fit_gp9c, ll_gp9c); gc()
    }
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

    if (length(loo_list) > 1) {
      cat("\n── Confronto LOO gp8 / gp9b / gp9c / gp10 ─────────────────────\n")
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

  org_draws <- fit_gp10$draws("eta_org_B_out", format = "matrix")
  P_draws   <- fit_gp10$draws("eta_P_B_out",   format = "matrix")
  N_draws   <- fit_gp10$draws("alpha_N_out",   format = "matrix")

  spatial_df_gp10 <- coords_km |>
    mutate(
      eta_org_mean = colMeans(org_draws),
      eta_org_sd   = apply(org_draws, 2, sd),
      alpha_N_mean = colMeans(N_draws),
      eta_P_mean   = colMeans(P_draws),
      eta_P_sd     = apply(P_draws, 2, sd)
    )

  rm(org_draws, P_draws, N_draws); gc()

  save(diag_list, smry_gp10, spatial_df_gp10, file = risultati_path_gp10)
  cat("\nRisultati gp10 salvati in:", risultati_path_gp10, "\n")

} else {
  cat("\nCarico risultati gp10 da:", risultati_path_gp10, "\n")
  load(risultati_path_gp10)
}


# ── 11. TRACE PLOTS ──────────────────────────────────────────────────────────

if (requireNamespace("bayesplot", quietly = TRUE)) {
  params_key <- c("gamma_SOC[1]", "gamma_N[1]",
                  "gamma_SOC[3]", "gamma_N[2]",
                  "gamma_N_B[1]", "gamma_N_B[2]",   # OnFarm, N_Natural
                  "sigma_GP_org", "rho_org",
                  "sigma_W_SOC",  "sigma_W_N")
  draws_key <- fit_gp10$draws(variables = params_key)
  print(mcmc_trace(draws_key) + theme_minimal())
  rm(draws_key); gc()
}


# ── 12. PLOT: effetti within ──────────────────────────────────────────────────

gSOC_df <- smry_gp10 |>
  filter(startsWith(variable, "gamma_SOC")) |>
  mutate(risposta = "logSOC", predittore = X_W_cols, stimato = TRUE)

gN_df <- smry_gp10 |>
  filter(startsWith(variable, "gamma_N[")) |>
  mutate(risposta = "logN", predittore = X_W_N_cols, stimato = TRUE)

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
    title    = "Effetti within-field: gp10 (sparsificato)",
    subtitle = paste0("CI 90%.  b_slope = ", round(b_hat, 3)),
    x = "Predittore within-field", y = "Coefficiente (log scale, standardizzato)",
    colour = NULL, shape = NULL
  ) +
  theme_minimal() + theme(legend.position = "bottom")
print(p_gamma)


# ── 13. PLOT SPAZIALE ────────────────────────────────────────────────────────

if (exists("spatial_df_gp10") && nrow(spatial_df_gp10) > 0) {

  p_org <- ggplot(spatial_df_gp10,
                  aes(x = x_km, y = y_km, fill = eta_org_mean)) +
    geom_point(shape = 21, size = 4, stroke = 0.4) +
    scale_fill_viridis_c(option = "D", name = "Fattore\norganico") +
    labs(title    = "Fattore organico between (gp10, fully sparse)",
         subtitle = sprintf("b_slope = %.3f | rho_org ≈ %.1f km", b_hat,
                            smry_gp10$median[smry_gp10$variable == "rho_org"]),
         x = "Est (km)", y = "Nord (km)") +
    theme_minimal()
  print(p_org)
}

cat("\n── Fine script 13 (gp10) ────────────────────────────────────────────\n")
