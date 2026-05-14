# =============================================================================
# 7. bayesian_gp7.R
#
# Rispetto a scripts/6. bayesian_gp3.R, unica modifica strutturale:
# il trend lineare sulle coordinate (beta_N_x, beta_N_y) viene sostituito
# da covariati di management between per logN: gamma_N_B * X_B.
#
# Motivazione:
# gp3 ha mostrato beta_N_x = -0.106 [CI -0.96, +0.77] e
# beta_N_y = -0.040 [CI -0.94, +0.85]: nessun gradiente direzionale.
# sigma_B_N = 2.10 (identico al caso i.i.d. prima del GP in gp2).
# Ne' il GP Matern (rho_N = 20.9 km, degenere) ne' il trend lineare
# catturano la variazione between di logN. La spiegazione piu' plausibile
# sono fattori di management (OnFarm, Irrigate, Fertilised, N_Natural),
# gia' usati come predittori between per logP. Si aggiunge gamma_N_B.
#
# RISULTATI run gp7 (modello finale):
# sigma_B_N = 0.94 — scende da 2.10 (gp3) ma resta sostanziale.
# I covariati spiegano ~80% della varianza between-field di logN (1-(0.94/2.10)^2),
# quasi interamente attraverso N_Natural:
#   gamma_N_B[1] OnFarm:     +0.364, CI [-0.284, +1.00]  — zero
#   gamma_N_B[2] Irrigate:   -0.228, CI [-0.745, +0.285] — zero
#   gamma_N_B[3] Fertilised: -0.238, CI [-0.809, +0.339] — zero
#   gamma_N_B[4] N_Natural:  -2.26,  CI [-2.67, -1.83]   — robusto
# Il sigma_B_N residuo = 0.94 e' irriducibile. Non si escludono OnFarm/Irrigate/
# Fertilised: la prior normal(0,1) li shrinka gia', e la simmetria con gamma_P_B
# e' preferibile. Struttura a 2 cluster (SW vs resto) valutata e rigettata:
# N_Natural come covariata cattura gia' la differenza di media tra i blocchi.
# =============================================================================


# ── 0. SETUP ─────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr",
                   repos = c("https://mc-stan.org/r-packages/", "https://cloud.r-project.org"))
}
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

if (requireNamespace("posterior", quietly = TRUE)) library(posterior)
if (requireNamespace("bayesplot", quietly = TRUE)) library(bayesplot)
if (requireNamespace("lavaan",    quietly = TRUE)) library(lavaan)


# ── 1. DATI ──────────────────────────────────────────────────────────────────

dati <- readRDS(here("data", "dati.rds"))
dati <- dati |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom)
  ) |>
  mutate(across(c(logBottom, PH, Texture1, Texture2, BulkDensity), ~ c(scale(.x))))

fit4rr_path <- here("stan", "fit4rr_freq.rds")
if (file.exists(fit4rr_path)) {
  fit4rr <- readRDS(fit4rr_path)
} else {
  message("fit4rr non trovato: esegui 3. mimic.R per generarlo.")
}

crop_raw <- readRDS(here("data", "crop.rds"))
coords_field <- crop_raw |>
  group_by(Field) |>
  summarise(Lat = mean(Lat, na.rm = TRUE), Long = mean(Long, na.rm = TRUE)) |>
  arrange(Field)
rm(crop_raw); gc()


# ── 2. PREPARAZIONE DATI STAN ────────────────────────────────────────────────

field_levels <- sort(unique(dati$Field))
J <- length(field_levels)
N <- nrow(dati)
cat("N osservazioni:", N, " | J field:", J, "\n")

dati <- dati |>
  mutate(field_int = as.integer(factor(Field, levels = field_levels)))

X_W_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity")
X_W <- as.matrix(dati[, X_W_cols])

X_B_cols <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")
X_B <- dati |>
  distinct(field_int, across(all_of(X_B_cols))) |>
  arrange(field_int) |>
  select(all_of(X_B_cols)) |>
  as.matrix()

stopifnot(nrow(X_B) == J)


# ── 3. COORDINATE E MATRICE DISTANZE ─────────────────────────────────────────
# Coordinate usate solo per la matrice distanze (GP organico).
# Non piu' standardizzate per trend lineare (rimosso in gp7).

ref_lat  <- mean(coords_field$Lat, na.rm = TRUE)
ref_long <- mean(coords_field$Long, na.rm = TRUE)

coords_km <- coords_field |>
  arrange(Field) |>
  mutate(
    x_km = (Long - ref_long) * 111.32 * cos(ref_lat * pi / 180),
    y_km = (Lat  - ref_lat)  * 110.54
  ) |>
  select(Field, x_km, y_km)

dist_mat <- as.matrix(dist(coords_km[, c("x_km", "y_km")]))

cat("Distanze inter-field (km): min =", round(min(dist_mat[dist_mat > 0]), 2),
    " mediana =", round(median(dist_mat[dist_mat > 0]), 1),
    " max =", round(max(dist_mat), 1), "\n")


# ── 4. LAMBDA_N FISSO ────────────────────────────────────────────────────────

if (exists("fit4rr")) {
  lambda_N_fixed <- parameterEstimates(fit4rr) |>
    filter(label == "a", op == "=~") |>
    slice(1) |>
    pull(est)
} else {
  lambda_N_fixed <- 0.6355
}
cat(sprintf("lambda_N within fisso: %.4f\n", lambda_N_fixed))


# ── 5. DIMENSIONE K_W ────────────────────────────────────────────────────────

K_W <- length(X_W_cols)
cat(sprintf("K_W = %d predittori within\n", K_W))


# ── 6. LISTA DATI STAN ────────────────────────────────────────────────────────

stan_data <- list(
  N          = N,
  J          = J,
  K_W        = K_W,
  K_B        = length(X_B_cols),
  field_id   = dati$field_int,
  logSOC     = dati$logSOC,
  logN       = dati$logN,
  logP       = dati$logP,
  X_W        = X_W,
  X_B        = X_B,
  dist_mat   = dist_mat,
  lambda_N   = lambda_N_fixed
)

rm(X_W, X_B, dist_mat); gc()


# ── 7. COMPILAZIONE MODELLO ───────────────────────────────────────────────────

stan_file <- here("stan", "m4rr_gp7.stan")
cat("\nCompilazione del modello Stan (usa cache se gia' compilato)...\n")
mod <- cmdstan_model(stan_file, compile = TRUE)
cat("OK. CmdStan version:", cmdstan_version(), "\n")


# ── 8. TEST RUN ───────────────────────────────────────────────────────────────

# fit_test <- mod$sample(
#   data = stan_data, seed = 42, chains = 1, parallel_chains = 1,
#   iter_warmup = 200, iter_sampling = 50, adapt_delta = 0.90, refresh = 50
# )
# print(fit_test$summary(c("gamma_org","gamma_N_B","gamma_P_B",
#                           "sigma_GP_org","rho_org","sigma_B_N","sigma_B_SOC")))
# cat("Divergenze test:", sum(fit_test$sampler_diagnostics()[,,"divergent__"]), "\n")
# rm(fit_test); gc()


# ── 9. MCMC ───────────────────────────────────────────────────────────────────

fit_path <- here("stan", "fit_m4rr_gp7.rds")

if (!file.exists(fit_path)) {

  cat("\nAvvio MCMC (4 catene x 3000 warmup + 5000 sampling)...\n")
  cat("Aggiornamenti ogni 200 iterazioni. Stima: 10-30 min.\n\n")

  fit <- mod$sample(
    data            = stan_data,
    seed            = 2024,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 3000,
    iter_sampling   = 5000,
    adapt_delta     = 0.97,
    max_treedepth   = 12,
    refresh         = 200,
    show_messages   = TRUE,
    output_dir      = here("stan")
  )

  fit$save_object(fit_path)
  cat("\nFit salvato in:", fit_path, "\n")

  csv_files  <- list.files(here("stan"), pattern = "\\.csv$",  full.names = TRUE)
  json_files <- list.files(here("stan"), pattern = "\\.json$", full.names = TRUE)
  invisible(file.remove(c(csv_files, json_files)))
  cat("CSV e JSON delle catene eliminati.\n")

} else {
  cat("\nCarico fit salvato da:", fit_path, "\n")
  fit <- readRDS(fit_path)
}

rm(stan_data); gc()


# ── 10. ANALISI ───────────────────────────────────────────────────────────────

risultati_path <- here("stan", "risultati_bayesiano7.RData")

if (!file.exists(risultati_path)) {

  # ── 10a. Diagnostiche ────────────────────────────────────────────────────
  cat("\nCalcolo diagnostiche...\n")

  smry_full <- fit$summary()
  np        <- fit$sampler_diagnostics()

  diag_list <- list(
    bad_rhat = smry_full |> filter(rhat > 1.05) |> arrange(desc(rhat)),
    bad_ess  = smry_full |> filter(ess_bulk < 400) |> arrange(ess_bulk),
    n_divs   = sum(np[, , "divergent__"]),
    max_td   = max(np[, , "treedepth__"]),
    ebfmi    = apply(np[, , "energy__"], 2, function(e) mean(diff(e)^2) / var(e))
  )

  cat("\n── Rhat > 1.05 ──\n");    print(diag_list$bad_rhat)
  cat("\n── ESS bulk < 400 ──\n"); print(head(diag_list$bad_ess, 20))
  cat("Divergenze:", diag_list$n_divs, "| Max treedepth:", diag_list$max_td, "\n")
  cat("E-BFMI per catena:", round(diag_list$ebfmi, 3), "\n")

  rm(smry_full, np); gc()

  # ── 10b. Summary parametri strutturali ───────────────────────────────────
  cat("\nCalcolo summary parametri...\n")

  cat(sprintf("\nlambda_N within fisso: %.4f (MLE M4rr, SE=0.059)\n", lambda_N_fixed))

  vars_key <- c(
    "gamma_org",       # [1]=logBottom [2]=Texture1 [3]=Texture2 [4]=BulkDensity
    "gamma_P_B",       # [1]=OnFarm [2]=Irrigate [3]=Fertilised [4]=N_Natural (*)
    "gamma_N_B",       # [1]=OnFarm [2]=Irrigate [3]=Fertilised [4]=N_Natural (*)
    # (*) [4] N_Natural: confounding geografico per P; per N potrebbe essere
    #     interpretabile (differenza di management reale) ma monitorare.
    "psi_W_org", "psi_W_P", "theta_W_SOC", "theta_W_N",
    "sigma_B_SOC", "sigma_B_N",
    "psi_B_org",
    "sigma_P_between",
    "sigma_GP_org", "rho_org"
  )

  smry <- fit$summary(vars_key) |>
    select(variable, mean, median, sd, q5, q95, rhat, ess_bulk) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  cat("\n── Parametri strutturali ──\n")
  print(smry, n = 35)

  # ── 10c. Factor scores ────────────────────────────────────────────────────
  cat("\nEstrazione factor scores between...\n")

  org_draws  <- fit$draws("eta_org_B_out", format = "matrix")
  P_draws    <- fit$draws("eta_P_B_out",   format = "matrix")
  N_draws    <- fit$draws("alpha_N_out",   format = "matrix")

  spatial_df <- coords_km |>
    mutate(
      eta_org_mean = colMeans(org_draws),
      eta_org_sd   = apply(org_draws, 2, sd),
      alpha_N_mean = colMeans(N_draws),
      alpha_N_sd   = apply(N_draws,   2, sd),
      eta_P_mean   = colMeans(P_draws),
      eta_P_sd     = apply(P_draws,   2, sd)
    )

  rm(org_draws, P_draws, N_draws); gc()

  # ── 10d. Confronto frequentista-bayesiano ─────────────────────────────────
  confronto_df <- NULL
  if (exists("fit4rr") && requireNamespace("lavaan", quietly = TRUE)) {
    library(lavaan)
    freq <- parameterEstimates(fit4rr) |>
      filter(op == "~") |>
      select(lhs, rhs, est_freq = est, se_freq = se, z_freq = z)

    gorg_draws <- fit$draws("gamma_org", format = "matrix")
    gPB_draws  <- fit$draws("gamma_P_B", format = "matrix")
    gNB_draws  <- fit$draws("gamma_N_B", format = "matrix")

    bayes <- bind_rows(
      tibble(lhs = "fertil_org_W",
             rhs = X_W_cols,
             est_bay = colMeans(gorg_draws),
             sd_bay  = apply(gorg_draws, 2, sd),
             q5_bay  = apply(gorg_draws, 2, quantile, 0.05),
             q95_bay = apply(gorg_draws, 2, quantile, 0.95)),
      tibble(lhs = "fertil_P_B",
             rhs = c("OnFarm","Irrigate","Fertilised","N_Natural"),
             est_bay = colMeans(gPB_draws),
             sd_bay  = apply(gPB_draws, 2, sd),
             q5_bay  = apply(gPB_draws, 2, quantile, 0.05),
             q95_bay = apply(gPB_draws, 2, quantile, 0.95)),
      tibble(lhs = "fertil_N_B",
             rhs = c("OnFarm","Irrigate","Fertilised","N_Natural"),
             est_bay = colMeans(gNB_draws),
             sd_bay  = apply(gNB_draws, 2, sd),
             q5_bay  = apply(gNB_draws, 2, quantile, 0.05),
             q95_bay = apply(gNB_draws, 2, quantile, 0.95))
    )
    confronto_df <- bayes |> left_join(freq, by = c("lhs","rhs"))
    rm(gorg_draws, gPB_draws, gNB_draws, freq, bayes); gc()

    cat("\n── Confronto frequentista vs bayesiano ──\n")
    print(confronto_df |> mutate(across(where(is.numeric), ~ round(.x, 3))))
  } else {
    cat("\nfit4rr non disponibile: esegui 3. mimic.R per il confronto.\n")
  }

  # ── Salva ────────────────────────────────────────────────────────────────
  save(diag_list, smry, spatial_df, confronto_df, file = risultati_path)
  cat("\nRisultati salvati in:", risultati_path, "\n")

} else {
  cat("\nCarico risultati salvati da:", risultati_path, "\n")
  load(risultati_path)
  cat("Oggetti caricati: diag_list, smry, spatial_df, confronto_df\n")
}


# ── 11. TRACE PLOTS ───────────────────────────────────────────────────────────

if (requireNamespace("bayesplot", quietly = TRUE)) {
  params_key <- c("gamma_org[1]", "gamma_org[2]",
                  "gamma_org[3]", "gamma_org[4]",
                  "gamma_N_B[1]", "gamma_N_B[2]",
                  "gamma_N_B[3]", "gamma_N_B[4]",
                  "sigma_GP_org", "rho_org",
                  "sigma_B_N",    "sigma_P_between")
  draws_key <- fit$draws(variables = params_key)
  print(mcmc_trace(draws_key) + theme_minimal())
  rm(draws_key); gc()
}


# ── 12. PLOT SPAZIALE ─────────────────────────────────────────────────────────

if (exists("spatial_df") && nrow(spatial_df) > 0) {

  p_org <- ggplot(spatial_df, aes(x = x_km, y = y_km, fill = eta_org_mean)) +
    geom_point(shape = 21, size = 4, stroke = 0.4) +
    scale_fill_viridis_c(option = "D", name = "Fattore\norganico") +
    labs(title = "Fattore organico between (SOC, GP Matern 3/2)",
         x = "Est (km)", y = "Nord (km)") +
    theme_minimal()

  p_N <- ggplot(spatial_df, aes(x = x_km, y = y_km, fill = alpha_N_mean)) +
    geom_point(shape = 21, size = 4, stroke = 0.4) +
    scale_fill_viridis_c(option = "A", name = "alpha_N") +
    labs(title = "Intercetta N between (management + i.i.d., media posteriore)",
         subtitle = "X_B * gamma_N_B + sigma_B_N * raw",
         x = "Est (km)", y = "Nord (km)") +
    theme_minimal()

  p_P <- ggplot(spatial_df, aes(x = x_km, y = y_km, fill = eta_P_mean)) +
    geom_point(shape = 21, size = 4, stroke = 0.4) +
    scale_fill_viridis_c(option = "C", name = "Fattore P") +
    labs(title = "Fattore fosforo between (media posteriore)",
         x = "Est (km)", y = "Nord (km)") +
    theme_minimal()

  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    print(p_org + p_N + p_P)
  } else {
    print(p_org); print(p_N); print(p_P)
  }
}
