# ── 0. SETUP ────────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr",
                   repos = c("https://mc-stan.org/r-packages/", "https://cloud.r-project.org"))
}
library(cmdstanr)

# Fix PATH per Rtools45 su Windows
if (.Platform$OS.type == "windows") {
  rtools_path <- Sys.getenv("RTOOLS45_HOME", unset = "C:/rtools45")
  Sys.setenv(
    PATH = paste(file.path(rtools_path, "ucrt64/bin"),
                 file.path(rtools_path, "usr/bin"),
                 Sys.getenv("PATH"), sep = ";"),
    RTOOLS44_HOME = rtools_path
  )
}

# Prima esecuzione: decommenta una volta sola
# check_cmdstan_toolchain(fix = TRUE)
# install_cmdstan(cores = 4)

if (requireNamespace("posterior",  quietly = TRUE)) library(posterior)
if (requireNamespace("bayesplot",  quietly = TRUE)) library(bayesplot)


# ── 1. DATI ──────────────────────────────────────────────────────────────────────

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

# Modello frequentista (per confronto in sezione 9d)
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


# ── 2. PREPARAZIONE DATI STAN ───────────────────────────────────────────────────

field_levels <- sort(unique(dati$Field))
J <- length(field_levels)
N <- nrow(dati)
cat("N osservazioni:", N, " | J field:", J, "\n")

dati <- dati |>
  mutate(field_int = as.integer(factor(Field, levels = field_levels)))

X_W_cols <- c("logBottom", "PH", "Texture1", "Texture2", "BulkDensity")
X_W <- as.matrix(dati[, X_W_cols])

X_B_cols <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")
X_B <- dati |>
  distinct(field_int, across(all_of(X_B_cols))) |>
  arrange(field_int) |>
  select(all_of(X_B_cols)) |>
  as.matrix()

stopifnot(nrow(X_B) == J)


# ── 3. COORDINATE E MATRICE DISTANZE ────────────────────────────────────────────

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


# ── 4. HORSESHOE HYPERPARAMETRI ──────────────────────────────────────────────────
# tau0 = (p0 / (K_W - p0)) * sigma_noise / sqrt(N)
# p0=3 (logBottom, Texture2, BulkDensity attesi non-nulli)

p0         <- 3
K_W        <- length(X_W_cols)
sigma_appr <- 0.57          # sqrt(psi_W_org + theta_W_SOC) da M4rr
tau0       <- (p0 / (K_W - p0)) * sigma_appr / sqrt(N)
cat(sprintf("tau0 horseshoe: %.4f\n", tau0))

slab_scale <- 2.0
slab_df    <- 4.0


# ── 5. LISTA DATI STAN ───────────────────────────────────────────────────────────

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
  tau0       = tau0,
  slab_scale = slab_scale,
  slab_df    = slab_df
)

rm(X_W, X_B, dist_mat); gc()


# ── 6. COMPILAZIONE MODELLO ──────────────────────────────────────────────────────

stan_file <- here("stan", "m4rr_gp.stan")
cat("\nCompilazione del modello Stan (usa cache se già compilato)...\n")
mod <- cmdstan_model(stan_file, compile = TRUE)
cat("OK. CmdStan version:", cmdstan_version(), "\n")


# ── 7. TEST RUN ──────────────────────────────────────────────────────────────────
# Decommenta per verificare che il modello giri prima del run completo.

# fit_test <- mod$sample(
#   data = stan_data, seed = 42, chains = 1, parallel_chains = 1,
#   iter_warmup = 200, iter_sampling = 50, adapt_delta = 0.90, refresh = 50
# )
# print(fit_test$summary(c("lambda_N","tau_hs","sigma_GP_org","rho_org","sigma_GP_P","rho_P")))
# cat("Divergenze test:", sum(fit_test$sampler_diagnostics()[,,"divergent__"]), "\n")
# rm(fit_test); gc()


# ── 8. MCMC: run o carica da file ────────────────────────────────────────────────

fit_path <- here("stan", "fit_m4rr_gp.rds")

if (!file.exists(fit_path)) {

  cat("\nAvvio MCMC (4 catene x 3000 warmup + 5000 sampling)...\n")
  cat("Aggiornamenti ogni 200 iterazioni. Stima: 20-60 min.\n\n")

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

} else {
  cat("\nCarico fit salvato da:", fit_path, "\n")
  fit <- readRDS(fit_path)
}

rm(stan_data); gc()


# ── 9. ANALISI: calcola o carica da file ─────────────────────────────────────────

risultati_path <- here("stan", "risultati_bayesiano.RData")

if (!file.exists(risultati_path)) {

  # ── 9a. Diagnostiche ─────────────────────────────────────────────────────────
  cat("\nCalcolo diagnostiche...\n")

  smry_full <- fit$summary()
  np        <- fit$sampler_diagnostics()

  diag_list <- list(
    bad_rhat = smry_full |> filter(rhat > 1.05) |> arrange(desc(rhat)),
    bad_ess  = smry_full |> filter(ess_bulk < 400) |> arrange(ess_bulk),
    n_divs   = sum(np[, , "divergent__"]),
    max_td   = max(np[, , "treedepth__"]),
    ebfmi    = apply(np[, , "energy__"], 2, function(e) {
      mean(diff(e)^2) / var(e)  # E-BFMI per catena
    })
  )

  cat("\n── Rhat > 1.05 ──\n");    print(diag_list$bad_rhat)
  cat("\n── ESS bulk < 400 ──\n"); print(head(diag_list$bad_ess, 20))
  cat("Divergenze:", diag_list$n_divs, "| Max treedepth:", diag_list$max_td, "\n")
  cat("E-BFMI per catena:", round(diag_list$ebfmi, 3), "\n")

  rm(smry_full, np); gc()

  # ── 9b. Summary parametri strutturali ───────────────────────────────────────
  cat("\nCalcolo summary parametri...\n")

  vars_key <- c(
    "lambda_N",
    "gamma_org",    # [1]=logBottom [2]=PH [3]=Texture1 [4]=Texture2 [5]=BulkDensity
    "gamma_P_B",    # [1]=OnFarm [2]=Irrigate [3]=Fertilised [4]=N_Natural
    "psi_W_org", "psi_W_P", "theta_W_SOC", "theta_W_N",
    "sigma_B_SOC", "sigma_B_N",
    "psi_B_org", "psi_B_P",
    "sigma_GP_org", "rho_org",
    "sigma_GP_P",   "rho_P",
    "tau_hs"
  )

  smry <- fit$summary(vars_key) |>
    select(variable, mean, median, sd, q5, q95, rhat, ess_bulk) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  cat("\n── Parametri strutturali ──\n")
  print(smry)

  # ── 9c. Factor scores spaziali ──────────────────────────────────────────────
  cat("\nEstrazione factor scores between...\n")

  org_draws <- fit$draws("eta_org_B_out", format = "matrix")
  P_draws   <- fit$draws("eta_P_B_out",   format = "matrix")

  spatial_df <- coords_km |>
    mutate(
      eta_org_mean = colMeans(org_draws),
      eta_org_sd   = apply(org_draws, 2, sd),
      eta_P_mean   = colMeans(P_draws),
      eta_P_sd     = apply(P_draws,   2, sd)
    )

  rm(org_draws, P_draws); gc()

  # ── 9d. Confronto frequentista-bayesiano ─────────────────────────────────────
  # Richiede fit4rr da 3. mimic.R — salta se non disponibile
  confronto_df <- NULL
  if (exists("fit4rr") && requireNamespace("lavaan", quietly = TRUE)) {
    library(lavaan)
    freq <- parameterEstimates(fit4rr) |>
      filter(op == "~") |>
      select(lhs, rhs, est_freq = est, se_freq = se, z_freq = z)

    gorg_draws <- fit$draws("gamma_org", format = "matrix")
    gPB_draws  <- fit$draws("gamma_P_B", format = "matrix")

    bayes <- bind_rows(
      tibble(lhs = "fertil_org_W",
             rhs = c("logBottom","PH","Texture1","Texture2","BulkDensity"),
             est_bay = colMeans(gorg_draws),
             sd_bay  = apply(gorg_draws, 2, sd),
             q5_bay  = apply(gorg_draws, 2, quantile, 0.05),
             q95_bay = apply(gorg_draws, 2, quantile, 0.95)),
      tibble(lhs = "fertil_P_B",
             rhs = c("OnFarm","Irrigate","Fertilised","N_Natural"),
             est_bay = colMeans(gPB_draws),
             sd_bay  = apply(gPB_draws, 2, sd),
             q5_bay  = apply(gPB_draws, 2, quantile, 0.05),
             q95_bay = apply(gPB_draws, 2, quantile, 0.95))
    )
    confronto_df <- bayes |> left_join(freq, by = c("lhs","rhs"))
    rm(gorg_draws, gPB_draws, freq, bayes); gc()

    cat("\n── Confronto frequentista vs bayesiano ──\n")
    print(confronto_df |> mutate(across(where(is.numeric), ~ round(.x, 3))))
  } else {
    cat("\nfit4rr non disponibile: esegui 3. mimic.R per il confronto.\n")
  }

  # ── Salva tutti i risultati ──────────────────────────────────────────────────
  save(diag_list, smry, spatial_df, confronto_df, file = risultati_path)
  cat("\nRisultati salvati in:", risultati_path, "\n")

} else {
  cat("\nCarico risultati salvati da:", risultati_path, "\n")
  load(risultati_path)
  cat("Oggetti caricati: diag_list, smry, spatial_df, confronto_df\n")
}


# ── 10. TRACE PLOTS ─────────────────────────────────────────────────────────────
# Richiede il fit in memoria. Decommenta per visualizzare.

# if (requireNamespace("bayesplot", quietly = TRUE)) {
#   params_key <- c("lambda_N", "gamma_org[1]", "gamma_org[2]", "gamma_org[3]",
#                   "gamma_org[4]", "gamma_org[5]",
#                   "sigma_GP_org", "rho_org", "sigma_GP_P", "rho_P", "tau_hs")
#   draws_key <- fit$draws(variables = params_key)
#   print(mcmc_trace(draws_key) + theme_minimal())
#   rm(draws_key); gc()
# }


# ── 11. PLOT SPAZIALE FACTOR SCORES ─────────────────────────────────────────────

if (exists("spatial_df") && nrow(spatial_df) > 0) {
  p_org <- ggplot(spatial_df, aes(x = x_km, y = y_km, fill = eta_org_mean)) +
    geom_point(shape = 21, size = 4, stroke = 0.4) +
    scale_fill_viridis_c(option = "D", name = "Fattore\norganico") +
    labs(title = "Fattore organico between (media posteriore)",
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
    print(p_org + p_P)
  } else {
    print(p_org); print(p_P)
  }
}
