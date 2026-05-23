# =============================================================================
# 10. gp8_random_slope.R
#
# Estensione di gp7: random slope implicita per il profilo verticale di SOC.
#
# MOTIVAZIONE (da script 9):
#   Per logSOC, il random slope su logBottom migliora l'AIC di 12.3 punti
#   (LRT p=0.0003), ma il modello frequentista è singolare (Corr=1.0).
#   La singolarità indica che la varianza della slope è interamente spiegata
#   dall'intercetta: slope_j ≈ b * eta_org_B_j.
#
# STRUTTURA DEL MODELLO (gp8):
#   Come gp7, con una sola modifica: eta_org_B[j] entra nella media di logSOC
#   come termine depth-dipendente invece di semplice intercetta:
#
#     mu_SOC[n] = eta_org_B[j] * (1 + b * logBottom_n) + alpha_SOC_res[j] + gamma_org' * X_W[n]
#
#   dove alpha_SOC_res[j] = sigma_B_SOC * raw_SOC[j] (residuo i.i.d. di SOC).
#   logN è invariato: la modifica riguarda solo come eta_org_B entra in SOC.
#   b è fissato al valore stimato nella SEZIONE 1 (no bilineare in Stan).
#
# SEZIONE 1 — stima di b tramite LMM bayesiano (brms, senza GP):
#   logSOC ~ logBottom + Texture1 + Texture2 + BulkDensity +
#            OnFarm + Irrigate + Fertilised + N_Natural + (1 + logBottom | Field)
#   Prior LKJ(2) sulla correlazione intercetta-slope.
#   b = rho * sigma_slope / sigma_int  (coefficiente di regressione slope su int.)
#
# SEZIONE 2 — Stan gp8:
#   Stesso setup di gp7; b passato come dato (data); file stan/m4rr_gp8.stan.
#
# RISULTATI run gp8:
#   b_hat (mediana posteriore brms) = 0.467, CI 90% [0.21, 0.70] — chiaramente != 0.
#   Interpretazione: campi con alto SOC medio (eta_org_B > 0) hanno profilo di
#   decadimento con la profondità meno ripido. Biologicamente plausibile (SOM
#   più uniformemente distribuita nel profilo in suoli organicamente ricchi).
#
#   Diagnostica Stan (4 catene × 5000 sampling + 3000 warmup):
#     Divergenze: 0 | Rhat > 1.05: 0 | ESS bulk min: ~3500
#     E-BFMI: 0.68-0.71 (ok, > 0.3) | Max treedepth: 10
#
#   Confronto LOO-CV gp7 vs gp8:
#     gp7: elpd_loo = -520.3 (SE 35.3)
#     gp8: elpd_loo = -509.0 (SE 34.9)
#     DELTA_ELPD = +11.3 (SE 5.6) - rapporto 2.0 SE: miglioramento predittivo reale.
#     p_loo = 122 vs 114: +8 parametri effettivi (struttura random slope).
#
#   Parametri chiave modificati rispetto a gp7:
#     theta_W_SOC:  0.302 -> 0.261 (-14%): meno varianza within SOC inspiegata
#     sigma_B_SOC:  0.156 -> 0.100 (-36%): meno residuo between SOC i.i.d.
#     psi_B_org:    0.076 -> 0.133 (+75%): nugget GP riassorbe varianza liberata
#     sigma_GP_org: 0.234 -> 0.182 (-22%): ampiezza GP lievemente piu' piccola
#     rho_org:      5.12  -> 5.21  km     — invariato (struttura spaziale stabile)
#     Tutti gli altri parametri (gamma_org, gamma_N_B, gamma_P_B, N e P): stabili.
#
#   Conclusione: gp8 e' il modello finale. Migliora gp7 catturando la
#   proporzionalita' tra intercetta e slope del profilo verticale di SOC.
# =============================================================================


# ── 0. SETUP ─────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)

if (!requireNamespace("brms",     quietly = TRUE)) install.packages("brms")
if (!requireNamespace("cmdstanr", quietly = TRUE))
  install.packages("cmdstanr",
                   repos = c("https://mc-stan.org/r-packages/",
                             "https://cloud.r-project.org"))
library(brms)
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
if (requireNamespace("lavaan",     quietly = TRUE)) library(lavaan)


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

fit4rr_path <- here("stan", "fit4rr_freq.rds")
if (file.exists(fit4rr_path)) {
  fit4rr <- readRDS(fit4rr_path)
} else {
  message("fit4rr non trovato: esegui 3. mimic.R per generarlo.")
}

lambda_N_fixed <- if (exists("fit4rr")) {
  parameterEstimates(fit4rr) |>
    filter(label == "a", op == "=~") |>
    slice(1) |>
    pull(est)
} else { 0.6355 }
cat(sprintf("lambda_N within fisso: %.4f\n", lambda_N_fixed))


# ══════════════════════════════════════════════════════════════════════════════
# SEZIONE 1 — LMM BAYESIANO PER STIMA DI b  (brms, senza GP)
# ══════════════════════════════════════════════════════════════════════════════
#
# Modello: logSOC ~ logBottom + Texture1 + Texture2 + BulkDensity +
#                   OnFarm + Irrigate + Fertilised + N_Natural + (1 + logBottom | Field)
#
# Non include PH (coerente con M4rr e gp7).
# Non include la struttura spaziale (GP): sostituita dal random intercept i.i.d.
# Questo introduce un errore di approssimazione: l'intercetta casuale in brms
# cattura tutta la variazione between-Field (GP + nugget + residuo), mentre in
# gp8 solo il GP e il nugget compongono eta_org_B. L'approssimazione è
# accettabile perché sigma_B_SOC << sigma_GP_org (da gp7: 0.16 vs > 1), quindi
# l'intercetta i.i.d. in brms ≈ eta_org_B in gp7.
#
# Prior:
#   LKJ(2): favorisce |rho| < 1, evita la singolarità frequentista
#   normal(0, 1) per i fixed effects
#   exponential(1) per le SD (half-prior su scala log)
#
# b = rho * sigma_slope / sigma_int  (regressione OLS slope su intercetta)
# Intervallo credibile al 90% per decidere se b è distinguibile da zero.

brms_path <- here("stan", "fit_lmm_b_slope.rds")

if (!file.exists(brms_path)) {

  cat("\n── Sezione 1: LMM bayesiano per stima di b ──────────────────────\n")

  prior_lmm <- c(
    prior(normal(0, 1),    class = b),
    prior(normal(0, 1),    class = Intercept),
    prior(exponential(1),  class = sd),
    prior(exponential(1),  class = sigma),
    prior(lkj(2),          class = cor)     # tira rho lontano da ±1
  )

  fit_lmm <- brm(
    logSOC ~ logBottom + Texture1 + Texture2 + BulkDensity +
             OnFarm + Irrigate + Fertilised + N_Natural +
             (1 + logBottom | Field),
    data    = dati,
    prior   = prior_lmm,
    chains  = 4, cores = 4,
    iter    = 4000, warmup = 2000,
    seed    = 2024,
    backend = "cmdstanr",
    file    = brms_path
  )

  saveRDS(fit_lmm, brms_path)
  cat("LMM salvato in:", brms_path, "\n")

} else {
  cat("\nCarico LMM bayesiano da:", brms_path, "\n")
  fit_lmm <- readRDS(brms_path)
}


# ── 1a. Estrazione e calcolo di b ─────────────────────────────────────────────

cat("\n── Summary LMM (varianze e correlazione) ──────────────────────────\n")
print(summary(fit_lmm))

# Estrai draws di sigma_int, sigma_slope, rho
draws_lmm <- as_draws_df(fit_lmm)

# Nomi standard brms: sd_Field__Intercept, sd_Field__logBottom, cor_Field__Intercept__logBottom
sigma_int  <- draws_lmm[["sd_Field__Intercept"]]
sigma_slp  <- draws_lmm[["sd_Field__logBottom"]]
rho        <- draws_lmm[["cor_Field__Intercept__logBottom"]]

# b = coefficiente di regressione slope ~ intercetta nel modello bivariato
b_draws <- rho * sigma_slp / sigma_int

cat("\n── Posterior di b (slope_j ≈ b * intercetta_j) ─────────────────────\n")
cat(sprintf("  Media:   %.3f\n", mean(b_draws)))
cat(sprintf("  Mediana: %.3f\n", median(b_draws)))
cat(sprintf("  SD:      %.3f\n", sd(b_draws)))
cat(sprintf("  CI 90%%: [%.3f, %.3f]\n",
            quantile(b_draws, 0.05), quantile(b_draws, 0.95)))
cat(sprintf("  P(b > 0): %.3f\n", mean(b_draws > 0)))

# Scegli b_hat come mediana (più robusta della media con distribuzioni asimmetriche)
b_hat <- median(b_draws)
cat(sprintf("\nb_hat (mediana posteriore) = %.4f  — usato in gp8 come dato fisso\n",
            b_hat))

# Plot della distribuzione posteriore di b
df_b <- tibble(b = b_draws)
p_b <- ggplot(df_b, aes(x = b)) +
  geom_histogram(bins = 60, fill = "steelblue", alpha = 0.8, colour = "white") +
  geom_vline(xintercept = b_hat, colour = "firebrick", linewidth = 1,
             linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray40", linewidth = 0.7) +
  labs(
    title    = "Posterior di b (slope proporzionale all'intercetta casuale)",
    subtitle = sprintf("Mediana = %.3f  |  CI 90%% [%.3f, %.3f]",
                       b_hat,
                       quantile(b_draws, 0.05),
                       quantile(b_draws, 0.95)),
    x = "b", y = "Conteggio"
  ) +
  theme_minimal()
print(p_b)

# Plot diagnostico: scatter (random intercept, random slope) da brms
re_lmm <- ranef(fit_lmm)$Field[, "Estimate", ] |>
  as.data.frame() |>
  rownames_to_column("Field") |>
  rename(u_int = Intercept, u_slp = logBottom)

p_re <- ggplot(re_lmm, aes(x = u_int, y = u_slp)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "firebrick", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray50") +
  labs(
    title    = "Random effects LMM bayesiano: intercetta vs slope",
    subtitle = sprintf("Pendenza OLS = %.3f  (≈ b_hat per il modello lineare)", b_hat),
    x = "Random intercetta (Field)", y = "Random slope di logBottom (Field)"
  ) +
  theme_minimal()
print(p_re)


# ══════════════════════════════════════════════════════════════════════════════
# SEZIONE 2 — STAN GP8  (gp7 + modifica random slope implicita)
# ══════════════════════════════════════════════════════════════════════════════

crop_raw     <- readRDS(here("data", "crop.rds"))
coords_field <- crop_raw |>
  group_by(Field) |>
  summarise(Lat  = mean(Lat,  na.rm = TRUE),
            Long = mean(Long, na.rm = TRUE)) |>
  arrange(Field)
rm(crop_raw); gc()

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

ref_lat  <- mean(coords_field$Lat,  na.rm = TRUE)
ref_long <- mean(coords_field$Long, na.rm = TRUE)

coords_km <- coords_field |>
  arrange(Field) |>
  mutate(
    x_km = (Long - ref_long) * 111.32 * cos(ref_lat * pi / 180),
    y_km = (Lat  - ref_lat)  * 110.54
  )

dist_mat <- as.matrix(dist(coords_km[, c("x_km", "y_km")]))

stan_data_gp8 <- list(
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
  lambda_N = lambda_N_fixed,
  b_slope  = b_hat           # fissato dalla Sezione 1
)

cat(sprintf("\nb_slope passato a Stan: %.4f\n", b_hat))
rm(X_W, X_B, dist_mat); gc()


# ── 2a. Compilazione ──────────────────────────────────────────────────────────

stan_file_gp8 <- here("stan", "m4rr_gp8.stan")
cat("\nCompilazione gp8...\n")
mod_gp8 <- cmdstan_model(stan_file_gp8, compile = TRUE)
cat("OK. CmdStan version:", cmdstan_version(), "\n")


# ── 2b. MCMC ─────────────────────────────────────────────────────────────────

fit_path_gp8 <- here("stan", "fit_m4rr_gp8.rds")

if (!file.exists(fit_path_gp8)) {

  cat("\nAvvio MCMC gp8 (4 catene × 5000 sampling + 3000 warmup)...\n")

  fit_gp8 <- mod_gp8$sample(
    data            = stan_data_gp8,
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

  fit_gp8$save_object(fit_path_gp8)
  cat("Fit salvato in:", fit_path_gp8, "\n")

  csv_files  <- list.files(here("stan"), pattern = "\\.csv$",  full.names = TRUE)
  json_files <- list.files(here("stan"), pattern = "\\.json$", full.names = TRUE)
  invisible(file.remove(c(csv_files, json_files)))

} else {
  cat("\nCarico gp8 da:", fit_path_gp8, "\n")
  fit_gp8 <- readRDS(fit_path_gp8)
}

rm(stan_data_gp8); gc()


# ── 2c. Diagnostica ──────────────────────────────────────────────────────────

risultati_path_gp8 <- here("stan", "risultati_gp8.RData")

if (!file.exists(risultati_path_gp8)) {

  smry_full <- fit_gp8$summary()
  np        <- fit_gp8$sampler_diagnostics()

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


  # ── 2d. Summary parametri e confronto con gp7 ───────────────────────────────
  cat("\n── Parametri strutturali gp8 ───────────────────────────────────────\n")
  cat(sprintf("b_slope fisso: %.4f\n", b_hat))

  vars_key <- c(
    "gamma_org",
    "gamma_P_B", "gamma_N_B",
    "psi_W_org", "psi_W_P", "theta_W_SOC", "theta_W_N",
    "sigma_B_SOC", "sigma_B_N",
    "psi_B_org", "sigma_P_between",
    "sigma_GP_org", "rho_org"
  )

  smry_gp8 <- fit_gp8$summary(vars_key) |>
    select(variable, mean, median, sd, q5, q95, rhat, ess_bulk) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  cat("\n── Parametri gp8 ──\n")
  print(smry_gp8, n = 35)

  # Confronto gp7 vs gp8: carica risultati gp7 se disponibili
  gp7_path <- here("stan", "risultati_bayesiano7.RData")
  if (file.exists(gp7_path)) {
    load(gp7_path)   # carica: diag_list, smry (di gp7), spatial_df, confronto_df
    smry_gp7 <- smry  # rinomina per chiarezza

    cat("\n── Confronto gp7 vs gp8 (parametri comuni) ────────────────────────\n")
    confronto_gp7_gp8 <- smry_gp7 |>
      select(variable, mean_gp7 = mean, sd_gp7 = sd) |>
      inner_join(
        smry_gp8 |> select(variable, mean_gp8 = mean, sd_gp8 = sd),
        by = "variable"
      ) |>
      mutate(delta = round(mean_gp8 - mean_gp7, 3))
    print(as.data.frame(confronto_gp7_gp8))
  }


  # ── 2e. ELPD (confronto predittivo) ────────────────────────────────────────
  # Se il pacchetto loo è disponibile, confronta gp7 e gp8 con LOO-CV.
  # log_lik deve essere nella generated_quantities di entrambi i modelli.

  if (requireNamespace("loo", quietly = TRUE)) {
    library(loo)
    ll_gp8 <- fit_gp8$draws("log_lik", format = "matrix")
    loo_gp8 <- loo(ll_gp8)
    cat("\n── LOO-CV gp8 ────────────────────────────────────────────────────\n")
    print(loo_gp8)

    gp7_fit_path <- here("stan", "fit_m4rr_gp7.rds")
    if (file.exists(gp7_fit_path)) {
      fit_gp7 <- readRDS(gp7_fit_path)
      ll_gp7  <- fit_gp7$draws("log_lik", format = "matrix")
      loo_gp7 <- loo(ll_gp7)
      cat("\n── LOO-CV gp7 ────────────────────────────────────────────────────\n")
      print(loo_gp7)
      cat("\n── Confronto LOO gp7 vs gp8 ──────────────────────────────────────\n")
      print(loo_compare(loo_gp7, loo_gp8))
      rm(fit_gp7, ll_gp7); gc()
    }
    rm(ll_gp8); gc()
  } else {
    cat("\nInstalla 'loo' per il confronto predittivo.\n")
  }


  # ── 2f. Factor scores spaziali ──────────────────────────────────────────────
  org_draws <- fit_gp8$draws("eta_org_B_out", format = "matrix")
  P_draws   <- fit_gp8$draws("eta_P_B_out",   format = "matrix")
  N_draws   <- fit_gp8$draws("alpha_N_out",   format = "matrix")

  spatial_df_gp8 <- coords_km |>
    mutate(
      eta_org_mean = colMeans(org_draws),
      eta_org_sd   = apply(org_draws, 2, sd),
      alpha_N_mean = colMeans(N_draws),
      eta_P_mean   = colMeans(P_draws),
      eta_P_sd     = apply(P_draws, 2, sd)
    )

  rm(org_draws, P_draws, N_draws); gc()

  save(diag_list, smry_gp8, spatial_df_gp8, file = risultati_path_gp8)
  cat("\nRisultati gp8 salvati in:", risultati_path_gp8, "\n")

} else {
  cat("\nCarico risultati gp8 da:", risultati_path_gp8, "\n")
  load(risultati_path_gp8)
}


# ── 2g. Trace plots ───────────────────────────────────────────────────────────

if (requireNamespace("bayesplot", quietly = TRUE)) {
  params_key <- c("gamma_org[1]", "gamma_org[2]",
                  "gamma_org[3]", "gamma_org[4]",
                  "sigma_GP_org", "rho_org",
                  "sigma_B_SOC",  "sigma_B_N",
                  "sigma_P_between")
  draws_key <- fit_gp8$draws(variables = params_key)
  print(mcmc_trace(draws_key) + theme_minimal())
  rm(draws_key); gc()
}


# ── 2h. Plot spaziale factor scores ──────────────────────────────────────────

if (exists("spatial_df_gp8") && nrow(spatial_df_gp8) > 0) {

  p_org <- ggplot(spatial_df_gp8,
                  aes(x = x_km, y = y_km, fill = eta_org_mean)) +
    geom_point(shape = 21, size = 4, stroke = 0.4) +
    scale_fill_viridis_c(option = "D", name = "Fattore\norganico") +
    labs(title    = "Fattore organico between (gp8, GP Matern 3/2)",
         subtitle = sprintf("b_slope = %.3f incluso nel profilo verticale", b_hat),
         x = "Est (km)", y = "Nord (km)") +
    theme_minimal()

  p_P <- ggplot(spatial_df_gp8,
                aes(x = x_km, y = y_km, fill = eta_P_mean)) +
    geom_point(shape = 21, size = 4, stroke = 0.4) +
    scale_fill_viridis_c(option = "C", name = "Fattore P") +
    labs(title = "Fattore fosforo between (gp8)",
         x = "Est (km)", y = "Nord (km)") +
    theme_minimal()

  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    print(p_org + p_P)
  } else {
    print(p_org); print(p_P)
  }
}

cat("\n── Fine script 10 ───────────────────────────────────────────────────\n")
