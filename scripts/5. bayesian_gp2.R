# =============================================================================
# 5. bayesian_gp2.R
#
# Seconda versione del modello bayesiano MIMIC 2F con GP.
# Rispetto a scripts/4. bayesian_gp.R, incorpora tre correzioni motivate
# dai risultati del primo run (stan/risultati_bayesiano.RData):
#
# [1] lambda_N_B libero [TENTATIVO FALLITO, soprasseduto da correzione 4]
#     Il primo run ha restituito sigma_B_N = 2.03 >> sigma_B_SOC = 0.16,
#     segnalando che il loading within (lambda_N = 0.636) non funziona al
#     livello between. Tentativo: campionare lambda_N_B come parametro libero.
#     Risultato: stesso Rhat=1.53 / ESS=7 del run con lambda_N libero.
#     Causa: lambda_N_B * eta_org_B[j] e' un accoppiamento bilineare identico
#     a quello che aveva bloccato NUTS in precedenza.
#     Soluzione finale (correzione 4): lambda_N_B rimosso; alpha_N[j] i.i.d.
#
# [2] GP per P rimosso, sostituito da random intercept
#     Il primo run ha restituito rho_P (mediana) = 18 km > distanza massima
#     tra field (20.4 km): il GP per P degenerava verso una covarianza
#     costante. Sostituito con sigma_P_between * z_P_B (NCP, i.i.d.).
#
# [3] PH rimosso da X_W
#     Il horseshoe nel primo run ha shrunkato PH a -0.032 (CI attraversa zero),
#     confermando la restrizione frequentista. Nel modello bayesiano non
#     c'e' il vincolo tecnico dell'LRT, quindi si rimuove direttamente.
#     K_W passa da 5 a 4: logBottom, Texture1, Texture2, BulkDensity.
#
# [4] Horseshoe rimosso da gamma_org; prior normal(0,1) per tutti (modifica 5)
#     Con K_W=4 e tutti e quattro i predittori meccanisticamente plausibili,
#     il horseshoe era sovrabbondante. La collinearita' BulkDensity-Texture2
#     (r=0.552) portava il horseshoe a shrinkare BulkDensity (-0.107, CI
#     [-0.251, +0.005]) molto piu' del frequentista (z=-4.01). Prior normal(0,1)
#     allinea tutti i gamma_org allo stesso prior. Rimossi tau0, slab_scale,
#     slab_df da stan_data.
#
# NOTA — N_Natural come predittore di P between (gamma_P_B[4]):
#     Stima bayesiana: -2.58 (CI [-3.03, -2.12]); frequentista: -0.309.
#     Discrepanza spuria: N_Natural=0 conta solo 5 field, tutti nel cluster SW
#     (distanza media intra-gruppo 2.3 km). La differenza grezza di logP tra
#     gruppi e' 0.11 log-unita' (~12%), non 2.58. Il random intercept i.i.d.
#     non assorbe la struttura spaziale del cluster, e l'effetto geografico
#     "SW ha P piu' alto" viene attribuito a N_Natural. Con J=5 controlli tutti
#     in un angolo, gamma_P_B[4] non e' interpretabile causalmente.
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

# Modello frequentista (per confronto in sezione 10d)
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

# PH rimosso da X_W (modifica 3): K_W = 4
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


# ── 4. LAMBDA_N FISSO (loading within) ───────────────────────────────────────
# lambda_N e' il loading within di logN sul fattore organico.
# Rimane fisso al MLE di M4rr per evitare l'accoppiamento bilineare
# lambda_N * eta_org_B che nel primo run causava ESS~7 e Rhat~1.53.
# lambda_N_B (between) e' stato rimosso (correzione 4): stesso problema.

if (exists("fit4rr")) {
  lambda_N_fixed <- parameterEstimates(fit4rr) |>
    filter(label == "a", op == "=~") |>
    slice(1) |>
    pull(est)
} else {
  lambda_N_fixed <- 0.6355  # MLE da M4rr (lavaan MLR), SE = 0.059
}
cat(sprintf("lambda_N within fisso: %.4f\n", lambda_N_fixed))


# ── 5. DIMENSIONE K_W (modifica 5: horseshoe rimosso, sezione semplificata) ───
# tau0, slab_scale, slab_df non piu' necessari: horseshoe sostituito da
# normal(0,1) per tutti i gamma_org (vedi modifica 5 nell'header e nel .stan).

K_W <- length(X_W_cols)   # 4: logBottom, Texture1, Texture2, BulkDensity
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
  # tau0, slab_scale, slab_df rimossi (modifica 5: horseshoe rimosso)
  lambda_N   = lambda_N_fixed    # loading within fisso (vedi sezione 4)
)

rm(X_W, X_B, dist_mat); gc()


# ── 7. COMPILAZIONE MODELLO ───────────────────────────────────────────────────

stan_file <- here("stan", "m4rr_gp2.stan")
cat("\nCompilazione del modello Stan (usa cache se gia' compilato)...\n")
mod <- cmdstan_model(stan_file, compile = TRUE)
cat("OK. CmdStan version:", cmdstan_version(), "\n")


# ── 8. TEST RUN ───────────────────────────────────────────────────────────────
# Decommenta per verificare che il modello giri prima del run completo.

# fit_test <- mod$sample(
#   data = stan_data, seed = 42, chains = 1, parallel_chains = 1,
#   iter_warmup = 200, iter_sampling = 50, adapt_delta = 0.90, refresh = 50
# )
# print(fit_test$summary(c("gamma_org","sigma_GP_org","rho_org",
#                           "sigma_P_between","sigma_B_N","sigma_B_SOC")))
# cat("Divergenze test:", sum(fit_test$sampler_diagnostics()[,,"divergent__"]), "\n")
# rm(fit_test); gc()


# ── 9. MCMC: run o carica da file ─────────────────────────────────────────────

fit_path <- here("stan", "fit_m4rr_gp2.rds")

if (!file.exists(fit_path)) {

  cat("\nAvvio MCMC (4 catene x 3000 warmup + 5000 sampling)...\n")
  cat("Aggiornamenti ogni 200 iterazioni. Stima: 15-50 min.\n\n")

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

  # Elimina CSV e JSON delle catene (ridondanti dopo save_object)
  csv_files  <- list.files(here("stan"), pattern = "\\.csv$",  full.names = TRUE)
  json_files <- list.files(here("stan"), pattern = "\\.json$", full.names = TRUE)
  invisible(file.remove(c(csv_files, json_files)))
  cat("CSV e JSON delle catene eliminati.\n")

} else {
  cat("\nCarico fit salvato da:", fit_path, "\n")
  fit <- readRDS(fit_path)
}

rm(stan_data); gc()


# ── 10. ANALISI: calcola o carica da file ─────────────────────────────────────

risultati_path <- here("stan", "risultati_bayesiano2.RData")

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

  # lambda_N (within) e' fisso; lambda_N_B rimosso (correzione 4)
  cat(sprintf("\nlambda_N within fisso: %.4f (MLE M4rr, SE=0.059)\n", lambda_N_fixed))

  vars_key <- c(
    # lambda_N_B rimosso (correzione 4); tau_hs rimosso (modifica 5)
    "gamma_org",    # [1]=logBottom [2]=Texture1 [3]=Texture2 [4]=BulkDensity
    "gamma_P_B",    # [1]=OnFarm [2]=Irrigate [3]=Fertilised [4]=N_Natural (*)
    # (*) gamma_P_B[4] (N_Natural) non interpretabile: confounding geografico,
    #     vedi NOTA nell'header dello script e nel file .stan
    "psi_W_org", "psi_W_P", "theta_W_SOC", "theta_W_N",
    "sigma_B_SOC", "sigma_B_N",          # sigma_B_N = residuo i.i.d. dopo GP N
    "psi_B_org", "psi_B_N",              # nugget dei due GP
    "sigma_P_between",
    "sigma_GP_org", "rho_org",           # GP SOC between
    "sigma_GP_N",   "rho_N"              # GP N between (Moran's I = 0.65, p<.0001)
  )

  smry <- fit$summary(vars_key) |>
    select(variable, mean, median, sd, q5, q95, rhat, ess_bulk) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  cat("\n── Parametri strutturali ──\n")
  print(smry, n = 30)

  # ── 10c. Factor scores spaziali ──────────────────────────────────────────
  cat("\nEstrazione factor scores between...\n")

  org_draws <- fit$draws("eta_org_B_out", format = "matrix")
  N_draws   <- fit$draws("eta_N_B_out",   format = "matrix")
  P_draws   <- fit$draws("eta_P_B_out",   format = "matrix")

  spatial_df <- coords_km |>
    mutate(
      eta_org_mean = colMeans(org_draws),
      eta_org_sd   = apply(org_draws, 2, sd),
      eta_N_mean   = colMeans(N_draws),
      eta_N_sd     = apply(N_draws,   2, sd),
      eta_P_mean   = colMeans(P_draws),
      eta_P_sd     = apply(P_draws,   2, sd)
    )

  rm(org_draws, N_draws, P_draws); gc()

  # ── 10d. Confronto frequentista-bayesiano ─────────────────────────────────
  # Nota: il confronto riguarda solo gamma_org (4 predittori, PH escluso) e
  # gamma_P_B. Il frequentista ha PH con est=0; non compare nel join bayesiano.
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
             rhs = X_W_cols,          # 4 predittori (no PH)
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

  # ── Salva ────────────────────────────────────────────────────────────────
  save(diag_list, smry, spatial_df, confronto_df, file = risultati_path)
  cat("\nRisultati salvati in:", risultati_path, "\n")

} else {
  cat("\nCarico risultati salvati da:", risultati_path, "\n")
  load(risultati_path)
  cat("Oggetti caricati: diag_list, smry, spatial_df, confronto_df\n")
}


# ── 11. TRACE PLOTS ───────────────────────────────────────────────────────────
# Decommenta per visualizzare. lambda_N_B e' stato rimosso (correzione 4).

if (requireNamespace("bayesplot", quietly = TRUE)) {
 params_key <- c("gamma_org[1]", "gamma_org[2]",
                 "gamma_org[3]", "gamma_org[4]",
                 "sigma_GP_org", "rho_org",
                 "sigma_GP_N",   "rho_N",
                 "sigma_P_between")
 # tau_hs rimosso: horseshoe non piu' nel modello (modifica 5)
 draws_key <- fit$draws(variables = params_key)
 print(mcmc_trace(draws_key) + theme_minimal())
 rm(draws_key); gc()
}


# ── 12. PLOT SPAZIALE FACTOR SCORES ──────────────────────────────────────────

if (exists("spatial_df") && nrow(spatial_df) > 0) {
  p_org <- ggplot(spatial_df, aes(x = x_km, y = y_km, fill = eta_org_mean)) +
    geom_point(shape = 21, size = 4, stroke = 0.4) +
    scale_fill_viridis_c(option = "D", name = "Fattore\norganico") +
    labs(title = "Fattore organico between (SOC, media posteriore)",
         x = "Est (km)", y = "Nord (km)") +
    theme_minimal()

  p_N <- ggplot(spatial_df, aes(x = x_km, y = y_km, fill = eta_N_mean)) +
    geom_point(shape = 21, size = 4, stroke = 0.4) +
    scale_fill_viridis_c(option = "A", name = "GP N") +
    labs(title = "Fattore N between (GP separato, media posteriore)",
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
