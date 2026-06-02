# =============================================================================
# 14_confronto_modelli.R  —  Confronto modelli A / B / M-SP (LOO + parametri)
#
# Stima e confronta tre versioni semplificate del modello M-SP:
#
#   Opzione A (m_final_A.stan):
#     Stessa struttura di M-SP. Rimossi X_B (management) e Texture1.
#     K_W = 3: logBottom (col 1), Texture2 (col 2), BulkDensity (col 3).
#     Struttura proporzionale per tutte e 3 le risposte.
#
#   Opzione B (m_final_B.stan):
#     Struttura response-specific fedele a projpred:
#       SOC: slope prop. + {logBottom, Texture2} fissi
#       N:   solo random intercept + {logBottom, Texture2, BulkDensity} fissi (no eta_N)
#       P:   solo slope prop., nessun effetto fisso
#
#   Opzione C (projpred projected posteriors):
#     Già calcolati in script 09 — confronto numerico da cache.
#
# Output:
#   stan/fit_final_A.rds              → fit Stan modello A
#   stan/fit_final_B.rds              → fit Stan modello B
#   output/cache/loo_all_models.rds   → lista LOO: M-SP, A, B
#   output/tables/tab_14_params.csv   → tabella confronto parametri
#   output/tables/tab_14_loo.csv      → tabella confronto LOO
#
# MCMC: 4 catene × 2000 warmup + 5000 sampling, adapt_delta=0.97, seed=2024
# (identici a M-SP per confronto LOO diretto)
# Dipende da: stan/fit_msp.rds, output/cache/proj_posteriors.rds, data/dati.rds
# =============================================================================


# ── 0. SETUP ──────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
library(cmdstanr)
library(posterior)
library(loo)
source(here("scripts", "00_utilities.R"))
setup_rtools()

dir.create(here("output", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "cache"),  recursive = TRUE, showWarnings = FALSE)
tab_dir   <- here("output", "tables")
cache_dir <- here("output", "cache")

cat("CmdStan version:", cmdstan_version(), "\n")

# Helper: summary con mediana e CI90%
# Usa la summary predefinita di posterior (q5, q95) e rinomina.
smry_draws <- function(fit, vars) {
  fit$summary(variables = vars) |>
    select(variable, median, q5, q95) |>
    rename(q05 = q5) |>
    mutate(across(c(median, q05, q95), ~round(.x, 3)))
}

# Helper: formatta una riga come "mediana [q05, q95]" o "—" se mancante
fmt_row <- function(smry_df, var_name) {
  row <- smry_df |> filter(variable == var_name)
  if (nrow(row) == 0) return("—")
  sprintf("%.3f [%.3f, %.3f]", row$median[1], row$q05[1], row$q95[1])
}

# Helper: estrae mediana e CI90% da una matrice di draws proiettati (projpred)
fmt_proj <- function(proj_obj, col_name) {
  if (is.null(proj_obj)) return("—")
  tryCatch({
    m <- as.matrix(proj_obj)
    if (!col_name %in% colnames(m)) return("—")
    x <- m[, col_name]
    sprintf("%.3f [%.3f, %.3f]", median(x), quantile(x, .05), quantile(x, .95))
  }, error = function(e) "—")
}


# ── 1. DATI ───────────────────────────────────────────────────────────────────

cat("\nPreparazione dati...\n")
dati <- readRDS(here("data", "dati.rds")) |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom)
  ) |>
  mutate(across(c(logBottom, Texture1, Texture2, BulkDensity),
                ~ c(scale(.x)))) |>
  mutate(Field = factor(Field))

field_levels <- sort(unique(as.integer(as.character(dati$Field))))
J <- length(field_levels)
N <- nrow(dati)

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))
cat(sprintf("N = %d | J = %d\n", N, J))

# Dati modello A: X_W = {logBottom[col1], Texture2[col2], BulkDensity[col3]}
# logBottom in col 1 è obbligatorio: m_final_A.stan usa X_W[n,1] per la slope proporzionale
X_W_A <- as.matrix(dati_int[, c("logBottom", "Texture2", "BulkDensity")])
stan_data_A <- list(
  N        = N,  J = J,  K_W = 3L,
  field_id = dati_int$field_int,
  logSOC   = dati_int$logSOC,
  logN     = dati_int$logN,
  logP     = dati_int$logP,
  X_W      = X_W_A
)
rm(X_W_A)

# Dati modello B: vettori separati (nomi espliciti, zero ambiguità)
stan_data_B <- list(
  N             = N,  J = J,
  field_id      = dati_int$field_int,
  logSOC        = dati_int$logSOC,
  logN          = dati_int$logN,
  logP          = dati_int$logP,
  logBottom_s   = dati_int$logBottom,
  Texture2_s    = dati_int$Texture2,
  BulkDensity_s = dati_int$BulkDensity
)
gc()


# ── 2. MODELLO A — FIT ────────────────────────────────────────────────────────

cat("\n═══ MODELLO A ═══════════════════════════════════════════════════\n")
fit_A_path <- here("stan", "fit_final_A.rds")

if (!file.exists(fit_A_path)) {
  cat("Compilazione m_final_A.stan...\n")
  mod_A <- cmdstan_model(here("stan", "m_final_A.stan"), compile = TRUE)
  cat("Compilazione OK.\n\n")
  fit_A <- mod_A$sample(
    data            = stan_data_A,
    seed            = 2024,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 3000,
    iter_sampling   = 5000,
    adapt_delta     = 0.97,
    max_treedepth   = 11,
    refresh         = 300,
    output_dir      = here("stan")
  )
  fit_A$save_object(fit_A_path)
  cat("\nFit A salvato in:", fit_A_path, "\n")
  invisible(file.remove(list.files(here("stan"), pattern = "\\.csv$", full.names = TRUE)))
} else {
  cat("Carico fit A da:", fit_A_path, "\n")
  fit_A <- readRDS(fit_A_path)
}


# ── 3. MODELLO A — DIAGNOSTICA ────────────────────────────────────────────────

cat("\n--- Diagnostica modello A ---\n")
np_A     <- fit_A$sampler_diagnostics()
n_divs_A <- sum(np_A[, , "divergent__"])
max_td_A <- max(np_A[, , "treedepth__"])
smry_A_full  <- fit_A$summary()
bad_rhat_A   <- sum(smry_A_full$rhat > 1.05, na.rm = TRUE)
bad_ess_A    <- sum(smry_A_full$ess_bulk < 400, na.rm = TRUE)
cat(sprintf("  Divergenze: %d | Max treedepth: %d | Rhat>1.05: %d | ESS<400: %d\n",
            n_divs_A, max_td_A, bad_rhat_A, bad_ess_A))
if (n_divs_A > 5) warning("ATTENZIONE: molte divergenze in modello A")

cat("\nParametri strutturali A:\n")
smry_A <- smry_draws(fit_A, c(
  "alpha_SOC", "psi_SOC", "eta_SOC", "b_SOC", "sigma_SOC",
  "alpha_N",   "psi_N",   "eta_N",   "b_N",   "sigma_N",
  "alpha_P",   "psi_P",   "eta_P",   "b_P",   "sigma_P",
  "gamma_SOC[1]", "gamma_SOC[2]", "gamma_SOC[3]",   # logBottom, Texture2, BulkDensity
  "gamma_N[1]",   "gamma_N[2]",   "gamma_N[3]",
  "gamma_P[1]",   "gamma_P[2]",   "gamma_P[3]"
))
print(smry_A, n = 30)
rm(smry_A_full); gc()


# ── 4. MODELLO B — FIT ────────────────────────────────────────────────────────

cat("\n═══ MODELLO B ═══════════════════════════════════════════════════\n")
fit_B_path <- here("stan", "fit_final_B.rds")

if (!file.exists(fit_B_path)) {
  cat("Compilazione m_final_B.stan...\n")
  mod_B <- cmdstan_model(here("stan", "m_final_B.stan"), compile = TRUE)
  cat("Compilazione OK.\n\n")
  fit_B <- mod_B$sample(
    data            = stan_data_B,
    seed            = 2024,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 3000,
    iter_sampling   = 5000,
    adapt_delta     = 0.97,
    max_treedepth   = 11,
    refresh         = 300,
    output_dir      = here("stan")
  )
  fit_B$save_object(fit_B_path)
  cat("\nFit B salvato in:", fit_B_path, "\n")
  invisible(file.remove(list.files(here("stan"), pattern = "\\.csv$", full.names = TRUE)))
} else {
  cat("Carico fit B da:", fit_B_path, "\n")
  fit_B <- readRDS(fit_B_path)
}


# ── 5. MODELLO B — DIAGNOSTICA ────────────────────────────────────────────────

cat("\n--- Diagnostica modello B ---\n")
np_B     <- fit_B$sampler_diagnostics()
n_divs_B <- sum(np_B[, , "divergent__"])
max_td_B <- max(np_B[, , "treedepth__"])
smry_B_full  <- fit_B$summary()
bad_rhat_B   <- sum(smry_B_full$rhat > 1.05, na.rm = TRUE)
bad_ess_B    <- sum(smry_B_full$ess_bulk < 400, na.rm = TRUE)
cat(sprintf("  Divergenze: %d | Max treedepth: %d | Rhat>1.05: %d | ESS<400: %d\n",
            n_divs_B, max_td_B, bad_rhat_B, bad_ess_B))
if (n_divs_B > 5) warning("ATTENZIONE: molte divergenze in modello B")

cat("\nParametri strutturali B:\n")
# NOTA: eta_N non esiste in modello B — non richiederlo
smry_B <- smry_draws(fit_B, c(
  "alpha_SOC", "psi_SOC", "eta_SOC", "b_SOC", "sigma_SOC",
  "alpha_N",   "psi_N",               "sigma_N",    # no eta_N, no b_N
  "alpha_P",   "psi_P",   "eta_P",   "b_P",   "sigma_P",
  "gamma_SOC[1]", "gamma_SOC[2]",               # logBottom, Texture2
  "gamma_N[1]",   "gamma_N[2]",   "gamma_N[3]", # logBottom, Texture2, BulkDensity
  "gamma_P_logB"                                # logBottom fisso per P (scalare)
))
print(smry_B, n = 20)
rm(smry_B_full); gc()


# ── 6. M-SP: SUMMARY + LOO (caricato una sola volta) ─────────────────────────

cat("\n═══ M-SP (reference) ════════════════════════════════════════════\n")
loo_cache_path <- file.path(cache_dir, "loo_all_models.rds")

cat("Carico M-SP...\n")
fit_msp <- readRDS(here("stan", "fit_msp.rds"))

# Summary M-SP (sempre necessaria per la tabella di confronto)
# Indici gamma M-SP: [1]=logBottom, [2]=Texture1, [3]=Texture2, [4]=BulkDensity
smry_msp <- smry_draws(fit_msp, c(
  "alpha_SOC", "psi_SOC", "eta_SOC", "b_SOC", "sigma_SOC",
  "alpha_N",   "psi_N",   "eta_N",   "b_N",   "sigma_N",
  "alpha_P",   "psi_P",   "eta_P",   "b_P",   "sigma_P",
  "gamma_SOC[1]", "gamma_SOC[3]", "gamma_SOC[4]",  # logBottom, Texture2, BulkDensity
  "gamma_N[1]",   "gamma_N[3]",   "gamma_N[4]",
  "gamma_P[1]",   "gamma_P[3]",   "gamma_P[4]"
))
cat("Summary M-SP estratta.\n")

# LOO M-SP (solo se non c'è la cache)
if (!file.exists(loo_cache_path)) {
  cat("Calcolo LOO M-SP...\n")
  ll_msp  <- fit_msp$draws("log_lik", format = "matrix")
  loo_msp <- loo(ll_msp, cores = 4)
  rm(ll_msp); gc()
  cat("LOO M-SP calcolato.\n")
}

rm(fit_msp); gc()


# ── 7. LOO ────────────────────────────────────────────────────────────────────

cat("\n═══ LOO-CV ══════════════════════════════════════════════════════\n")

if (!file.exists(loo_cache_path)) {
  cat("Calcolo LOO modello A...\n")
  ll_A   <- fit_A$draws("log_lik", format = "matrix")
  loo_A  <- loo(ll_A, cores = 4)
  rm(ll_A); gc()

  cat("Calcolo LOO modello B...\n")
  ll_B   <- fit_B$draws("log_lik", format = "matrix")
  loo_B  <- loo(ll_B, cores = 4)
  rm(ll_B); gc()

  loo_list <- list(`M-SP` = loo_msp, `A (ridotto)` = loo_A, `B (resp-spec)` = loo_B)
  saveRDS(loo_list, loo_cache_path)
  cat("Cache LOO salvata:", loo_cache_path, "\n")
  rm(loo_msp, loo_A, loo_B); gc()
} else {
  cat("Carico LOO da cache:", loo_cache_path, "\n")
}

loo_list <- readRDS(loo_cache_path)

cat("\n--- LOO individuali ---\n")
for (nm in names(loo_list)) {
  cat(sprintf("\n%s:\n", nm)); print(loo_list[[nm]])
}

cat("\n--- Confronto LOO (migliore in cima) ---\n")
loo_comp <- loo_compare(loo_list)
print(loo_comp)

# Salva tabella LOO
loo_tab <- as.data.frame(loo_comp) |>
  tibble::rownames_to_column("modello") |>
  mutate(across(where(is.numeric), ~round(.x, 2)))
write.csv(loo_tab, file.path(tab_dir, "tab_14_loo.csv"), row.names = FALSE)
cat("  Salvato: output/tables/tab_14_loo.csv\n")


# ── 8. TABELLA CONFRONTO PARAMETRI ────────────────────────────────────────────

cat("\n═══ TABELLA CONFRONTO PARAMETRI ═════════════════════════════════\n")
cat("Legenda indici:\n")
cat("  M-SP:  gamma_r[1]=logBottom, [2]=Texture1, [3]=Texture2, [4]=BulkDensity\n")
cat("  A:     gamma_r[1]=logBottom, [2]=Texture2, [3]=BulkDensity\n")
cat("  B SOC: gamma_SOC[1]=logBottom, [2]=Texture2\n")
cat("  B N:   gamma_N[1]=logBottom, [2]=Texture2, [3]=BulkDensity\n")
cat("  B P:   gamma_P_logB (scalare) + slope proporzionale\n\n")

# Projected posteriors da cache
proj_cache <- file.path(cache_dir, "proj_posteriors.rds")
proj_list  <- if (file.exists(proj_cache)) {
  cat("Carico projected posteriors da:", proj_cache, "\n")
  readRDS(proj_cache)
} else {
  cat("AVVISO: proj_posteriors.rds non trovato — colonna C sarà vuota\n")
  list(SOC = NULL, N = NULL, P = NULL)
}

# Funzione che costruisce le righe per una risposta
build_resp_rows <- function(
    resp,
    msp_logB, msp_tex2, msp_bd,   # nomi parametro gamma in M-SP
    A_logB,   A_tex2,   A_bd,     # nomi parametro gamma in A
    B_logB,   B_tex2,   B_bd,     # nomi parametro gamma in B (NA se assente)
    proj_obj                       # projected posteriors per questa risposta
) {
  r <- resp
  rows <- tribble(
    ~Parametro,                 ~`M-SP`,                            ~`A (ridotto)`,                    ~`B (resp-spec)`,                                                              ~`C (projpred)`,
    paste0("alpha_",  r),       fmt_row(smry_msp, paste0("alpha_",  r)), fmt_row(smry_A, paste0("alpha_",  r)), fmt_row(smry_B, paste0("alpha_",  r)),                                     fmt_proj(proj_obj, "(Intercept)"),
    paste0("psi_",    r),       fmt_row(smry_msp, paste0("psi_",    r)), fmt_row(smry_A, paste0("psi_",    r)), fmt_row(smry_B, paste0("psi_",    r)),                                     "—",
    paste0("eta_",    r),       fmt_row(smry_msp, paste0("eta_",    r)), fmt_row(smry_A, paste0("eta_",    r)), if (r == "N") "— (no eta in B)" else fmt_row(smry_B, paste0("eta_", r)),  "—",
    paste0("sigma_",  r),       fmt_row(smry_msp, paste0("sigma_",  r)), fmt_row(smry_A, paste0("sigma_",  r)), fmt_row(smry_B, paste0("sigma_",  r)),                                     "—",
    paste0("γ_", r, " logBottom"), fmt_row(smry_msp, msp_logB),          fmt_row(smry_A, A_logB),              if (is.na(B_logB)) "—" else fmt_row(smry_B, B_logB),                       fmt_proj(proj_obj, "logBottom"),
    paste0("γ_", r, " Texture2"), fmt_row(smry_msp, msp_tex2),           fmt_row(smry_A, A_tex2),              if (is.na(B_tex2)) "—" else fmt_row(smry_B, B_tex2),                       fmt_proj(proj_obj, "Texture2"),
    paste0("γ_", r, " BulkDensity"), fmt_row(smry_msp, msp_bd),         fmt_row(smry_A, A_bd),               if (is.na(B_bd))   "—" else fmt_row(smry_B, B_bd),                         fmt_proj(proj_obj, "BulkDensity")
  )
  rows$Risposta <- resp
  rows
}

tab_SOC <- build_resp_rows("SOC",
  msp_logB = "gamma_SOC[1]", msp_tex2 = "gamma_SOC[3]", msp_bd = "gamma_SOC[4]",
  A_logB   = "gamma_SOC[1]", A_tex2   = "gamma_SOC[2]", A_bd   = "gamma_SOC[3]",
  B_logB   = "gamma_SOC[1]", B_tex2   = "gamma_SOC[2]", B_bd   = NA,
  proj_obj = proj_list$SOC
)

tab_N <- build_resp_rows("N",
  msp_logB = "gamma_N[1]",   msp_tex2 = "gamma_N[3]",   msp_bd = "gamma_N[4]",
  A_logB   = "gamma_N[1]",   A_tex2   = "gamma_N[2]",   A_bd   = "gamma_N[3]",
  B_logB   = "gamma_N[1]",   B_tex2   = "gamma_N[2]",   B_bd   = "gamma_N[3]",
  proj_obj = proj_list$N
)

tab_P <- build_resp_rows("P",
  msp_logB = "gamma_P[1]",      msp_tex2 = "gamma_P[3]",   msp_bd = "gamma_P[4]",
  A_logB   = "gamma_P[1]",      A_tex2   = "gamma_P[2]",   A_bd   = "gamma_P[3]",
  B_logB   = "gamma_P_logB",    B_tex2   = NA,             B_bd   = NA,
  proj_obj = proj_list$P
)

tab_all <- bind_rows(tab_SOC, tab_N, tab_P) |>
  select(Risposta, Parametro, `M-SP`, `A (ridotto)`, `B (resp-spec)`, `C (projpred)`)

cat("\n--- SOC ---\n")
print(tab_SOC |> select(-Risposta) |> as.data.frame(), row.names = FALSE)
cat("\n--- N ---\n")
print(tab_N   |> select(-Risposta) |> as.data.frame(), row.names = FALSE)
cat("\n--- P ---\n")
print(tab_P   |> select(-Risposta) |> as.data.frame(), row.names = FALSE)

write.csv(tab_all, file.path(tab_dir, "tab_14_params.csv"), row.names = FALSE)
cat("\n  Salvato: output/tables/tab_14_params.csv\n")


# ── 9. SINTESI DIAGNOSTICA ────────────────────────────────────────────────────

cat("\n═══ SINTESI DIAGNOSTICA ══════════════════════════════════════════\n")
cat(sprintf("%-14s | %10s | %12s | %9s | %8s\n",
            "Modello", "Divergenze", "MaxTreedepth", "Rhat>1.05", "ESS<400"))
cat(strrep("-", 62), "\n")
cat(sprintf("%-14s | %10d | %12d | %9d | %8d\n",
            "A (ridotto)",  n_divs_A, max_td_A, bad_rhat_A, bad_ess_A))
cat(sprintf("%-14s | %10d | %12d | %9d | %8d\n",
            "B (resp-spec)", n_divs_B, max_td_B, bad_rhat_B, bad_ess_B))

cat("\n── Fine script 10 ──────────────────────────────────────────────\n")
