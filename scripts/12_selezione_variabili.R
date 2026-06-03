# =============================================================================
# 12_selezione_variabili.R  —  Selezione variabili con projpred (reference: M-SP-RIRS-MVRE)
#
# COME FUNZIONA projpred (NON fitta modelli Stan per ogni subset):
#   1. Reference model: usa le draws da M-SP-RIRS-MVRE (script 10).
#      Fornisce la distribuzione predittiva posteriore come "target".
#   2. Proiezione (economica): per ogni sottmodello candidato, minimizza
#      la KL-divergenza rispetto alle predizioni del reference.
#      Per Gaussian = regressione pesata (lm/lmer), NON Stan.
#   3. LOO approssimato: confronta i sottmodelli con il reference tramite
#      PSIS-LOO, cerca il più piccolo indistinguibile dal reference.
#
# REFERENCE MODEL: M-SP-RIRS-MVRE (fit_msp_rirs_mvre.rds) via draws cmdstanr.
#   Per ogni risposta (SOC, N, P) separatamente, costruisce la matrice
#   mu (draws × osservazioni) dal fit già salvato.
#
# STRUTTURA MVRE: V[6,J] con righe [SOC_int, SOC_slope, N_int, N_slope, P_int, P_slope]
#   mu_r[i,j] = alpha_r + V[2r-1,j] + V[2r,j]*logBottom_i + gamma_r*X_W + beta_r*X_B
#
# X_W (K_W=5): logBottom[1], Texture1[2], Texture2[3], BulkDensity[4], PH[5]
#
# Dipende da: stan/fit_msp_rirs_mvre.rds, scripts/00_utilities.R
# =============================================================================


# ── 0. SETUP ──────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
source(here("scripts", "00_utilities.R"))

if (!requireNamespace("projpred", quietly = TRUE)) install.packages("projpred")
library(projpred)

setup_rtools()

dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "cache"),   recursive = TRUE, showWarnings = FALSE)
fig_dir   <- here("output", "figures")
cache_dir <- here("output", "cache")

save_fig <- function(fname, p, w = 14, h = 9, u = "cm") {
  ggplot2::ggsave(file.path(fig_dir, fname), plot = p,
                  width = w, height = h, units = u, device = "pdf")
  cat(sprintf("  [fig] Salvato: %s\n", fname))
}

projpred_cache_path <- file.path(cache_dir, "projpred_varsel_mvre.rds")


# ── 1. DATI (identici a script 09) ────────────────────────────────────────────

dati <- readRDS(here("data", "dati.rds")) |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom)
  ) |>
  mutate(across(c(logBottom, Texture1, Texture2, BulkDensity, PH),
                ~ c(scale(.x)))) |>
  mutate(Field = factor(Field))

field_levels <- sort(unique(as.integer(as.character(dati$Field))))
J <- length(field_levels)
N <- nrow(dati)

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))

X_W_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")
X_W      <- as.matrix(dati_int[, X_W_cols])

X_B_cols <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")
X_B <- dati_int |>
  distinct(field_int, across(all_of(X_B_cols))) |>
  arrange(field_int) |>
  select(all_of(X_B_cols)) |>
  as.matrix()

field_id <- dati_int$field_int

cat(sprintf("N = %d | J = %d | K_W = %d\n", N, J, length(X_W_cols)))


# ── 2. CARICA FIT M-SP-RIRS-MVRE ─────────────────────────────────────────────

fit_path <- here("stan", "fit_msp_rirs_mvre.rds")
cat("Carico M-SP-RIRS-MVRE da:", fit_path, "\n")
fit_ref <- readRDS(fit_path)


# ── 3. CALCOLO MATRICE MU DAL MODELLO M-SP-RIRS-MVRE ─────────────────────────
# Per ogni risposta r, costruisce la matrice mu (D_thin × N) dalle draws.
#
# Struttura M-SP-RIRS-MVRE: V[6,J] con transformed parameter
#   righe: [SOC_int=1, SOC_slope=2, N_int=3, N_slope=4, P_int=5, P_slope=6]
#   mu_r[d,n] = alpha_r[d]
#             + V[2r-1, field_id[n]]           ← random intercept
#             + V[2r,   field_id[n]] * X_W[n,1] ← random slope su logBottom
#             + X_W[n,] %*% gamma_r[d,]
#             + X_B[field_id[n],] %*% beta_r[d,]

draws_raw  <- fit_ref$draws(format = "matrix")
draws_full <- matrix(as.double(draws_raw),
                     nrow     = nrow(draws_raw),
                     ncol     = ncol(draws_raw),
                     dimnames = dimnames(draws_raw))
D_total    <- nrow(draws_full)
thin_idx   <- seq(1, D_total, by = ceiling(D_total / 500))
draws_t    <- draws_full[thin_idx, ]
D          <- nrow(draws_t)
cat(sprintf("Draw totali: %d | Draw thinned: %d\n", D_total, D))
rm(draws_full); gc()

# Funzione: estrae matrice D × K per parametri vettoriali
extract_mat <- function(draws, prefix, K) {
  cols <- paste0(prefix, "[", 1:K, "]")
  draws[, cols, drop = FALSE]
}

# Funzione: estrae matrice D × J da V[row, j] (MVRE: matrice 6×J in Stan)
extract_V_row <- function(draws, row) {
  cols <- paste0("V[", row, ",", 1:J, "]")
  draws[, cols, drop = FALSE]
}

# Funzione: calcola mu_r per una risposta (M-SP-RIRS-MVRE)
# Mappa risposta → righe di V[6,J]:
#   SOC: int=V[1,j], slope=V[2,j]
#   N:   int=V[3,j], slope=V[4,j]
#   P:   int=V[5,j], slope=V[6,j]
compute_mu <- function(r) {
  cat(sprintf("  Calcolo mu_%s...\n", r))

  row_int <- c(SOC = 1L, N = 3L, P = 5L)[r]
  row_slo <- c(SOC = 2L, N = 4L, P = 6L)[r]

  alpha   <- draws_t[, paste0("alpha_", r)]               # D
  sigma   <- draws_t[, paste0("sigma_", r)]               # D
  u_int   <- extract_V_row(draws_t, row_int)              # D × J (random intercepts)
  u_slope <- extract_V_row(draws_t, row_slo)              # D × J (random slopes)
  gamma   <- extract_mat(draws_t, paste0("gamma_", r), length(X_W_cols))  # D × K_W
  beta    <- extract_mat(draws_t, paste0("beta_", r),  length(X_B_cols))  # D × K_B

  # Espandi X_B al livello osservazione
  X_B_n <- X_B[field_id, ]  # N × K_B

  # Parte casuale: u_int[field_id] + u_slope[field_id] * logBottom
  rand_int   <- u_int[,   field_id]                              # D × N
  rand_slope <- u_slope[, field_id] * outer(rep(1, D), X_W[, 1]) # D × N  (col 1 = logBottom)

  # Effetti fissi within: gamma %*% t(X_W)
  fix_within  <- tcrossprod(gamma, X_W)   # D × N

  # Effetti fissi between: beta %*% t(X_B_n)
  fix_between <- tcrossprod(beta, X_B_n)  # D × N

  # Intercetta globale
  intercept <- outer(alpha, rep(1, N))  # D × N

  mu <- intercept + rand_int + rand_slope + fix_within + fix_between  # D × N

  list(mu = mu, sigma = sigma)
}

cat("Calcolo matrici mu per tutte le risposte (M-SP-RIRS-MVRE)...\n")
res_SOC <- compute_mu("SOC")
res_N   <- compute_mu("N")
res_P   <- compute_mu("P")
rm(draws_t); gc()


# ── 4. DATA FRAME PER I SUBMODELLI DI PROIEZIONE ─────────────────────────────
# I submodelli di proiezione (lm/lmer interni a projpred) hanno bisogno di
# un data frame con le covariate e la variabile di raggruppamento per Field.

# X_B espanso a livello osservazione
X_B_expanded <- as.data.frame(X_B[field_id, ])
names(X_B_expanded) <- X_B_cols

data_pp <- bind_cols(
  as.data.frame(X_W),            # logBottom, Texture1, Texture2, BulkDensity, PH
  X_B_expanded,                  # OnFarm, Irrigate, Fertilised, N_Natural
  data.frame(Field = dati_int$Field)  # raggruppamento
)

# Le risposte saranno passate esplicitamente a init_refmodel


# ── 5. COSTRUZIONE REFERENCE MODEL E VARSEL ───────────────────────────────────
#
# init_refmodel() di projpred accetta:
#   mu  = matrice D × N (draws delle medie posteriori)
#   dis = vettore D (draws di sigma, la "dispersione" per Gaussian)
#
# La formula del submodello definisce il search space.
# logBottom è incluso sempre (mandatory) perché è la covariata di profondità
# principale. Le altre 7 covariate sono opzionali.

# Formula del search space (massima per i submodelli).
# Usa (logBottom | Field) per riflettere la struttura del M-SP:
#   z_nu_r[j] * (psi_r + eta_r * logBottom)  ←  random intercept + random slope
# In lme4, (logBottom | Field) = intercetta casuale + slope casuale per logBottom
# con correlazione libera. Non è la struttura proporzionale esatta (che ha un solo
# parametro eta invece di una matrice di covarianza 2×2), ma è la più vicina
# disponibile nei submodelli di proiezione.
# logBottom è incluso anche nei fissi perché nel M-SP c'è sia gamma_r[1]
# (fisso, uguale per tutti) sia z_nu_r[j]*eta_r (casuale, per campo).
formula_full <- logY ~ logBottom + Texture1 + Texture2 + BulkDensity + PH +
                       OnFarm + Irrigate + Fertilised + N_Natural +
                       (logBottom | Field)

run_projpred <- function(mu_draws, sigma_draws, y_vec, risposta,
                         data_submod, formula_full) {

  cat(sprintf("\n═══ projpred per log%s ══════════════════════════════════\n",
              risposta))

  # Aggiunge risposta al data frame dei submodelli
  data_sub <- data_submod |>
    mutate(logY = y_vec)

  # ref_predfun: deve restituire matrice N × D (osservazioni × draws)
  # mu_draws è D × N internamente → trasponiamo
  ref_pf <- function(fit, newdata = NULL) t(fit$mu)  # N × D

  # extract_model_data: restituisce y, weights e offset per le osservazioni
  extr_fn <- function(fit, newdata = NULL, ...) {
    nd <- if (is.null(newdata)) data_sub else newdata
    list(
      y       = nd$logY,
      weights = rep(1, nrow(nd)),
      offset  = rep(0, nrow(nd))
    )
  }

  # Crea reference model dal M-SP (API projpred 2.10.0)
  refmod <- tryCatch(
    suppressWarnings(
      init_refmodel(
        object             = list(mu = mu_draws),  # D × N portato in fit$mu
        data               = data_sub,
        formula            = formula_full,
        family             = gaussian(),
        ref_predfun        = ref_pf,
        dis                = sigma_draws,          # D (sigma posteriori)
        extract_model_data = extr_fn
      )
    ),
    error = function(e) {
      cat("  init_refmodel fallito:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(refmod)) {
    cat("  Impossibile costruire il reference model.\n")
    cat("  Versione installata:", as.character(packageVersion("projpred")), "\n")
    return(invisible(NULL))
  }

  # Variable selection con forward search e PSIS-LOO interno (veloce).
  # varsel() usa già PSIS-LOO per confrontare i submodelli con il reference:
  # è sufficiente per la selezione. cv_varsel aggiunge un outer-CV che con
  # N=220 obs e lmer interni richiede ore — non necessario qui.
  cat("  Avvio varsel (forward search, PSIS-LOO interno)...\n")
  vs <- tryCatch(
    varsel(
      refmod,
      method     = "forward",
      nterms_max = 8,
      verbose    = FALSE,
      seed       = 2024
    ),
    error = function(e) {
      cat("  varsel fallito:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(vs)) return(invisible(NULL))

  # Dimensione suggerita
  size_sug <- tryCatch(suggest_size(vs, alpha = 0.1), error = function(e) NA)
  terms_ord <- solution_terms(vs)

  cat(sprintf("  Ordine selezione: %s\n", paste(terms_ord, collapse = " → ")))
  cat(sprintf("  Dimensione suggerita (α=0.10): %s\n",
              ifelse(is.na(size_sug), "NA", as.character(size_sug))))

  if (!is.na(size_sug) && size_sug > 0) {
    cat(sprintf("  Variabili selezionate: %s\n",
                paste(head(terms_ord, size_sug), collapse = ", ")))
  } else {
    cat("  Nessuna variabile fissa necessaria oltre al random intercept\n")
  }

  # Plot ELPD path
  p <- tryCatch(
    plot(vs, stats = "elpd", deltas = TRUE) +
      ggtitle(sprintf("ELPD path — log%s (reference = M-SP-RIRS-MVRE)", risposta)) +
      theme_minimal(),
    error = function(e) NULL
  )
  if (!is.null(p)) print(p)

  invisible(vs)
}

# ── Cache: carica se esiste, altrimenti calcola e salva ──────────────────────
if (file.exists(projpred_cache_path)) {
  cat("\nCarico projpred cache:", projpred_cache_path, "\n")
  vs_cache <- readRDS(projpred_cache_path)
  vs_SOC   <- vs_cache$SOC
  vs_N     <- vs_cache$N
  vs_P     <- vs_cache$P
  rm(vs_cache)
} else {
  cat("\nCache projpred non trovata — avvio calcolo (può richiedere 15-30 min)...\n")

  vs_SOC <- run_projpred(res_SOC$mu, res_SOC$sigma,
                         dati_int$logSOC, "SOC", data_pp, formula_full)

  vs_N   <- run_projpred(res_N$mu,   res_N$sigma,
                         dati_int$logN,   "N",   data_pp, formula_full)

  vs_P   <- run_projpred(res_P$mu,   res_P$sigma,
                         dati_int$logP,   "P",   data_pp, formula_full)

  cat("\nSalvo projpred cache...\n")
  saveRDS(list(SOC = vs_SOC, N = vs_N, P = vs_P), projpred_cache_path)
  cat("  Cache salvata:", projpred_cache_path, "\n")
}


# ── 6. RIEPILOGO ─────────────────────────────────────────────────────────────

cat("\n═══ RIEPILOGO SELEZIONE VARIABILI ═══════════════════════════════\n")
cat(sprintf("%-6s | %-40s | %s\n", "Resp.", "Ordine selezione (prime 5)", "n suggerito"))
cat(strrep("-", 65), "\n")

for (nm in c("SOC", "N", "P")) {
  vs_obj <- get(paste0("vs_", nm))
  if (!is.null(vs_obj)) {
    terms  <- solution_terms(vs_obj)
    n_sug  <- tryCatch(suggest_size(vs_obj, alpha = 0.1), error = function(e) NA)
    top5   <- paste(head(terms, 5), collapse = " → ")
    cat(sprintf("log%-3s | %-40s | %s\n", nm, top5,
                ifelse(is.na(n_sug), "NA", as.character(n_sug))))
  } else {
    cat(sprintf("log%-3s | [varsel non riuscito]\n", nm))
  }
}

cat("\n  Interpretazione:\n")
cat("  - n=0: solo (1|Field) + intercetta è sufficiente\n")
cat("  - n=k: le prime k variabili nell'ordine di selezione sono necessarie\n")
cat("  - 'necessario' = aggiungere altre variabili non migliora la predizione\n")
cat("    rispetto al reference (M-SP) in modo statisticamente rilevante\n")

# ── Salva figure projpred ─────────────────────────────────────────────────────
cat("\nSalvo figure projpred...\n")
risposta_nomi <- list(SOC = vs_SOC, N = vs_N, P = vs_P)
pp_plots <- list()

for (nm in names(risposta_nomi)) {
  vs_obj <- risposta_nomi[[nm]]
  if (!is.null(vs_obj)) {
    p <- tryCatch(
      plot(vs_obj, stats = "elpd", deltas = TRUE) +
        ggtitle(sprintf("ELPD path — log%s (reference = M-SP-RIRS)", nm)) +
        labs(subtitle = "ΔELPD rispetto al reference model M-SP-RIRS-MVRE",
             x = "Numero di predittori inclusi", y = "ΔELPD") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold")),
      error = function(e) { cat("  Errore plot", nm, ":", conditionMessage(e), "\n"); NULL }
    )
    pp_plots[[nm]] <- p
    if (!is.null(p)) save_fig(sprintf("fig_08_projpred_%s.pdf", tolower(nm)), p, w = 12, h = 8)
  }
}

# Salva anche summary strutturato come RDS per uso in script 08
projpred_summary <- lapply(names(risposta_nomi), function(nm) {
  vs_obj <- risposta_nomi[[nm]]
  if (is.null(vs_obj)) return(NULL)
  list(
    risposta   = nm,
    ordine     = tryCatch(solution_terms(vs_obj), error = function(e) character(0)),
    n_suggerito = tryCatch(suggest_size(vs_obj, alpha = 0.1), error = function(e) NA_integer_)
  )
})
names(projpred_summary) <- names(risposta_nomi)
saveRDS(projpred_summary, file.path(cache_dir, "projpred_summary_mvre.rds"))
cat("  Sommario projpred salvato in output/cache/projpred_summary_mvre.rds\n")

cat("\n── Fine script 12 ──────────────────────────────────────────────\n")
