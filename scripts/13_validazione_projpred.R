# =============================================================================
# 13_validazione_projpred.R  —  Validazione: posteriors proiettati da projpred vs M-SP-RIRS-MVRE
#
# Per ogni risposta r, costruisce il "modello finale" come il sottosistema
# selezionato da projpred (script 12) e ne verifica la coerenza con M-SP-RIRS-MVRE.
#
# Strategia:
#   1. Tenta project() di projpred per ottenere posteriori proiettati
#      (nessun MCMC: proiezione diretta dal reference model M-SP-RIRS-MVRE)
#   2. Fit frequentista con lme4 (sanity check, stessa formula)
#   3. Confronta coefficienti fissi: MVRE gamma vs Projected vs lme4
#   4. Verifica che gamma di MVRE per le variabili NON selezionate siano ~ 0
#
# X_W (K_W=5): logBottom[1], Texture1[2], Texture2[3], BulkDensity[4], PH[5]
#
# Output:
#   output/cache/proj_posteriors_mvre.rds       -> draws proiettate
#   output/cache/lme4_final_mvre.rds            -> fit lme4
#   output/figures/fig_13_proj_vs_msp.pdf       -> confronto coeff. fissi
#   output/figures/fig_14_nonselected_gamma.pdf -> verifica gamma ~ 0
#   output/tables/tab_13_final_comparison.csv   -> tabella riassuntiva
#
# Dipende da: stan/fit_msp_rirs_mvre.rds, output/cache/projpred_varsel_mvre.rds
# =============================================================================


# ── 0. SETUP ──────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
library(lme4)
library(posterior)
source(here("scripts", "00_utilities.R"))

if (!requireNamespace("projpred", quietly = TRUE)) install.packages("projpred")
library(projpred)

setup_rtools()

dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "cache"),   recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables"),  recursive = TRUE, showWarnings = FALSE)
fig_dir   <- here("output", "figures")
cache_dir <- here("output", "cache")
tab_dir   <- here("output", "tables")

save_fig <- function(fname, p, w = 18, h = 12, u = "cm") {
  ggplot2::ggsave(file.path(fig_dir, fname), plot = p,
                  width = w, height = h, units = u, device = "pdf")
  cat(sprintf("  [fig] Salvato: %s\n", fname))
}


# ── 1. DATI (identici a script 11) ────────────────────────────────────────────

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

cat(sprintf("N = %d | J = %d\n", N, J))


# ── 2. CARICA PROJPRED CACHE ──────────────────────────────────────────────────

projpred_cache_path <- file.path(cache_dir, "projpred_varsel_mvre.rds")
cat("Carico projpred cache da:", projpred_cache_path, "\n")
vs_cache <- readRDS(projpred_cache_path)
vs_SOC   <- vs_cache$SOC
vs_N     <- vs_cache$N
vs_P     <- vs_cache$P
rm(vs_cache)

n_SOC <- tryCatch(suggest_size(vs_SOC, alpha = 0.1), error = function(e) 4L)
n_N   <- tryCatch(suggest_size(vs_N,   alpha = 0.1), error = function(e) 4L)
n_P   <- tryCatch(suggest_size(vs_P,   alpha = 0.1), error = function(e) 3L)

terms_SOC <- head(solution_terms(vs_SOC), n_SOC)
terms_N   <- head(solution_terms(vs_N),   n_N)
terms_P   <- head(solution_terms(vs_P),   n_P)

cat(sprintf("\nModelli selezionati da projpred:\n"))
cat(sprintf("  SOC (n=%d): %s\n", n_SOC, paste(terms_SOC, collapse = ", ")))
cat(sprintf("  N   (n=%d): %s\n", n_N,   paste(terms_N,   collapse = ", ")))
cat(sprintf("  P   (n=%d): %s\n", n_P,   paste(terms_P,   collapse = ", ")))


# ── 3. PROIEZIONE CON project() (tenta, fallback a solo lme4) ─────────────────

proj_cache_path <- file.path(cache_dir, "proj_posteriors_mvre.rds")

extract_proj_summary <- function(proj_obj, resp_name) {
  if (is.null(proj_obj)) return(NULL)
  tryCatch({
    prj_draws <- as.matrix(proj_obj)
    cat(sprintf("  Colonne projected %s: %s\n", resp_name,
                paste(head(colnames(prj_draws), 10), collapse = ", ")))
    # Estrae solo i coefficienti fissi (escludi sigma, Sigma, intercetta)
    skip_pat <- "^sigma|^Sigma|^\\(Intercept\\)"
    fixed_cols <- colnames(prj_draws)[!grepl(skip_pat, colnames(prj_draws))]
    if (length(fixed_cols) == 0) return(NULL)
    bind_rows(lapply(fixed_cols, function(nm) {
      tibble(
        parameter = nm,
        median    = median(prj_draws[, nm]),
        q05       = quantile(prj_draws[, nm], 0.05),
        q95       = quantile(prj_draws[, nm], 0.95),
        risposta  = resp_name,
        fonte     = "Projected (projpred)"
      )
    }))
  }, error = function(e) {
    cat(sprintf("  Errore estrazione projected %s: %s\n", resp_name, conditionMessage(e)))
    NULL
  })
}

if (file.exists(proj_cache_path)) {
  cat("\nCarico posteriori proiettati da cache...\n")
  proj_list <- readRDS(proj_cache_path)
} else {
  cat("\nCalcolo proiezioni projpred (alcuni minuti)...\n")
  project_safe <- function(vs_obj, nterms, resp_name) {
    cat(sprintf("  Proiezione %s (n=%d)...\n", resp_name, nterms))
    tryCatch(
      project(vs_obj, nterms = nterms, seed = 2024, verbose = FALSE),
      error = function(e) {
        cat(sprintf("  project() fallito per %s: %s\n", resp_name, conditionMessage(e)))
        NULL
      }
    )
  }
  proj_SOC <- project_safe(vs_SOC, n_SOC, "SOC")
  proj_N   <- project_safe(vs_N,   n_N,   "N")
  proj_P   <- project_safe(vs_P,   n_P,   "P")
  proj_list <- list(SOC = proj_SOC, N = proj_N, P = proj_P)
  saveRDS(proj_list, proj_cache_path)
  cat("  Proiezioni salvate in cache:", proj_cache_path, "\n")
}

proj_summary <- bind_rows(
  extract_proj_summary(proj_list$SOC, "SOC"),
  extract_proj_summary(proj_list$N,   "N"),
  extract_proj_summary(proj_list$P,   "P")
)
cat(sprintf("Posterior proiettati estratti: %d righe\n", nrow(proj_summary)))


# ── 4. FIT LMER SANITY CHECK ──────────────────────────────────────────────────

lme4_cache_path <- file.path(cache_dir, "lme4_final_mvre.rds")

if (file.exists(lme4_cache_path)) {
  cat("\nCarico lme4 da cache.\n")
  lme4_fits <- readRDS(lme4_cache_path)
} else {
  cat("\nFit lme4 per i modelli finali proiettati...\n")

  # SOC: logBottom + Texture2 + random slope per logBottom
  m_soc <- lmer(logSOC ~ logBottom + Texture2 + (logBottom | Field),
                data = dati_int, REML = TRUE)
  cat("  SOC: convergenza =", isSingular(m_soc) == FALSE, "\n")

  # N: logBottom + Texture2 + BulkDensity + random intercept
  m_n <- lmer(logN ~ logBottom + Texture2 + BulkDensity + (1 | Field),
              data = dati_int, REML = TRUE)
  cat("  N:   convergenza =", isSingular(m_n) == FALSE, "\n")

  # P: logBottom + random slope per logBottom (nessun fixed effect oltre logBottom)
  m_p <- lmer(logP ~ logBottom + (logBottom | Field),
              data = dati_int, REML = TRUE)
  cat("  P:   convergenza =", isSingular(m_p) == FALSE, "\n")

  lme4_fits <- list(SOC = m_soc, N = m_n, P = m_p)
  saveRDS(lme4_fits, lme4_cache_path)
  cat("  Cache lme4 salvata.\n")
}

lme4_summary <- bind_rows(lapply(c("SOC", "N", "P"), function(r) {
  m <- lme4_fits[[r]]
  coefs <- fixef(m)
  se    <- sqrt(diag(as.matrix(vcov(m))))
  coefs_no_int <- coefs[names(coefs) != "(Intercept)"]
  se_no_int    <- se[names(se) != "(Intercept)"]
  if (length(coefs_no_int) == 0) return(NULL)
  tibble(
    parameter = names(coefs_no_int),
    median    = unname(coefs_no_int),
    q05       = unname(coefs_no_int) - 1.645 * unname(se_no_int),
    q95       = unname(coefs_no_int) + 1.645 * unname(se_no_int),
    risposta  = r,
    fonte     = "lme4 (REML)"
  )
}))

cat("\nRiepilogo lme4:\n")
for (r in c("SOC", "N", "P")) {
  cat(sprintf("\n  log%s:\n", r))
  print(summary(lme4_fits[[r]])$coefficients |> round(3))
  print(VarCorr(lme4_fits[[r]]))
}


# ── 5. CARICA M-SP PER CONFRONTO ─────────────────────────────────────────────

cat("\nCarico M-SP-RIRS-MVRE per confronto gamma...\n")
fit_msp <- readRDS(here("stan", "fit_msp_rirs_mvre.rds"))

draws_raw <- fit_msp$draws(format = "matrix")
draws_msp <- matrix(as.double(draws_raw), nrow = nrow(draws_raw),
                    ncol = ncol(draws_raw), dimnames = dimnames(draws_raw))
rm(draws_raw); gc()
cat(sprintf("  Draw M-SP-RIRS: %d\n", nrow(draws_msp)))

X_W_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")
X_B_cols <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")

# Riassunto posteriori gamma_r e beta_r
summarize_param <- function(v) {
  tibble(median = median(v), q05 = quantile(v, .05), q95 = quantile(v, .95))
}

msp_summary <- bind_rows(lapply(c("SOC", "N", "P"), function(r) {
  bind_rows(
    # gamma (within-field)
    lapply(seq_along(X_W_cols), function(k) {
      v <- draws_msp[, paste0("gamma_", r, "[", k, "]")]
      summarize_param(v) |>
        mutate(parameter = X_W_cols[k], risposta = r,
               fonte = "M-SP-RIRS-MVRE (reference)", tipo = "within-field")
    }),
    # beta (between-field)
    lapply(seq_along(X_B_cols), function(k) {
      v <- draws_msp[, paste0("beta_", r, "[", k, "]")]
      summarize_param(v) |>
        mutate(parameter = X_B_cols[k], risposta = r,
               fonte = "M-SP-RIRS-MVRE (reference)", tipo = "management")
    })
  )
}))

# tau_alpha, tau_beta, rho (struttura random bivariata M-SP-RIRS)
struct_summary <- bind_rows(lapply(c("SOC", "N", "P"), function(r) {
  ta  <- draws_msp[, paste0("tau_alpha_", r)]
  tb  <- draws_msp[, paste0("tau_beta_",  r)]
  rho <- draws_msp[, paste0("rho_",       r)]
  bind_rows(
    summarize_param(ta)  |> mutate(parameter = paste0("tau_alpha_", r), risposta = r,
                                   fonte = "M-SP-RIRS-MVRE (reference)", tipo = "random_struct"),
    summarize_param(tb)  |> mutate(parameter = paste0("tau_beta_",  r), risposta = r,
                                   fonte = "M-SP-RIRS-MVRE (reference)", tipo = "random_struct"),
    summarize_param(rho) |> mutate(parameter = paste0("rho_",       r), risposta = r,
                                   fonte = "M-SP-RIRS-MVRE (reference)", tipo = "random_struct")
  )
}))

cat("\n═══ M-SP gamma (within-field) ═══════════════════════════════\n")
print(msp_summary |> filter(tipo == "within-field") |>
  select(risposta, parameter, median, q05, q95) |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  as.data.frame())

cat("\n═══ M-SP beta (management) ══════════════════════════════════\n")
cat("(nessuno di questi dovrebbe essere rilevante - projpred ha escluso tutto il management)\n")
print(msp_summary |> filter(tipo == "management") |>
  select(risposta, parameter, median, q05, q95) |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  as.data.frame())

cat("\n═══ M-SP-RIRS tau_alpha, tau_beta, rho (struttura bivariata) ═══\n")
print(struct_summary |>
  select(parameter, median, q05, q95) |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  as.data.frame())


# ── 6. FIGURA 1: confronto coeff. fissi selezionati ───────────────────────────

# Variabili selezionate per ogni risposta (solo fixed effects)
selected_fixed <- list(
  SOC = c("logBottom", "Texture2"),
  N   = c("logBottom", "Texture2", "BulkDensity"),
  P   = character(0)  # nessun fixed effect per logP
)

msp_sel <- msp_summary |>
  filter(tipo == "within-field") |>
  filter(
    (risposta == "SOC" & parameter %in% selected_fixed$SOC) |
    (risposta == "N"   & parameter %in% selected_fixed$N)
  )

compare_df <- bind_rows(
  msp_sel,
  proj_summary |> filter(!is.na(parameter)),
  lme4_summary
) |>
  filter(risposta %in% c("SOC", "N"))  # P non ha fixed effects

if (nrow(compare_df) > 0) {
  p_compare <- compare_df |>
    ggplot(aes(x = parameter, y = median, color = fonte,
               ymin = q05, ymax = q95, shape = fonte)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(position = position_dodge(width = 0.5), size = 2.8) +
    geom_errorbar(position = position_dodge(width = 0.5), width = 0.2,
                  na.rm = TRUE) +
    facet_wrap(~ risposta, scales = "free_y", ncol = 2) +
    scale_color_brewer(palette = "Set1", name = "Stima") +
    scale_shape_manual(name = "Stima",
                       values = c("M-SP-RIRS-MVRE (reference)" = 16,
                                  "Projected (projpred)" = 17,
                                  "lme4 (REML)" = 15)) +
    labs(
      title    = "Confronto coefficienti fissi: M-SP vs Projected vs lme4",
      subtitle = "Solo covariate selezionate da projpred | Mediana e CI90%",
      x = "Covariata", y = "Coefficiente (scala standardizzata)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold"),
      legend.position  = "bottom",
      axis.text.x      = element_text(angle = 20, hjust = 1)
    )
  save_fig("fig_13_proj_vs_msp.pdf", p_compare, w = 20, h = 12)
} else {
  cat("  Nessun dato per fig_13 (projected posteriors non disponibili).\n")
}


# ── 7. FIGURA 2: gamma non selezionati in M-SP (dovrebbero essere ~ 0) ────────

nonsel_fixed <- list(
  SOC = setdiff(X_W_cols, selected_fixed$SOC),  # Texture1, BulkDensity
  N   = setdiff(X_W_cols, selected_fixed$N),     # Texture1
  P   = X_W_cols                                 # tutti (nessun fixed selezionato)
)

msp_nonsel <- msp_summary |>
  filter(tipo == "within-field") |>
  filter(
    (risposta == "SOC" & parameter %in% nonsel_fixed$SOC) |
    (risposta == "N"   & parameter %in% nonsel_fixed$N)   |
    (risposta == "P"   & parameter %in% nonsel_fixed$P)
  )

management_df <- msp_summary |> filter(tipo == "management")

nonsel_all <- bind_rows(
  msp_nonsel |> mutate(gruppo = "Within-field (non sel.)"),
  management_df |> mutate(gruppo = "Management (nessuno sel.)")
)

p_nonsel <- nonsel_all |>
  ggplot(aes(x = parameter, y = median, ymin = q05, ymax = q95, color = risposta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(position = position_dodge(width = 0.4), size = 2.5) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.2) +
  facet_wrap(~ gruppo, scales = "free_x", ncol = 1) +
  scale_color_brewer(palette = "Set2", name = "Risposta") +
  labs(
    title    = "Covariate NON selezionate da projpred: gamma in M-SP",
    subtitle = "Atteso: CI90% copre zero (conferma selezione) | Mediana e CI90%",
    x = "Covariata", y = "gamma o beta (scala standardizzata)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 30, hjust = 1)
  )

save_fig("fig_14_nonselected_gamma.pdf", p_nonsel, w = 20, h = 16)


# ── 8. TABELLA RIEPILOGATIVA FINALE ───────────────────────────────────────────

tab_final <- bind_rows(
  msp_summary |>
    mutate(ci90 = sprintf("[%.3f, %.3f]", q05, q95)) |>
    select(risposta, tipo, parameter, median, ci90, fonte),
  struct_summary |>
    mutate(ci90 = sprintf("[%.3f, %.3f]", q05, q95)) |>
    select(risposta, tipo, parameter, median, ci90, fonte),
  lme4_summary |>
    mutate(ci90 = ifelse(is.na(q05), "-", sprintf("[%.3f, %.3f]", q05, q95)),
           tipo = "lme4_final") |>
    select(risposta, tipo, parameter, median, ci90, fonte)
)

write.csv(tab_final, file.path(tab_dir, "tab_13_final_comparison.csv"),
          row.names = FALSE)
cat("\n  Tabella salvata: output/tables/tab_13_final_comparison.csv\n")


# ── 9. DIAGNOSI DIVERGENZE ────────────────────────────────────────────────────

cat("\n═══ DIAGNOSI DIVERGENZE projpred vs M-SP ═══════════════════\n")
cat("Confronto coefficienti fissi selezionati: M-SP gamma vs lme4\n\n")

for (r in c("SOC", "N")) {
  cat(sprintf("--- log%s ---\n", r))
  msp_r  <- msp_sel |> filter(risposta == r)
  lme4_r <- lme4_summary |> filter(risposta == r)
  merged <- inner_join(msp_r, lme4_r, by = "parameter",
                        suffix = c("_msp", "_lme4"))
  if (nrow(merged) == 0) {
    cat("  Nessun parametro comune trovato.\n")
    next
  }
  merged <- merged |>
    mutate(
      diff_abs  = abs(median_msp - median_lme4),
      divergenza = ifelse(diff_abs > 0.2, "ATTENZIONE", "OK")
    )
  print(merged |>
    select(parameter, median_msp, median_lme4, diff_abs, divergenza) |>
    mutate(across(where(is.numeric), ~ round(.x, 3))) |>
    as.data.frame())
  cat("\n")
}

cat("Nota: differenze > 0.2 (scala standardizzata) indicate come 'ATTENZIONE'.\n")
cat("Atteso: coerenza tra M-SP e lme4 per le covariate selezionate.\n")
cat("M-SP usa la struttura proporzionale (eta) - lme4 usa (logBottom|Field) come approssimazione.\n")

cat("\n── Fine script 13 ──────────────────────────────────────────────\n")
