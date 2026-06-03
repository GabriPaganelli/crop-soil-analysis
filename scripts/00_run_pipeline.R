# =============================================================================
# 00_run_pipeline.R  —  Master script: esegue l'intera pipeline downstream
#
# Esegue in ordine gli script 10–19 dopo che i modelli principali (07-09) sono
# stati fittati. Il modello finale è M-SP-RIRS-MVRE (fit_msp_rirs_mvre.rds).
#
# PREREQUISITI:
#   stan/fit_msp_rirs_mvre.rds   — modello finale
#   stan/fit_msp_rirs.rds        — confronto M-SP-RIRS
#   stan/fit_msp.rds             — confronto M-SP
#   stan/fit_mri.rds             — confronto M-RI
#   data/dati.rds                — dati preprocessati
#
# STRUTTURA:
#   1. Validazione prerequisiti
#   2. Esecuzione sequenziale degli script con gestione errori
#   3. Riepilogo finale
#
# TEMPI ATTESI (prima esecuzione, senza cache):
#   Script 10 (MVRE LOO):    ~5 min (solo LOO, fit già salvato)
#   Script 12 (projpred):    ~20-40 min
#   Script 13 (validazione): ~10 min
#   Script 14 (A/B):         ~2h (due MCMC da zero)
#   Script 15 (sensitivity): ~1h (un MCMC da zero)
#   Script 16 (MV LOO):      ~5 min (solo LOO, fit già salvato)
#   Script 17 (frequentista): ~5 min
#   Script 18 (figure):      ~10 min
#   Script 19 (report):      ~3 min
# =============================================================================

library(here)

# ── 0. SETUP ──────────────────────────────────────────────────────────────────

source(here("scripts", "00_utilities.R"))
setup_rtools()

# ── 1. VALIDAZIONE PREREQUISITI ───────────────────────────────────────────────

validate_prerequisites <- function() {
  cat("═══ VALIDAZIONE PREREQUISITI ═══════════════════════════════════\n")

  required_files <- list(
    fit_mvre  = here("stan", "fit_msp_rirs_mvre.rds"),
    fit_rirs  = here("stan", "fit_msp_rirs.rds"),
    fit_msp   = here("stan", "fit_msp.rds"),
    fit_mri   = here("stan", "fit_mri.rds"),
    dati      = here("data", "dati.rds"),
    stan_mvre = here("stan", "m_sp_rirs_mvre.stan"),
    stan_A    = here("stan", "m_sp_rirs_mvre_A.stan"),
    stan_B    = here("stan", "m_sp_rirs_mvre_B.stan")
  )

  for (nm in names(required_files)) {
    if (!file.exists(required_files[[nm]])) {
      stop(sprintf("File mancante: %s\n  Percorso: %s", nm, required_files[[nm]]))
    }
    cat(sprintf("  OK %-12s: %s\n", nm, basename(required_files[[nm]])))
  }

  # Verifica parametri MVRE accessibili
  cat("  Verifica parametri MVRE...\n")
  fit_test <- tryCatch(readRDS(here("stan", "fit_msp_rirs_mvre.rds")),
                       error = function(e) stop("Impossibile caricare fit MVRE: ", conditionMessage(e)))
  pnames <- tryCatch(fit_test$metadata()$stan_variables,
                     error = function(e) character(0))
  needed <- c("V", "tau_alpha_SOC", "rho_SOC", "rho_int_SOC_N")
  missing_p <- setdiff(needed, pnames)
  if (length(missing_p) > 0) {
    stop("Parametri mancanti nel fit MVRE: ", paste(missing_p, collapse = ", "),
         "\n  Il fit potrebbe essere corrotto o da una versione diversa dello Stan file.")
  }
  rm(fit_test); gc()
  cat("  OK parametri MVRE verificati\n")
  cat("═══ VALIDAZIONE OK ══════════════════════════════════════════════\n\n")
  invisible(TRUE)
}

# ── 2. RUNNER ─────────────────────────────────────────────────────────────────

run_safe <- function(script_path, label) {
  cat(sprintf("\n%s\n=== %s ===\n%s\n", strrep("=", 60), label, strrep("=", 60)))
  t0 <- proc.time()["elapsed"]
  result <- withCallingHandlers(
    tryCatch({
      source(script_path, local = new.env(parent = globalenv()))
      "OK"
    }, error = function(e) {
      msg <- sprintf("ERRORE: %s", conditionMessage(e))
      cat("\n!!!", msg, "\n")
      msg
    }),
    warning = function(w) {
      invokeRestart("muffleWarning")
    }
  )
  elapsed <- proc.time()["elapsed"] - t0
  cat(sprintf("\n[%s] %.0f sec\n", result, elapsed))
  invisible(result)
}

# ── 3. PIPELINE ───────────────────────────────────────────────────────────────

tryCatch(validate_prerequisites(), error = function(e) {
  cat("\n!!! VALIDAZIONE FALLITA:", conditionMessage(e), "\n")
  cat("Correggere i prerequisiti prima di eseguire la pipeline.\n")
  stop(conditionMessage(e))
})

pipeline <- c(
  "10 Modello finale MVRE"   = "10_model_final_mvre.R",
  "11 Landuse variante"      = "11_model_msp_landuse.R",
  "12 Projpred (MVRE)"       = "12_selezione_variabili.R",
  "13 Validazione projpred"  = "13_validazione_projpred.R",
  "14 Confronto A/B"         = "14_confronto_modelli.R",
  "15 Sensitivity Pareto"    = "15_sensitivity_pareto.R",
  "16 MV residui"            = "16_robustezza_msp_mv.R",
  "17 Frequentista corexp"   = "17_frequentista_corexp.R",
  "18 Figure principali"     = "18_figure_principali.R",
  "19 Figure report"         = "19_figure_report.R"
)

results <- setNames(vector("list", length(pipeline)), names(pipeline))

for (lbl in names(pipeline)) {
  script_file <- here("scripts", pipeline[[lbl]])
  if (!file.exists(script_file)) {
    cat(sprintf("\n  SKIP %s: file non trovato (%s)\n", lbl, pipeline[[lbl]]))
    results[[lbl]] <- "SKIP (file non trovato)"
    next
  }
  results[[lbl]] <- run_safe(script_file, lbl)
}

# ── 4. RIEPILOGO ──────────────────────────────────────────────────────────────

cat("\n\n")
cat(strrep("═", 60), "\n")
cat("=== RIEPILOGO PIPELINE ===\n")
cat(strrep("═", 60), "\n")
for (lbl in names(results)) {
  status <- results[[lbl]]
  icon   <- if (startsWith(status, "OK")) "✓" else if (startsWith(status, "SKIP")) "─" else "✗"
  cat(sprintf("  %s %-30s: %s\n", icon, lbl, status))
}

n_ok   <- sum(startsWith(unlist(results), "OK"))
n_err  <- sum(startsWith(unlist(results), "ERRORE"))
n_skip <- sum(startsWith(unlist(results), "SKIP"))

cat(strrep("─", 60), "\n")
cat(sprintf("  Completati: %d | Errori: %d | Saltati: %d / %d totali\n",
            n_ok, n_err, n_skip, length(results)))

if (n_err == 0) {
  cat("\n  Tutti gli script completati. Compila il PDF con:\n")
  cat("    cd report && pdflatex report_skeleton.tex\n")
} else {
  cat("\n  Correggere gli errori riportati prima di compilare il PDF.\n")
}
