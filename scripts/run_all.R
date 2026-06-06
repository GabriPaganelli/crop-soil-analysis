# =============================================================================
# run_all.R  —  Master script: esegue l'intera pipeline (07–22)
#
# Esegue in ordine gli script 07–22. I fit (07-09) sono cachati: se esistono
# già vengono solo caricati. Il modello finale è M-SP-RIRS-MVRE (script 10).
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
#   Script 07-09 (M-RI/SP/RIRS MCMC): ~10 min ciascuno se non cachati
#   Script 10 (MVRE LOO):             ~5 min (fit già in cache)
#   Script 12 (projpred):             ~20-35 min
#   Script 13 (validazione):          ~5 min
#   Script 14 (A/B MCMC):             ~10 min ciascuno se non cachati
#   Script 15 (sensitivity MCMC):     ~10 min se non cachato
#   Script 16 (MVRE-FULL LOO):        ~3 min (fit già in cache)
#   Script 17 (frequentista):         ~2 min
#   Script 18-19 (figure + report):   ~10 min
#   Script 20 (spatial confounding):  ~2 min
#   Script 21 (GP robustezza MCMC):   ~30 min se non cachato
#   Script 22 (Lat/Long MCMC):        ~10 min se non cachato
# =============================================================================

library(here)

# ── 0. SETUP ──────────────────────────────────────────────────────────────────

source(here("scripts", "00_utilities.R"))
setup_rtools()

# ── 1. VALIDAZIONE PREREQUISITI ───────────────────────────────────────────────

validate_prerequisites <- function() {
  cat("═══ VALIDAZIONE PREREQUISITI ═══════════════════════════════════\n")

  # Solo i prerequisiti che non vengono prodotti dalla pipeline stessa
  required_files <- list(
    dati      = here("data", "dati.rds"),
    stan_07   = here("stan", "m4rr_v2_no_gp_mu.stan"),
    stan_08   = here("stan", "m4rr_v2_ri_slope_mu.stan"),
    stan_09   = here("stan", "m_sp_rirs.stan"),
    stan_mvre = here("stan", "m_sp_rirs_mvre.stan"),
    stan_A    = here("stan", "m_sp_rirs_mvre_A.stan"),
    stan_B    = here("stan", "m_sp_rirs_mvre_B.stan"),
    stan_full = here("stan", "m_sp_rirs_mvre_full.stan"),
    stan_gp   = here("stan", "m_sp_rirs_mvre_gp.stan")
  )

  for (nm in names(required_files)) {
    if (!file.exists(required_files[[nm]])) {
      stop(sprintf("File mancante: %s\n  Percorso: %s", nm, required_files[[nm]]))
    }
    cat(sprintf("  OK %-12s: %s\n", nm, basename(required_files[[nm]])))
  }

  # Se il fit MVRE esiste già, verifica che i parametri siano accessibili
  if (file.exists(here("stan", "fit_msp_rirs_mvre.rds"))) {
    cat("  Verifica parametri MVRE (fit già in cache)...\n")
    fit_test <- tryCatch(readRDS(here("stan", "fit_msp_rirs_mvre.rds")),
                         error = function(e) stop("Impossibile caricare fit MVRE: ", conditionMessage(e)))
    pnames <- tryCatch(fit_test$metadata()$stan_variables, error = function(e) character(0))
    needed <- c("V", "tau_alpha_SOC", "rho_SOC", "rho_int_SOC_N")
    missing_p <- setdiff(needed, pnames)
    if (length(missing_p) > 0) {
      stop("Parametri mancanti nel fit MVRE: ", paste(missing_p, collapse = ", "),
           "\n  Il fit potrebbe essere corrotto o da una versione diversa dello Stan file.")
    }
    rm(fit_test); gc()
    cat("  OK parametri MVRE verificati\n")
  } else {
    cat("  fit_msp_rirs_mvre.rds non trovato — verrà creato dallo script 10\n")
  }
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
  "07 M-RI (confronto)"      = "07_model_ri.R",
  "08 M-SP (confronto)"      = "08_model_msp.R",
  "09 M-SP-RIRS (confronto)" = "09_model_msp_rirs.R",
  "10 Modello finale MVRE"   = "10_model_final_mvre.R",
  "11 Landuse variante"      = "11_model_msp_landuse.R",
  "12 Projpred (MVRE)"       = "12_selezione_variabili.R",
  "13 Validazione projpred"  = "13_validazione_projpred.R",
  "14 Confronto A/B"         = "14_confronto_modelli.R",
  "15 Sensitivity Pareto"    = "15_sensitivity_pareto.R",
  "16 MV residui"            = "16_robustezza_msp_mv.R",
  "17 Frequentista corexp"   = "17_frequentista_corexp.R",
  "18 Figure principali"     = "18_figure_principali.R",
  "19 Figure report"         = "19_figure_report.R",
  "20 Spatial confounding"   = "20_spatial_confounding.R",
  "21 GP robustezza"         = "21_robustezza_gp.R",
  "22 Lat/Long robustezza"   = "22_robustezza_latlong.R"
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
  cat("    cd report && pdflatex report.tex\n")
} else {
  cat("\n  Correggere gli errori riportati prima di compilare il PDF.\n")
}
