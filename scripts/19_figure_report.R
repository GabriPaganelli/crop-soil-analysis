# =============================================================================
# 19_figure_report.R  —  Export figure per script 13/14 + copia in report/images/
#
# Produce:
#   output/figures/fig_15_loo_ab_comparison.pdf    — LOO: A / B / M-SP-RIRS
#   output/figures/fig_16_sensitivity_rho.pdf      — rho_r/tau_beta_r: full vs no-infl
#   output/figures/fig_17_projpred_panel.pdf       — selection path × 3 risposte
#
# Alla fine, copia tutte le figure necessarie al report in report/images/.
#
# Dipende da:
#   output/tables/tab_13_loo.csv           (LOO comparison A/B/M-SP-RIRS)
#   output/tables/tab_14_sensitivity.csv   (sensitivity analysis su M-SP-RIRS)
#   output/cache/projpred_varsel_rirs.rds  (varsel objects SOC/N/P)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(patchwork)
})

dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("report", "images"),  recursive = TRUE, showWarnings = FALSE)

fig_dir <- here("output", "figures")
img_dir <- here("report", "images")

resp_colors <- c("SOC" = "#2166AC", "N" = "#1A9850", "P" = "#D73027")


# ── 1. LOO: A / B / M-SP ─────────────────────────────────────────────────────

cat("=== Fig 15: LOO comparison A / B / M-SP-RIRS-MVRE ===\n")

loo_path <- here("output", "tables", "tab_14_loo.csv")
if (!file.exists(loo_path)) {
  cat("  SKIP: tab_14_loo.csv non trovato (eseguire script 14 prima).\n")
  p_loo_ab <- NULL
} else {

loo_ab <- read.csv(loo_path) |>
  mutate(
    modello  = factor(modello, levels = rev(c("B (resp-spec)", "A (ridotto)", "M-SP-RIRS-MVRE"))),
    ci_lo    = elpd_diff - se_diff,
    ci_hi    = elpd_diff + se_diff
  )

p_loo_ab <- ggplot(loo_ab, aes(x = elpd_diff, y = modello)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 1.1, colour = "steelblue") +
  geom_point(size = 3.5, colour = "steelblue") +
  geom_text(aes(label = sprintf("%.1f ± %.1f", elpd_diff, se_diff)),
            hjust = -0.15, vjust = 0.4, size = 3.2, colour = "grey30") +
  scale_x_continuous(
    limits = c(-15, 15),
    name   = expression(Delta * "ELPD rispetto al modello migliore")
  ) +
  labs(
    title    = "Confronto LOO-CV: M-SP-RIRS-MVRE e versioni semplificate",
    subtitle = "Differenze tra modelli A, B e M-SP-RIRS-MVRE | ΔELPD rispetto al migliore",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

save_fig("fig_15_loo_ab_comparison.pdf", p_loo_ab, w = 14, h = 7)

} # end if loo_path exists


# ── 2. SENSITIVITY: rho_r / tau_beta_r ───────────────────────────────────────

cat("\n=== Fig 16: Sensitivity rho_r / tau_beta_r ===\n")

sens_path <- here("output", "tables", "tab_15_sensitivity.csv")
if (!file.exists(sens_path)) {
  cat("  SKIP: tab_15_sensitivity.csv non trovato (eseguire script 15 prima).\n")
} else {

sens_raw <- read.csv(sens_path)

sens_key <- sens_raw |>
  filter(grepl("^(rho_|tau_beta_)", variable)) |>
  mutate(
    tipo_par = if_else(grepl("^rho_", variable), "rho[r]", "tau[beta[r]]"),
    risposta = factor(
      case_when(grepl("SOC", variable) ~ "SOC", grepl("_N$", variable) ~ "N",
                grepl("_P$", variable) ~ "P"),
      levels = c("SOC", "N", "P")
    )
  ) |>
  pivot_longer(cols = c(med_full, q05_full, q95_full, med_ni, q05_ni, q95_ni),
               names_to = "stat", values_to = "val") |>
  mutate(
    dataset    = if_else(grepl("_ni$|ni", stat), "Senza influenti", "Completo"),
    stat_clean = case_when(grepl("^med",stat)~"med",grepl("q05",stat)~"q05",grepl("q95",stat)~"q95")
  ) |>
  pivot_wider(id_cols = c(risposta, tipo_par, dataset),
              names_from = stat_clean, values_from = val) |>
  mutate(dataset = factor(dataset, levels = c("Completo", "Senza influenti")))

p_sens <- ggplot(sens_key, aes(x = med, y = dataset, colour = risposta)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = resp_colors, name = "Risposta") +
  facet_grid(risposta ~ tipo_par, scales = "free_x", labeller = label_parsed) +
  labs(
    title    = expression("Robustezza di " * rho[r] * " e " * tau[beta[r]] * " (Pareto k ≥ 0.7)"),
    subtitle = "Confronto M-SP-RIRS-MVRE completo vs senza osservazioni influenti",
    x = "Valore (mediana e CI 90%)", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none",
        panel.grid.major.y = element_blank(), strip.text = element_text(face = "bold"))

save_fig("fig_16_sensitivity_rho.pdf", p_sens, w = 16, h = 12)

} # end if sens_path exists


# ── 3. PROJPRED PANEL (3 risposte) ────────────────────────────────────────────

cat("\n=== Fig 17: Projpred selection panel ===\n")

varsel_path <- here("output", "cache", "projpred_varsel_mvre.rds")

if (file.exists(varsel_path)) {
  # Carica pacchetti projpred (solo se già installato, non necessario rieseguire)
  has_projpred <- requireNamespace("projpred", quietly = TRUE)
  if (!has_projpred) {
    cat("  AVVISO: projpred non disponibile. Salto fig_17.\n")
  } else {
    suppressPackageStartupMessages(library(projpred))
    cat("  Carico varsel objects...\n")
    varsel_list <- tryCatch(readRDS(varsel_path), error = function(e) {
      cat("  Errore caricamento varsel:", conditionMessage(e), "\n"); NULL
    })

    if (!is.null(varsel_list)) {
      resp_names <- c("SOC", "N", "P")
      pp_plots   <- list()

      for (nm in resp_names) {
        vs_obj <- varsel_list[[nm]]
        if (!is.null(vs_obj)) {
          p <- tryCatch(
            plot(vs_obj, stats = "elpd", deltas = TRUE) +
              ggtitle(sprintf("log%s", nm)) +
              labs(x = "Numero di predittori", y = "ΔELPD") +
              theme_minimal(base_size = 11) +
              theme(
                plot.title  = element_text(face = "bold", colour = resp_colors[nm]),
                axis.title  = element_text(size = 9),
                axis.text   = element_text(size = 8)
              ),
            error = function(e) {
              cat(sprintf("  Errore plot %s: %s\n", nm, conditionMessage(e))); NULL
            }
          )
          pp_plots[[nm]] <- p
        }
      }

      # Combina con patchwork (solo i pannelli riusciti)
      pp_ok <- Filter(Negate(is.null), pp_plots)
      if (length(pp_ok) >= 2) {
        p_panel <- wrap_plots(pp_ok, ncol = 3) +
          plot_annotation(
            title    = "Selezione variabili (projpred forward search, PSIS-LOO)",
            subtitle = "ΔELPD rispetto al reference model M-SP-RIRS. Linea tratteggiata = soglia α=0.10",
            theme    = theme(
              plot.title    = element_text(face = "bold", size = 12),
              plot.subtitle = element_text(size = 10, colour = "grey40")
            )
          )
        save_fig("fig_17_projpred_panel.pdf", p_panel, w = 22, h = 8)
      } else {
        cat("  AVVISO: meno di 2 pannelli prodotti, fig_17 non salvata.\n")
        cat("  I 3 grafici individuali in fig_08_projpred_*.pdf rimangono validi.\n")
      }
    }
  }
} else {
  cat("  AVVISO: projpred_varsel.rds non trovato. Salto fig_17.\n")
  cat("  Usa fig_08_projpred_soc/n/p.pdf come alternativa individuale.\n")
}


# ── 4. COPIA FIGURE IN report/images/ ─────────────────────────────────────────

cat("\n=== Copia figure in report/images/ ===\n")

# Mappa: nome_sorgente (in output/figures/) → nome_destinazione (in report/images/)
# Per le 3 figure da script 12, la sorgente è già in fig_dir.
fig_map <- list(
  "fig_01_aic_forma_funzionale.pdf"     = "fig01_aic.pdf",
  "fig_03_spaghetti_soc_n.pdf"          = "fig02_spaghetti.pdf",
  "fig_02_scatter_intercetta_slope.pdf" = "fig03_scatter.pdf",
  "fig_04_posterior_rho.pdf"            = "fig04_rho.pdf",
  "fig_05_forest_gamma.pdf"             = "fig05_gamma.pdf",
  "fig_17_projpred_panel.pdf"           = "fig06_projpred.pdf",
  "fig_09_loo_comparison.pdf"           = "fig07_loo4.pdf",
  "fig_07_ppc.pdf"                      = "fig08_ppc.pdf",
  "fig_15_loo_ab_comparison.pdf"        = "fig09_loo_ab.pdf",
  "fig_18_cross_corr.pdf"               = "fig10_cross_corr.pdf",
  "fig_16_sensitivity_rho.pdf"          = "fig11_sensitivity.pdf",
  "fig_13_proj_vs_msp.pdf"              = "fig12_proj_vs_msp.pdf",
  "fig_14_nonselected_gamma.pdf"        = "fig13_nonselected_gamma.pdf",
  "fig_texture_triangle.pdf"            = "fig_texture_triangle.pdf",
  "fig_06_forest_beta.pdf"              = "appfig01_beta.pdf",
  "fig_10_trace_key.pdf"                = "appfig02_trace.pdf"
)

for (src_name in names(fig_map)) {
  src  <- file.path(fig_dir, src_name)
  dst  <- file.path(img_dir, fig_map[[src_name]])
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    cat(sprintf("  Copiato: %s → images/%s\n", src_name, fig_map[[src_name]]))
  } else {
    cat(sprintf("  MANCANTE: %s  (non copiato)\n", src_name))
  }
}

cat(sprintf("\nContenuto report/images/ (%d file):\n",
            length(list.files(img_dir))))
for (f in sort(list.files(img_dir))) cat(sprintf("  %s\n", f))

cat("\n── Fine script 19 ──────────────────────────────────────────────\n")
