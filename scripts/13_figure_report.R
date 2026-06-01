# =============================================================================
# 13_figure_report.R  —  Export figure per script 10/11 + copia in report/images/
#
# Produce:
#   output/figures/fig_15_loo_ab_comparison.pdf   — LOO: A / B / M-SP
#   output/figures/fig_16_sensitivity_eta.pdf     — η_r: full vs no-influential
#   output/figures/fig_17_projpred_panel.pdf      — selection path × 3 risposte
#
# Alla fine, copia tutte le figure necessarie al report in report/images/
# (con nomi puliti: fig01_aic.pdf … appfig02_trace.pdf)
#
# Dipende da:
#   output/tables/tab_10_loo.csv         (LOO comparison A/B/M-SP)
#   output/tables/tab_11_sensitivity.csv (sensitivity analysis)
#   output/cache/projpred_varsel.rds     (varsel objects SOC/N/P)
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

# Helper salvataggio
save_fig <- function(fname, p, w = 16, h = 9, u = "cm") {
  ggplot2::ggsave(file.path(fig_dir, fname), plot = p,
                  width = w, height = h, units = u, device = "pdf")
  cat(sprintf("  [fig] Salvato: output/figures/%s\n", fname))
}

resp_colors <- c("SOC" = "#2166AC", "N" = "#1A9850", "P" = "#D73027")


# ── 1. LOO: A / B / M-SP ─────────────────────────────────────────────────────

cat("=== Fig 15: LOO comparison A / B / M-SP ===\n")

loo_ab <- read.csv(here("output", "tables", "tab_10_loo.csv")) |>
  mutate(
    modello  = factor(modello, levels = rev(c("B (resp-spec)", "A (ridotto)", "M-SP"))),
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
    title    = "Confronto LOO-CV: M-SP e versioni semplificate",
    subtitle = "Differenze non significative (|ΔELPD/SE| < 0.5 per tutti i confronti)",
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


# ── 2. SENSITIVITY: η_r full vs no-influential ────────────────────────────────

cat("\n=== Fig 16: Sensitivity η_r ===\n")

sens_raw <- read.csv(here("output", "tables", "tab_11_sensitivity.csv"))

eta_df <- sens_raw |>
  filter(grepl("^eta_", variable)) |>
  mutate(
    risposta = sub("eta_", "", variable),
    risposta = factor(risposta, levels = c("SOC", "N", "P"))
  ) |>
  pivot_longer(
    cols      = c(med_full, q05_full, q95_full, med_ni, q05_ni, q95_ni),
    names_to  = "stat",
    values_to = "val"
  ) |>
  mutate(
    tipo = if_else(grepl("_ni$|_ni$|ni", stat), "Senza influenti", "Modello completo"),
    stat_clean = case_when(
      grepl("^med", stat) ~ "med",
      grepl("q05", stat)  ~ "q05",
      grepl("q95", stat)  ~ "q95"
    )
  ) |>
  pivot_wider(
    id_cols   = c(risposta, tipo),
    names_from  = stat_clean,
    values_from = val
  ) |>
  mutate(
    tipo = factor(tipo, levels = c("Modello completo", "Senza influenti")),
    col  = resp_colors[as.character(risposta)]
  )

p_sens <- ggplot(eta_df, aes(x = med, y = tipo, colour = risposta)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 1.2, position = position_dodge(0.0)) +
  geom_point(size = 3.5, position = position_dodge(0.0)) +
  scale_colour_manual(values = resp_colors, name = "Risposta") +
  facet_wrap(~risposta, ncol = 3, scales = "fixed") +
  scale_x_continuous(name = expression(eta[r] ~ "(mediana e CI 90%)")) +
  labs(
    title    = expression("Robustezza di " * eta[r] * " alle osservazioni influenti (Pareto k ≥ 0.7)"),
    subtitle = "10 osservazioni rimosse su 220 (4.5%). La linea tratteggiata indica η = 0.",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 10, colour = "grey40"),
    legend.position    = "none",
    panel.grid.major.y = element_blank(),
    strip.text         = element_text(face = "bold")
  )

save_fig("fig_16_sensitivity_eta.pdf", p_sens, w = 16, h = 7)


# ── 3. PROJPRED PANEL (3 risposte) ────────────────────────────────────────────

cat("\n=== Fig 17: Projpred selection panel ===\n")

varsel_path <- here("output", "cache", "projpred_varsel.rds")

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
            subtitle = "ΔELPD rispetto al reference model M-SP. Linea tratteggiata = soglia α=0.10",
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
  "fig_01_aic_forma_funzionale.pdf"  = "fig01_aic.pdf",
  "fig_03_spaghetti_soc_n.pdf"       = "fig02_spaghetti.pdf",
  "fig_02_scatter_intercetta_slope.pdf" = "fig03_scatter.pdf",
  "fig_04_posterior_eta.pdf"         = "fig04_eta.pdf",
  "fig_05_forest_gamma.pdf"          = "fig05_gamma.pdf",
  "fig_17_projpred_panel.pdf"        = "fig06_projpred.pdf",   # da questo script
  "fig_09_loo_comparison.pdf"        = "fig07_loo4.pdf",
  "fig_07_ppc.pdf"                   = "fig08_ppc.pdf",
  "fig_15_loo_ab_comparison.pdf"     = "fig09_loo_ab.pdf",     # da questo script
  "fig_16_sensitivity_eta.pdf"       = "fig10_sensitivity.pdf",# da questo script
  "fig_06_forest_beta.pdf"           = "appfig01_beta.pdf",
  "fig_10_trace_key.pdf"             = "appfig02_trace.pdf"
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

cat("\n── Fine script 12 ──────────────────────────────────────────────\n")
