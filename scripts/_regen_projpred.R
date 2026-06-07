suppressPackageStartupMessages({
  library(tidyverse); library(here); library(patchwork); library(projpred)
})
source(here("scripts", "00_utilities.R"))

fig_dir <- here("output", "figures")
img_dir <- here("report", "images")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)

resp_colors <- c("SOC" = "#2166AC", "N" = "#1A9850", "P" = "#D73027")

varsel_list <- readRDS(here("output", "cache", "projpred_varsel_mvre.rds"))
resp_names  <- c("SOC", "N", "P")
pp_plots    <- list()

for (nm in resp_names) {
  vs_obj <- varsel_list[[nm]]
  if (!is.null(vs_obj)) {
    p <- tryCatch(
      plot(vs_obj, stats = "elpd", deltas = TRUE) +
        ggtitle(sprintf("log%s", nm)) +
        labs(x = "N. predittori", y = expression(Delta * "ELPD"),
             subtitle = NULL, caption = NULL) +
        theme_minimal(base_size = 13) +
        theme(
          plot.title       = element_text(face = "bold", size = 14,
                                          colour = resp_colors[nm]),
          plot.subtitle    = element_blank(),
          plot.caption     = element_blank(),
          axis.title       = element_text(size = 12),
          axis.text.x      = element_text(size = 11, angle = 40, hjust = 1),
          axis.text.y      = element_text(size = 11),
          panel.grid.minor = element_blank()
        ),
      error = function(e) { cat(sprintf("Errore plot %s: %s\n", nm, conditionMessage(e))); NULL }
    )
    pp_plots[[nm]] <- p
  }
}

pp_ok <- Filter(Negate(is.null), pp_plots)
p_panel <- wrap_plots(pp_ok, ncol = 3) +
  plot_annotation(
    title = "Selezione variabili (projpred forward search, PSIS-LOO)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

save_fig("fig_17_projpred_panel.pdf", p_panel, w = 24, h = 12)
cat("OK: fig_17_projpred_panel.pdf\n")

file.copy(file.path(fig_dir, "fig_17_projpred_panel.pdf"),
          file.path(img_dir, "fig06_projpred.pdf"), overwrite = TRUE)
cat("Copiato in report/images/fig06_projpred.pdf\n")
