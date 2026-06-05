# ============================================================
#  EDA — Traiettorie per covariata
#
#  Versioni alternative di eda_03_traiettorie:
#  linee mediane separate per gruppi di covariata.
#
#  Output (EDA/):
#    eda_traj_onfarm.pdf      OnFarm vs Off-farm (2 gruppi)
#    eda_traj_irrigate.pdf    Irrigato vs Non irrigato (2 gruppi)
#    eda_traj_fertilised.pdf  Fertilizzato vs Non (2 gruppi)
#    eda_traj_gestione.pdf    7 classi di gestione
#    eda_traj_gestione3.pdf   3 macro-gruppi (Natural / Farm / Off-farm)
# ============================================================

library(tidyverse)
library(patchwork)

dati_raw <- readRDS(here::here("data", "crop.rds")) |>
  mutate(
    logSOC   = log(PercSOC),
    logN     = log(PercTotNitro),
    logP     = log(PercTotPhos),
    Gestione = factor(Landuse,
                      levels = c(4, 1, 2, 3, 7, 6, 5),
                      labels = c("Natural", "Farm NI", "Farm I",
                                 "Off-farm NINF", "Off-farm NIF",
                                 "Off-farm INF",  "Off-farm IF")),
    MacroGest = case_when(
      Landuse == 4             ~ "Natural",
      Landuse %in% c(1, 2)    ~ "Farm",
      TRUE                    ~ "Off-farm"
    ) |> factor(levels = c("Natural", "Farm", "Off-farm"))
  )

col_trio  <- c("SOC" = "#C0392B", "N" = "#2471A3", "P" = "#1E8449")
depth_num <- c(20, 30, 40, 50, 60, 80)

# ── dati long (tutte le covariate incluse) ────────────────────
dati_long <- dati_raw |>
  pivot_longer(c(logSOC, logN, logP),
               names_to  = "Variabile",
               values_to = "logVal") |>
  mutate(
    Variabile = recode(Variabile, logSOC = "SOC", logN = "N", logP = "P"),
    Variabile = factor(Variabile, levels = c("SOC", "N", "P"))
  ) |>
  arrange(Variabile, Field, Bottom)

x_labels <- c("SOC" = "log(SOC %)", "N" = "log(N %)", "P" = "log(P %)")

th_traj <- theme_minimal(base_size = 11) +
  theme(
    plot.background    = element_rect(fill = "#FAFAF8", color = NA),
    panel.background   = element_rect(fill = "#F4F2ED", color = NA),
    panel.grid.major.y = element_line(color = "#E5E2DB", linewidth = 0.3,
                                      linetype = "dashed"),
    panel.grid.major.x = element_line(color = "#EEEBE4", linewidth = 0.25),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(size = 9, color = "#555555"),
    axis.title         = element_text(size = 10, color = "#333333"),
    plot.title         = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title       = element_text(size = 9, face = "bold"),
    legend.text        = element_text(size = 8.5),
    plot.margin        = margin(10, 12, 8, 10)
  )

# ── funzione pannello ─────────────────────────────────────────
make_cov_panel <- function(var, cov_col, pal,
                           show_y_title = FALSE, show_legend = FALSE) {

  d   <- filter(dati_long, Variabile == var)
  med <- d |>
    group_by(Bottom, .data[[cov_col]]) |>
    summarise(med = median(logVal, na.rm = TRUE), .groups = "drop") |>
    arrange(.data[[cov_col]], Bottom)

  ggplot(d, aes(x = logVal, y = Bottom, group = Field)) +
    geom_path(color = "#BBBBBB", linewidth = 0.45, alpha = 0.4) +
    geom_path(data = med,
              aes(x = med, y = Bottom,
                  color = .data[[cov_col]],
                  group = .data[[cov_col]]),
              linewidth = 1.3, lineend = "round",
              inherit.aes = FALSE) +
    geom_point(data = med,
               aes(x = med, y = Bottom,
                   color = .data[[cov_col]],
                   group = .data[[cov_col]]),
               size = 2.2, shape = 21, fill = "#FAFAF8", stroke = 1.1,
               inherit.aes = FALSE) +
    scale_y_reverse(
      breaks = depth_num,
      labels = paste0(depth_num, " cm"),
      expand = expansion(mult = c(0.04, 0.06))
    ) +
    scale_color_manual(values = pal, name = cov_col) +
    labs(title = var,
         x     = x_labels[var],
         y     = if (show_y_title) "Profondità (cm)" else NULL) +
    th_traj +
    theme(
      legend.position = if (show_legend) "right" else "none",
      axis.title.y    = if (show_y_title) element_text() else element_blank(),
      axis.text.y     = element_text(size = 9, color = "#555555")
    )
}

# ── helper per assemblare figura 3-pannelli ───────────────────
make_fig <- function(cov_col, pal, title) {
  p1 <- make_cov_panel("SOC", cov_col, pal, show_y_title = TRUE)
  p2 <- make_cov_panel("N",   cov_col, pal)
  p3 <- make_cov_panel("P",   cov_col, pal, show_legend = TRUE)

  (p1 + p2 + p3) +
    plot_layout(ncol = 3, widths = c(1, 1, 1.28)) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title      = element_text(size = 14, face = "bold", hjust = 0.5,
                                       color = "#1A1A1A"),
        plot.background = element_rect(fill = "#FAFAF8", color = NA)
      )
    )
}

save_fig <- function(fname, p, w = 30, h = 18, u = "cm") {
  path <- here::here("EDA", fname)
  ggsave(path, plot = p, width = w, height = h, units = u,
         device = cairo_pdf, bg = "#FAFAF8")
  message("✓ ", path)
}

# ══════════════════════════════════════════════════════════════
#  1. OnFarm
# ══════════════════════════════════════════════════════════════
dati_long <- dati_long |>
  mutate(OnFarm_f = factor(OnFarm, levels = c(0, 1),
                           labels = c("Off-farm", "On-farm")))

pal_onfarm2 <- c("Off-farm" = "#5B87C0", "On-farm" = "#C97D30")

fig_onfarm <- make_fig("OnFarm_f", pal_onfarm2,
                       "Profili verticali per gestione On-farm vs Off-farm")
save_fig("eda_traj_onfarm.pdf", fig_onfarm)

# ══════════════════════════════════════════════════════════════
#  2. Irrigate
# ══════════════════════════════════════════════════════════════
dati_long <- dati_long |>
  mutate(Irrigate_f = factor(Irrigate, levels = c(0, 1),
                             labels = c("Non irrigato", "Irrigato")))

pal_irr <- c("Non irrigato" = "#888888", "Irrigato" = "#2471A3")

fig_irr <- make_fig("Irrigate_f", pal_irr,
                    "Profili verticali per irrigazione")
save_fig("eda_traj_irrigate.pdf", fig_irr)

# ══════════════════════════════════════════════════════════════
#  3. Fertilised
# ══════════════════════════════════════════════════════════════
dati_long <- dati_long |>
  mutate(Fertilised_f = factor(Fertilised, levels = c(0, 1),
                               labels = c("Non fertilizzato", "Fertilizzato")))

pal_fert <- c("Non fertilizzato" = "#888888", "Fertilizzato" = "#1E8449")

fig_fert <- make_fig("Fertilised_f", pal_fert,
                     "Profili verticali per fertilizzazione")
save_fig("eda_traj_fertilised.pdf", fig_fert)

# ══════════════════════════════════════════════════════════════
#  4. Gestione (7 classi)
# ══════════════════════════════════════════════════════════════
dati_long <- dati_long |>
  mutate(Gestione = factor(Landuse,
                           levels = c(4, 1, 2, 3, 7, 6, 5),
                           labels = c("Natural", "Farm NI", "Farm I",
                                      "Off-farm NINF", "Off-farm NIF",
                                      "Off-farm INF",  "Off-farm IF")))

pal_gest7 <- c(
  "Natural"       = "#4E9B6F",
  "Farm NI"       = "#C97D30",
  "Farm I"        = "#E8B84B",
  "Off-farm NINF" = "#5B87C0",
  "Off-farm NIF"  = "#8CB4D8",
  "Off-farm INF"  = "#9B6BAD",
  "Off-farm IF"   = "#C0392B"
)

fig_gest7 <- make_fig("Gestione", pal_gest7,
                      "Profili verticali per classe di gestione (7 classi)")
save_fig("eda_traj_gestione.pdf", fig_gest7, w = 32, h = 18)

# ══════════════════════════════════════════════════════════════
#  5. Macro-gruppo (Natural / Farm / Off-farm)
# ══════════════════════════════════════════════════════════════
dati_long <- dati_long |>
  mutate(MacroGest = case_when(
    Landuse == 4          ~ "Natural",
    Landuse %in% c(1, 2) ~ "Farm",
    TRUE                  ~ "Off-farm"
  ) |> factor(levels = c("Natural", "Farm", "Off-farm")))

pal_macro <- c("Natural" = "#4E9B6F", "Farm" = "#C97D30", "Off-farm" = "#5B87C0")

fig_macro <- make_fig("MacroGest", pal_macro,
                      "Profili verticali per macro-gruppo di gestione")
save_fig("eda_traj_gestione3.pdf", fig_macro)

message("\n✓ Tutte le figure covariate generate in EDA/")
