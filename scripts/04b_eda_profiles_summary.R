# =============================================================================
# 04b_eda_profiles_summary.R  —  EDA: soil chemical-physical properties
#
#  Output (output/figures/):
#    eda_01_distributions.pdf   istogramma + densità (scala originale e log)
#    eda_02_management.pdf        boxplot + beeswarm per classe di gestione
#    eda_03_trajectories.pdf     spaghetti plot profondità × campo × variabile
#    eda_04_trajectories_orig.pdf spaghetti plot scala originale (%)
#
#  Pacchetti: tidyverse, patchwork, ggbeeswarm
# ============================================================

library(tidyverse)
library(patchwork)
library(ggbeeswarm)

# ── 0. Dati ──────────────────────────────────────────────────

dati <- readRDS(here::here("data", "crop_full.rds")) |>
  mutate(
    logSOC   = log(PercSOC),
    logN     = log(PercTotNitro),
    logP     = log(PercTotPhos),
    Gestione = factor(Landuse,
                      levels = c(4, 1, 2, 3, 7, 6, 5),
                      labels = c("Natural", "Farm NI", "Farm I",
                                 "Off-farm NINF", "Off-farm NIF",
                                 "Off-farm INF",  "Off-farm IF"))
  )

col_soc  <- "#C0392B"
col_n    <- "#2471A3"
col_p    <- "#1E8449"
col_trio <- c("SOC" = col_soc, "N" = col_n, "P" = col_p)

pal_gest <- c(
  "Natural"       = "#4E9B6F",
  "Farm NI"       = "#C97D30",
  "Farm I"        = "#E8B84B",
  "Off-farm NINF" = "#5B87C0",
  "Off-farm NIF"  = "#8CB4D8",
  "Off-farm INF"  = "#9B6BAD",
  "Off-farm IF"   = "#C9A0DC"
)

th <- theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "#FAFAF8", color = NA),
    panel.background = element_rect(fill = "#F4F2ED", color = NA),
    panel.grid.major = element_line(color = "#E5E2DB", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 11, color = "#1A1A1A"),
    plot.subtitle    = element_text(size = 8.5, color = "#555555",
                                    margin = margin(b = 5)),
    axis.text        = element_text(size = 9, color = "#444444"),
    axis.title       = element_text(size = 9.5, color = "#333333"),
    strip.text       = element_text(face = "bold", size = 10,
                                    color = "#1A1A1A"),
    plot.margin      = margin(10, 14, 8, 10)
  )

save_fig <- function(fname, p, w = 28, h = 20, u = "cm") {
  path <- here::here("output", "figures", fname)
  ggsave(path, plot = p, width = w, height = h, units = u,
         device = cairo_pdf, bg = "#FAFAF8")
  message("✓ ", path)
}

# long format (usato in entrambe le figure)
dati_long <- dati |>
  pivot_longer(c(logSOC, logN, logP),
               names_to  = "Risposta",
               values_to = "logVal") |>
  mutate(
    Val      = exp(logVal),
    Risposta = recode(Risposta, logSOC = "SOC", logN = "N", logP = "P"),
    Risposta = factor(Risposta, levels = c("SOC", "N", "P"))
  )


# ══════════════════════════════════════════════════════════════
#  FIGURA 1 — Distribuzioni: scala originale e log
#  Riga A: istogramma + densità su scala originale
#  Riga B: istogramma + densità su scala log
# ══════════════════════════════════════════════════════════════

p_raw <- dati_long |>
  ggplot(aes(x = Val, fill = Risposta, color = Risposta)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 25, alpha = 0.35,
                 color = NA, position = "identity") +
  geom_density(linewidth = 0.9, alpha = 0) +
  scale_fill_manual(values  = col_trio, guide = "none") +
  scale_color_manual(values = col_trio, guide = "none") +
  facet_wrap(~ Risposta, scales = "free", ncol = 3) +
  labs(
    title = "A  Scala originale (%)",
    x = "Concentrazione (%)", y = "Densità"
  ) +
  th

p_log <- dati_long |>
  ggplot(aes(x = logVal, fill = Risposta, color = Risposta)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 25, alpha = 0.35,
                 color = NA, position = "identity") +
  geom_density(linewidth = 0.9, alpha = 0) +
  geom_rug(alpha = 0.25, linewidth = 0.3) +
  scale_fill_manual(values  = col_trio, guide = "none") +
  scale_color_manual(values = col_trio, guide = "none") +
  facet_wrap(~ Risposta, scales = "free_y", ncol = 3) +
  labs(
    title = "B  Scala logaritmica",
    x = "log(Concentrazione)", y = "Densità"
  ) +
  th

fig1 <- (p_raw / p_log) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title   = "Distribuzione delle variabili risposta",
    caption = "N = 220 osservazioni · SOC = carbonio organico · N = azoto totale · P = fosforo totale",
    theme   = theme(
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.caption    = element_text(size = 7.5, color = "#999999", hjust = 0.5),
      plot.background = element_rect(fill = "#FAFAF8", color = NA)
    )
  )

save_fig("eda_01_distributions.pdf", fig1, w = 28, h = 22)


# ══════════════════════════════════════════════════════════════
#  FIGURA 2 — Differenze tra classi di gestione
#  Boxplot + beeswarm, un pannello per risposta
# ══════════════════════════════════════════════════════════════

p_gest <- dati_long |>
  ggplot(aes(x = Gestione, y = logVal,
             color = Gestione, fill = Gestione)) +
  geom_boxplot(
    alpha         = 0.15,
    outlier.shape = NA,
    linewidth     = 0.6,
    width         = 0.5
  ) +
  ggbeeswarm::geom_beeswarm(
    size     = 1.4,
    alpha    = 0.55,
    cex      = 0.65,
    priority = "random"
  ) +
  scale_color_manual(values = pal_gest, guide = "none") +
  scale_fill_manual(values  = pal_gest, guide = "none") +
  facet_wrap(~ Risposta, scales = "free_y", ncol = 1) +
  labs(
    title    = "Distribuzione per classe di gestione",
    subtitle = "Boxplot + osservazioni individuali · tutte le profondità aggregate",
    caption  = "N = 220 · J = 40 campi · log-scala",
    x        = NULL,
    y        = "log(Concentrazione)"
  ) +
  th +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 9),
    strip.text      = element_text(face = "bold", size = 11, color = "#1A1A1A"),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 9, color = "#555555", hjust = 0.5),
    plot.caption    = element_text(size = 7.5, color = "#999999", hjust = 0.5),
    plot.background = element_rect(fill = "#FAFAF8", color = NA)
  )

save_fig("eda_02_management.pdf", p_gest, w = 22, h = 30)


# ══════════════════════════════════════════════════════════════
#  FIGURA 3 — Traiettorie per campo × profondità
#  Spaghetti plot: 3 pannelli (SOC / N / P), 40 linee per pannello
#  Asse Y = profondità cm (invertito, orientamento soil-science)
#  Asse X = log(concentrazione %)
#  Colore = quartile intercetta OLS del campo (Q1 blu → Q4 rosso)
#  Linea nera spessa = mediana globale per profondità
# ══════════════════════════════════════════════════════════════

depth_num <- c(20, 30, 40, 50, 60, 80)

# Long format con log-valori originali
dati_traj <- dati |>
  pivot_longer(c(logSOC, logN, logP),
               names_to  = "Variabile",
               values_to = "logVal") |>
  mutate(
    Variabile = recode(Variabile, logSOC = "SOC", logN = "N", logP = "P"),
    Variabile = factor(Variabile, levels = c("SOC", "N", "P"))
  )

# Mediana + banda IQR per profondità × variabile
mediana_glob <- dati_traj |>
  group_by(Variabile, Bottom) |>
  summarise(med = median(logVal, na.rm = TRUE),
            q25 = quantile(logVal, 0.25, na.rm = TRUE),
            q75 = quantile(logVal, 0.75, na.rm = TRUE),
            .groups = "drop") |>
  arrange(Variabile, Bottom)

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
    plot.margin        = margin(10, 12, 8, 10)
  )

make_traj_panel <- function(var, show_y_title = FALSE) {
  d   <- filter(dati_traj,   Variabile == var) |> arrange(Field, Bottom)
  med <- filter(mediana_glob, Variabile == var) |> arrange(Bottom)
  col <- col_trio[var]

  ggplot(d, aes(x = logVal, y = Bottom, group = Field)) +
    geom_path(color = "#999999", linewidth = 0.5, alpha = 0.45) +
    geom_point(color = "#999999", size = 1.0, alpha = 0.4) +
    geom_ribbon(data = med, aes(xmin = q25, xmax = q75, y = Bottom, group = 1),
                fill = col, alpha = 0.18, inherit.aes = FALSE) +
    geom_path(data = med, aes(x = med, y = Bottom, group = 1),
              color = col, linewidth = 1.4, lineend = "round",
              inherit.aes = FALSE) +
    geom_point(data = med, aes(x = med, y = Bottom, group = 1),
               color = col, size = 2.8, shape = 21,
               fill = "#FAFAF8", stroke = 1.3,
               inherit.aes = FALSE) +
    scale_y_reverse(
      breaks = depth_num,
      labels = paste0(depth_num, " cm"),
      expand = expansion(mult = c(0.04, 0.06))
    ) +
    labs(
      title = var,
      x     = x_labels[var],
      y     = if (show_y_title) "Profondità (cm)" else NULL
    ) +
    th_traj +
    theme(
      axis.title.y = if (show_y_title) element_text() else element_blank(),
      axis.text.y  = element_text(size = 9, color = "#555555")
    )
}

p3_SOC <- make_traj_panel("SOC", show_y_title = TRUE)
p3_N   <- make_traj_panel("N")
p3_P   <- make_traj_panel("P")

fig3 <- (p3_SOC + p3_N + p3_P) +
  plot_layout(ncol = 3) +
  plot_annotation(
    title    = "Profili verticali di SOC, N e P nei 40 campi",
    theme    = theme(
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5,
                                     color = "#1A1A1A"),
      plot.subtitle   = element_text(size = 9, color = "#555555", hjust = 0.5),
      plot.caption    = element_text(size = 7.5, color = "#999999", hjust = 0.5),
      plot.background = element_rect(fill = "#FAFAF8", color = NA)
    )
  )

save_fig("eda_03_trajectories.pdf", fig3, w = 28, h = 18)


# ══════════════════════════════════════════════════════════════
#  FIGURA 4 — Traiettorie su scala originale (%)
#  Identica alla figura 3 ma asse X = concentrazione %
# ══════════════════════════════════════════════════════════════

dati_traj_orig <- dati |>
  pivot_longer(c(PercSOC, PercTotNitro, PercTotPhos),
               names_to  = "Variabile",
               values_to = "Val") |>
  mutate(
    Variabile = recode(Variabile,
                       PercSOC = "SOC", PercTotNitro = "N", PercTotPhos = "P"),
    Variabile = factor(Variabile, levels = c("SOC", "N", "P"))
  )

mediana_glob_orig <- dati_traj_orig |>
  group_by(Variabile, Bottom) |>
  summarise(med = median(Val, na.rm = TRUE),
            q25 = quantile(Val, 0.25, na.rm = TRUE),
            q75 = quantile(Val, 0.75, na.rm = TRUE),
            .groups = "drop") |>
  arrange(Variabile, Bottom)

x_labels_orig <- c("SOC" = "SOC (%)", "N" = "N (%)", "P" = "P (%)")

make_traj_panel_orig <- function(var, show_y_title = FALSE) {
  d   <- filter(dati_traj_orig,   Variabile == var) |> arrange(Field, Bottom)
  med <- filter(mediana_glob_orig, Variabile == var) |> arrange(Bottom)
  col <- col_trio[var]

  ggplot(d, aes(x = Val, y = Bottom, group = Field)) +
    geom_path(color = "#999999", linewidth = 0.5, alpha = 0.45) +
    geom_point(color = "#999999", size = 1.0, alpha = 0.4) +
    geom_ribbon(data = med, aes(xmin = q25, xmax = q75, y = Bottom, group = 1),
                fill = col, alpha = 0.18, inherit.aes = FALSE) +
    geom_path(data = med, aes(x = med, y = Bottom, group = 1),
              color = col, linewidth = 1.4, lineend = "round",
              inherit.aes = FALSE) +
    geom_point(data = med, aes(x = med, y = Bottom, group = 1),
               color = col, size = 2.8, shape = 21,
               fill = "#FAFAF8", stroke = 1.3,
               inherit.aes = FALSE) +
    scale_y_reverse(
      breaks = depth_num,
      labels = paste0(depth_num, " cm"),
      expand = expansion(mult = c(0.04, 0.06))
    ) +
    labs(
      title = var,
      x     = x_labels_orig[var],
      y     = if (show_y_title) "Profondità (cm)" else NULL
    ) +
    th_traj +
    theme(
      axis.title.y = if (show_y_title) element_text() else element_blank(),
      axis.text.y  = element_text(size = 9, color = "#555555")
    )
}

p4_SOC <- make_traj_panel_orig("SOC", show_y_title = TRUE)
p4_N   <- make_traj_panel_orig("N")
p4_P   <- make_traj_panel_orig("P")

fig4 <- (p4_SOC + p4_N + p4_P) +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Profili verticali di SOC, N e P nei 40 campi (scala originale)",
    theme = theme(
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5,
                                     color = "#1A1A1A"),
      plot.background = element_rect(fill = "#FAFAF8", color = NA)
    )
  )

save_fig("eda_04_trajectories_orig.pdf", fig4, w = 28, h = 18)


message("\n✓ EDA completata — output in output/figures/:")
message("  eda_01_distributions.pdf")
message("  eda_02_management.pdf")
message("  eda_03_trajectories.pdf")
message("  eda_04_trajectories_orig.pdf")
