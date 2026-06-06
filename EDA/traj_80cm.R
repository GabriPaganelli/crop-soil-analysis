# Traiettorie — solo campi con osservazione a 80 cm

library(tidyverse)
library(patchwork)

dati <- readRDS(here::here("data", "crop.rds")) |>
  mutate(
    logSOC = log(PercSOC),
    logN   = log(PercTotNitro),
    logP   = log(PercTotPhos)
  )

# Campi con almeno un'osservazione a Bottom == 80
campi_80 <- dati |> filter(Bottom == 80) |> pull(Field) |> unique()
dati_80  <- filter(dati, Field %in% campi_80)

message(sprintf("Campi con osservazione a 80 cm: %d su %d totali",
                length(campi_80), n_distinct(dati$Field)))

depth_num <- c(20, 30, 40, 50, 60, 80)
col_trio  <- c("SOC" = "#C0392B", "N" = "#2471A3", "P" = "#1E8449")

dati_traj <- dati_80 |>
  pivot_longer(c(logSOC, logN, logP),
               names_to  = "Variabile",
               values_to = "logVal") |>
  mutate(
    Variabile = recode(Variabile, logSOC = "SOC", logN = "N", logP = "P"),
    Variabile = factor(Variabile, levels = c("SOC", "N", "P"))
  )

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
  d   <- filter(dati_traj,    Variabile == var) |> arrange(Field, Bottom)
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

p_SOC <- make_traj_panel("SOC", show_y_title = TRUE)
p_N   <- make_traj_panel("N")
p_P   <- make_traj_panel("P")

n_campi <- length(campi_80)

fig <- (p_SOC + p_N + p_P) +
  plot_layout(ncol = 3) +
  plot_annotation(
    title    = sprintf("Profili verticali di SOC, N e P — %d campi con profilo completo (fino a 80 cm)", n_campi),
    theme    = theme(
      plot.title      = element_text(size = 13, face = "bold", hjust = 0.5,
                                     color = "#1A1A1A"),
      plot.background = element_rect(fill = "#FAFAF8", color = NA)
    )
  )

out <- here::here("EDA", "traj_80cm.pdf")
ggsave(out, fig, width = 28, height = 18, units = "cm")
message("Salvato: ", out)
