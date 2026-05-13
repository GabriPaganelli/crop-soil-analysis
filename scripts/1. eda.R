library(tidyverse)
library(GGally)
library(patchwork)
library(ggpubr)
library(ggtern)
library(soiltexture)
library(here)

# ── 1. CARICAMENTO ────────────────────────────────────────────────────────────
crop = readRDS(here("data", "crop.rds"))
source(here('scripts', 'utilities.R'))

target_vars <- c("PercSOC", "PercTotNitro", "PercTotPhos")

# ── 3. ANALISI ESPLORATIVA UNIVARIATA ─────────────────────────────────────────
analizza_dataset(crop)
grafico_distribuzioni(crop[, !names(crop) %in% c("Field", "Plot", "Lat", "Long")])

# ── Triangolo di texture ──────────────────────────────────────────────────────
plot_texture_triangle(crop, color_var = "Landuse")
plot_texture_triangle(crop, color_var = "Landuse", version = "ttplot")

# Correlazioni e scatter di ogni target vs tutti i predittori
plot_x_vs_y(crop, target_vars)


# ── 4. PROFILI VERTICALI (depth profiles) ─────────────────────────────────────
# Porta i crop in formato long (usato anche nelle sezioni successive)
crop_long <- crop |>
  select(Field, Bottom, Lat, Long, all_of(target_vars), Landuse) |>
  pivot_longer(cols = all_of(target_vars),
               names_to = "target",
               values_to = "valore")

# Per Field (40 livelli → viridis_d automatico): mostra variabilità individuale
plot_depth_profiles(crop, target_vars, color_var = "Field", ncol = 1)

# Per Landuse (7 livelli → Set1): mostra separazione per classe gestionale
plot_depth_profiles(crop, target_vars, color_var = "Landuse", ncol = 1)


# ── 5. STRUTTURA GERARCHICA ───────────────────────────────────────────────────
# Varianza tra Field vs dentro Field (guidata da profondità)

# Boxplot per Landuse (7 categorie, leggibile) x profondità
ggplot(crop_long,
       aes(x = Landuse, y = valore, fill = as.factor(Bottom))) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_wrap(~ target, scales = "free_y") +
  labs(title = "Distribuzione per Landuse e profondità",
       x = "Landuse", y = "Valore", fill = "Profondità (cm)") +
  theme_minimal()

# Rapporto sd_between / sd_within per quantificare la gerarchia
crop_long |>
  group_by(target, Field) |>
  summarise(media = mean(valore), sd_within = sd(valore), .groups = "drop") |>
  group_by(target) |>
  summarise(
    sd_between      = sd(media),
    sd_within_media = mean(sd_within),
    rapporto        = sd_between / sd_within_media
  )
# PercSOC: rapporto < 1 → variabilità guidata dalla profondità più che dal field
# PercTotNitro, PercTotPhos: rapporto ≈ 1 → profondità e field contano in modo simile

# ── 6. CORRELAZIONE TRA LE TRE VARIABILI TARGET ───────────────────────────────

# Colorato per Landuse
ggpairs(
  crop,
  columns  = target_vars,
  aes(colour = Landuse, alpha = 0.6),
  upper = list(continuous = wrap("cor", size = 3)),
  lower = list(continuous = wrap("points", size = 0.8)),
  diag  = list(continuous = wrap("densityDiag", alpha = 0.5))
) +
  labs(title = "Correlazioni tra variabili target — per Landuse") +
  theme_minimal()

# Colorato per Fertilised
ggpairs(
  crop,
  columns  = target_vars,
  aes(colour = Fertilised, alpha = 0.6),
  upper = list(continuous = wrap("cor", size = 3)),
  lower = list(continuous = wrap("points", size = 0.8)),
  diag  = list(continuous = wrap("densityDiag", alpha = 0.5))
) +
  labs(title = "Correlazioni tra variabili target — per Fertilised") +
  theme_minimal()

# Colorato per Texture
ggpairs(
  crop,
  columns  = target_vars,
  aes(colour = Texture, alpha = 0.6),
  upper = list(continuous = wrap("cor", size = 3)),
  lower = list(continuous = wrap("points", size = 0.8)),
  diag  = list(continuous = wrap("densityDiag", alpha = 0.5))
) +
  labs(title = "Correlazioni tra variabili target — per Texture") +
  theme_minimal()

# ── 7. EFFETTO SPAZIALE ───────────────────────────────────────────────────────

# Un plot per variabile (scala colore separata, altrimenti dominata da PercSOC)
plots_spaziali <- map(target_vars, ~ {
  crop_long |>
    filter(target == .x) |>
    ggplot(aes(x = Long, y = Lat, colour = valore)) +
    geom_point(alpha = 0.7, size = 3) +
    scale_colour_viridis_c(option = "magma", direction = -1) +
    labs(title = .x, x = "Longitudine", y = "Latitudine", colour = "Valore") +
    theme_minimal()
})
wrap_plots(plots_spaziali, nrow = 1)

# Correlazione numerica con le coordinate
crop |>
  summarise(across(all_of(target_vars), list(
    cor_Lat  = ~ cor(.x, Lat,  use = "complete.obs"),
    cor_Long = ~ cor(.x, Long, use = "complete.obs")
  ))) |>
  pivot_longer(everything(),
               names_to  = c("variabile", "coord"),
               names_sep = "_cor_",
               values_to = "correlazione")


# 8. Plot 1:5 appartengono alla stessa farm?
crop |>
  mutate(
    Field_Num   = as.numeric(gsub("\\D", "", as.character(Field))),
    Field_Group = as.factor((Field_Num - 1) %/% 5)
  ) |>
  ggplot(aes(x = Long, y = Lat, color = Field_Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  scale_color_viridis_d(option = "turbo") +
  labs(
    title    = "Distribuzione Lat/Long raggruppata per Field",
    subtitle = "I colori cambiano ogni 5 campi (es. 1-5, 6-10...)",
    x        = "Longitudine",
    y        = "Latitudine",
    color    = "Blocco Field"
  ) +
  theme(panel.grid.minor = element_blank())


# ── 9. LINEARITÀ DI LOG(TARGET) VS BOTTOM ────────────────────────────────────
# Per ogni field × target: R² del modello lineare e test F lineare vs quadratico.
# Se R² è alto e il termine quadratico non è significativo (p ≥ 0.05),
# la relazione log-lineare con la profondità è supportata.

lin_results <- crop |>
  pivot_longer(cols = all_of(target_vars), names_to = "target", values_to = "value") |>
  group_by(Field, target) |>
  summarise(
    r2 = {
      y  <- log(value)
      x  <- Bottom
      ok <- is.finite(y) & !is.na(x)
      if (sum(ok) < 4) NA_real_
      else summary(lm(y[ok] ~ x[ok]))$r.squared
    },
    p_quad = {
      y  <- log(value)
      x  <- Bottom
      ok <- is.finite(y) & !is.na(x)
      if (sum(ok) < 5) NA_real_
      else anova(lm(y[ok] ~ x[ok]),
                 lm(y[ok] ~ x[ok] + I(x[ok]^2)))$`Pr(>F)`[2]
    },
    .groups = "drop"
  ) |>
  mutate(sig_quad = p_quad < 0.05)

# Etichette leggibili per i grafici
target_labels <- c(
  PercSOC      = "SOC (%)",
  PercTotNitro = "Azoto totale (%)",
  PercTotPhos  = "Fosforo totale (%)"
)

# ── Pannello 1: R² ────────────────────────────────────────────────────────────
make_r2_plot <- function(tgt) {
  lin_results |>
    filter(target == tgt) |>
    mutate(Field = reorder(Field, r2)) |>
    ggplot(aes(x = r2, y = Field)) +
    geom_point(size = 2, colour = "#2166ac") +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "gray50") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.7, 0.75, 1)) +
    labs(title = target_labels[tgt], x = "R²", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
}

wrap_plots(map(target_vars, make_r2_plot), nrow = 1) +
  plot_annotation(
    title    = "R² di lm(log(target) ~ Bottom) per field",
    subtitle = "Linea tratteggiata = soglia orientativa 0.7",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

# ── Pannello 2: test F (lineare vs quadratico) ────────────────────────────────
make_ftest_plot <- function(tgt) {
  lin_results |>
    filter(target == tgt) |>
    mutate(Field = reorder(Field, p_quad)) |>
    ggplot(aes(x = p_quad, y = Field, colour = sig_quad)) +
    geom_point(size = 2) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = "gray50") +
    scale_colour_manual(
      values = c("TRUE" = "#d73027", "FALSE" = "#1a9850"),
      labels = c("TRUE" = "p < 0.05  (quadratico significativo)",
                 "FALSE" = "p ≥ 0.05  (lineare OK)"),
      name   = NULL
    ) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(title = target_labels[tgt], x = "p-value termine quadratico", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title     = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

wrap_plots(map(target_vars, make_ftest_plot), nrow = 1) +
  plot_annotation(
    title    = "Test F: modello lineare vs quadratico per field",
    subtitle = "Rosso = il termine quadratico migliora significativamente il fit",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )
