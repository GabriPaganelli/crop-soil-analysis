library(tidyverse)
library(GGally)
library(patchwork)
library(ggpubr)
library(ggtern)
library(soiltexture)
library(corrplot)
library(here)

crop <- readRDS(here("data", "crop.rds"))
dati <- readRDS(here("data", "dati.rds"))
source(here("scripts", "utilities.R"))

target_vars <- c("PercSOC", "PercTotNitro", "PercTotPhos")

# Formato long delle target — usato in più sezioni del file
crop_long <- crop |>
  select(Field, Bottom, Lat, Long, all_of(target_vars), Landuse) |>
  pivot_longer(cols      = all_of(target_vars),
               names_to  = "target",
               values_to = "valore")


# ── ANALISI UNIVARIATA ────────────────────────────────────────────────────────
analizza_dataset(crop)
grafico_distribuzioni(crop[, !names(crop) %in% c("Field", "Plot", "Lat", "Long")])


# ── TESSITURA ─────────────────────────────────────────────────────────────────

# Triangolo USDA (versione ggtern e versione TT.plot)
plot_texture_triangle(crop, color_var = "Landuse")
plot_texture_triangle(crop, color_var = "Landuse", version = "ttplot")

# Le coordinate ILR (Texture1, Texture2) sostituiscono le 3 percentuali in dati.
# Qui le visualizziamo nel contesto delle classi USDA per dare un'interpretazione
# geometrica alle due componenti.
ggplot(dati, aes(x = Texture1, y = Texture2, colour = crop$Texture)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title   = "Coordinate ILR della tessitura per classe USDA",
       x       = "Texture1  (finezza: argilla+limo vs sabbia)",
       y       = "Texture2  (argilla vs limo nella frazione fine)",
       colour  = "Classe USDA") +
  theme_minimal()

# Correlazione di Texture1 e Texture2 con le variabili target
plot_x_vs_y(dati[, c("Texture1", "Texture2", target_vars)], target_vars)


# ── TARGET VS PREDITTORI ──────────────────────────────────────────────────────
plot_x_vs_y(crop, target_vars)


# ── PROFILI VERTICALI ─────────────────────────────────────────────────────────
# crop_long è usato anche nelle sezioni successive
crop_long <- crop |>
  select(Field, Bottom, Lat, Long, all_of(target_vars), Landuse) |>
  pivot_longer(cols      = all_of(target_vars),
               names_to  = "target",
               values_to = "valore")

# Per Field (40 livelli): mostra la variabilità individuale
plot_depth_profiles(crop, target_vars, color_var = "Field", ncol = 1)

# Per Landuse (7 livelli): mostra la separazione per classe gestionale
plot_depth_profiles(crop, target_vars, color_var = "Landuse", ncol = 1)


# ── STRUTTURA GERARCHICA ──────────────────────────────────────────────────────
# Varianza tra Field vs varianza dentro Field (guidata dalla profondità)

ggplot(crop_long,
       aes(x = Landuse, y = valore, fill = as.factor(Bottom))) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_wrap(~ target, scales = "free_y") +
  labs(title = "Distribuzione per Landuse e profondità",
       x = "Landuse", y = "Valore", fill = "Profondità (cm)") +
  theme_minimal()

# Rapporto sd_between / sd_within per quantificare la gerarchia
# PercSOC: rapporto < 1 → variabilità guidata dalla profondità più che dal field
# PercTotNitro, PercTotPhos: rapporto ≈ 1 → profondità e field contano in modo simile
crop_long |>
  group_by(target, Field) |>
  summarise(media = mean(valore), sd_within = sd(valore), .groups = "drop") |>
  group_by(target) |>
  summarise(
    sd_between      = sd(media),
    sd_within_media = mean(sd_within),
    rapporto        = sd_between / sd_within_media
  )


# ── CORRELAZIONI TRA LE VARIABILI TARGET ──────────────────────────────────────

ggpairs(
  crop,
  columns = target_vars,
  aes(colour = Landuse, alpha = 0.6),
  upper = list(continuous = wrap("cor", size = 3)),
  lower = list(continuous = wrap("points", size = 0.8)),
  diag  = list(continuous = wrap("densityDiag", alpha = 0.5))
) +
  labs(title = "Correlazioni tra variabili target — per Landuse") +
  theme_minimal()

ggpairs(
  crop,
  columns = target_vars,
  aes(colour = Fertilised, alpha = 0.6),
  upper = list(continuous = wrap("cor", size = 3)),
  lower = list(continuous = wrap("points", size = 0.8)),
  diag  = list(continuous = wrap("densityDiag", alpha = 0.5))
) +
  labs(title = "Correlazioni tra variabili target — per Fertilised") +
  theme_minimal()

ggpairs(
  crop,
  columns = target_vars,
  aes(colour = Texture, alpha = 0.6),
  upper = list(continuous = wrap("cor", size = 3)),
  lower = list(continuous = wrap("points", size = 0.8)),
  diag  = list(continuous = wrap("densityDiag", alpha = 0.5))
) +
  labs(title = "Correlazioni tra variabili target — per Texture") +
  theme_minimal()

# Correlazioni stratificate per Landuse e profondità
# Verifica se la struttura di correlazione tra SOC, N e P cambia
# a seconda della classe di gestione o della profondità di campionamento.

cor_per_landuse <- crop |>
  group_split(Landuse) |>
  setNames(levels(crop$Landuse)) |>
  lapply(\(d) cor(d[, target_vars], use = "complete.obs"))

op <- par(mfrow = c(2, 4), mar = c(0, 0, 2, 0))
for (nome in names(cor_per_landuse)) {
  corrplot(cor_per_landuse[[nome]],
           method      = "color",
           type        = "upper",
           addCoef.col = "black",
           number.cex  = 0.9,
           tl.col      = "black",
           tl.cex      = 0.8,
           cl.lim      = c(-1, 1),
           col         = colorRampPalette(c("tomato", "white", "steelblue"))(100),
           title       = nome,
           mar         = c(0, 0, 1.5, 0))
}
par(op)

cor_per_depth <- crop |>
  group_split(Bottom) |>
  setNames(paste0(sort(unique(crop$Bottom)), " cm")) |>
  lapply(\(d) cor(d[, target_vars], use = "complete.obs"))

op <- par(mfrow = c(2, 3), mar = c(0, 0, 2, 0))
for (nome in names(cor_per_depth)) {
  corrplot(cor_per_depth[[nome]],
           method      = "color",
           type        = "upper",
           addCoef.col = "black",
           number.cex  = 0.9,
           tl.col      = "black",
           tl.cex      = 0.8,
           cl.lim      = c(-1, 1),
           col         = colorRampPalette(c("tomato", "white", "steelblue"))(100),
           title       = nome,
           mar         = c(0, 0, 1.5, 0))
}
par(op)

# Tabelle numeriche e grafici di variazione
tab_cor_landuse <- crop |>
  group_by(Landuse) |>
  summarise(
    `SOC-N` = cor(PercSOC, PercTotNitro, use = "complete.obs"),
    `SOC-P` = cor(PercSOC, PercTotPhos,  use = "complete.obs"),
    `N-P`   = cor(PercTotNitro, PercTotPhos, use = "complete.obs"),
    n       = n(),
    .groups = "drop"
  )
print(tab_cor_landuse)

tab_cor_depth <- crop |>
  group_by(Bottom) |>
  summarise(
    `SOC-N` = cor(PercSOC, PercTotNitro, use = "complete.obs"),
    `SOC-P` = cor(PercSOC, PercTotPhos,  use = "complete.obs"),
    `N-P`   = cor(PercTotNitro, PercTotPhos, use = "complete.obs"),
    n       = n(),
    .groups = "drop"
  )
print(tab_cor_depth)

tab_cor_depth |>
  pivot_longer(cols      = c(`SOC-N`, `SOC-P`, `N-P`),
               names_to  = "coppia",
               values_to = "correlazione") |>
  ggplot(aes(x = Bottom, y = correlazione, color = coppia, group = coppia)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(breaks = sort(unique(crop$Bottom))) +
  labs(title = "Correlazioni tra SOC, N e P al variare della profondità",
       x     = "Profondità (cm)",
       y     = "Correlazione di Pearson",
       color = "Coppia") +
  theme_minimal()

tab_cor_landuse |>
  pivot_longer(cols      = c(`SOC-N`, `SOC-P`, `N-P`),
               names_to  = "coppia",
               values_to = "correlazione") |>
  ggplot(aes(x = Landuse, y = correlazione, fill = coppia)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Correlazioni tra SOC, N e P per classe di gestione",
       x     = "Landuse",
       y     = "Correlazione di Pearson",
       fill  = "Coppia") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


# ── EFFETTO SPAZIALE ──────────────────────────────────────────────────────────

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

# Correlazione numerica con le coordinate geografiche
crop |>
  summarise(across(all_of(target_vars), list(
    cor_Lat  = ~ cor(.x, Lat,  use = "complete.obs"),
    cor_Long = ~ cor(.x, Long, use = "complete.obs")
  ))) |>
  pivot_longer(everything(),
               names_to  = c("variabile", "coord"),
               names_sep = "_cor_",
               values_to = "correlazione")

# I field 1-5 appartengono alla stessa farm?
crop |>
  mutate(
    Field_Num   = as.numeric(gsub("\\D", "", as.character(Field))),
    Field_Group = as.factor((Field_Num - 1) %/% 5)
  ) |>
  ggplot(aes(x = Long, y = Lat, color = Field_Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_viridis_d(option = "turbo") +
  labs(title    = "Distribuzione Lat/Long raggruppata per Field",
       subtitle = "I colori cambiano ogni 5 campi (es. 1-5, 6-10...)",
       x        = "Longitudine",
       y        = "Latitudine",
       color    = "Blocco Field") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())


# ── LINEARITÀ DI LOG(TARGET) VS BOTTOM ────────────────────────────────────────
# Per ogni field × target: R² del modello lineare e test F lineare vs quadratico.
# R² alto + termine quadratico non significativo (p ≥ 0.05) → relazione log-lineare
# con la profondità supportata.

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

target_labels <- c(
  PercSOC      = "SOC (%)",
  PercTotNitro = "Azoto totale (%)",
  PercTotPhos  = "Fosforo totale (%)"
)

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
      plot.title      = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

wrap_plots(map(target_vars, make_ftest_plot), nrow = 1) +
  plot_annotation(
    title    = "Test F: modello lineare vs quadratico per field",
    subtitle = "Rosso = il termine quadratico migliora significativamente il fit",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )
