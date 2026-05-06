library(tidyverse)
library(GGally)
library(patchwork)
library(here)

# ── 1. CARICAMENTO ────────────────────────────────────────────────────────────
crop = readRDS(here("data", "crop.rds"))
source(here('scripts', 'utilities.R'))

target_vars <- c("PercSOC", "PercTotNitro", "PercTotPhos")

# ── 3. ANALISI ESPLORATIVA UNIVARIATA ─────────────────────────────────────────
analizza_dataset(crop)
grafico_distribuzioni(crop[,-c(2,3,18,19)])

# Plotta i tuoi crop direttamente sul triangolo USDA
soiltexture::TT.plot(class.sys = "USDA.TT", tri.data = data.frame(
  SAND = crop$PercSand,
  SILT = crop$PercSilt,
  CLAY = crop$PercClay
))

# Correlazioni e scatter di ogni target vs tutti i predittori
plot_x_vs_y(crop, target_vars)


# ── 4. PROFILI VERTICALI (depth profiles) ─────────────────────────────────────
# Ogni osservazione è un campionamento a una certa profondità (Bottom)
# dello stesso plot: ha senso visualizzare le traiettorie per plot

# Porta i crop in formato long
crop_long <- crop |>
  select(Field, Bottom, all_of(target_vars)) |>
  pivot_longer(cols = all_of(target_vars),
               names_to = "target",
               values_to = "valore")

# Grafico
ggplot(crop_long, aes(x = Bottom, y = valore,
                      group = Field, color = Field)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ target, scales = "free_y", ncol = 1,
             labeller = labeller(target = c(
               PercSOC      = "SOC (%)",
               PercTotNitro = "Azoto totale (%)",
               PercTotPhos  = "Fosforo totale (%)"
             ))) +
  scale_x_continuous(breaks = c(20, 30, 40, 50, 60, 80)) +
  scale_color_viridis_d(option = "turbo") +
  labs(title = "Andamento delle target vars per profondità",
       x = "Bottom (profondità, cm)",
       y = "Valore (%)",
       color = "Field") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold"))


# ── 5. STRUTTURA GERARCHICA ───────────────────────────────────────────────────
# Varianza tra Field vs dentro Field (guidata da profondità)

# Boxplot per Landuse (7 categorie, leggibile) x profondità
ggplot(crop_long_target,
       aes(x = Landuse, y = valore, fill = as.factor(Bottom))) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_wrap(~ variabile, scales = "free_y") +
  labs(title = "Distribuzione per Landuse e profondità",
       x = "Landuse", y = "Valore", fill = "Profondità (cm)") +
  theme_minimal()

# Rapporto sd_between / sd_within per quantificare la gerarchia
crop_long_target |>
  group_by(variabile, Field) |>
  summarise(media = mean(valore), sd_within = sd(valore), .groups = "drop") |>
  group_by(variabile) |>
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

# ── 7. EFFETTO SPAZIALE ───────────────────────────────────────────────────────

# Un plot per variabile (scala colore separata, altrimenti dominata da PercSOC)
plots_spaziali <- map(target_vars, ~ {
  crop_long_target |>
    filter(variabile == .x) |>
    ggplot(aes(x = Lat, y = Long, size = valore, colour = valore)) +
    geom_point(alpha = 0.6) +
    scale_size_continuous(range = c(1, 8)) +
    scale_colour_viridis_c(option = "magma", direction = -1) +
    labs(title = .x, x = "Latitudine", y = "Longitudine",
         size = "Valore", colour = "Valore") +
    theme_minimal()
})

wrap_plots(plots_spaziali, ncol = 1)

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