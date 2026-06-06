# ============================================================
#  EDA — Grafici motivazione modello
#
#  Output (EDA/):
#    eda_mod_04_mappa.pdf       mappa spaziale campi (colorata per SOC medio)
#    eda_mod_05_blup_cross.pdf  correlazione BLUP intercette (LMM univariati)
# ============================================================

library(tidyverse)
library(patchwork)
library(lme4)

dati <- readRDS(here::here("data", "crop.rds")) |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom),
    MacroGest = factor(case_when(
      Landuse == 4          ~ "Natural",
      Landuse %in% c(1, 2) ~ "Farm",
      TRUE                  ~ "Off-farm"
    ), levels = c("Natural", "Farm", "Off-farm"))
  )

pal_macro <- c("Natural" = "#4E9B6F", "Farm" = "#C97D30", "Off-farm" = "#5B87C0")

th <- theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "#FAFAF8", color = NA),
    panel.background = element_rect(fill = "#F4F2ED", color = NA),
    panel.grid.major = element_line(color = "#E5E2DB", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text        = element_text(size = 9, color = "#444444"),
    axis.title       = element_text(size = 9.5, color = "#333333"),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.margin      = margin(10, 14, 8, 10)
  )

save_fig <- function(fname, p, w = 28, h = 18, u = "cm") {
  path <- here::here("EDA", fname)
  ggsave(path, plot = p, width = w, height = h, units = u,
         device = cairo_pdf, bg = "#FAFAF8")
  message("✓ ", path)
}


# ══════════════════════════════════════════════════════════════
#  FIGURA 4 — Mappa spaziale dei campi
# ══════════════════════════════════════════════════════════════

field_map <- dati |>
  group_by(Field, Lat, Long, MacroGest) |>
  summarise(SOC_medio = mean(logSOC, na.rm = TRUE), .groups = "drop")

fig4 <- ggplot(field_map, aes(x = Long, y = Lat)) +
  geom_point(aes(fill = SOC_medio, shape = MacroGest),
             size = 4.0, color = "#444444", stroke = 0.5) +
  scale_fill_viridis_c(option = "viridis", direction = 1,
                       name = "log(SOC %)\nmedio") +
  scale_shape_manual(values = c("Natural" = 21, "Farm" = 22, "Off-farm" = 24),
                     name = NULL) +
  labs(title = "Distribuzione spaziale dei 40 campi",
       x = "Longitudine", y = "Latitudine") +
  th +
  theme(legend.position = "right")

save_fig("eda_mod_04_mappa.pdf", fig4, w = 28, h = 18)


# ══════════════════════════════════════════════════════════════
#  FIGURA 5 — Correlazione intercette casuali (ranef da LMM)
#  Modello completo: stesse covariate del modello finale MVRE
#  (within: logBottom, Texture1, Texture2, BulkDensity, PH)
#  (between: OnFarm, Irrigate, Fertilised, N_Natural)
# ══════════════════════════════════════════════════════════════

# dati.rds ha già Texture1, Texture2, log-risposte
dati_lmm <- readRDS(here::here("data", "dati.rds")) |>
  mutate(
    logSOC       = log(PercSOC),
    logN         = log(PercTotNitro),
    logP         = log(PercTotPhos),
    logBottom_c  = log(Bottom) - mean(log(Bottom)),
    BulkDensity_s = scale(BulkDensity)[,1],
    PH_s          = scale(PH)[,1],
    MacroGest     = factor(case_when(
      N_Natural == 0 ~ "Natural",
      OnFarm == 1    ~ "Farm",
      TRUE           ~ "Off-farm"
    ), levels = c("Natural", "Farm", "Off-farm"))
  )

form <- . ~ logBottom_c + Texture1 + Texture2 + BulkDensity_s + PH_s +
            OnFarm + Irrigate + Fertilised + N_Natural + (1 | Field)

lmm_SOC <- lmer(update(form, logSOC ~ .), data = dati_lmm, REML = TRUE)
lmm_N   <- lmer(update(form, logN   ~ .), data = dati_lmm, REML = TRUE)
lmm_P   <- lmer(update(form, logP   ~ .), data = dati_lmm, REML = TRUE)

blup_wide <- tibble(
  Field = rownames(ranef(lmm_SOC)$Field),
  SOC   = ranef(lmm_SOC)$Field[[1]],
  N     = ranef(lmm_N)$Field[[1]],
  P     = ranef(lmm_P)$Field[[1]]
) |>
  left_join(dati_lmm |> distinct(Field, MacroGest) |> mutate(Field = as.character(Field)),
            by = "Field")

make_blup_scatter <- function(xvar, yvar) {
  r <- cor(blup_wide[[xvar]], blup_wide[[yvar]])
  ggplot(blup_wide, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_hline(yintercept = 0, color = "#CCCCCC", linewidth = 0.4, linetype = "dashed") +
    geom_vline(xintercept = 0, color = "#CCCCCC", linewidth = 0.4, linetype = "dashed") +
    geom_point(size = 2.6, alpha = 0.8, color = "#555555") +
    geom_smooth(method = "lm", se = TRUE,
                color = "#333333", fill = "#CCCCCC",
                linewidth = 0.7, alpha = 0.2) +
    annotate("text", x = -Inf, y = Inf,
             label = paste0("r = ", round(r, 2)),
             hjust = -0.2, vjust = 1.6,
             size = 3.8, color = "#333333", fontface = "italic") +
    labs(x = paste0("Intercetta casuale — ", xvar),
         y = paste0("Intercetta casuale — ", yvar),
         title = paste0(xvar, " vs ", yvar)) +
    th +
    theme(plot.title = element_text(size = 11))
}

p_b_soc_n <- make_blup_scatter("SOC", "N")
p_b_soc_p <- make_blup_scatter("SOC", "P")
p_b_n_p   <- make_blup_scatter("N",   "P")

fig5 <- (p_b_soc_n + p_b_soc_p + p_b_n_p) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "#FAFAF8", color = NA))
  )

save_fig("eda_mod_05_blup_cross.pdf", fig5, w = 28, h = 14)

message("\n✓ Figure motivazione modello generate in EDA/:")
message("  eda_mod_04_mappa.pdf")
message("  eda_mod_05_blup_cross.pdf")
