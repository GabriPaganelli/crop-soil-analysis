# =============================================================================
# 03_eda_profili_np.R  —  EDA: analisi approfondita profili N e P
#
# Motivazione: dallo script 04 emergono due pattern:
#   - N: pochi campi erratici sembrano mascherare un segnale nei campi "ricchi"
#   - P: erraticità più diffusa
#
# Questo script:
#   1. Stratifica per livello (tercili) e calcola Corr dentro ogni strato
#   2. Identifica i campi "influenti" (outlier nella scatter intercetta-slope)
#   3. Verifica se i campi erratici condividono caratteristiche (management)
#   4. Confronta profili N e P stratificati per livello
#
# Dipende da: data/dati.rds, scripts/00_utilities.R
# =============================================================================


library(tidyverse)
library(here)
source(here("scripts", "00_utilities.R"))


# ── 1. DATI ───────────────────────────────────────────────────────────────────

dati_raw <- readRDS(here("data", "dati.rds")) |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom)
  ) |>
  mutate(Field = factor(Field))

mu_logB <- mean(log(dati_raw$Bottom))
dati <- dati_raw |> mutate(logBottom_c = logBottom - mu_logB)

# Caratteristiche di campo (una riga per campo)
field_chars <- dati |>
  distinct(Field, OnFarm, Irrigate, Fertilised, N_Natural) |>
  mutate(Field = as.character(Field))


# ── 2. OLS PER CAMPO ─────────────────────────────────────────────────────────
# fit_field() definita in 00_utilities.R

ols_N <- fit_field(dati, "logN")   |> mutate(Field = as.character(Field)) |> left_join(field_chars, by = "Field")
ols_P   <- fit_field(dati, "logP")   |> mutate(Field = as.character(Field)) |> left_join(field_chars, by = "Field")
ols_SOC <- fit_field(dati, "logSOC") |> mutate(Field = as.character(Field)) |> left_join(field_chars, by = "Field")


# ── 3. ANALISI PER STRATO: correlazione nei tercili di intercetta ──────────────

analizza_strati <- function(ols_df, risposta) {
  ols_s <- ols_df |>
    mutate(strato = ntile(int, 3),
           strato_lbl = case_when(
             strato == 1 ~ "Basso (T1)",
             strato == 2 ~ "Medio (T2)",
             strato == 3 ~ "Alto  (T3)"
           ))

  cat(sprintf("\n═══ %s — Corr(int, slope) per strato di livello ═══\n", risposta))
  cat(sprintf("%-12s | %8s %8s | %8s\n", "Strato", "Pearson", "Spearman", "n"))
  cat(strrep("-", 42), "\n")

  ols_s |>
    group_by(strato_lbl) |>
    summarise(
      r_p = cor(int, slope, method = "pearson"),
      r_s = cor(int, slope, method = "spearman"),
      n   = n(),
      .groups = "drop"
    ) |>
    arrange(strato_lbl) |>
    pwalk(function(strato_lbl, r_p, r_s, n) {
      cat(sprintf("%-12s | %8.3f %8.3f | %8d\n", strato_lbl, r_p, r_s, n))
    })

  # Correlazione totale per confronto
  cat(sprintf("%-12s | %8.3f %8.3f | %8d\n", "TOTALE",
              cor(ols_s$int, ols_s$slope, method = "pearson"),
              cor(ols_s$int, ols_s$slope, method = "spearman"),
              nrow(ols_s)))

  ols_s
}

ols_N_s   <- analizza_strati(ols_N,   "logN")
ols_P_s   <- analizza_strati(ols_P,   "logP")
ols_SOC_s <- analizza_strati(ols_SOC, "logSOC")


# ── 4. IDENTIFICAZIONE CAMPI INFLUENTI (N e P) ────────────────────────────────

identifica_outlier <- function(ols_df, risposta) {
  lm_fit <- lm(slope ~ int, data = ols_df)
  ols_aug <- ols_df |>
    mutate(
      fitted_slope  = fitted(lm_fit),
      resid         = residuals(lm_fit),
      std_resid     = rstudent(lm_fit),
      cooksd        = cooks.distance(lm_fit),
      leverage      = hatvalues(lm_fit),
      is_outlier    = abs(std_resid) > 2.0,
      is_influential = cooksd > (4 / nrow(ols_df))
    )

  cat(sprintf("\n═══ %s — Campi influenti o outlier ═══════════════════\n", risposta))
  flagged <- ols_aug |>
    filter(is_outlier | is_influential) |>
    arrange(desc(abs(std_resid))) |>
    select(Field, int, slope, std_resid, cooksd,
           OnFarm, Irrigate, Fertilised, N_Natural)

  if (nrow(flagged) == 0) {
    cat("  Nessun campo con |std_resid| > 2 o Cook's D > 4/n\n")
  } else {
    print(flagged, n = 20)
  }

  ols_aug
}

ols_N_aug   <- identifica_outlier(ols_N,   "logN")
ols_P_aug   <- identifica_outlier(ols_P,   "logP")


# ── 5. SCATTER ANNOTATO: intercetta vs slope con campi etichettati ─────────────

scatter_annotato <- function(ols_aug, risposta, colore = "steelblue") {

  # Campi da etichettare: outlier/influenti
  da_etichettare <- ols_aug |> filter(is_outlier | is_influential)

  # Aggiungi strato come shape
  ols_aug2 <- ols_aug |>
    mutate(strato_lbl = case_when(
      ntile(int, 3) == 1 ~ "Basso",
      ntile(int, 3) == 2 ~ "Medio",
      ntile(int, 3) == 3 ~ "Alto"
    ))

  r_totale  <- cor(ols_aug$int, ols_aug$slope, method = "pearson")
  r_alto    <- with(ols_aug2 |> filter(strato_lbl == "Alto"),
                    cor(int, slope, method = "pearson"))
  r_basso   <- with(ols_aug2 |> filter(strato_lbl == "Basso"),
                    cor(int, slope, method = "pearson"))

  p <- ggplot(ols_aug2, aes(int, slope)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(aes(color = strato_lbl, size = cooksd * 100),
               alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE,
                color = "grey40", fill = "grey80", alpha = 0.2,
                linetype = "dashed") +
    # Linea di regressione solo per strato Alto
    geom_smooth(data = ols_aug2 |> filter(strato_lbl == "Alto"),
                method = "lm", se = FALSE,
                color = "tomato", linewidth = 1) +
    ggrepel::geom_text_repel(
      data = da_etichettare,
      aes(label = Field),
      size = 3, color = "grey30", max.overlaps = 20
    ) +
    scale_color_manual(values = c("Basso" = "#440154", "Medio" = "#21908C",
                                  "Alto"  = "#FDE725"),
                       name = "Livello campo") +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r tot = %+.3f\nr alto = %+.3f\nr basso = %+.3f",
                             r_totale, r_alto, r_basso),
             hjust = 1.1, vjust = 1.5, size = 3.5, color = "grey20") +
    labs(
      title    = sprintf("log%s — scatter intercetta vs slope (annotato)", risposta),
      subtitle = "Dimensione punto ∝ Cook's D | Etichette = campi influenti",
      x        = sprintf("Intercetta per campo (log%s a profondità media)", risposta),
      y        = "Slope (decadimento con log-profondità)"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  print(p)
  invisible(p)
}

# Installa ggrepel se necessario
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")

scatter_annotato(ols_N_aug,   "N",   "forestgreen")
scatter_annotato(ols_P_aug,   "P",   "tomato")
scatter_annotato(ols_SOC_aug <- identifica_outlier(ols_SOC, "logSOC"), "SOC")


# ── 6. PROFILI VERTICALI STRATIFICATI: N e P ──────────────────────────────────
# Spaghetti con colorazione per strato, per vedere se il pattern emerge visivamente

spaghetti_stratificato <- function(ols_s_df, dati_df, y_var, titolo) {

  level_map <- ols_s_df |> select(Field, strato_lbl)
  dati_r <- dati_df |>
    select(Field, logBottom, y = all_of(y_var)) |>
    mutate(Field = as.character(Field)) |>
    left_join(level_map, by = "Field") |>
    mutate(strato_lbl = factor(strato_lbl,
                               levels = c("Basso (T1)", "Medio (T2)", "Alto  (T3)")))

  ggplot(dati_r, aes(x = logBottom, y = y,
                     group = Field, color = strato_lbl)) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    geom_point(size = 1, alpha = 0.5) +
    # Media per strato
    stat_summary(aes(group = strato_lbl), fun = mean,
                 geom = "line", linewidth = 2, linetype = "dashed") +
    scale_color_manual(values = c("Basso (T1)" = "#440154",
                                  "Medio (T2)" = "#21908C",
                                  "Alto  (T3)" = "#FDE725"),
                       name = "Strato") +
    facet_wrap(~ strato_lbl, ncol = 3) +
    labs(
      title    = titolo,
      subtitle = "Tratteggio = media di strato | Se slope-intercept correlati, i panel T3 sono più piatti",
      x        = "log(profondità in cm)",
      y        = paste0("log(", y_var, ")")
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold"))
}

print(spaghetti_stratificato(ols_N_s, dati, "logN",
  "Profili logN per strato di livello — pattern intercetta-slope?"))
print(spaghetti_stratificato(ols_P_s, dati, "logP",
  "Profili logP per strato di livello — erraticità confronto"))


# ── 7. CARATTERISTICHE DEI CAMPI ERRATICI ────────────────────────────────────
# Per N: i campi outlier nella scatter differiscono nel management?

cat("\n═══ MANAGEMENT DEI CAMPI INFLUENTI (N) ═════════════════════════\n")
cat("Campi influenti vs resto:\n")

ols_N_aug |>
  mutate(gruppo = ifelse(is_outlier | is_influential, "Influente", "Normale")) |>
  group_by(gruppo) |>
  summarise(
    n          = n(),
    OnFarm_pct    = mean(OnFarm)    * 100,
    Irrigate_pct  = mean(Irrigate)  * 100,
    Fertilised_pct= mean(Fertilised)* 100,
    N_Natural_pct = mean(N_Natural) * 100,
    int_medio     = mean(int),
    slope_medio   = mean(slope),
    .groups = "drop"
  ) |>
  print()

cat("\n═══ MANAGEMENT DEI CAMPI INFLUENTI (P) ═════════════════════════\n")
ols_P_aug |>
  mutate(gruppo = ifelse(is_outlier | is_influential, "Influente", "Normale")) |>
  group_by(gruppo) |>
  summarise(
    n             = n(),
    OnFarm_pct    = mean(OnFarm)    * 100,
    Irrigate_pct  = mean(Irrigate)  * 100,
    Fertilised_pct= mean(Fertilised)* 100,
    N_Natural_pct = mean(N_Natural) * 100,
    int_medio     = mean(int),
    slope_medio   = mean(slope),
    .groups = "drop"
  ) |>
  print()


# ── 8. SLOPE NOISE: coefficiente di variazione della slope ────────────────────
# Quanto è rumorosa la stima OLS della slope per campo?
# Se slope_se è grande rispetto a slope, la stima è inaffidabile.

cat("\n═══ QUALITÀ STIMA SLOPE (slope / slope_se) ══════════════════════\n")
for (nm in c("SOC", "N", "P")) {
  ols_df <- get(paste0("ols_", nm))
  snr <- abs(ols_df$slope) / ols_df$slope_se
  cat(sprintf("log%s: SNR mediana = %.2f | % SNR < 1 = %.0f%%\n",
              nm,
              median(snr, na.rm = TRUE),
              mean(snr < 1, na.rm = TRUE) * 100))
}

cat("\n  SNR < 1 = slope stimata è più piccola del suo errore standard\n")
cat("  (normale con pochi punti per campo, ma segnala alta incertezza)\n")

cat("\n── Fine script 21b ──────────────────────────────────────────────\n")
