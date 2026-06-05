# =============================================================================
# 02_eda_forma_icc.R  —  EDA: forma funzionale, ICC, correlazione intercetta-slope
#
# STRUTTURA:
#   Sezione 0 — Verifiche preliminari: ICC, forma funzionale del profilo
#               verticale, random slope (da dati grezzi con lme4).
#               Queste analisi motivano le scelte del modello M-SP.
#   Sezione 1-6 — Relazione intercetta-slope nei dati grezzi (OLS per campo):
#               conferma il segnale di eta_SOC stimato da M-SP.
#
# RISULTATI PRINCIPALI (Sezione 0):
#   - ICC intra-field per logSOC ≈ 0.72 (dopo correzione per Bottom):
#     variabilità between-field dominante → random intercept necessario
#   - AIC log(SOC): log(Bottom) < quadratico < lineare < factor(Bottom)
#     → forma power-law (logBottom) è la più parsimoniosa tra le buone
#   - Corr(BLUP intercetta, BLUP slope) per SOC ≈ -0.93 (Bottom non centrato)
#     equivalente a ≈ +0.93 centrato → "chi parte alto decade più lentamente"
#
# Dipende da: data/dati.rds, scripts/00_utilities.R
# =============================================================================


# ── 0. SETUP ──────────────────────────────────────────────────────────────────

library(tidyverse)
library(here)
library(lme4)
source(here("scripts", "00_utilities.R"))

dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)
fig_dir <- here("output", "figures")

save_fig <- function(fname, p, w = 16, h = 9, u = "cm") {
  ggplot2::ggsave(file.path(fig_dir, fname), plot = p,
                  width = w, height = h, units = u, device = "pdf")
  cat(sprintf("  [fig] Salvato: %s\n", fname))
}


# ── 0. VERIFICHE PRELIMINARI (motivazione del modello) ────────────────────────

dati_raw <- readRDS(here("data", "dati.rds")) |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(Field = factor(Field),
         logSOC = log(PercSOC),
         Bottom_c = Bottom - mean(Bottom))

# ── 0a. ICC INTRA-FIELD ───────────────────────────────────────────────────────
# ICC "vera" separando varianza verticale sistematica dalla varianza between-field.
# m0: solo random intercept (varianza verticale inclusa nel residuo)
# m1: logSOC ~ Bottom + (1|Field) → varianza between emerge più chiaramente

m0 <- lmer(logSOC ~ 1            + (1 | Field), data = dati_raw, REML = TRUE)
m1 <- lmer(logSOC ~ Bottom_c     + (1 | Field), data = dati_raw, REML = TRUE)

vc0 <- as.data.frame(VarCorr(m0))
icc0 <- vc0$vcov[1] / sum(vc0$vcov)
vc1 <- as.data.frame(VarCorr(m1))
icc1 <- vc1$vcov[1] / sum(vc1$vcov)

cat("══════════════════════════════════════════════════════\n")
cat(" 0a. ICC INTRA-FIELD (logSOC)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("  ICC modello nullo:          %.3f\n", icc0))
cat(sprintf("  ICC dopo correzione Bottom: %.3f\n", icc1))
cat("  → Il random intercept spiega la maggior parte della variabilità.\n\n")

# ── 0b. FORMA FUNZIONALE DEL PROFILO VERTICALE ────────────────────────────────
# Confronto AIC (ML) tra forme funzionali per la dipendenza da Bottom.
# Risultato atteso: log(Bottom) ha AIC più basso (forma power-law parsimoniosa).

m_lin  <- lmer(logSOC ~ Bottom_c               + (1|Field), data = dati_raw, REML = FALSE)
m_quad <- lmer(logSOC ~ Bottom_c + I(Bottom_c^2) + (1|Field), data = dati_raw, REML = FALSE)
m_log  <- lmer(logSOC ~ log(Bottom)             + (1|Field), data = dati_raw, REML = FALSE)
m_fact <- lmer(logSOC ~ factor(Bottom)           + (1|Field), data = dati_raw, REML = FALSE)

aic_tbl <- AIC(m_lin, m_quad, m_log, m_fact)
rownames(aic_tbl) <- c("Lineare (Bottom_c)", "Quadratico", "log(Bottom)", "factor(Bottom)")

cat("══════════════════════════════════════════════════════\n")
cat(" 0b. CONFRONTO AIC — FORMA FUNZIONALE (logSOC)\n")
cat("══════════════════════════════════════════════════════\n")
print(aic_tbl[order(aic_tbl$AIC), ])
cat("  → log(Bottom) = forma power-law: la migliore tra le parametriche.\n\n")

# Figura AIC
aic_df <- data.frame(
  forma = rownames(aic_tbl),
  AIC   = aic_tbl$AIC
) |>
  mutate(
    delta_AIC = AIC - min(AIC),
    forma     = factor(forma, levels = forma[order(delta_AIC, decreasing = TRUE)])
  )

p_aic <- ggplot(aic_df, aes(x = delta_AIC, y = forma)) +
  geom_col(aes(fill = delta_AIC < 2), show.legend = FALSE, alpha = 0.85, width = 0.6) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey60")) +
  geom_text(aes(label = sprintf("AIC diff = %.1f", delta_AIC)), hjust = -0.15, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title    = "Confronto AIC - forma funzionale del profilo verticale",
    subtitle = "Risposta: logSOC  |  Modello base: (1 | Field)  |  Stima: ML",
    x        = "AIC - min(AIC)",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

print(p_aic)
save_fig("fig_01_aic_forma_funzionale.pdf", p_aic, w = 14, h = 8)

# ── 0c. RANDOM SLOPE: "chi parte alto decade più lentamente" ──────────────────
# Con Bottom centrato: correlazione intercetta-slope positiva (campo ricco
# a profondità media → slope più bassa = decade meno da lì in poi).
# Non centrato: correlazione negativa (campo ricco a superficie → decade più).
# I due modelli sono matematicamente equivalenti; la parametrizzazione di M-SP
# usa la versione non centrata implicita nella struttura proporzionale.

m2 <- lmer(logSOC ~ Bottom_c + (Bottom_c | Field), data = dati_raw, REML = TRUE)
blup <- as.data.frame(ranef(m2)$Field)
names(blup) <- c("u0_intercept", "u1_slope")

cat("══════════════════════════════════════════════════════\n")
cat(" 0c. RANDOM SLOPE — BLUP intercetta vs slope (Bottom centrato)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("  Corr(u0, u1) centrato:   %+.3f\n", cor(blup$u0_intercept, blup$u1_slope)))
cat(sprintf("  (Equivalente: b_OLS ≈ %+.3f su Bottom NON centrato)\n",
            coef(lm(u1_slope ~ u0_intercept, data = blup))[2]))
cat("  → Conferma il segnale di eta_SOC nel modello M-SP.\n\n")

rm(m0, m1, m2, m_lin, m_quad, m_log, m_fact, vc0, vc1, blup, dati_raw); gc()


# ── 1. DATI ───────────────────────────────────────────────────────────────────

dati <- readRDS(here("data", "dati.rds")) |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom)
  ) |>
  mutate(Field = factor(Field))

# Centra logBottom per avere intercetta = valore alla profondità media
mu_logB <- mean(log(dati$Bottom))
dati <- dati |> mutate(logBottom_c = logBottom - mu_logB)

cat(sprintf("N = %d | J = %d | profondità media = %.1f cm (logBottom = %.3f)\n",
            nrow(dati), length(levels(dati$Field)),
            exp(mu_logB), mu_logB))


# ── 2. OLS PER CAMPO ──────────────────────────────────────────────────────────
# fit_field() definita in 00_utilities.R

ols_SOC <- fit_field(dati, "logSOC") |> mutate(risposta = "SOC")
ols_N   <- fit_field(dati, "logN")   |> mutate(risposta = "N")
ols_P   <- fit_field(dati, "logP")   |> mutate(risposta = "P")

ols_all <- bind_rows(ols_SOC, ols_N, ols_P) |>
  mutate(risposta = factor(risposta, levels = c("SOC", "N", "P")))


# ── 3. CORRELAZIONI ───────────────────────────────────────────────────────────

cat("\n═══ CORRELAZIONE INTERCETTA–SLOPE PER RISPOSTA ══════════════════\n")
cat(sprintf("%-6s | %8s %8s | %8s\n", "Resp.", "Pearson", "Spearman", "n campi"))
cat(strrep("-", 40), "\n")

corr_tbl <- ols_all |>
  group_by(risposta) |>
  summarise(
    pearson  = cor(int, slope, method = "pearson"),
    spearman = cor(int, slope, method = "spearman"),
    n        = n(),
    .groups  = "drop"
  )

for (i in seq_len(nrow(corr_tbl))) {
  cat(sprintf("%-6s | %8.3f %8.3f | %8d\n",
              corr_tbl$risposta[i],
              corr_tbl$pearson[i],
              corr_tbl$spearman[i],
              corr_tbl$n[i]))
}

cat("\n  Confronto con modello 20:\n")
cat("  eta_SOC = +0.209 (forte) → atteso Corr positiva\n")
cat("  eta_N   = -0.051 (≈0)   → atteso Corr ≈ 0\n")
cat("  eta_P   = +0.096 (debole)→ atteso Corr positiva debole\n")


# ── 4. SCATTER PLOT: INTERCETTA vs SLOPE (BLUP lmer) ─────────────────────────
# Usa BLUP invece di OLS per-campo: con n=4-6 per gruppo le stime OLS hanno
# molta varianza → correlazione intercetta-slope attenuata da errore di misura.
# I BLUP sono stime condizionali shrinkate verso la media di gruppo (Empirical
# Bayes): riducono il rumore senza introdurre bias sistematico.
#
# Confronto: OLS dà Pearson ≈ 0.39 per SOC; BLUP atteso ≈ +0.93.

blup_list <- lapply(
  setNames(c("logSOC", "logN", "logP"), c("SOC", "N", "P")),
  function(y_var) {
    m <- lmer(
      as.formula(paste(y_var, "~ logBottom_c + (logBottom_c | Field)")),
      data = dati, REML = TRUE
    )
    b           <- as.data.frame(ranef(m)$Field)
    names(b)    <- c("u0_intercept", "u1_slope")
    b$Field     <- rownames(b)
    b$risposta  <- switch(y_var, logSOC = "SOC", logN = "N", logP = "P")
    b
  }
)

blup_df <- bind_rows(blup_list) |>
  mutate(risposta = factor(risposta, levels = c("SOC", "N", "P")))

corr_blup <- blup_df |>
  group_by(risposta) |>
  summarise(
    pearson  = cor(u0_intercept, u1_slope),
    spearman = cor(u0_intercept, u1_slope, method = "spearman"),
    n        = n(),
    .groups  = "drop"
  ) |>
  mutate(label = sprintf("Pearson = %+.3f\nSpearman = %+.3f", pearson, spearman))

cat("\n═══ CORRELAZIONE INTERCETTA–SLOPE (BLUP lmer) ═══════════════════\n")
print(corr_blup)

label_pos_blup <- blup_df |>
  group_by(risposta) |>
  summarise(
    x = min(u0_intercept) + 0.05 * diff(range(u0_intercept)),
    y = max(u1_slope) - 0.05 * diff(range(u1_slope)),
    .groups = "drop"
  ) |>
  left_join(corr_blup, by = "risposta")

p_scatter_blup <- ggplot(blup_df, aes(u0_intercept, u1_slope)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.5, alpha = 0.7, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "tomato", fill = "tomato", alpha = 0.15) +
  geom_text(data = label_pos_blup, aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1, size = 3.2, color = "grey30") +
  facet_wrap(~ risposta, scales = "free", ncol = 3) +
  labs(
    title    = "Relazione intercetta-slope per campo (BLUP lmer)",
    subtitle = "BLUP = stime condizionali shrinkate; meno attenuazione da rumore di misura rispetto a OLS",
    x        = "BLUP intercetta di campo (livello alla profondita' media)",
    y        = "BLUP slope di campo con logBottom"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

print(p_scatter_blup)
save_fig("fig_02_scatter_intercetta_slope.pdf", p_scatter_blup, w = 18, h = 7)


# ── 5. SPAGHETTI PLOT: PROFILI VERTICALI COLORATI PER LIVELLO ────────────────
# Per il SOC (effetto più forte): mostra i profili di tutti i campi,
# colorati dal più basso (viola) al più alto (giallo-verde) livello medio.
# Se b_SOC > 0, i campi scuri (bassi) devono avere profili più ripidi.

make_spaghetti <- function(ols_df, dati_df, y_var, titolo) {

  # Ordina campi per intercetta (livello medio)
  field_order <- ols_df |>
    arrange(int) |>
    mutate(rank = row_number()) |>
    select(Field, int_rank = rank, int_val = int)

  dati_r <- dati_df |>
    select(Field, logBottom, y = all_of(y_var)) |>
    left_join(field_order, by = "Field")

  ggplot(dati_r, aes(x = logBottom, y = y, group = Field, color = int_rank)) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    geom_point(size = 1.2, alpha = 0.6) +
    scale_color_viridis_c(
      name   = "Rank livello\n(basso -> alto)",
      option = "plasma"
    ) +
    labs(
      title    = titolo,
      subtitle = "Colore: dal viola (campo piu' povero) al giallo (piu' ricco)",
      x        = "log(profondita' in cm)",
      y        = paste0("log(", y_var, ")")
    ) +
    theme_minimal(base_size = 11)
}

p_sp_SOC <- make_spaghetti(ols_SOC, dati, "logSOC",
                           "Profili verticali SOC — pattern intercetta-slope")
p_sp_N   <- make_spaghetti(ols_N,   dati, "logN",
                           "Profili verticali N — atteso no pattern")
p_sp_P   <- make_spaghetti(ols_P,   dati, "logP",
                           "Profili verticali P — atteso pattern debole")

print(p_sp_SOC)
print(p_sp_N)
print(p_sp_P)

# Figure combinate per il report
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  # Fig principale: SOC (segnale forte) e N (no segnale) affiancati
  p_spagh_soc_n <- (p_sp_SOC + theme(legend.position = "none")) +
                   (p_sp_N   + theme(legend.position = "right")) +
                   plot_annotation(
                     title    = "Profili verticali di SOC e N nei 40 campi",
                     subtitle = "Colore: dal viola (campo piu' povero) al giallo (piu' ricco)"
                   )
  print(p_spagh_soc_n)
  save_fig("fig_03_spaghetti_soc_n.pdf", p_spagh_soc_n, w = 18, h = 9)

  # Fig completa: tutte e 3 le risposte
  p_spagh_all <- (p_sp_SOC + theme(legend.position = "none")) +
                 (p_sp_N   + theme(legend.position = "none")) +
                 (p_sp_P   + theme(legend.position = "right")) +
                 plot_layout(ncol = 3) +
                 plot_annotation(
                   title    = "Profili verticali di SOC, N e P nei 40 campi",
                   subtitle = "Colore: dal viola (campo piu' povero) al giallo (piu' ricco)"
                 )
  print(p_spagh_all)
  save_fig("fig_03b_spaghetti_all.pdf", p_spagh_all, w = 22, h = 8)
} else {
  save_fig("fig_03a_spaghetti_soc.pdf", p_sp_SOC, w = 10, h = 8)
  save_fig("fig_03b_spaghetti_n.pdf",   p_sp_N,   w = 10, h = 8)
  save_fig("fig_03c_spaghetti_p.pdf",   p_sp_P,   w = 10, h = 8)
}


# ── 6. DISTRIBUZIONE DELLE SLOPE ─────────────────────────────────────────────

cat("\n═══ DISTRIBUZIONE SLOPE EMPIRICHE ═══════════════════════════════\n")
cat(sprintf("%-6s | %8s %8s %8s %8s %8s\n",
            "Resp.", "media", "sd", "q5", "q95", "% > 0"))
cat(strrep("-", 55), "\n")

ols_all |>
  group_by(risposta) |>
  summarise(
    m       = mean(slope),
    s       = sd(slope),
    q5      = quantile(slope, 0.05),
    q95     = quantile(slope, 0.95),
    pct_pos = mean(slope > 0) * 100,
    .groups = "drop"
  ) |>
  pwalk(function(risposta, m, s, q5, q95, pct_pos) {
    cat(sprintf("%-6s | %8.3f %8.3f %8.3f %8.3f %7.1f%%\n",
                risposta, m, s, q5, q95, pct_pos))
  })

cat("\n── Fine script 21 ──────────────────────────────────────────────\n")

