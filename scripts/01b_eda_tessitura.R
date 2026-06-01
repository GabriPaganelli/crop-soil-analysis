# =============================================================================
# 01b_eda_tessitura.R  —  EDA composizionale tessitura (ILR, triangolo)
#
# Motivazione: PercClay, PercSilt, PercSand sommano sempre a 100% → vivono su
# un simplesso 2D, non in R³. Usarle grezze in un modello causa due problemi:
#   1. dipendenza lineare perfetta (rango 2, non 3)
#   2. correlazioni spurie per costruzione (Pearson, 1897): se Clay sale,
#      Sand e Silt devono scendere per definizione, non per motivi biologici
#
# Questo script mostra:
#   1. CLR + PCA: interpretazione visiva (biplot, loadings)
#   2. ILR: 2 coordinate ortogonali usate nei modelli (Texture1, Texture2)
#   3. Confronto ILR / CLR+PCA: stessa informazione, assi ruotate
#   4. Triangolo tessitura USDA colorato per campo
#
# Dipende da: data/dati.rds (prodotto da 01_preprocessing.R)
# =============================================================================

library(tidyverse)
library(compositions)
library(soiltexture)
library(patchwork)
library(here)
source(here("scripts", "00_utilities.R"))


# ── 1. DATI ───────────────────────────────────────────────────────────────────

dati <- readRDS(here("data", "dati.rds"))


# ── 2. CLR + PCA: interpretazione ────────────────────────────────────────────
# La CLR proietta la composizione in R³ con vincolo di somma zero.
# La PCA estrae 2 componenti reali (la terza ha varianza nulla per costruzione).
# Lettura tipica dei loadings:
#   PC1 ≈ contrasto argilla vs limo nella frazione fine
#   PC2 ≈ gradiente finezza (argilla+limo vs sabbia)

gran     <- acomp(dati[, c("PercClay", "PercSilt", "PercSand")])
clr_gran <- clr(gran)
pca_gran <- princomp(clr_gran)

cat("══════════════════════════════════════════════════\n")
cat(" VARIANZA SPIEGATA — CLR+PCA\n")
cat("══════════════════════════════════════════════════\n")
summary(pca_gran)
cat("\nLoadings:\n")
print(pca_gran$loadings)

biplot(pca_gran, main = "Biplot CLR+PCA — tessitura")


# ── 3. ILR: coordinate per i modelli ─────────────────────────────────────────
# L'ILR produce esattamente D-1 = 2 coordinate ortogonali in R² senza vincoli.
# Si comportano come qualsiasi predittore continuo in mixed effects, ecc.
#
#   Texture1 (ILR1) ≈ contrasto argilla vs limo nella frazione fine
#   Texture2 (ILR2) ≈ gradiente finezza (argilla+limo vs sabbia)

ilr_gran <- ilr(gran)

cat("\n══════════════════════════════════════════════════\n")
cat(" STATISTICHE ILR\n")
cat("══════════════════════════════════════════════════\n")
cat(sprintf("Texture1: media = %.3f | sd = %.3f | range [%.3f, %.3f]\n",
            mean(ilr_gran[, 1]), sd(ilr_gran[, 1]),
            min(ilr_gran[, 1]), max(ilr_gran[, 1])))
cat(sprintf("Texture2: media = %.3f | sd = %.3f | range [%.3f, %.3f]\n",
            mean(ilr_gran[, 2]), sd(ilr_gran[, 2]),
            min(ilr_gran[, 2]), max(ilr_gran[, 2])))
cat(sprintf("Correlazione Texture1-Texture2: %.3f\n",
            cor(ilr_gran[, 1], ilr_gran[, 2])))
cat("  → Se |r| ≈ 0, le coordinate sono ortogonali come atteso dall'ILR\n")


# ── 4. CONFRONTO ILR / CLR+PCA ───────────────────────────────────────────────
# CLR+PCA è una rotazione data-driven dello spazio ILR: il pattern visivo è
# identico nei due scatter, cambiano solo orientamento e scala degli assi.

mat_tex <- TT.points.in.classes(
  tri.data  = data.frame(SAND = dati$PercSand, SILT = dati$PercSilt,
                         CLAY = dati$PercClay),
  class.sys = "USDA.TT"
)
texture_tmp <- apply(mat_tex, 1, function(x) {
  nm <- names(x)[x > 0]
  if (length(nm)) nm[1] else NA_character_
})

df_confronto <- data.frame(
  ILR1    = ilr_gran[, 1],
  ILR2    = ilr_gran[, 2],
  PC1     = pca_gran$scores[, 1],
  PC2     = pca_gran$scores[, 2],
  Texture = texture_tmp
)

p_ilr <- ggplot(df_confronto, aes(x = ILR1, y = ILR2, colour = Texture)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title    = "Spazio ILR",
       subtitle = "Coordinate usate nei modelli (Texture1, Texture2)",
       x        = "ILR1  (Texture1)",
       y        = "ILR2  (Texture2)") +
  theme_minimal() +
  theme(legend.position = "none")

p_pca <- ggplot(df_confronto, aes(x = PC1, y = PC2, colour = Texture)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title    = "Spazio CLR+PCA",
       subtitle = "Rotazione data-driven dello stesso spazio",
       x        = "PC1",
       y        = "PC2") +
  theme_minimal()

print(
  p_ilr + p_pca +
    plot_annotation(
      title    = "ILR e CLR+PCA codificano la stessa informazione",
      subtitle = "Il pattern è identico — le assi sono ruotate",
      theme    = theme(plot.title = element_text(face = "bold", size = 13))
    )
)


# ── 5. TRIANGOLO TESSITURA USDA ───────────────────────────────────────────────
# Mostra la distribuzione delle osservazioni sul simplesso.
# plot_texture_triangle() definita in 00_utilities.R.

cat("\n══════════════════════════════════════════════════\n")
cat(" DISTRIBUZIONE CLASSI TESSITURA USDA\n")
cat("══════════════════════════════════════════════════\n")
print(
  dati |>
    mutate(texture = texture_tmp) |>
    count(texture, sort = TRUE)
)

plot_texture_triangle(dati, version = "ttplot")

cat("\n── Fine script 01b ─────────────────────────────────────────────\n")
