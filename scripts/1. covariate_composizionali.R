library(tidyverse)
library(compositions)
library(soiltexture)
library(patchwork)
library(here)

crop <- readRDS(here("data", "crop.rds"))
dati <- readRDS(here("data", "dati.rds"))


# ── DATI COMPOSIZIONALI ───────────────────────────────────────────────────────
# PercClay, PercSilt, PercSand sommano sempre a 100%: vivono su un simplesso
# 2D, non in R³. Usarle grezze in un modello causa due problemi:
#   1. dipendenza lineare perfetta (rango 2, non 3)
#   2. correlazioni spurie per costruzione (Pearson, 1897): se Clay sale,
#      Sand e Silt devono scendere per definizione, non per motivi biologici
#
# Strategia:
#   • CLR + PCA  →  interpretazione visiva (biplot, loadings)
#   • ILR        →  2 coordinate ortogonali in R² per i modelli

gran <- acomp(dati[, c("PercClay", "PercSilt", "PercSand")])


# ── CLR + PCA: interpretazione ────────────────────────────────────────────────
# La CLR proietta la composizione in R³ con vincolo di somma zero.
# La PCA estrae 2 componenti reali (la terza ha varianza nulla per costruzione).
# Lettura tipica dei loadings:
#   PC1 ≈ gradiente finezza  (argilla+limo vs sabbia)
#   PC2 ≈ contrasto argilla vs limo nella frazione fine

clr_gran <- clr(gran)
pca_gran <- princomp(clr_gran)

summary(pca_gran)    # varianza spiegata per componente
pca_gran$loadings    # conferma l'interpretazione
biplot(pca_gran)


# ── ILR: coordinate per i modelli ─────────────────────────────────────────────
# L'ILR produce esattamente D-1 = 2 coordinate ortogonali in R² senza vincoli.
# Si comportano come qualsiasi predittore continuo in mixed effects, regressione, ecc.
#
#   Texture1 ≈ gradiente finezza (argilla+limo vs sabbia)
#   Texture2 ≈ contrasto argilla vs limo nella frazione fine

ilr_gran <- ilr(gran)


# ── CLR+PCA e ILR: stessa informazione, assi ruotate ─────────────────────────
# CLR+PCA è una rotazione data-driven dello spazio ILR: il pattern visivo è
# identico nei due scatter, cambiano solo orientamento e scala degli assi.

# Classe texture USDA per colorare i punti (calcolata al volo, non salvata in dati)
mat_tex     <- TT.points.in.classes(
  tri.data  = data.frame(SAND = dati$PercSand, SILT = dati$PercSilt, CLAY = dati$PercClay),
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
       subtitle = "Coordinate usate nei modelli",
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

p_ilr + p_pca +
  plot_annotation(
    title    = "ILR e CLR+PCA codificano la stessa informazione",
    subtitle = "Il pattern è identico — le assi sono ruotate",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )


# ── Aggiornamento di dati ─────────────────────────────────────────────────────
# Sostituisce le 3 percentuali composizionali con le 2 coordinate ILR.

dati <- dati |>
  mutate(
    Texture1 = ilr_gran[, 1],
    Texture2 = ilr_gran[, 2]
  ) |>
  select(-PercClay, -PercSilt, -PercSand)

saveRDS(dati, here("data", "dati.rds"))
