# =============================================================================
# 02_eda_tessitura.R  —  EDA composizionale tessitura (ILR, triangolo)
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
# dati.rds non ha PercClay/Silt/Sand (rimosse dopo calcolo ILR).
# Usiamo crop.rds (dataset completo) per le analisi di tessitura grezza.

dati <- readRDS(here("data", "dati.rds"))  # per Texture1/Texture2/altri dati
crop <- readRDS(here("data", "crop.rds"))  # per PercClay/PercSilt/PercSand


# ── 2. CLR + PCA: interpretazione ────────────────────────────────────────────
# La CLR proietta la composizione in R³ con vincolo di somma zero.
# La PCA estrae 2 componenti reali (la terza ha varianza nulla per costruzione).
# Lettura tipica dei loadings:
#   PC1 ≈ contrasto argilla vs limo nella frazione fine
#   PC2 ≈ gradiente finezza (argilla+limo vs sabbia)

gran     <- acomp(crop[, c("PercClay", "PercSilt", "PercSand")])
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
  tri.data  = data.frame(SAND = crop$PercSand, SILT = crop$PercSilt,
                         CLAY = crop$PercClay),
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
  crop |>
    mutate(texture = texture_tmp) |>
    count(texture, sort = TRUE)
)

# ── Classe tessiturale USDA per ogni osservazione ────────────────────────────
mat_tex_full <- TT.points.in.classes(
  tri.data  = data.frame(SAND = crop$PercSand, SILT = crop$PercSilt,
                         CLAY = crop$PercClay),
  class.sys = "USDA.TT"
)
texture_class <- apply(mat_tex_full, 1, function(x) {
  nm <- names(x)[x > 0]; if (length(nm)) nm[1] else NA_character_
})
crop$TextureClass <- factor(texture_class)

# ── Triangolo USDA (ggtern, punti colorati per classe) ────────────────────────
tex_ttplot <- data.frame(SAND = crop$PercSand, SILT = crop$PercSilt, CLAY = crop$PercClay)

p_triangle <- plot_texture_triangle(crop, color_var = "TextureClass",
                                    show_legend = FALSE, version = "ggtern")

# ── Palette Landuse (triangolo) e TextureClass (ILR) ─────────────────────────
land_lvls <- as.character(sort(unique(crop$Landuse)))
pal_land  <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#A65628", "#F781BF")[seq_along(land_lvls)],
  land_lvls
)

tex_lvls <- levels(crop$TextureClass)
base_pal <- c("#2166ac","#4393c3","#1a9850","#74c476","#fdae61",
              "#d73027","#762a83","#e08214","#a6761d","#666666",
              "#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")
tex_pal  <- setNames(base_pal[seq_len(min(length(tex_lvls), length(base_pal)))], tex_lvls)

# ── ILR biplot: colorato per TextureClass ────────────────────────────────────
df_ilr <- data.frame(
  Texture1 = ilr_gran[, 1],
  Texture2 = ilr_gran[, 2],
  Classe   = crop$TextureClass
)

p_ilr_report <- ggplot(df_ilr, aes(x = Texture1, y = Texture2, colour = Classe)) +
  geom_point(alpha = 0.75, size = 2.2) +
  scale_colour_manual(values = tex_pal, name = "Classe USDA") +
  labs(
    title = NULL,
    x     = expression("Texture"[1] ~ "(argilla / limo)"),
    y     = expression("Texture"[2] ~ "(finezza / sabbia)")
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

# ── Figura combinata per il report ────────────────────────────────────────────
dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)

# Salva il triangolo da solo
ggplot2::ggsave(here("output", "figures", "fig_texture_triangle.pdf"),
                plot = p_triangle, width = 14, height = 14, units = "cm",
                device = cairo_pdf)
cat("Salvato: output/figures/fig_texture_triangle.pdf\n")

# Salva il pannello combinato: TT.plot (Landuse) + ILR biplot (Landuse)
col_land <- pal_land[as.character(crop$Landuse)]

p_tri_ttplot <- wrap_elements(full = ~{
  TT.plot(
    class.sys      = "USDA.TT",
    tri.data       = tex_ttplot,
    bg             = "white", frame.bg.col = "white",
    class.p.bg.col = colori_trasparenti,
    class.line.col = "grey60",
    class.lab.show = "abr", class.lab.col = "#3A2408", cex.lab = 0.75,
    pch = 16, col = col_land, cex = 0.7, cex.axis = 0.6,
    arrows.show = FALSE
  )
  legend("topright", legend = land_lvls,
         col = pal_land, pch = 16, bty = "n", cex = 0.75, title = "Landuse")
})

p_tex_combined <- p_tri_ttplot + p_ilr_report +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(
    title = "Tessitura del suolo per categoria d'uso",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )
ggplot2::ggsave(here("output", "figures", "fig_texture_combined.pdf"),
                plot = p_tex_combined, width = 26, height = 13, units = "cm",
                device = cairo_pdf)
cat("Salvato: output/figures/fig_texture_combined.pdf\n")

cat("\n── Fine script 02 ─────────────────────────────────────────────\n")
