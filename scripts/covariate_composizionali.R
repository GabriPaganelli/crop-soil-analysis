# ── 8. STRUTTURA DI COVARIANZA STRATIFICATA ───────────────────────────────────
# L'idea è capire se la correlazione tra SOC, N e P cambia
# a seconda della classe di gestione o della profondità

library(corrplot)

# --- 8a. Correlazioni per classe di gestione (Landuse) ---

cor_per_landuse <- crop |>
  group_split(Landuse) |>
  setNames(levels(crop$Landuse)) |>
  lapply(\(d) cor(d[, target_vars], use = "complete.obs"))

par(mfrow = c(2, 4), mar = c(0, 0, 2, 0))
for (nome in names(cor_per_landuse)) {
  corrplot(cor_per_landuse[[nome]],
           method     = "color",
           type       = "upper",
           addCoef.col = "black",
           number.cex = 0.9,
           tl.col     = "black",
           tl.cex     = 0.8,
           cl.lim     = c(-1, 1),
           col        = colorRampPalette(c("tomato", "white", "steelblue"))(100),
           title      = nome,
           mar        = c(0, 0, 1.5, 0))
}
par(mfrow = c(1, 1))

# --- 8b. Correlazioni per profondità (Bottom) ---

cor_per_depth <- crop |>
  group_split(Bottom) |>
  setNames(paste0(sort(unique(crop$Bottom)), " cm")) |>
  lapply(\(d) cor(d[, target_vars], use = "complete.obs"))

par(mfrow = c(2, 3), mar = c(0, 0, 2, 0))
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
par(mfrow = c(1, 1))

# --- 8c. Tabelle numeriche e grafici di variazione ---

# Correlazioni pairwise per Landuse
tab_cor_landuse <- crop |>
  group_by(Landuse) |>
  summarise(
    `SOC-N` = cor(PercSOC, PercTotNitro, use = "complete.obs"),
    `SOC-P` = cor(PercSOC, PercTotPhos,  use = "complete.obs"),
    `N-P`   = cor(PercTotNitro, PercTotPhos, use = "complete.obs"),
    n = n(),
    .groups = "drop"
  )
print(tab_cor_landuse)

# Correlazioni pairwise per profondità
tab_cor_depth <- crop |>
  group_by(Bottom) |>
  summarise(
    `SOC-N` = cor(PercSOC, PercTotNitro, use = "complete.obs"),
    `SOC-P` = cor(PercSOC, PercTotPhos,  use = "complete.obs"),
    `N-P`   = cor(PercTotNitro, PercTotPhos, use = "complete.obs"),
    n = n(),
    .groups = "drop"
  )
print(tab_cor_depth)

# Come cambiano le correlazioni con la profondità?
tab_cor_depth |>
  pivot_longer(cols = c(`SOC-N`, `SOC-P`, `N-P`),
               names_to  = "coppia",
               values_to = "correlazione") |>
  ggplot(aes(x = Bottom, y = correlazione, color = coppia, group = coppia)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(breaks = sort(unique(crop$Bottom))) +
  labs(title  = "Come cambiano le correlazioni tra SOC, N e P con la profondità",
       x      = "Profondità (cm)",
       y      = "Correlazione di Pearson",
       color  = "Coppia") +
  theme_minimal()

# Come cambiano per classe di gestione?
tab_cor_landuse |>
  pivot_longer(cols = c(`SOC-N`, `SOC-P`, `N-P`),
               names_to  = "coppia",
               values_to = "correlazione") |>
  ggplot(aes(x = Landuse, y = correlazione, fill = coppia)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Correlazioni tra SOC, N e P per classe di gestione",
       x = "Landuse", y = "Correlazione di Pearson", fill = "Coppia") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


## Covariate composizionali

library(compositions)

# Oggetto composizionale
gran <- acomp(crop[, c("PercClay", "PercSilt", "PercSand")])

# Opzione A: CLR + PCA
clr_gran <- clr(gran)
pca_gran <- princomp(clr_gran)
summary(pca_gran)# quanta varianza spiega PC1?
pca_gran$loadings
biplot(pca_gran)           # interpretazione delle componenti

# Aggiunge i punteggi al dataset
crop$PC1_gran <- pca_gran$scores[, 1]
crop$PC2_gran <- pca_gran$scores[, 2]

# Opzione B: ILR diretto (equivalente, più pulito)
ilr_gran <- ilr(gran)
crop$ILR1 <- ilr_gran[, 1]
crop$ILR2 <- ilr_gran[, 2]

plot(ilr_gran)
