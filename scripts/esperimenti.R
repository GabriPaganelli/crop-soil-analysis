op <- par(mfrow = c(1, 2))

hist(crop$PercSand,          breaks = 20)
hist(log10(crop$PercSand),   breaks = 20, main = "log10(PercSand)")

hist(crop$PercTotNitro,      breaks = 20)
hist(log10(crop$PercTotNitro), breaks = 20, main = "log10(PercTotNitro)")

hist(crop$PercSOC,           breaks = 20)
hist(log10(crop$PercSOC),    breaks = 20, main = "log10(PercSOC)")

hist(crop$PercTotPhos,       breaks = 20)
hist(log10(crop$PercTotPhos), breaks = 20, main = "log10(PercTotPhos)")

hist(crop$PercClay,          breaks = 20)
hist(log10(crop$PercClay),   breaks = 20, main = "log10(PercClay)")

hist(crop$PercSilt,          breaks = 20)
hist(log10(crop$PercSilt),   breaks = 20, main = "log10(PercSilt)")

hist(crop$PH,                breaks = 20)
hist(log10(crop$PH),         breaks = 20, main = "log10(PH)")

par(op)


# Bimodalità di BulkDensity
plot_density_y_by_x(crop, 'BulkDensity', numeric_x_bins = 'quartile')

# pca

library(FactoMineR)
library(factoextra)

# Variabili numeriche per la PCA
num_vars <- c("Bottom", "PH", "PercClay", "PercSilt", "PercSand",
              "PercTotNitro", "PercSOC", "PercTotPhos", "BulkDensity")

# Variabili qualitative supplementari
quali_vars <- c("Landuse", "OnFarm", "Irrigate", "Fertilised", "N_Natural", "Texture")

# Indici delle colonne quali supplementari nel subset
df_pca <- crop[, c(num_vars, quali_vars)]
quali_idx <- (length(num_vars) + 1):ncol(df_pca)

# PCA
res_pca <- PCA(df_pca,
               scale.unit  = TRUE,
               quali.sup   = quali_idx,
               graph       = FALSE)

# Scree plot
fviz_eig(res_pca, addlabels = TRUE, ncp = 10)

# Biplot individui + variabili
fviz_pca_biplot(res_pca,
                habillage   = "Texture",   # colora per Landuse, cambia a piacere
                addEllipses = TRUE,
                ellipse.level = 0.95,
                label       = "var",
                repel       = TRUE,
                title       = "PCA – Biplot (Dim1 × Dim2)")

# Contributo delle variabili
fviz_contrib(res_pca, choice = "var", axes = 1, title = "Contributo variabili – PC1")
fviz_contrib(res_pca, choice = "var", axes = 2, title = "Contributo variabili – PC2")

# Cos² variabili
fviz_cos2(res_pca, choice = "var", axes = 1:2)

# Summary numerico
summary(res_pca, nbelements = length(num_vars))



