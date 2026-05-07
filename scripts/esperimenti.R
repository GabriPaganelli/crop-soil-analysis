par(mfrow=c(1,2))

hist(crop$PercSand, breaks = 20)
hist(crop$PercSand %>% log(base = 10), breaks = 20)

hist(crop$PercTotNitro, breaks = 20)
hist(crop$PercTotNitro %>% log(base = 10), breaks = 20)

hist(crop$PercSOC, breaks = 20)
hist(crop$PercSOC %>% log(base = 10), breaks = 20)

hist(crop$PercTotPhos, breaks = 20)
hist(crop$PercTotPhos %>% log(base = 10), breaks = 20)

hist(crop$PercClay, breaks = 20)
hist(crop$PercClay %>% log(base = 10), breaks = 20)

hist(crop$PercSilt, breaks = 20)
hist(crop$PercSilt %>% log(base = 10), breaks = 20)

hist(crop$PH, breaks = 20)
hist(crop$PH %>% log(base = 10), breaks = 20)


# Bimodalità di BulkDensity
plot_x_vs_y(crop, 'BulkDensity')

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



