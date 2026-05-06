# ── 1. CARICAMENTO ────────────────────────────────────────────────────────────
setwd("C:/Users/gabri/Desktop/Università/Statistica iterazione/Lavoro di gruppo")

library(soiltexture)

crop <- read.csv("crop.csv", sep = ";")

# ── 2. PREPROCESSING ──────────────────────────────────────────────────────────
crop <- crop |>
  select(-OwnId) |>
  mutate(across(c(Landuse, Field, Plot, OnFarm, Irrigate,
                  Fertilised, N_Natural, Class), as.factor))

# Prepara il dataframe con i nomi colonna richiesti dal pacchetto
# (Assicurati di sostituire df$sand, df$silt, df$clay con i nomi esatti delle tue colonne)
tex_data <- data.frame(
  SAND = crop$PercSand,
  SILT = crop$PercSilt,
  CLAY = crop$PercClay
)

# Classifica i punti usando il sistema USDA
classes <- TT.points.in.classes(tri.data = tex_data, class.sys = "USDA.TT")

# Estrae il nome della classe per ogni riga (se un punto è sul bordo, prende la prima classe)
crop$texture_class <- apply(classes, 1, function(x) names(x)[which(x > 0)][1])

head(crop)

# Vettore di conversione per le 12 classi USDA
usda_lookup <- c(
  "Cl"     = "clay",
  "SiCl"   = "silty clay",
  "SaCl"   = "sandy clay",
  "ClLo"   = "clay loam",
  "SiClLo" = "silty clay loam",
  "SaClLo" = "sandy clay loam",
  "Lo"     = "loam",
  "SiLo"   = "silt loam",
  "SaLo"   = "sandy loam",
  "Si"     = "silt",
  "LoSa"   = "loamy sand",
  "Sa"     = "sand"
)

# Applica la conversione alla colonna creata in precedenza
crop$texture <- usda_lookup[crop$texture_class]

crop$Class = NULL
crop$texture_class = NULL
colnames(crop) = c(colnames(crop)[1:18],'Texture')

head(crop)
saveRDS(crop, file = "crop.rds")
