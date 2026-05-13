# ── 1. CARICAMENTO ────────────────────────────────────────────────────────────

library(tidyverse)
library(soiltexture)
library(here)

crop <- read.csv(here("data", "crop.csv"), sep = ";")

# ── 2. TIPI ───────────────────────────────────────────────────────────────────
crop <- crop |>
  select(-OwnId) |>
  mutate(across(c(Landuse, Field, Plot, OnFarm, Irrigate,
                  Fertilised, N_Natural, Class), as.factor))

# ── 3. CLASSIFICAZIONE TEXTURE USDA ───────────────────────────────────────────
# I punti sul bordo tra due classi vengono assegnati alla prima classe adiacente.
tex_data <- data.frame(
  SAND = crop$PercSand,
  SILT = crop$PercSilt,
  CLAY = crop$PercClay
)

classes <- TT.points.in.classes(tri.data = tex_data, class.sys = "USDA.TT")

crop$texture_class <- apply(classes, 1, function(x) names(x)[which(x > 0)][1])

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

crop$texture <- usda_lookup[crop$texture_class]
crop$Class <- NULL
crop$texture_class <- NULL
crop <- rename(crop, Texture = texture)

# ── 4. OUTPUT ─────────────────────────────────────────────────────────────────
# crop  → dataset completo (tutte le variabili originali + Texture)
# dati  → subset analitico (identificatori di campo + variabili numeriche)

dati <- crop |>
  select(Field, OnFarm, Irrigate, Fertilised, N_Natural,
         Bottom, PH, PercClay, PercSilt, PercSand,
         PercTotNitro, PercSOC, PercTotPhos, BulkDensity)

saveRDS(crop, here("data", "crop.rds"))
saveRDS(dati, here("data", "dati.rds"))
