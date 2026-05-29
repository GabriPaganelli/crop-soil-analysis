# =============================================================================
# 01_preprocessing.R  —  Preprocessing del dataset crop
#
# Input:  data/crop.csv
# Output: data/crop.rds  (dataset completo con Texture USDA e ILR)
#         data/dati.rds  (subset analitico usato dai modelli)
#
# Operazioni:
#   1. Carica crop.csv, sistema i tipi
#   2. Classifica la texture secondo USDA (soiltexture)
#   3. Calcola Texture1, Texture2 via compositions::ilr()
#   4. Salva crop.rds e dati.rds
#
# ILR: compositions::ilr(acomp(Clay, Silt, Sand)), base di Helmert standard.
#   Texture1 = ILR1: contrasto argilla vs limo nella frazione fine
#   Texture2 = ILR2: gradiente finezza (argilla+limo vs sabbia)
#   Vedi 01b_eda_texture.R per la motivazione e l'analisi visiva.
# =============================================================================

library(tidyverse)
library(soiltexture)
library(here)


# ── 1. CARICAMENTO ────────────────────────────────────────────────────────────

crop <- read.csv(here("data", "crop.csv"), sep = ";")


# ── 2. TIPI ───────────────────────────────────────────────────────────────────

crop <- crop |>
  select(-OwnId) |>
  mutate(across(c(Landuse, Field, Plot, OnFarm, Irrigate,
                  Fertilised, N_Natural, Class), as.factor))


# ── 3. CLASSIFICAZIONE TEXTURE USDA ───────────────────────────────────────────
# Punti sul bordo tra due classi → prima classe adiacente.

tex_data <- data.frame(
  SAND = crop$PercSand,
  SILT = crop$PercSilt,
  CLAY = crop$PercClay
)

classes <- TT.points.in.classes(tri.data = tex_data, class.sys = "USDA.TT")

crop$texture_class <- apply(classes, 1, function(x) names(x)[which(x > 0)][1])

usda_lookup <- c(
  "Cl"     = "clay",            "SiCl"   = "silty clay",
  "SaCl"   = "sandy clay",      "ClLo"   = "clay loam",
  "SiClLo" = "silty clay loam", "SaClLo" = "sandy clay loam",
  "Lo"     = "loam",            "SiLo"   = "silt loam",
  "SaLo"   = "sandy loam",      "Si"     = "silt",
  "LoSa"   = "loamy sand",      "Sa"     = "sand"
)

crop$texture <- usda_lookup[crop$texture_class]
crop$Class        <- NULL
crop$texture_class <- NULL
crop <- rename(crop, Texture = texture)


# ── 4. ILR TEXTURE (Texture1, Texture2) ──────────────────────────────────────
# Clay + Silt + Sand ≈ 100 → vincolo compositivo: le tre frazioni non sono
# indipendenti, includerle tutte crea collinearità perfetta.
# ILR le trasforma in 2 coordinate ortogonali prive di vincolo.

gran     <- compositions::acomp(cbind(crop$PercClay, crop$PercSilt, crop$PercSand))
ilr_gran <- compositions::ilr(gran)

crop <- crop |>
  mutate(
    Texture1 = ilr_gran[, 1],
    Texture2 = ilr_gran[, 2]
  )


# ── 5. OUTPUT ─────────────────────────────────────────────────────────────────
# crop  → dataset completo (variabili originali + Texture USDA + ILR)
# dati  → subset analitico usato dai modelli Stan (no Lat/Long, no Landuse)

dati <- crop |>
  select(Field, OnFarm, Irrigate, Fertilised, N_Natural,
         Bottom, PH, PercClay, PercSilt, PercSand,
         PercTotNitro, PercSOC, PercTotPhos, BulkDensity,
         Texture1, Texture2)

saveRDS(crop, here("data", "crop.rds"))
saveRDS(dati, here("data", "dati.rds"))

cat(sprintf("Salvati crop.rds (%d righe, %d col) e dati.rds (%d righe, %d col)\n",
            nrow(crop), ncol(crop), nrow(dati), ncol(dati)))
