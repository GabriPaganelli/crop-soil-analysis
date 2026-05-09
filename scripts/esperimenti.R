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


# ==============================================================================
# VERIFICA EMPIRICA DELLE OSSERVAZIONI ESPLORATIVE
# ==============================================================================
library(lme4)
library(dplyr)

target_vars <- c("PercSOC", "PercTotNitro", "PercTotPhos")

# Bottom centrato: riduce la correlazione intercept-slope nel mixed model
crop$Bottom_c <- crop$Bottom - mean(crop$Bottom)


# ── 1. PROFILI VERTICALI COLORATI PER LANDUSE ─────────────────────────────────
# Già fatto in eda.R colorando per Field. Qui vogliamo vedere se le traiettorie
# si separano visivamente per classe gestionale, sia in scala originale
# che in log: nel secondo caso, linee approssimativamente rette suggeriscono
# un decadimento esponenziale (log Y ~ lineare in Bottom).

plot_depth_profiles(crop, target_vars, color_var = "Landuse", ncol = 3)
plot_depth_profiles(crop, target_vars, color_var = "Landuse", log_y = TRUE, ncol = 3)


# ── 2. CONFONDIMENTO TEXTURE × LANDUSE (variabili continue) ───────────────────
# La tabella Texture × Landuse (già vista) mostra forte associazione.
# I boxplot delle componenti continue clay/silt/sand mostrano quanto siano
# distinguibili le landuse in termini di substrato pedologico.
# Se le distribuzioni si sovrappongono poco → texture e management sono quasi
# indistinguibili statisticamente.

vars_texture <- c("PercClay", "PercSilt", "PercSand", "PH", "BulkDensity")

plots_tex <- lapply(vars_texture, function(v) {
  ggplot(crop, aes(x = Landuse, y = .data[[v]], fill = Landuse)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, show.legend = FALSE) +
    labs(x = "Landuse", y = v, title = v) +
    theme_minimal(base_size = 10)
})
wrap_plots(plots_tex, ncol = 3) +
  plot_annotation(
    title    = "Distribuzione delle covariate fisiche per Landuse",
    subtitle = "Overlap ridotto → texture e management sono statisticamente confusi"
  )


# ── 3. STRUTTURA DEL DESIGN ────────────────────────────────────────────────────
# Verifica che Landuse sia determinata quasi deterministicamente dai 3 binari
# e che OnFarm = 1 implichi sempre Fertilised = 1.

crop |>
  distinct(Field, Landuse, Irrigate, Fertilised, OnFarm) |>
  arrange(Landuse)

# Cella (OnFarm=1, Fertilised=0) deve essere vuota
table(OnFarm = crop$OnFarm, Fertilised = crop$Fertilised)

# Struttura depth per Landuse: chi ha Bottom=80 e chi no
table(Landuse = crop$Landuse, Bottom = crop$Bottom)


# ── 4. ICC VERA E DECOMPOSIZIONE DELLA VARIANZA ────────────────────────────────
# La pseudo-ICC in eda.R confonde varianza verticale (sistematica) con rumore.
# Il modello nullo lme4 separa correttamente varianza tra field da varianza
# residua (entro field, lungo le profondità).
# m0: solo random intercept
# m1: aggiunto Bottom_c → atteso che sigma_field aumenti (emerge l'eterogeneità
#     tra field, prima "mascherata" dal trend verticale condiviso)

m0 <- lmer(log(PercSOC) ~ 1            + (1 | Field), data = crop, REML = TRUE)
m1 <- lmer(log(PercSOC) ~ Bottom_c     + (1 | Field), data = crop, REML = TRUE)

vc0 <- as.data.frame(VarCorr(m0))
icc <- vc0$vcov[1] / sum(vc0$vcov)
cat("ICC intra-Field (log SOC, modello nullo):", round(icc, 3), "\n\n")

cat("VarCorr m0 (nullo):\n");         print(VarCorr(m0))
cat("\nVarCorr m1 (+ Bottom_c):\n");  print(VarCorr(m1))


# ── 5. RANDOM SLOPE: "chi parte alto scende più in fretta" ────────────────────
# Con Bottom centrato la convergenza è stabile (il modello non centrato dava
# un warning di convergenza).
# ATTENZIONE alla centratura: Corr(u0j, u1j) cambia segno rispetto al modello
# non centrato. Con Bottom_c l'intercetta = SOC alla PROFONDITÀ MEDIA (~45 cm),
# non alla superficie:
#   - Corr > 0 (centrato):  field con più SOC a ~45 cm decadono MENO da lì in poi
#   - Corr < 0 (non centr.): field con più SOC stimato a 0 cm decadono più in fretta
# I due modelli sono matematicamente equivalenti; la correlazione del non centrato
# cattura meglio l'osservazione qualitativa "chi parte alto scende più in fretta".
# Il BLUP plot visualizza la relazione intercetta-slope per tutti i 40 field.

m2 <- lmer(log(PercSOC) ~ Bottom_c + (Bottom_c | Field), data = crop, REML = TRUE)
cat("\nVarCorr m2 (random slope su Bottom_c):\n"); print(VarCorr(m2))

blup <- as.data.frame(ranef(m2)$Field)
names(blup) <- c("u0_intercept", "u1_slope")

ggplot(blup, aes(x = u0_intercept, y = u1_slope)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title    = "BLUP: intercetta vs slope per Field (log SOC ~ Bottom_c)",
    subtitle = "Pendenza negativa = field ricchi in superficie decadono più rapidamente",
    x        = "Random intercept (u0j)",
    y        = "Random slope su Bottom_c (u1j)"
  ) +
  theme_minimal()


# ── 6. NON-LINEARITÀ DEL PROFILO VERTICALE ────────────────────────────────────
# Il modello con factor(Bottom) non impone nessuna forma funzionale:
# i coefficienti stimano il profilo medio in modo non parametrico.
# Li confrontiamo per AIC (ML, non REML) con le alternative lineari e quadratica.
# RISULTATO: il modello quadratico ha AIC più basso del fattoriale → la non-
# linearità è reale ma catturabile con una semplice parabola (parsimonia!).
# Il profilo fattoriale mostra la forma: forte caduta fino a ~50 cm, poi plateau.

m_lin  <- lmer(log(PercSOC) ~ Bottom_c               + (1|Field), data=crop, REML=FALSE)
m_quad <- lmer(log(PercSOC) ~ Bottom_c + I(Bottom_c^2) + (1|Field), data=crop, REML=FALSE)
m_log  <- lmer(log(PercSOC) ~ log(Bottom)             + (1|Field), data=crop, REML=FALSE)
m_fact <- lmer(log(PercSOC) ~ factor(Bottom)           + (1|Field), data=crop, REML=FALSE)
# m_log:  log(Y) ~ log(Bottom) → Y = A·Bottom^β  (power-law, non esponenziale)
# m_lin:  log(Y) ~ Bottom      → Y = A·e^{β·Bottom} (esponenziale)
# Sono forme diverse: l'esponenziale ha tasso di decadimento costante per cm,
# il power-law rallenta più in fretta vicino alla superficie.

print(AIC(m_lin, m_quad, m_log, m_fact))

# Visualizza i coefficienti del modello fattoriale come profilo medio
coef_fact <- fixef(m_fact)
depth_levels <- c(20, 30, 40, 50, 60, 80)
# primo livello è l'intercetta (Bottom=20), i restanti sono differenze rispetto a 20
profilo <- data.frame(
  Bottom = depth_levels,
  log_SOC_medio = c(coef_fact[1], coef_fact[1] + coef_fact[-1])
)
ggplot(profilo, aes(x = Bottom, y = log_SOC_medio)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue") +
  labs(
    title    = "Profilo medio di log(SOC) per profondità (fixed effects, factor(Bottom))",
    subtitle = "Nessuna forma funzionale imposta: evidenzia plateau in profondità",
    x        = "Bottom (cm)",
    y        = "log(SOC) stimato"
  ) +
  theme_minimal()


# ── 7. ETEROSCHEDASTICITÀ LUNGO LA PROFONDITÀ ─────────────────────────────────
# Se la varianza residua di log(SOC) decresce con la profondità, questo
# giustifica strutture di varianza esplicite (varIdent, varExp in nlme).
# Usiamo i residui di m1 (trend verticale rimosso, random intercept condizionato).

df_res <- data.frame(
  Bottom  = crop$Bottom,
  residuo = residuals(m1)
)

ggplot(df_res, aes(x = factor(Bottom), y = abs(residuo))) +
  geom_boxplot(fill = "steelblue", alpha = 0.65, outlier.size = 0.8) +
  labs(
    title    = "|Residui| di log(SOC) ~ Bottom_c + (1|Field) per profondità",
    subtitle = "Dispersione decrescente con la profondità → eteroschedasticità",
    x        = "Bottom (cm)",
    y        = "|Residuo condizionato|"
  ) +
  theme_minimal()
