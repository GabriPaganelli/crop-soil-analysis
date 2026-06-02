# ═══════════════════════════════════════════════════════════
# STEP 0 — Preparazione dati
# ═══════════════════════════════════════════════════════════

rm(list = ls())
library(nlme)
library(tidyverse)
library(here)

crop <- readRDS(here("data", "dati.rds"))

crop <- crop %>%
  mutate(
    logSOC = log(PercSOC),
    logN   = log(PercTotNitro),
    logP   = log(PercTotPhos),
    Field  = factor(Field),
    OnFarm = factor(OnFarm),
    Irrigate = factor(Irrigate),
    Fertilised = factor(Fertilised),
    N_Natural = factor(N_Natural)
  ) %>%
  mutate(across(c(PH, Texture1, Texture2, BulkDensity),
                ~ c(scale(.x))))

# Bottom resta in scala originale (cm)

# Formato long
crop_long <- crop %>%
  mutate(obs_id = row_number()) %>%
  pivot_longer(
    cols = c(logSOC, logN, logP),
    names_to = "response",
    values_to = "y"
  ) %>%
  mutate(response = factor(response))


# ═══════════════════════════════════════════════════════════
# STEP 1a — Trivariato base: R piena (corSymm) + G correlata
# ═══════════════════════════════════════════════════════════

fit_step1 <- lme(
  y ~ response:(Texture1 + Texture2 + Irrigate + Fertilised +
                  Bottom + BulkDensity + N_Natural + OnFarm + PH) - 1 + response,
  random = ~ 0 + response | Field,
  correlation = corSymm(form = ~ 1 | Field/obs_id),
  weights = varIdent(form = ~ 1 | response),
  data = crop_long,
  method = "REML"
)

cat("══════ STEP 1: TRIVARIATO BASE (R piena + G) ══════\n")
cat("AIC:", AIC(fit_step1), "\n\n")
summary(fit_step1)


# ═══════════════════════════════════════════════════════════
# STEP 1b — Visualizzazione correlazione residua R
# ═══════════════════════════════════════════════════════════

# Estrai correlazione residua da corSymm
cs <- corMatrix(fit_step1$modelStruct$corStruct)[[1]]
nomi <- c("logN", "logP", "logSOC")
rownames(cs) <- colnames(cs) <- nomi

cs_long <- cs %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "cor")

p_R <- ggplot(cs_long, aes(var1, var2, fill = cor)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = round(cor, 3)), size = 6, fontface = "bold") +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Correlazione residua R (within-field)",
       subtitle = "Valori prossimi a zero: R può essere diagonale",
       x = "", y = "", fill = "r") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_R)

# Correlazione G (random effects) per confronto
re <- ranef(fit_step1)
re_df <- as.data.frame(unclass(re))
re_cor <- cor(re_df)

re_long <- re_cor %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "cor")

p_G <- ggplot(re_long, aes(var1, var2, fill = cor)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = round(cor, 3)), size = 6, fontface = "bold") +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Correlazione G (between-field)",
       subtitle = "SOC-N correlati: la covarianza tra risposte è a livello di campo",
       x = "", y = "", fill = "r") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_G)


# ═══════════════════════════════════════════════════════════
# STEP 2 — Semplificazione: togli corSymm (R diagonale)
#
# Motivazione: le correlazioni residue sono trascurabili
# (0.005, 0.009, 0.112). Togliendo corSymm risparmiamo
# 3 parametri senza perdita di informazione.
# ═══════════════════════════════════════════════════════════

fit_step2 <- lme(
  y ~ response:(Texture1 + Texture2 + Irrigate + Fertilised +
                  Bottom + BulkDensity + N_Natural + OnFarm + PH) - 1 + response,
  random = ~ 0 + response | Field,
  weights = varIdent(form = ~ 1 | response),
  data = crop_long,
  method = "REML"
)

cat("══════ STEP 2: R DIAGONALE (senza corSymm) ══════\n")
cat("AIC step 1 (R piena):    ", AIC(fit_step1), "\n")
cat("AIC step 2 (R diagonale):", AIC(fit_step2), "\n")
cat("Differenza AIC:", AIC(fit_step1) - AIC(fit_step2), 
    "(negativo = R diagonale migliore)\n\n")
summary(fit_step2)

# ═══════════════════════════════════════════════════════════
# STEP 3 — Test corExp per ogni risposta (modelli separati)
#
# Scopo: verificare per quale risposta l'autocorrelazione
# verticale è significativa. Usiamo modelli univariati
# perché nel formato long non è possibile specificare
# corExp risposta-specifico.
# ═══════════════════════════════════════════════════════════

# --- logSOC: con e senza corExp ---
fit_SOC_base <- lme(
  logSOC ~ Texture1 + Texture2 + Irrigate + Fertilised +
    Bottom + BulkDensity + N_Natural + OnFarm + PH,
  random = ~ 1 | Field, data = crop, method = "REML"
)

fit_SOC_exp <- lme(
  logSOC ~ Texture1 + Texture2 + Irrigate + Fertilised +
    Bottom + BulkDensity + N_Natural + OnFarm + PH,
  random = ~ 1 | Field,
  correlation = corExp(form = ~ Bottom | Field),
  data = crop, method = "REML"
)

# --- logN: con e senza corExp ---
fit_N_base <- lme(
  logN ~ Texture1 + Texture2 + Irrigate + Fertilised +
    Bottom + BulkDensity + N_Natural + OnFarm + PH,
  random = ~ 1 | Field, data = crop, method = "REML"
)

fit_N_exp <- lme(
  logN ~ Texture1 + Texture2 + Irrigate + Fertilised +
    Bottom + BulkDensity + N_Natural + OnFarm + PH,
  random = ~ 1 | Field,
  correlation = corExp(form = ~ Bottom | Field),
  data = crop, method = "REML"
)

# --- logP: con e senza corExp ---
fit_P_base <- lme(
  logP ~ Texture1 + Texture2 + Irrigate + Fertilised +
    Bottom + BulkDensity + N_Natural + OnFarm + PH,
  random = ~ 1 | Field, data = crop, method = "REML"
)

fit_P_exp <- lme(
  logP ~ Texture1 + Texture2 + Irrigate + Fertilised +
    Bottom + BulkDensity + N_Natural + OnFarm + PH,
  random = ~ 1 | Field,
  correlation = corExp(form = ~ Bottom | Field),
  data = crop, method = "REML"
)

# --- Confronto LRT ---
cat("══════ STEP 3: TEST corExp PER RISPOSTA ══════\n")

cat("\n--- logSOC ---\n")
print(anova(fit_SOC_base, fit_SOC_exp))

cat("\n--- logN ---\n")
print(anova(fit_N_base, fit_N_exp))

cat("\n--- logP ---\n")
print(anova(fit_P_base, fit_P_exp))

cat("\n══════ RANGE STIMATI ══════\n")
cat("logSOC:", coef(fit_SOC_exp$modelStruct$corStruct, unconstrained = FALSE), "cm\n")
cat("logN:  ", coef(fit_N_exp$modelStruct$corStruct, unconstrained = FALSE), "cm\n")
cat("logP:  ", coef(fit_P_exp$modelStruct$corStruct, unconstrained = FALSE), "cm\n")


# ═══════════════════════════════════════════════════════════
# STEP 4 — corExp selettivo per logSOC nel modello trivariato
#
# Scopo: integrare la correlazione esponenziale nel modello
# trivariato mantenendo G correlata, ma applicando corExp
# solo a logSOC (dove è significativo).
#
# Tecnica: coordinate bidimensionali artificiali. La prima
# dimensione è la profondità (reale per SOC, gonfiata per
# N e P). La seconda dimensione separa le risposte con
# distanza enorme. corExp opera sulla distanza euclidea 2D:
# tra profondità diverse di SOC la distanza è reale (~10 cm),
# tra qualsiasi altra coppia è enorme (~5000) → α ≈ 0.
# ═══════════════════════════════════════════════════════════

crop_long <- crop_long %>%
  mutate(
    coord_bottom = case_when(
      response == "logSOC" ~ Bottom,
      response == "logN"   ~ Bottom * 100,
      response == "logP"   ~ Bottom * 100
    ),
    coord_resp = case_when(
      response == "logSOC" ~ 0,
      response == "logN"   ~ 5000,
      response == "logP"   ~ 10000
    )
  )

fit_step4 <- lme(
  y ~ response:(Texture1 + Texture2 + Irrigate + Fertilised +
                  Bottom + BulkDensity + N_Natural + OnFarm + PH) - 1 + response,
  random = ~ 0 + response | Field,
  correlation = corExp(value = 12, form = ~ coord_bottom + coord_resp | Field,
                       fixed = FALSE),
  weights = varIdent(form = ~ 1 | response),
  data = crop_long,
  method = "REML",
  control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim")
)

cat("══════ STEP 4: TRIVARIATO CON corExp SELETTIVO ══════\n")
cat("AIC step 2 (senza corExp):     ", AIC(fit_step2), "\n")
cat("AIC step 4 (corExp selettivo): ", AIC(fit_step4), "\n")
cat("Differenza:", AIC(fit_step2) - AIC(fit_step4), "\n\n")

cat("Range corExp stimato (cm):\n")
cat(coef(fit_step4$modelStruct$corStruct, unconstrained = FALSE), "\n\n")

summary(fit_step4)


# ═══════════════════════════════════════════════════════════
# STEP 5 — Test interazioni Bottom × gestione
#
# Scopo: verificare se il gradiente verticale dei nutrienti
# cambia con il tipo di gestione. Es: l'irrigazione modifica
# la velocità con cui SOC decade in profondità?
#
# 4 interazioni (Bottom × Irrigate, Fertilised, OnFarm,
# N_Natural) × 3 risposte = 12 parametri aggiuntivi.
# Test con LRT (richiede stima ML, non REML).
# ═══════════════════════════════════════════════════════════

fit_step5_no <- lme(
  y ~ response:(Texture1 + Texture2 + Irrigate + Fertilised +
                  Bottom + BulkDensity + N_Natural + OnFarm + PH) - 1 + response,
  random = ~ 0 + response | Field,
  correlation = corExp(value = 12, form = ~ coord_bottom + coord_resp | Field,
                       fixed = FALSE),
  weights = varIdent(form = ~ 1 | response),
  data = crop_long,
  method = "ML",
  control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim")
)

fit_step5_inter <- lme(
  y ~ response:(Texture1 + Texture2 + Irrigate + Fertilised +
                  Bottom + BulkDensity + N_Natural + OnFarm + PH +
                  Bottom:Irrigate + Bottom:Fertilised +
                  Bottom:OnFarm + Bottom:N_Natural) - 1 + response,
  random = ~ 0 + response | Field,
  correlation = corExp(value = 12, form = ~ coord_bottom + coord_resp | Field,
                       fixed = FALSE),
  weights = varIdent(form = ~ 1 | response),
  data = crop_long,
  method = "ML",
  control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim")
)

cat("══════ STEP 5: TEST INTERAZIONI Bottom × gestione ══════\n\n")
cat("--- LRT globale (12 parametri aggiuntivi) ---\n")
print(anova(fit_step5_no, fit_step5_inter))

cat("\n--- Coefficienti delle interazioni ---\n")
tt <- summary(fit_step5_inter)$tTable
print(round(tt[grep("Bottom:", rownames(tt)), c("Value", "Std.Error", "t-value", "p-value")], 4))


# Test: converge anche con nlminb se gli dai più iterazioni?
fit_test <- lme(
  y ~ response:(Texture1 + Texture2 + Irrigate + Fertilised +
                  Bottom + BulkDensity + N_Natural + OnFarm + PH) - 1 + response,
  random = ~ 0 + response | Field,
  correlation = corExp(value = 12, form = ~ coord_bottom + coord_resp | Field,
                       fixed = FALSE),
  weights = varIdent(form = ~ 1 | response),
  data = crop_long,
  method = "REML",
  control = lmeControl(maxIter = 200, msMaxIter = 200)  # nlminb con più iterazioni
)
# errore convergenza


# Test stabilità: optim con valori iniziali diversi
fit_test2 <- lme(
  y ~ response:(Texture1 + Texture2 + Irrigate + Fertilised +
                  Bottom + BulkDensity + N_Natural + OnFarm + PH) - 1 + response,
  random = ~ 0 + response | Field,
  correlation = corExp(value = 5, form = ~ coord_bottom + coord_resp | Field,
                       fixed = FALSE),
  weights = varIdent(form = ~ 1 | response),
  data = crop_long,
  method = "REML",
  control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim")
)

cat("Range con value iniziale 12:", coef(fit_step4$modelStruct$corStruct, unconstrained = FALSE), "\n")
cat("Range con value iniziale 5: ", coef(fit_test2$modelStruct$corStruct, unconstrained = FALSE), "\n")


# ═══════════════════════════════════════════════════════════
# STEP 6 — Analisi di confondimento: gestione vs proprietà fisiche
#
# Scopo: verificare se l'effetto apparente della gestione
# è confuso con le proprietà fisiche del suolo (Texture,
# BulkDensity). I siti on-farm potrebbero avere più SOC
# non perché la gestione è migliore, ma perché occupano
# suoli con tessitura e densità intrinsecamente favorevoli.
#
# Approccio: confrontare i coefficienti di gestione nel
# modello completo (con covariate fisiche) vs modello
# ridotto (senza covariate fisiche).
# ═══════════════════════════════════════════════════════════

# Modello SENZA covariate fisiche (Texture1, Texture2, BulkDensity)
fit_step6_no_phys <- lme(
  y ~ response:(Irrigate + Fertilised + Bottom +
                  N_Natural + OnFarm + PH) - 1 + response,
  random = ~ 0 + response | Field,
  correlation = corExp(value = 12, form = ~ coord_bottom + coord_resp | Field,
                       fixed = FALSE),
  weights = varIdent(form = ~ 1 | response),
  data = crop_long,
  method = "REML",
  control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim")
)

cat("══════ STEP 6: ANALISI DI CONFONDIMENTO ══════\n\n")

# Tabella confondimento tra campi
cat("--- Proprietà fisiche per gruppo di gestione ---\n")
crop %>%
  group_by(OnFarm, Irrigate, Fertilised) %>%
  summarise(
    Texture1_mean = round(mean(Texture1), 3),
    Texture2_mean = round(mean(Texture2), 3),
    BD_mean = round(mean(BulkDensity), 3),
    n_field = n_distinct(Field),
    .groups = "drop"
  ) %>%
  print()

# Confronto coefficienti gestione
cat("\n--- Coefficienti gestione: con vs senza covariate fisiche ---\n")
tt_full <- summary(fit_step4)$tTable
tt_no   <- summary(fit_step6_no_phys)$tTable

for (resp in c("logSOC", "logN", "logP")) {
  cat(sprintf("\n%s:\n", resp))
  cat(sprintf("%-15s  %8s %8s  %8s %8s\n", "", "β_full", "p_full", "β_no", "p_no"))
  for (v in c("Irrigate1", "Fertilised1", "OnFarm1", "N_Natural1")) {
    nm <- paste0("response", resp, ":", v)
    if (nm %in% rownames(tt_full) & nm %in% rownames(tt_no)) {
      cat(sprintf("%-15s  %8.4f %8.4f  %8.4f %8.4f\n",
                  v,
                  tt_full[nm, 1], tt_full[nm, 5],
                  tt_no[nm, 1], tt_no[nm, 5]))
    }
  }
}


# ═══════════════════════════════════════════════════════════
# STEP 7 - MODELLO FINALE: senza PH, logBottom come fisso, Bottom cm per corExp
#
# Modifiche rispetto a Step 4:
#   - Rimosso PH (non significativo per nessuna risposta)
#   - logBottom standardizzato come effetto fisso (power-law)
#   - Bottom in cm come coordinata per corExp (range in cm)
# ═══════════════════════════════════════════════════════════

# Aggiungi logBottom standardizzato al dataset
crop$logBottom <- c(scale(log(crop$Bottom)))

# Ricostruisci formato long con logBottom
crop_long <- crop %>%
  mutate(obs_id = row_number()) %>%
  pivot_longer(
    cols = c(logSOC, logN, logP),
    names_to = "response",
    values_to = "y"
  ) %>%
  mutate(response = factor(response))

# Coordinate 2D: Bottom in cm per SOC, gonfiato per N e P
crop_long <- crop_long %>%
  mutate(
    coord_bottom = case_when(
      response == "logSOC" ~ Bottom,
      response == "logN"   ~ Bottom * 100,
      response == "logP"   ~ Bottom * 100
    ),
    coord_resp = case_when(
      response == "logSOC" ~ 0,
      response == "logN"   ~ 5000,
      response == "logP"   ~ 10000
    )
  )

# Modello finale
fit_finale <- lme(
  y ~ response:(Texture1 + Texture2 + Irrigate + Fertilised +
                  logBottom + BulkDensity + N_Natural + OnFarm) - 1 + response,
  random = ~ 0 + response | Field,
  correlation = corExp(value = 12, form = ~ coord_bottom + coord_resp | Field,
                       fixed = FALSE),
  weights = varIdent(form = ~ 1 | response),
  data = crop_long,
  method = "REML",
  control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim")
)

cat("══════ MODELLO FINALE ══════\n")
cat("AIC Step 4 (con PH, Bottom lineare):", AIC(fit_step4), "\n")
cat("AIC finale (senza PH, logBottom):   ", AIC(fit_finale), "\n")
cat("Differenza:", AIC(fit_step4) - AIC(fit_finale), "\n\n")

cat("Range corExp (cm):\n")
cat(coef(fit_finale$modelStruct$corStruct, unconstrained = FALSE), "\n\n")

summary(fit_finale)


# ═══════════════════════════════════════════════════════════
# TEST: modello finale senza N_Natural
# ═══════════════════════════════════════════════════════════

fit_finale_noNat <- lme(
  y ~ response:(Texture2 + Irrigate + Fertilised +
                  logBottom + BulkDensity + OnFarm) - 1 + response,
  random = ~ 0 + response | Field,
  correlation = corExp(value = 12, form = ~ coord_bottom + coord_resp | Field,
                       fixed = FALSE),
  weights = varIdent(form = ~ 1 | response),
  data = crop_long,
  method = "REML",
  control = lmeControl(maxIter = 200, msMaxIter = 200, opt = "optim")
)

cat("══════ CONFRONTO ══════\n")
cat("AIC con N_Natural:   ", AIC(fit_finale), "\n")
cat("AIC senza N_Natural: ", AIC(fit_finale_noNat), "\n")
cat("Differenza:", AIC(fit_finale) - AIC(fit_finale_noNat), "\n\n")

# Confronto coefficienti
tt_con  <- summary(fit_finale)$tTable
tt_senza <- summary(fit_finale_noNat)$tTable

cat("--- Coefficienti gestione: con vs senza N_Natural ---\n")
gestione <- grep("Irrigate|Fertilised|OnFarm", rownames(tt_senza), value = TRUE)

cat(sprintf("%-35s  %8s %8s  %8s %8s\n", "", "β_con", "p_con", "β_senza", "p_senza"))
for (v in gestione) {
  v_con <- v
  if (v_con %in% rownames(tt_con)) {
    cat(sprintf("%-35s  %8.4f %8.4f  %8.4f %8.4f\n",
                v,
                tt_con[v_con, 1], tt_con[v_con, 5],
                tt_senza[v, 1], tt_senza[v, 5]))
  }
}

summary(fit_finale_noNat)
# questo è il modello finale da tenere


# ==============================================================================
# COMMENTO E INTERPRETAZIONE DEI COEFFICIENTI (MODELLO FINALE MULTIVARIATO REML)
# ==============================================================================
# Il modello analizza l'effetto del management e delle proprietà fisiche del 
# suolo su tre risposte simultanee in scala logaritmica: logN, logP e logSOC.
# Poiché le risposte sono logaritmiche, la variazione percentuale esatta della
# variabile originaria per un incremento unitario del predittore è pari a:
# %_Variazione = (exp(Beta) - 1) * 100
# ------------------------------------------------------------------------------

# 1. IMPATTO DEL MANAGEMENT (GESTIONE DEL SUOLO)
# ------------------------------------------------------------------------------
# A. IRRIGAZIONE (Irrigate1):
# - Esibisce un comportamento fortemente asimmetrico e specializzato.
# - Ha un fortissimo impatto positivo e altamente significativo sul Fosforo 
#   (logP: Beta = +0.4440, p = 0.0049). L'irrigazione incrementa il Fosforo 
#   disponibile di circa il 56% (exp(0.444)-1), stimolando la solubilizzazione.
# - Ha un impatto positivo e marginalmente significativo sulla Sostanza Organica
#   (logSOC: Beta = +0.2608, p = 0.0614, ~ +30%), indicando che la disponibilità
#   idrica costante favorisce l'accumulo di biomassa e residui nel tempo.
# - L'effetto sull'Azoto (logN) tende al negativo ma NON è significativo (p = 0.25).

# B. FERTILIZZAZIONE (Fertilised1):
# - Mostra un effetto negativo e marginalmente significativo sul Fosforo 
#   (logP: Beta = -0.3919, p = 0.0560). Questo calo di circa il 32% (exp(-0.39)-1)
#   può essere legato a un effetto "diluizione": la fertilizzazione spinge la 
#   crescita vegetativa che consuma prontamente il fosforo libero nel suolo.
# - Gli effetti su Azoto e SOC sono del tutto non significativi (p > 0.37),
#   mostrando che la fertilizzazione chimica standard da sola non sposta i livelli
#   basali di lungo periodo di queste due componenti.

# C. CONDUZIONE AZIENDALE (OnFarm1 vs Stazione Sperimentale):
# - Nei campi reali gestiti direttamente dai contadini (OnFarm), i livelli di
#   nutrienti chimici sono nettamente e significativamente più alti:
#   Azoto (logN: Beta = +0.3397, p = 0.0207, ovvero ~ +40%) e soprattutto 
#   Fosforo (logP: Beta = +0.7184, p = 0.0058, ovvero ~ +105%).
# - La gestione OnFarm non mostra invece differenze significative sulla 
#   struttura profonda della Sostanza Organica (logSOC: p = 0.6511).


# 2. PROPRIETÀ STRUTTURALI E FISICHE DEL SUOLO (FATTORI DI CONTROLLO CRUCIALI)
# ------------------------------------------------------------------------------
# D. TESSITURA (Texture2 - Seconda Componente Principale Ortogonale):
# - Ha un fortissimo impatto negativo e altamente significativo su Azoto 
#   (Beta = -0.2435, p = 0.0000) e SOC (Beta = -0.4382, p = 0.0000). Spostarsi
#   lungo l'asse positivo di questa componente riduce drasticamente la capacità 
#   strutturale del terreno di trattenere matrice organica e azoto.

# E. PROFONDITÀ (logBottom):
# - Fortissimo effetto negativo su Azoto (Beta = -0.2135, p = 0.0000) e SOC 
#   (Beta = -0.3485, p = 0.0000). La fertilità bio-organica crolla verticalmente
#   man mano che si scende nel profilo del suolo.

# F. COMPATTAZIONE (BulkDensity - Densità Apparente):
# - È un fattore limitante severo e altamente significativo.
# - All'aumentare della compattazione, l'Azoto cala vistosamente (Beta = -0.1715, 
#   p = 0.0013) così come la Sostanza Organica (Beta = -0.2438, p = 0.0042). Un 
#   suolo troppo denso ostacola l'aerazione e l'attività microbica di umidificazione.


# ==============================================================================
# RACCOMANDAZIONI AGRONOMICHE STRATEGICHE PER IL CONTADINO (MASSIMIZZAZIONE)
# ==============================================================================
# 1) IRRIGAZIONE CONTROLLATA/A GOCCIA: Sfruttare l'irrigazione per sbloccare il 
#    Fosforo (+56%) e accumulare SOC (+30%), ma evitare turni irrigui eccessivi 
#    a pioggia per scongiurare il rischio latente di dilavamento dell'Azoto.
# 2) GESTIONE DELLA COMPATTAZIONE: Poiché la Bulk Density deprime significativamente
#    N e SOC, è tassativo ridurre il calpestio da macchinari pesanti e adottare 
#    pratiche di minima lavorazione (minimum tillage) per mantenere il suolo areato.
# 3) RIVISIONE APPORTI CHIMICI: Integrare le fertilizzazioni standard (che mostrano
#    un effetto di diluizione/blocco sul Fosforo) con ammendanti organici (compost,
#    sovesci) capaci di contrastare il deficit strutturale dei suoli svantaggiati (Texture2).
# 4) PROTEZIONE DELLO STRATO SUPERFICIALE: Poiché N e SOC si azzerano in profondità 
#    (logBottom), lo strato 0-20 cm è il più prezioso e va protetto da erosione e 
#    dilavamento con tecniche di pacciamatura o mantenimento dei residui colturali.
# ==============================================================================
