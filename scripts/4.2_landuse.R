# ════════════════════════════════════════════════════════════════════════════
# Script 4 — Modello Mixed Multivariato (LMM-MV) LANDUSE RELEVELED
# Analisi congiunta di log(SOC), log(N), log(P) con struttura di covarianza
# ════════════════════════════════════════════════════════════════════════════
#
# RAZIONALE
# I modelli univariati (script 2) assumono implicitamente che i residui di
# SOC, N e P siano indipendenti — assunzione biologicamente insostenibile
# dato il ciclo stoichiometrico C:N:P nei suoli. Il modello multivariato:
#
#   (a) stima la matrice G (between-Field):
#       un sito con elevato random intercept per SOC ha anche elevato
#       random intercept per N e P? → "fattore latente di fertilità del sito"
#
#   (b) stima la matrice R (residui entro-Field):
#       dipendenza stoichiometrica C:N:P che persiste dopo aver rimosso
#       gestione, profondità e covariate pedologiche
#
#   (c) stima gli effetti fissi di gestione congiuntamente sulle 3 risposte,
#       migliorando l'efficienza (errori standard più piccoli) quando
#       le risposte sono positivamente correlate
#
# STRATEGIA DI ANALISI
#   1. Feasibility check — quanti parametri vs quanti gruppi
#   2. Motivazione empirica — correlazione dei residui univariati
#   3. MCMCglmm (principale): MCMC con prior IW piatte ≈ REML multivariato
#   4. sommer (verifica frequentista pura): REML multivariato via Newton-Raphson
#   5. Confronto struttura random: intercept only vs intercept + slope
#
# RIFERIMENTI
#   Hadfield JD (2010) MCMCglmm Course Notes — J Stat Softw 33(2)
#   Covarrubias-Pazaran G (2016) sommer: solving mixed model equations — Crop Sci
#   Henderson & Quaas (1976) Multiple trait evaluation — J Dairy Sci 59:1102
#   Sorensen & Gianola (2002) Likelihood, Bayesian and MCMC Methods — Springer
#   Piepho et al. (2008) Mixed models for plant breeding — J Agronomy Crop Sci

library(MCMCglmm)
library(sommer)
library(lme4)
library(lmerTest)
library(tidyverse)
library(here)
library(patchwork)
library(dplyr)

# MCMCglmm carica MASS, che maschera dplyr::select, dplyr::filter, dplyr::rename.
# Riportiamo le funzioni dplyr in primo piano nel search path.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

dati <- readRDS(here("data", "dati.rds"))
crop <- readRDS(here("data", "crop.rds"))


dati$Bottom_c  <- dati$Bottom - mean(dati$Bottom)
dati$Bottom_c2 <- dati$Bottom_c^2   # pre-calcolo: evita I() nelle formule MCMCglmm

dati <- dati |>
  mutate(
    logSOC = log(PercSOC),
    logN   = log(PercTotNitro),
    logP   = log(PercTotPhos),
    Field  = factor(Field),
    Landuse = factor(crop$Landuse)
  ) |>
  mutate(
    Landuse = relevel(factor(Landuse), ref = "4")
  ) 

target_vars <- c("logSOC", "logN", "logP")
n_risposte  <- length(target_vars)
n_fields    <- length(unique(dati$Field))
n_obs       <- nrow(dati)



# ════════════════════════════════════════════════════════════════════════════
# 2. MODELLO MULTIVARIATO — MCMCglmm (principale)
# ════════════════════════════════════════════════════════════════════════════
# Approccio: formato long (una riga per risposta × osservazione)
# La colonna "trait" identifica la risposta. Gli effetti fissi sono specificati
# come interazione trait × predittore → stime separate per ciascuna risposta.
#
# Struttura random:
#   G = us(trait):Field  → matrice 3×3 non strutturata tra-Field
#                           (3 varianze + 3 covarianze = 6 parametri liberi)
#
# Struttura residua:
#   R = us(trait):units  → matrice 3×3 non strutturata entro-Field
#                           (dipendenza stoichiometrica C:N:P)
#
# Prior Inverse-Wishart:
#   V = diag * 0.1, nu = n_risposte + 2 = 5
#   Con nu piccolo e V piccola → prior quasi non informativa per le correlazioni
#   E[V | IW(V, nu)] = V / (nu - p - 1) → tende a 0 per nu → p+1
#   Questa scelta approssima la soluzione REML del modello LMM univariato.
#   Riferimento: Gelman & Hill (2007) cap. 13 — prior per matrici di covarianza

dati_long <- dati |>
  dplyr::select(Field, Landuse,
                Bottom_c, Bottom_c2, PH, BulkDensity, Texture1, Texture2,
                logSOC, logN, logP) |>
  pivot_longer(cols = all_of(target_vars),
               names_to  = "componente",       # "trait" è parola riservata in MCMCglmm
               values_to = "y") |>
  mutate(componente = factor(componente, levels = target_vars)) |>
  as.data.frame()   # MCMCglmm non accetta tibble

prior_int <- list(
  G = list(G1 = list(V = diag(n_risposte) * 0.1, nu = n_risposte + 2)),
  R = list(       V = diag(n_risposte) * 0.1, nu = n_risposte + 2)
)

formula_mcmc <- y ~ componente - 1 +
  componente:Landuse +
  componente:Bottom_c + componente:Bottom_c2 +
  componente:PH + componente:BulkDensity + componente:Texture1 + componente:Texture2

cat("Avvio MCMCglmm — Random intercept multivariato (G 3×3) ...\n")
cat("  100 000 iterazioni, burnin 20 000, thin 80 → 1 000 campioni posteriori\n\n")

set.seed(42)
m_mv_int <- MCMCglmm(
  formula_mcmc,
  random   = ~ us(componente):Field,
  rcov     = ~ us(componente):units,
  family   = "gaussian",   # formato long: 1 colonna y → 1 famiglia
  prior    = prior_int,
  data     = dati_long,
  nitt     = 100000, burnin = 20000, thin = 80,
  verbose  = FALSE
)

cat("Completato.\n\n")
summary(m_mv_int)


# ── 2a. Diagnostica MCMC ─────────────────────────────────────────────────
# Trace plots: le catene devono sembrare "rumore bianco" (mixing rapido)
# Effective sample size (ESS): almeno 200 per ciascun parametro
# Gelman-Rubin: < 1.1 (richiede 2+ catene, qui usiamo ESS come proxy)

par(mfrow = c(2, 3))
plot(m_mv_int$VCV[, 1:6], auto.layout = FALSE, main = "VCV — catene MCMC (G matrix)")
par(mfrow = c(1, 1))

cat("── Effective Sample Size (G matrix + R matrix) ──\n")
ess_G <- effectiveSize(m_mv_int$VCV[, grep("Field", colnames(m_mv_int$VCV))])
ess_R <- effectiveSize(m_mv_int$VCV[, grep("units", colnames(m_mv_int$VCV))])
cat("G matrix:\n"); print(round(ess_G, 0))
cat("R matrix:\n"); print(round(ess_R, 0))
cat("\n(ESS < 200 → aumentare nitt o ridurre thin)\n\n")

cat("── Autocorrelazione MCMC (lag = thin) ──\n")
ac_G <- autocorr.diag(m_mv_int$VCV[, 1:6])
print(round(ac_G, 3))
cat("\n(Autocorrelazione al lag 1 < 0.1 indica buon mixing)\n\n")


# ── 2b. Estrazione matrici G e R dalla catena MCMC ───────────────────────

# Nomi colonne VCV: "componentelogSOC:componentelogSOC.Field" = Var(logSOC | Field)
#                  "componentelogSOC:componentelogN.Field"   = Cov(logSOC, logN | Field)

vcv_G  <- m_mv_int$VCV[, grep("Field", colnames(m_mv_int$VCV))]
vcv_R  <- m_mv_int$VCV[, grep("units", colnames(m_mv_int$VCV))]
G_post <- colMeans(vcv_G)
R_post <- colMeans(vcv_R)

# Ricostruisce le matrici 3×3 (mediane posteriori)
G_mat <- matrix(G_post, n_risposte, n_risposte,
                dimnames = list(c("SOC","N","P"), c("SOC","N","P")))
R_mat <- matrix(R_post, n_risposte, n_risposte,
                dimnames = list(c("SOC","N","P"), c("SOC","N","P")))

# Converti in matrici di correlazione
cov2cor_safe <- function(M) {
  sd_vec <- sqrt(diag(M))
  M / outer(sd_vec, sd_vec)
}
G_cor <- cov2cor_safe(G_mat)
R_cor <- cov2cor_safe(R_mat)

cat("── Matrice G: correlazioni tra-Field (mediana posteriore) ──\n")
print(round(G_cor, 3))
cat("\nInterpretazione diagonale (ICC per risposta):\n")
for (i in seq_len(n_risposte)) {
  icc <- G_mat[i, i] / (G_mat[i, i] + R_mat[i, i])
  cat(sprintf("  ICC(%s) = %.3f  → %.0f%% variabilità tra Field\n",
              rownames(G_mat)[i], icc, 100 * icc))
}
cat("\nInterpretazione fuori diagonale (G_cor):\n")
cat("  Cor(SOC, N | Field)  =", round(G_cor["SOC", "N"], 3), "\n")
cat("  Cor(SOC, P | Field)  =", round(G_cor["SOC", "P"], 3), "\n")
cat("  Cor(N,   P | Field)  =", round(G_cor["N",   "P"], 3), "\n")
cat("  (Valori alti → fattore latente di fertilità del sito)\n\n")

cat("── Matrice R: correlazioni residue entro-Field (mediana posteriore) ──\n")
print(round(R_cor, 3))
cat("\nInterpretazione:\n")
cat("  Dipendenza stoichiometrica C:N:P che persiste entro ogni Field\n")
cat("  dopo aver controllato per gestione, profondità e covariate.\n\n")


# ── 2c. Intervalli credibili per le correlazioni G e R ───────────────────

# Funzione ausiliaria: calcola la catena di correlazioni dato cov e var
get_cor_chain <- function(vcv_mat, nm1, nm2, tag) {
  col_cov <- paste0("componente", nm1, ":componente", nm2, ".", tag)
  col_v1  <- paste0("componente", nm1, ":componente", nm1, ".", tag)
  col_v2  <- paste0("componente", nm2, ":componente", nm2, ".", tag)
  if (!col_cov %in% colnames(vcv_mat)) return(NULL)
  vcv_mat[, col_cov] / sqrt(vcv_mat[, col_v1] * vcv_mat[, col_v2])
}

print_cor_ic <- function(vcv_mat, tag, label) {
  cat(sprintf("── IC 95%% posteriore — correlazioni %s ──\n", label))
  pairs <- list(c("logSOC","logN"), c("logSOC","logP"), c("logN","logP"))
  for (p in pairs) {
    chain <- get_cor_chain(vcv_mat, p[1], p[2], tag)
    if (is.null(chain)) next
    ic <- quantile(chain, c(0.025, 0.25, 0.5, 0.75, 0.975))
    prob_pos <- mean(chain > 0)
    cat(sprintf("  Cor(%s, %s): median=%.3f  IC95%%=[%.3f, %.3f]  P(r>0)=%.3f\n",
                p[1], p[2], ic["50%"], ic["2.5%"], ic["97.5%"], prob_pos))
  }
  cat("\n")
}

print_cor_ic(vcv_G, "Field", "G (tra-Field)")
print_cor_ic(vcv_R, "units", "R (residua entro-Field)")


# ── 2d. Heatmap matrici G e R ────────────────────────────────────────────

heatmap_cor <- function(mat, titolo, subtitolo = NULL) {
  df <- as.data.frame(as.table(mat)) |>
    dplyr::rename(Riga = Var1, Colonna = Var2, r = Freq)
  ggplot(df, aes(x = Colonna, y = Riga, fill = r)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.2f", r)),
              size = 5, fontface = "bold") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d73027",
                         midpoint = 0, limits = c(-1, 1)) +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(hjust = 0.5)) +
    labs(title = titolo, subtitle = subtitolo,
         x = NULL, y = NULL, fill = "r")
}

p_G_heat <- heatmap_cor(
  G_cor,
  "Matrice G — Correlazioni tra-Field",
  "Variabilità tra siti dopo rimozione dei predittori"
)
p_R_heat <- heatmap_cor(
  R_cor,
  "Matrice R — Correlazioni residue entro-Field",
  "Dipendenza stoichiometrica C:N:P"
)
print(p_G_heat | p_R_heat)


# ── 2e. Effetti fissi — confronto multivariato vs univariati ─────────────

fe_names <- colnames(m_mv_int$Sol)

# 1. ESTRAZIONE DINAMICA DEI LIVELLI DI LANDUSE
# Troviamo tutti i coefficienti che contengono "Landuse"
landuse_cols <- grep("componente.*:Landuse", fe_names, value = TRUE)

# Estraiamo i suffissi dei livelli (es. se c'è "Landuse1" estraiamo "1") Rimuovendo la parte iniziale
# Questo garantisce che il codice funzioni sia che i livelli siano numeri (1,2,3..) sia che siano stringhe ("Forest","Crop"..)
landuse_levels <- unique(sub("componente.*:Landuse", "", landuse_cols))

cat("Livelli di Landuse rilevati nel modello (contrasti attivi):", paste(landuse_levels, collapse = ", "), "\n")
cat("Nota: Il livello mancante è il livello '4', usato come riferimento (baseline).\n\n")

# 2. COSTRUZIONE DEL DATAFRAME DEGLI EFFETTI FISSI
fe_df <- bind_rows(lapply(landuse_levels, function(lvl) {
  bind_rows(lapply(target_vars, function(tr) {
    
    # Costruiamo il pattern esatto generato da MCMCglmm per l'interazione
    pattern <- paste0("componente", tr, ":Landuse", lvl)
    
    # Controllo di sicurezza se il coefficiente esiste
    if (!pattern %in% fe_names) return(NULL)
    
    chain <- m_mv_int$Sol[, pattern]
    ic <- quantile(chain, c(0.025, 0.1, 0.5, 0.9, 0.975))
    
    data.frame(
      risposta   = tr,
      livello    = lvl,
      confronto  = paste0("Landuse ", lvl, " vs 4"), # Etichetta esplicita del contrasto
      lwr95      = ic["2.5%"], lwr80 = ic["10%"],
      median     = ic["50%"],
      upr80      = ic["90%"], upr95 = ic["97.5%"],
      prob_pos   = mean(chain > 0)
    )
  }))
})) |>
  mutate(
    risposta   = factor(risposta,  levels = target_vars,
                        labels = c("SOC", "Azoto", "Fosforo")),
    # Ordiniamo i livelli per fare in modo che il grafico appaia ordinato
    livello    = factor(livello, levels = sort(landuse_levels)),
    confronto  = factor(confronto, levels = paste0("Landuse ", sort(landuse_levels), " vs 4"))
  )

# 3. GRAFICO GGPLOT AGGIORNATO PER VARIABILE CATEGORIALE
ggplot(fe_df, aes(x = confronto, color = risposta)) +
  geom_linerange(aes(ymin = lwr95, ymax = upr95),
                 position = position_dodge(width = 0.6), linewidth = 0.5) +
  geom_linerange(aes(ymin = lwr80, ymax = upr80),
                 position = position_dodge(width = 0.6), linewidth = 1.5) +
  geom_point(aes(y = median),
             position = position_dodge(width = 0.6), size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title    = "Contrasti di Landuse — Modello MCMCglmm Multivariato",
    subtitle = "Variazione dei 6 livelli rispetto a Landuse 4 (baseline a 0) — Linea spessa: IC 80%, Sottile: IC 95%",
    x = NULL, y = "Effetto stimato (log scale, differenza vs ref)", color = "Risposta"
  ) +
  theme_minimal(base_size = 12)

# 4. TABELLA DELLE PROBABILITÀ POSTERIORI P(β > 0)
cat("── Probabilità posteriori P(β > 0) per i contrasti di Landuse (vs Riferimento = 4) ──\n")
fe_df |>
  mutate(
    `P(β>0)` = sprintf("%.3f", prob_pos),
    CI_95   = sprintf("[%+.3f, %+.3f]", lwr95, upr95),
    segnale = case_when(
      prob_pos > 0.95 | prob_pos < 0.05 ~ "forte",
      prob_pos > 0.80 | prob_pos < 0.20 ~ "moderato",
      TRUE                               ~ "debole"
    )
  ) |>
  select(confronto, risposta, median, CI_95, `P(β>0)`, segnale) |>
  arrange(confronto, risposta) |>
  print()


