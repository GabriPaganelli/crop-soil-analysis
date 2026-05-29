# Note di progetto — crop-soil-analysis

## Contesto

Studio dei profili verticali di sostanza organica del suolo (SOC), azoto totale (N) e fosforo
totale (P) in 40 campi agricoli. Per ogni campo sono disponibili misurazioni a più profondità
(fino a 6), per un totale di circa 220 osservazioni.

L'obiettivo principale è capire se e come la struttura verticale del suolo (la velocità con
cui i nutrienti decadono con la profondità) dipende dal livello medio del campo.

---

## Struttura dei dati

### Gerarchia

- **OwnId**: proprietario (anonimizzato)
- **Field**: campo — unità di analisi principale (J = 40)
- **Plot**: sub-campionamento all'interno del campo (aggregato a livello field per l'analisi)
- **Bottom**: profondità inferiore dell'intervallo di campionamento (cm)
  - Valori tipici: 10, 20, 30, 60, 90, 120 cm, ma non tutti i campi hanno tutte le profondità

### Variabili risposta

| Variabile | Unità | Trasformazione |
|-----------|-------|----------------|
| PercSOC | % | log(PercSOC) → logSOC |
| PercTotNitro | % | log(PercTotNitro) → logN |
| PercTotPhos | % | log(PercTotPhos) → logP |

La trasformazione log è giustificata empiricamente e porta a un profilo verticale
lineare in log-log (legge di potenza: nutriente ∝ profondità^β).

### Covariati within-field (variano per osservazione)

| Variabile | Descrizione | Trattamento |
|-----------|-------------|-------------|
| logBottom | log(profondità cm) | standardizzata (scala) |
| Texture1 | ILR1 della tessitura | calcolata in preprocessing; standardizzata |
| Texture2 | ILR2 della tessitura | calcolata in preprocessing; standardizzata |
| BulkDensity | densità apparente | standardizzata |

**ILR (Isometric Log-Ratio)**: trasformazione composizionale per trattare Sand/Silt/Clay
come variabili reali, eliminando il vincolo a somma 100%:
- Texture1 = √(2/3) · log(Sand / √(Silt·Clay))
- Texture2 = √(1/2) · log(Silt / Clay)

**Nota**: crop.csv non contiene Texture1/Texture2. Vengono calcolate da 01_preprocessing.R
partendo dalle percentuali grezze. I dati grezzi della tessitura (PercSand, PercSilt, PercClay)
sono vincolati: PercSand + PercSilt + PercClay ≈ 100% — l'ILR risolve questo problema.

### Covariati between-field (uno per campo, binari 0/1)

| Variabile | Descrizione |
|-----------|-------------|
| OnFarm | compostaggio/fertilizzazione organica in loco |
| Irrigate | campo irrigato |
| Fertilised | fertilizzazione minerale |
| N_Natural | apporto di azoto da fonti naturali |

### Anomalie e osservazioni sui dati

- **Profondità asimmetriche**: non tutti i campi hanno misurazioni a tutte le profondità.
  Questo rende l'analisi per-campo con lme4 rumorosa (SNR slope < 1 per molti campi).
- **Multicollinearità nelle X_W**: Texture1, Texture2, BulkDensity sono correlate.
  Irrilevante per il modello bayesiano con prior regolarizzanti, ma da tenere a mente
  nell'interpretazione dei gamma_r.
- **BulkDensity**: presenta outlier; standardizzazione applicata.
- **Valori mancanti**: sporadici su alcune profondità; gestiti implicitamente dal modello
  multilevel (ogni osservazione ha il suo logBottom).

---

## Cronologia delle scelte analitiche

### 1. Preprocessing e scelta della forma funzionale

La forma log-log per il profilo verticale (logSOC ~ logBottom) è stata scelta
confrontando quattro alternative con AIC (lmer, ML):

| Forma | AIC |
|-------|-----|
| log(Bottom) ← **scelta** | più basso |
| quadratica in Bottom | intermedio |
| lineare in Bottom | alto |
| factor(Bottom) | simile a log ma non parsimoniosa |

La forma log-log è equivalente a un decadimento a legge di potenza:
SOC ∝ Bottom^β, dove β = slope nella scala log-log.

### 2. Struttura random: ICC e necessità del random intercept

ICC intra-field per logSOC calcolato su due modelli lmer:
- **Modello nullo** (solo random intercept): ICC ≈ 0.72
- **Dopo correzione per Bottom**: ICC rimane alto

La variabilità tra campi domina → il random intercept è necessario e ben identificato.

### 3. Segnale intercetta-slope

Con lmer (logSOC ~ Bottom_c + (Bottom_c | Field), REML=TRUE):
- **Corr(BLUP intercetta, BLUP slope) ≈ +0.93** (con Bottom centrato)
- I campi con SOC alto alla profondità media decadono più lentamente

Questo segnale ha motivato il modello M-SP con slope proporzionale.

### 4. Da modelli fattoriali a modello bayesiano

Inizialmente esplorate strutture MIMIC a 2 fattori (non documentate nel codice finale).
Abbandonate per:
- Difficoltà di identificazione con 3 risposte correlate
- La struttura proporzionale del decadimento verticale (η_r) cattura il segnale chiave
  in modo più diretto e interpretabile

### 5. Esplorazione spaziale e abbandono del GP

Analisi spaziale in 07_eda_spatial_gp.R:
- **SOC**: debole autocorrelazione spaziale (Moran's I positivo, significativo)
- **N, P**: nessuna autocorrelazione spaziale rilevante
- Il GP spaziale (modelli M-GPS e M-GP) non migliora il LOO rispetto a M-SP

Conclusione: la struttura spaziale non è necessaria; l'eterogeneità tra campi è catturata
dal random intercept proporzionale.

---

## Risultati principali

### Parametro η (velocità di decadimento)

| Risposta | η (mediana) | CI90 | Interpretazione |
|----------|-------------|------|-----------------|
| SOC | +0.209 | [0.125, 0.302] | effetto forte, ben identificato |
| N | −0.051 | ≈ 0 (span 0) | profilo verticale uniforme tra campi |
| P | +0.096 | [0.019, 0.179] | effetto debole ma presente |

**I campi con SOC alto alla profondità media decadono meno con la profondità.**
Possibile meccanismo: accumulo di materia organica stabile negli orizzonti profondi
nei suoli più ricchi; oppure selezione di pratiche agronomiche che favoriscono
sia i livelli alti che la distribuzione uniforme.

### Confronto modelli (LOO-CV)

Ranking: **M-SP > M-RI > M-GPS > M-GP**

M-SP (slope proporzionale) è il migliore pur avendo meno parametri del GP.
La semplicità del modello proporzionale è una virtù, non un compromesso.

### Effetti del management (β_r)

In generale deboli dopo aver controllato per i covariati within-field.
Il management (irrigazione, fertilizzazione) non spiega grandi differenze nei profili
una volta che l'eterogeneità between-field è catturata dal random intercept.

---

## Decisioni prese (non ovvie)

1. **Tre risposte in parallelo** in un unico modello Stan anziché tre modelli separati:
   consente di condividere la struttura del codice Stan e di confrontare i parametri
   direttamente, ma ogni risposta ha i propri parametri (non pooling cross-risposta).

2. **Non-Centered Parameterization (NCP)**: `z_nu_r[j] ~ N(0,1)`,
   poi `nu_r[j] = z_nu_r[j] * (psi_r + eta_r * logBottom)`.
   Necessaria per evitare divergenze NUTS con prior deboli sui random effects.

3. **Intercetta globale α_r**: aggiunta dopo che senza di essa beta_r (management)
   esplodeva a ±2. Il random intercept ha media zero per costruzione (NCP);
   senza α_r, i beta cercano di assorbire la media globale della risposta.

4. **Texture come ILR invece di PercSand/PercSilt/PercClay separati**:
   le percentuali di tessitura sono compositional data (somma = 100%) e non possono
   essere usate come regressori indipendenti senza violare i vincoli geometrici.
   L'ILR proietta sul simplesso in R² senza vincoli.

5. **logBottom standardizzato in X_W**: diversamente dagli script EDA (dove logBottom
   è centrato ma non scalato), nei modelli Stan logBottom è standardizzato (z-score)
   per rendere gamma_r[1] comparabile in scala con gli altri gamma_r.
   Questo cambia la scala di η_r rispetto all'analisi OLS per campo.

6. **projpred sui soli effetti fissi**: la struttura proporzionale (η_r) è già nel
   reference model; i submodelli di proiezione usano (1|Field) come approssimazione
   accettabile per selezionare gamma_r e beta_r.
