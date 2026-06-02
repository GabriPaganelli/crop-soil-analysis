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
- Texture1 = √(1/2) · log(Clay / Silt)          ← contrasto argilla vs limo
- Texture2 = √(2/3) · log(√(Clay·Silt) / Sand)  ← gradiente finezza vs sabbia

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

Analisi spaziale in 04_eda_spaziale.R:
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

### Robustezza: M-SP-MV (script 10)

Estensione multivariata di M-SP: i residui (eps_SOC, eps_N, eps_P) alla stessa osservazione
sono modellati con MVNormal e matrice di correlazione Omega_3×3 (3 parametri rho aggiuntivi).
**Stan file**: `stan/m_msp_mv.stan` | **Fit**: `stan/fit_msp_mv.rds`

| Correlazione | Mediana | CI90% | Interpretazione |
|--------------|---------|-------|-----------------|
| rho_SOC_N | +0.067 | [−0.063, +0.197] | copre zero |
| rho_SOC_P | ≈ 0 | — | copre zero |
| rho_N_P | ≈ 0 | — | copre zero |

**eta_r invariante**: Δmediana < 0.005 per tutte e tre le risposte → la struttura proporzionale
non dipende dall'ipotesi di indipendenza dei residui.

**LOO**: M-SP vince su M-SP-MV di ΔELPD = +4.7 (SE=1.4, |z|=3.4) → aggiungere
correlazioni quasi-zero penalizza la predizione.

**Conclusione**: i tre nutrienti sono condizionalmente indipendenti dati i predittori di M-SP.

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

---

## Risultati script 10 — Confronto modelli A / B / M-SP (LOO-CV)

### LOO — tutti e tre statisticamente equivalenti

| Modello | elpd_loo | p_loo | elpd_diff | se_diff | ratio |
|---------|----------|-------|-----------|---------|-------|
| B (resp-spec) | −504.3 | 110.4 | 0.0 | — | — |
| A (ridotto) | −504.7 | 119.6 | −0.4 | 4.6 | 0.09 |
| M-SP | −506.7 | 126.7 | −2.3 | 7.2 | 0.32 |

**Soglia standard**: |elpd_diff / se_diff| ≥ 2. Tutti i rapporti ≪ 2 → nessun modello è
predittivamente superiore. B vince formalmente solo per parsimonia (meno parametri effettivi).

### Confronto parametri chiave (mediana [CI90%])

**SOC — risultato solido**

| Parametro | M-SP | A (ridotto) | B (resp-spec) | C (projpred) |
|-----------|------|-------------|---------------|--------------|
| eta_SOC | 0.209 [0.125, 0.302] | 0.214 [0.132, 0.306] | 0.207 [0.127, 0.300] | — |
| gamma logBottom | −0.378 | −0.383 | −0.432 | −0.447 |
| gamma Texture2 | −0.529 | −0.494 | −0.588 | −0.477 |
| gamma BulkDensity | −0.223 | −0.189 | — | — |

eta_SOC è **invariante** su tutti i modelli e sempre lontano da zero. Conclusione robusta.

**N — confermato piatto tra campi**

| Parametro | M-SP | A (ridotto) | B (resp-spec) |
|-----------|------|-------------|---------------|
| eta_N | −0.051 [−0.124, 0.021] | −0.023 [−0.096, 0.050] | — (no eta in B) |

Segno variabile, CI sempre spanno zero. Nessuna evidenza di slope proporzionale per N.

**P — evidenza moderata, dipendente dal modello**

| Parametro | M-SP | A (ridotto) | B (resp-spec) |
|-----------|------|-------------|---------------|
| eta_P | 0.096 [0.019, 0.179] | 0.105 [0.031, 0.185] | 0.054 [−0.016, 0.126] |
| psi_P | 0.392 | 0.474 | 0.568 |
| gamma_P_logB | — | — | −0.121 [−0.185, −0.056] |

In M-SP e A eta_P è chiaramente positivo. In B (che separa l'effetto fisso di logBottom),
eta_P si riduce e il CI include zero. Il motivo: in M-SP/A la componente fissa e casuale
della profondità per P non sono distinguibili, gonfiando eta_P. In B sono separate → il segnale
proporzionale per P è più debole di quello per SOC. **Conclusione: moderata evidenza per P.**

**Nota su alpha_r**: in M-SP, alpha_SOC = −0.240 con CI larghissimo [−0.527, 0.049]; in A/B
si stringe (−0.116/−0.126). Normale: in M-SP l'intercetta globale è aggiustata dai beta_r
(management); rimuovendo X_B, alpha_r assorbe quegli effetti e diventa più precisa.

---

## Risultati script 11 — Sensitivity analysis (Pareto k)

### Distribuzione dei k

| Categoria | N osservazioni |
|-----------|---------------|
| k < 0.5 (buono) | 189 |
| 0.5 ≤ k < 0.7 (moderato) | 21 |
| 0.7 ≤ k < 1.0 (alto) | 9 |
| k ≥ 1.0 (molto alto) | 1 |

**10 osservazioni rimosse** (k ≥ 0.7). Obs più influente: obs 146, Field 26, Bottom=20cm,
k=1.253. Campi con più influenti: Field 23 (3 obs), Field 26 (2 obs), Field 33 (2 obs).
Concentrazione sulle profondità intermedie-profonde (40–80cm), dove i profili sono più variabili.

### Diagnostica MCMC post-rimozione
Divergenze: 0 | Max treedepth: 9 | Rhat>1.05: 0 | ESS<400: 0 ✅

### Confronto parametri full vs no-influential

**Flag ATTENZIONE** (mediana no-infl fuori dal CI90% full):

| Parametro | Full | No-influential | Variazione |
|-----------|------|----------------|-----------|
| sigma_SOC | 0.531 [0.486, 0.582] | 0.440 [0.402, 0.484] | −17% |
| sigma_N | 0.345 [0.316, 0.378] | 0.266 [0.243, 0.292] | −23% |

Le obs influenti sono outlier di residuo (valori anomali di SOC e N), non leve nella struttura.
La riduzione di sigma è attesa e non preoccupante per l'inferenza sui parametri di interesse.

**Parametri eta_r — tutti OK o moderati:**

| Parametro | Full [CI90%] | No-infl [CI90%] | Flag |
|-----------|-------------|-----------------|------|
| eta_SOC | 0.209 [0.125, 0.302] | 0.160 [0.089, 0.239] | OK |
| eta_N | −0.051 [−0.124, 0.021] | 0.016 [−0.033, 0.061] | nota |
| eta_P | 0.096 [0.019, 0.179] | 0.119 [0.039, 0.206] | OK |

- **eta_SOC**: si riduce leggermente (0.209 → 0.160) ma CI rimane chiaramente positivo → robusto.
- **eta_N**: cambia segno (−0.051 → +0.016), entrambi vicini a zero → conclusione invariata (N piatto).
- **eta_P**: leggermente più forte senza influenti (0.096 → 0.119) → evidenza si consolida.

**Tutti i gamma** (Texture2, BulkDensity, logBottom) sono stabili (flag OK) tra full e no-infl.

### Sintesi robustezza

| Risultato | Robustezza |
|-----------|-----------|
| eta_SOC > 0 (campi con SOC alto decadono meno) | ✅ Molto robusto |
| eta_N ≈ 0 (N uniforme tra campi) | ✅ Molto robusto |
| eta_P > 0 (evidenza moderata) | ⚠️ Moderato — dipende dal modello, si consolida senza influenti |
| gamma Texture2 < 0 (più argilla/limo → più SOC/N/P) | ✅ Robusto |
| gamma BulkDensity < 0 (densità maggiore → meno nutrienti) | ✅ Robusto |
| Effetti management (beta_r) deboli | ✅ Confermato da LOO su A e B |

---

## Strategia per il report (Opzione C)

- **Modello principale**: M-SP (interpretazione biologica piena, beta management inclusi)
- **Model comparison**: A e B come evidenza di parsimonia (LOO equivalente con meno parametri)
- **Sensitivity**: script 11 come sezione di robustezza (le conclusioni chiave reggono)
- **Narrativa**: dalla domanda scientifica (il decadimento dipende dal livello medio del campo?)
  ai dati, al modello, ai risultati, alla loro robustezza
