# Note sul modello fattoriale latente per l'analisi del suolo

Dialogo strutturato su come modellare PercSOC, PercTotNitro e PercTotPhos come
espressioni di un fattore latente di fertilità, con covariate nel modello strutturale
e struttura gerarchica per Field. Le variabili di lavoro sono quelle in `data/dati.rds`.

---

## Indice

1. [Dal `factanal` al CFA](#topic-1)
2. [Il modello MIMIC](#topic-2)
3. [La struttura gerarchica](#topic-3)
4. [Multilevel MIMIC in R](#topic-4)
5. [Scelte specifiche sul dataset](#topic-5)
6. [Stima, fit e diagnostica](#topic-6)
7. [Interpretazione](#topic-7)
8. [Factor scores e diagnostica spaziale](#topic-8)
9. [Correlazione spaziale nei random effects](#topic-9)
10. [Estensione bayesiana (M4-Bayes)](#topic-10)

---

## Topic 1 — Dal `factanal` al CFA {#topic-1}

### 1.1 Cosa fa `factanal` (EFA)

In EFA lasci che siano i dati a decidere quanti fattori ci sono e come le variabili
si raggruppano su di essi. Il modello è:

$$\mathbf{y} = \boldsymbol{\Lambda}\boldsymbol{\eta} + \boldsymbol{\varepsilon}$$

dove **y** è il vettore delle variabili osservate, **Λ** è la matrice dei loadings,
**η** è il vettore dei fattori latenti e **ε** sono gli errori (varianze specifiche).
Il vincolo classico è che **η** ~ N(0, I) e **ε** ~ N(0, Ψ) con Ψ diagonale.

Il problema: con *p* variabili e *k* fattori, i parametri liberi in **Λ** sono *p×k*,
ma il modello è sottoidentificato senza ulteriori vincoli. `factanal` risolve questo
imponendo **Λ'Ψ⁻¹Λ** diagonale. Poi ruoti i loadings (varimax, oblimin…) per renderli
interpretabili — ma la rotazione non è parte del modello, è una post-elaborazione
arbitraria.

### 1.2 Cosa cambia nel CFA

Nel CFA si specifica *a priori* quale variabile carica su quale fattore. In questo caso
la specifica è immediata e teoricamente motivata:

```
PercSOC      → η (fertilità)
PercTotNitro → η (fertilità)
PercTotPhos  → η (fertilità)
```

Tutti e tre caricano su un unico fattore. Nessuna rotazione, nessuna ambiguità.

In `lavaan`:

```r
modello <- '
  fertilita =~ PercSOC + PercTotNitro + PercTotPhos
'
```

L'operatore `=~` si legge *"è misurato da"*.

### 1.3 Identificazione con 3 indicatori e 1 fattore

Il fattore latente non ha un'unità di misura naturale. Per stimarlo si fissa la scala
in uno dei due modi:

- **Opzione A — fixed loading (default `lavaan`):** fissi il loading del primo
  indicatore a 1. Il fattore prende la scala di `PercSOC`.
- **Opzione B — fixed variance:** fissi Var(η) = 1. Tutti e tre i loadings sono
  liberi. Più comoda quando si aggiungono le covariate (Topic 2).

Con 3 indicatori e 1 fattore il modello ha esattamente **0 gradi di libertà**
(è *esattamente identificato*): non è possibile valutare il fit, ma la struttura
a 1 fattore per (SOC, N, P) è teoricamente ben motivata, quindi questo è
pragmaticamente accettabile.

### 1.4 Una differenza pratica rispetto a `factanal`

`factanal` usa la matrice di correlazione (standardizza tutto internamente).
`lavaan` lavora sulla matrice di covarianza per default. Dato che tutte e tre le
risposte sono percentuali dello stesso ordine di grandezza, non è un problema critico.

---

> **Domanda:** non è prevista l'inclusione di covariate? ha senso partire senza
> covariate per esplorare, ma andranno incluse.

Le covariate entrano nel Topic 2 (MIMIC). Il CFA stand-alone serve come controllo
di plausibilità (i loadings sono tutti positivi? le varianze degli errori sono
ragionevoli?) più che come vero test, proprio perché il modello è esattamente
identificato. Il Topic 1 era necessario per costruire il vocabolario che serve per
parlare di MIMIC in modo preciso.

---

## Topic 2 — Il modello MIMIC {#topic-2}

### 2.1 La struttura

MIMIC = CFA + regressione sul fattore latente. Graficamente:

```
OnFarm ──┐
Irrigate ─┤
Fertilised┤
N_Natural ┤
Bottom ───┼──→  η (fertilità)  ──→  PercSOC      + ε₁
PH ───────┤                    ──→  PercTotNitro  + ε₂
Texture1 ─┤                    ──→  PercTotPhos   + ε₃
Texture2 ─┤
BulkDens ─┘
```

Le frecce a destra sono il **modello di misura**. Le frecce a sinistra sono il
**modello strutturale**: le covariate *causano* il fattore latente, non gli
indicatori direttamente.

In equazioni:

$$\eta = \gamma_1 \cdot \text{OnFarm} + \gamma_2 \cdot \text{Irrigate} + \ldots + \gamma_9 \cdot \text{BulkDensity} + \zeta$$

$$y_j = \lambda_j \cdot \eta + \varepsilon_j \quad (j = 1,2,3)$$

In `lavaan`:

```r
mimic <- '
  # Modello di misura
  fertilita =~ PercSOC + PercTotNitro + PercTotPhos

  # Modello strutturale
  fertilita ~ OnFarm + Irrigate + Fertilised + N_Natural +
              Bottom + PH + Texture1 + Texture2 + BulkDensity
'
fit <- sem(mimic, data = dati)
```

### 2.2 Perché non fare due step separati

L'alternativa sarebbe: (1) stimare il CFA, ricavare i factor scores η̂, (2) regredire
η̂ sulle covariate con `lm`. Il problema è che i factor scores sono **stime con
errore**, e questo errore non viene propagato nel secondo step: le varianze degli
errori standard nel secondo step sono sottostimate, e i p-value risultano troppo
ottimistici. Nel MIMIC tutto viene stimato in un unico passo con ML.

### 2.3 Parametri stimati

| Blocco | Parametri |
|---|---|
| Loadings λ (fixing λ₁=1) | 2 |
| Varianze errori ε₁, ε₂, ε₃ | 3 |
| Coefficienti strutturali γ₁…γ₉ | 9 |
| Varianza residua ζ del fattore | 1 |

### 2.4 Effetti diretti

Il MIMIC assume che le covariate influenzino gli indicatori **esclusivamente
attraverso il fattore latente**. Se ad esempio `Bottom` avesse un effetto diretto
su `PercTotPhos` indipendentemente dalla fertilità (il fosforo è poco mobile e
si concentra in superficie per ragioni fisiche), questo effetto verrebbe assorbito
nel residuo ε₃. Testabile a posteriori con i modification indices (Topic 6).

---

> **Nota (analogia autoencoder):** la struttura del MIMIC ricorda un autoencoder:
> encoder = indicatori → fattore latente, decoder = fattore latente → ricostruzione
> degli indicatori, con le covariate che condizionano la rappresentazione latente.
> Un VAE (Variational Autoencoder) è matematicamente molto vicino a un modello
> fattoriale con decoder non lineare. La differenza pratica: il MIMIC è lineare
> e dà inferenza statistica completa (p-value, IC); il VAE è non lineare ma senza
> inferenza formale.

---

## Topic 3 — La struttura gerarchica {#topic-3}

### 3.1 Perché ignorarla è un problema

Il dataset ha 40 Field × ~5-6 profondità ≈ 220 righe. Le righe non sono
indipendenti: le 6 misure dentro Field 1 condividono lo stesso suolo, la stessa
storia di gestione, la stessa posizione geografica. Trattarle come indipendenti
sottostima gli errori standard → p-value troppo piccoli → falsi positivi.

### 3.2 I due livelli di variabilità

| Variabile | Varia dentro Field? | Livello |
|---|---|---|
| `Bottom` | sì | L1 (within) |
| `PH`, `Texture1`, `Texture2`, `BulkDensity` | sì | L1 (within) |
| `OnFarm`, `Irrigate`, `Fertilised`, `N_Natural` | no | L2 (between) |
| `PercSOC`, `PercTotNitro`, `PercTotPhos` | sì | L1 (within) |

### 3.3 Il random intercept nel modello strutturale

$$\eta_{ij} = \underbrace{\gamma_1 \cdot \text{OnFarm}_j + \ldots}_{\text{between-Field}} + \underbrace{\gamma_5 \cdot \text{Bottom}_{ij} + \ldots}_{\text{within-Field}} + \underbrace{u_j}_{\text{random intercept Field}} + \zeta_{ij}$$

dove $u_j \sim N(0, \sigma^2_u)$. Il rapporto $\sigma^2_u / (\sigma^2_u + \sigma^2_\zeta)$
è l'**ICC** — già stimata negli esperimenti con `lme4` per log(SOC), dove era
sostanziale.

### 3.4 L'autocorrelazione lungo Bottom

Dopo aver tolto il trend di Bottom e il random intercept, i residui a 20 e 30 cm
dentro lo stesso Field sono probabilmente ancora più simili tra loro che i residui
a 20 e 80 cm. Strutture di correlazione residua:

- **AR(1):** correlazione tra residui adiacenti = ρ, a 2 passi = ρ², ecc.
- **Exponential (`corExp` in `nlme`):** correlazione = exp(−d/range) con d distanza
  reale in cm. Più flessibile perché le profondità non sono equispaziate.

Negli esperimenti è già emersa eteroschedasticità (varianza residua di log(SOC)
decrescente con la profondità) — ulteriore strato di complessità gestibile con
`nlme` (`varIdent`, `varExp`), non con `lavaan`.

### 3.5 Il modello completo ideale

1. Struttura fattoriale (3 indicatori → 1 fattore latente)
2. Regressione sul fattore latente (9 covariate, L1 e L2)
3. Random intercept per Field
4. Correlazione residua lungo Bottom
5. Eteroschedasticità dei residui per profondità

Nessun pacchetto R fa tutte e 5 le cose insieme in modo pulito (→ Topic 4).

---

## Topic 4 — Multilevel MIMIC in R {#topic-4}

### 4.1 Riformulazione come LMM multivariato

Sostituendo l'equazione strutturale dentro quella di misura:

$$\mathbf{y}_{ij} = \underbrace{\boldsymbol{\Lambda}\boldsymbol{\gamma}'\mathbf{x}_{ij}}_{\text{media fissa}} + \underbrace{\boldsymbol{\Lambda} u_j}_{\text{effetto random Field}} + \underbrace{\boldsymbol{\Lambda}\zeta_{ij} + \boldsymbol{\varepsilon}_{ij}}_{\text{residuo within}}$$

Questo è un **modello lineare misto multivariato con struttura fattoriale sulla
covarianza** (Roy & Lin, 2002). Le tre risposte condividono lo stesso random effect
Field, pesato dai loadings Λ.

### 4.2 Opzione A — `lavaan` multilevel

```r
modello_ml <- '
  level: 1
    eta_W =~ PercSOC + PercTotNitro + PercTotPhos
    eta_W ~ Bottom + PH + Texture1 + Texture2 + BulkDensity

  level: 2
    eta_B =~ PercSOC + PercTotNitro + PercTotPhos
    eta_B ~ OnFarm + Irrigate + Fertilised + N_Natural
'
fit_ml <- sem(modello_ml, data = dati, cluster = "Field")
```

`eta_W` cattura la variabilità within-Field (lungo profondità), `eta_B` quella
between-Field. I loadings dei due livelli vengono stimati separatamente.

**Pro:** sintassi familiare, indici di fit standard.  
**Contro:** solo random intercept, nessuna struttura di correlazione residua lungo
Bottom.

### 4.3 Opzione B — `galamm`

```r
library(galamm)

form <- list(
  PercSOC ~ Bottom + PH + Texture1 + Texture2 + BulkDensity +
            OnFarm + Irrigate + Fertilised + N_Natural + (1 | Field),
  PercTotNitro ~ Bottom + PH + Texture1 + Texture2 + BulkDensity +
                 OnFarm + Irrigate + Fertilised + N_Natural + (1 | Field),
  PercTotPhos ~ Bottom + PH + Texture1 + Texture2 + BulkDensity +
                OnFarm + Irrigate + Fertilised + N_Natural + (1 | Field)
)

load_matrix <- matrix(c(1, NA, NA), ncol = 1)
fit_galamm <- galamm(form, load = load_matrix, data = dati)
```

**Pro:** random effects e loadings stimati insieme, più vicino alla formulazione
di Roy & Lin.  
**Contro:** pacchetto giovane (2024), autocorrelazione lungo Bottom non supportata
nativamente.

### 4.4 Opzione C — stima ad hoc (Stan/EM)

La log-verosimiglianza marginale ha forma chiusa (tutto gaussiano). La covarianza
marginale per il Field $j$ ha struttura **compound symmetry + factor analytic**:

$$\mathbf{y}_j \sim N\left(\boldsymbol{\Lambda}\boldsymbol{\gamma}'\mathbf{X}_j,\; \sigma^2_u \boldsymbol{\Lambda}\boldsymbol{\Lambda}' \otimes \mathbf{1}\mathbf{1}' + \sigma^2_\zeta \boldsymbol{\Lambda}\boldsymbol{\Lambda}' \otimes \mathbf{I} + \boldsymbol{\Theta} \otimes \mathbf{I}\right)$$

**Stan/MCMC:** permette di aggiungere `corExp`, eteroschedasticità, e produce
incertezza calibrata su tutto. ~50-80 righe di Stan.  
**EM algoritmico:** più veloce, sfrutta la forma chiusa degli E-step (quello che
`galamm` fa sotto il cofano).

### 4.5 Raccomandazione pratica

| Obiettivo | Strumento |
|---|---|
| Prototipo rapido, verifica loadings | `lavaan` multilevel |
| Modello principale con random effects corretti | `galamm` |
| Correlazione residua lungo Bottom + eteroschedasticità | Stan o `nlme` custom |

Percorso suggerito: `lavaan` per prototipo → `galamm` per modello finale → Stan
solo se la diagnostica mostra autocorrelazione residua rilevante lungo Bottom.

---

## Topic 5 — Scelte specifiche sul dataset {#topic-5}

### 5.1 Log-trasformazione degli indicatori

SOC, N e P sono percentuali positive con distribuzione destra-asimmetrica.
La log-trasformazione è necessaria prima di usarli come indicatori:

```r
dati <- dati |>
  mutate(
    logSOC = log(PercSOC),
    logN   = log(PercTotNitro),
    logP   = log(PercTotPhos)
  )
```

I loadings λ_j saranno in unità log — adimensionali rispetto al fattore latente.

### 5.2 Forma funzionale di Bottom

Dagli esperimenti il modello quadratico batte quello lineare per AIC, ma il
profilo fattoriale mostra un plateau dopo ~50 cm. Scelta più parsimoniosa:
**log(Bottom)**:

- monotona decrescente
- un solo parametro
- cattura la caduta rapida in superficie e la decelerazione in profondità
- biologicamente motivata (decadimento a tasso decrescente)

```r
fertilita ~ log(Bottom) + PH + Texture1 + Texture2 + BulkDensity + ...
```

Alternativa: `Bottom + I(Bottom^2)` con Bottom centrato. Confrontare per AIC.

> **Nota:** Bottom come variabile continua risolve automaticamente il problema
> dell'assenza della modalità 80 cm per alcuni Landuse — quelle osservazioni
> semplicemente non esistono, nessun trattamento speciale necessario.

### 5.3 Collinearità tra le dummy di gestione

`OnFarm = 1` implica sempre `Fertilised = 1` — la cella (OnFarm=1, Fertilised=0)
è vuota. Non è collinearità perfetta:

- il coefficiente di `Fertilised` è identificato dal contrasto Offfarm-fertilizzato
  vs Offfarm-non-fertilizzato
- il coefficiente di `OnFarm` cattura **l'effetto aggiuntivo dell'essere on-farm,
  al netto della fertilizzazione** — lavorazione, rotazione, meccanizzazione, ecc.

Non è necessario scegliere tra le due dummy: il modello gira. Il punto è
interpretativo: `OnFarm` non misura "fertilizzazione" ma "gestione intensiva
complessiva".

### 5.4 Outlier di PercTotPhos

Il valore 0.92% è ~10× la mediana. Approccio: stimare il modello **con e senza**
e riportare entrambi come sensitivity analysis.

```r
dati |> filter(PercTotPhos > 0.5) |> select(Field, Bottom, PercTotPhos)
```

### 5.5 Centratura e scala delle covariate continue

Standardizzare tutte le covariate continue per rendere i γ confrontabili:

```r
dati <- dati |>
  mutate(across(c(PH, Texture1, Texture2, BulkDensity), scale))
```

`log(Bottom)` va centrato alla profondità media (log(45) ≈ 3.8). Le dummy
rimangono 0/1.

---

## Topic 6 — Stima, fit e diagnostica {#topic-6}

### 6.1 Lo stimatore

Con dati continui e approssimativamente normali dopo log-trasformazione:
**MLR** (ML con errori standard di Huber-White e statistica di Satorra-Bentler):

```r
fit <- sem(modello, data = dati, cluster = "Field", estimator = "MLR")
```

MLR è asintoticamente corretto anche con non-normalità moderata e dati cluster.

### 6.2 Indici di fit

| Indice | Soglia accettabile | Cosa misura |
|---|---|---|
| **χ²** (p-value) | p > 0.05 | test esatto, sensibile all'n |
| **RMSEA** | < 0.08 (ottimo < 0.05) | errore di approssimazione per df |
| **CFI** | > 0.95 | fit relativo al modello nullo |
| **SRMR** | < 0.08 | residui standardizzati medi |

Con n moderato (~220 righe, 40 Field), preferire RMSEA e CFI al χ².

```r
summary(fit, fit.measures = TRUE)
fitMeasures(fit, c("rmsea", "cfi", "srmr", "aic", "bic"))
```

### 6.3 Cosa guardare nei risultati

```r
summary(fit, standardized = TRUE)
```

- **Loadings** (`Latent Variables`): tutti positivi e significativi? Con
  log-trasformazione, λ_SOC e λ_N attesi alti (C:N ≈ 12), λ_P più variabile.
- **Coefficienti strutturali** (`Regressions`): i γ standardizzati (`Std.all`)
  sono direttamente confrontabili.
- **Varianze residue** (`Variances`): θ_j per ogni indicatore (quota non spiegata
  dal fattore) e ψ (varianza residua del fattore dopo le covariate).

### 6.4 Modification indices

```r
modindices(fit, sort = TRUE, maximum.number = 10)
```

Casi rilevanti:
- **Effetti diretti covariata → indicatore:** es. `Bottom → logP` (fosforo poco
  mobile, concentrato in superficie per ragioni fisiche).
- **Correlazioni tra errori:** es. `logSOC ~~ logN` (coppia C-N con dinamica
  specifica non catturata dal fattore unico).

Liberare percorsi solo se c'è giustificazione sostantiva — non data dredging.

### 6.5 Convergenza e problemi numerici

- **Heywood cases:** varianze residue negative o loadings > 1 → modello mal
  identificato o covariata perfettamente collineare.
- **Warning di convergenza:** provare `optim.method = "BFGS"` o scalare meglio
  le variabili (5.5 aiuta).
- **SE = NA o molto grandi:** verificare collinearità con `car::vif()`.

---

> **Domanda (emersa durante il Topic 6):** nei dati, P e SOC hanno correlazione
> bassa in generale, ma stratificando per alcune categoriali si vede alta
> correlazione per alcune modalità e zero per altre. Il modello ne tiene conto?

**Risposta:** no, il modello standard non ne tiene conto automaticamente —
e questo è un problema reale (measurement non-invariance). Vedere il Topic 7
per le soluzioni.

---

## Topic 7 — Interpretazione {#topic-7}

### 7.1 I loadings

I loadings standardizzati (`Std.all`) sono correlazioni tra l'indicatore e il
fattore latente:

- λ ≈ 0.9: l'indicatore è quasi interamente determinato dal fattore
- λ ≈ 0.5: metà della variabilità è spiegata dal fattore, metà è specifica
- λ < 0.3: l'indicatore misura male il fattore

**Varianza specifica** = 1 − λ²_std. Se la varianza specifica di logN è quasi
zero, N non porta informazione aggiuntiva rispetto a SOC. Se è alta per logP,
P ha una dinamica propria.

### 7.2 I coefficienti strutturali

I γ standardizzati rispondono a: *di quante σ cambia η quando la covariata
aumenta di una σ, tenendo fisse le altre?*

| Covariata | Effetto atteso su η |
|---|---|
| `log(Bottom)` | negativo (fertilità decade con profondità) |
| `Texture1` (finezza) | positivo (argilla protegge materia organica) |
| `PH` | negativo (pH basso rallenta decomposizione) |
| `BulkDensity` | negativo (suolo denso = meno materia organica) |
| `OnFarm` | ambiguo (gestione intensiva) |
| `Fertilised` | ambiguo (dipende dal tipo) |

Se un γ ha segno opposto all'atteso, potrebbe segnalare il confondimento
Texture × Landuse: i siti on-farm più argillosi nell'area potrebbero far
risultare `OnFarm` positivo su `fertil_organica`.

### 7.3 La varianza residua del fattore

$$R^2_\eta = 1 - \frac{\sigma^2_\zeta}{\text{Var}(\eta)}$$

Nel contesto multilevel:
- $\sigma^2_u$: variabilità tra Field non spiegata dalle dummy di gestione
- $\sigma^2_\zeta$: variabilità within-Field non spiegata da Bottom e covariate fisiche

Il rapporto indica dove si concentra l'incertezza residua.

### 7.4 Il problema della correlazione SOC-P eterogenea e le soluzioni

La correlazione tra logSOC e logP varia per gruppo — alta per Natural
(P fa parte del ciclo organico), bassa per siti fertilizzati (P aggiunto
esternamente). Il modello a 1 fattore stima una media ponderata di questa
correlazione: non descrive bene nessuno dei due gruppi.

**Tre strade:**

**A — Due fattori latenti (raccomandato)**

```r
modello_2f <- '
  fertil_organica =~ logSOC + logN
  fertil_P        =~ logP
  fertil_organica ~~ fertil_P

  fertil_organica ~ log(Bottom) + PH + Texture1 + Texture2 + BulkDensity + ...
  fertil_P        ~ log(Bottom) + PH + Texture1 + Texture2 + BulkDensity + ...
'
```

Con un solo indicatore per `fertil_P`, fissare λ_P = 1 e varianza residua di
logP = 0 (tratta logP come risposta diretta delle covariate con correlazione
residua con `fertil_organica`).

**B — Modello multi-gruppo**

```r
fit_mg <- sem(modello, data = dati, group = "OnFarm")
```

Testa l'invarianza progressivamente: prima i loadings, poi gli intercetti.

**C — Correlazione residua libera (soluzione minimalista)**

```r
modello_rescor <- '
  fertilita =~ logSOC + logN + logP
  fertilita ~ ...covariate...
  logSOC ~~ logP
'
```

Cattura la covarianza SOC-P non spiegata dal fattore, ma in modo omogeneo
su tutti i gruppi.

---

> **Domanda:** ha senso stimare sia con uno che con due fattori e confrontare?

**Risposta:** sì, ed è l'approccio corretto. I due modelli sono nested (il
modello a 1 fattore è un caso speciale di quello a 2 fattori con λ_P su
`fertil_organica` = 0), quindi il confronto Δχ² è valido. Il confronto dà:
>
> 1. **Test formale di fit** — AIC/BIC + Δχ²
> 2. **Risposta biologica** — se la covarianza tra i due fattori è alta per
>    Natural e bassa per Farm, stai separando due regimi del fosforo
> 3. **Robustezza** — se i γ sul fattore organico cambiano poco, le conclusioni
>    gestionali sono robuste rispetto alla specifica del modello

Schema di analisi:

```
1 fattore → fit scarso / MI grandi per logSOC~~logP?
     ↓ sì
2 fattori → fit migliore + interpretazione biologica
     ↓
covarianza tra fattori per gruppo → storia gestionale
```

### 7.5 Cosa riportare

1. **Fit del modello** — RMSEA, CFI, confronto 1 vs 2 fattori (AIC/BIC + Δχ²)
2. **Loadings standardizzati** — con IC al 95%
3. **Varianze specifiche** — ridondanza SOC-N, indipendenza parziale di P
4. **Coefficienti strutturali** — tabella γ standardizzati, SE, p-value
5. **R² del fattore** — quanto le covariate spiegano della fertilità latente
6. **Covarianza tra fattori per gruppo** (se modello a 2 fattori)
7. **Sensitivity analysis** — con e senza outlier di PercTotPhos

---

## Topic 8 — Factor scores e diagnostica spaziale {#topic-8}

### 8.1 Cosa sono i factor scores (BLUP)

Dopo aver stimato M4, ogni osservazione ha un valore "implicito" del fattore latente
condizionato ai dati osservati. In lavaan si estraggono con `lavPredict()`.

| Livello | Funzione | Output | Interpretazione |
|---------|----------|--------|-----------------|
| Between-Field | `lavPredict(fit4, level=2)` | 40 righe × 2 col | Fertilità media del sito, residua rispetto alle covariate di gestione |
| Within-Field | `lavPredict(fit4, level=1)` | 220 righe × 2 col | Deviazione dalla media di Field lungo il profilo verticale |

Nota: `lavPredict` calcola la media del posterior condizionato (Empirical Bayes /
stimatore MAP). Non fornisce automaticamente l'incertezza — per quella serve un
approccio bayesiano (vedi Topic 10).

### 8.2 Risultati empirici (M4, dati completi)

**Profilo verticale within-Field (eta^W, fertilità organica)**

| Profondità (cm) | Score medio | SD |
|---|---|---|
| 20 | +0.698 | 0.59 |
| 30 | +0.327 | 0.64 |
| 40 | −0.030 | 0.69 |
| 50 | −0.310 | 0.71 |
| 60 | −0.539 | 0.70 |
| 80 | −0.293 | 0.55 |

Chiaro gradiente decrescente con plateau in profondità (≥60 cm), coerente con la
letteratura sulla distribuzione verticale della materia organica.

**Field outlier (|z| > 2 su almeno un fattore)**

| Field | Landuse | z organico | z fosforo | Interpretazione |
|---|---|---|---|---|
| 20 | 4 | **−2.71** | −0.53 | Fertilità organica molto bassa, non spiegata dalle covariate |
| 18 | 4 | **−2.44** | −0.54 | Idem — entrambi in Landuse 4 → pattern sistematico |
| 32 | 7 | −0.78 | **−2.46** | Fosforo molto basso; Landuse 7 in generale ha z_P negativi |

La Landuse 7 mostra sistematicamente score di fosforo negativi (tutti e 5 i Field
hanno z_P < 0), suggerendo un regime del ciclo del fosforo non catturato dalle dummy.

### 8.3 Diagnostica spaziale

**Moran's I sui BLUP between-Field:**

| Fattore | I | p-value | Interpretazione |
|---|---|---|---|
| Fertilità organica | 0.121 | 0.018 | Autocorrelazione significativa ma debole |
| Fosforo | 0.245 | < 0.001 | Autocorrelazione forte — struttura spaziale reale |

**Regressione BLUP ~ Lat + Long:**

- Fertilità organica: R² = 0.02, p = 0.69 → nessuna tendenza geografica lineare
- Fosforo: R² = 0.17, p = 0.03 → tendenza geografica significativa

Conclusione: il fosforo ha struttura spaziale residua non catturata dalle covariate
(gestione + tessitura). Un GP spaziale sul between-Field è giustificato per il
fattore fosforo; per la fertilità organica il segnale è più debole.

---

## Topic 9 — Correlazione spaziale nei random effects {#topic-9}

### 9.1 Il problema

M4 assume $u_j \sim \mathcal{N}(0, \psi_B)$ i.i.d. Il Moran's I mostra che questa
ipotesi è violata, specialmente per il fosforo. Field vicini geograficamente hanno
random effects simili — questo produce errori standard sottostimati per i γ between.

### 9.2 Struttura generale

Si sostituisce la varianza scalare con una matrice di covarianza strutturata:

$$\mathbf{u} \sim \mathcal{N}(\mathbf{0},\; \mathbf{K}(\boldsymbol{\theta}))$$

dove $K_{jj'} = \text{Cov}(u_j, u_{j'})$ è funzione della distanza geografica $d_{jj'}$.

### 9.3 Scelte del kernel

**Matérn 3/2** (scelta di default — buon compromesso):
$$K(d;\,\ell,\sigma^2) = \sigma^2\!\left(1 + \frac{\sqrt{3}\,d}{\ell}\right)\exp\!\left(-\frac{\sqrt{3}\,d}{\ell}\right)$$

Parametri: $\ell$ (range di decorrelazione), $\sigma^2$ (varianza spaziale).

**Nugget** — sempre consigliato con pochi cluster (40 Field):
$$\mathbf{K} = \sigma^2_w\, \mathbf{H}_{\text{Matérn}}(d;\,\ell) + \sigma^2_\epsilon\, \mathbf{I}$$

Il nugget $\sigma^2_\epsilon$ cattura la variabilità idiosincratica non spaziale
di ogni Field. Senza nugget, Field vicini diventano quasi perfettamente correlati
→ matrice quasi singolare → instabilità numerica.

### 9.4 Implementazione in R

Per il MIMIC multilevel con GP spaziale, le opzioni non-bayesiane sono limitate:
lavaan e galamm non supportano strutture di correlazione spaziale sui random effects.

| Pacchetto | Approccio | Note |
|-----------|-----------|------|
| `glmmTMB` | ML con struttura `exp(dist)` | Solo singola risposta; no struttura fattoriale |
| `sdmTMB` | ML con Matérn via TMB | Più flessibile di glmmTMB, ancora solo single-response |
| `brms` | Bayes con `gp(Lat,Long)` | Integra GP + struttura fattoriale; unica opzione per MIMIC completo |
| `INLA` | Laplace approx + SPDE | Veloce, ma specifica manuale del modello fattoriale |

Per il MIMIC multilevel completo (struttura fattoriale + GP spaziale),
**brms è l'unica opzione pratica** che combina le due componenti in un singolo modello.

---

## Topic 10 — Estensione bayesiana (M4-Bayes) {#topic-10}

### 10.1 Struttura del modello

M4-Bayes è identico a M4 nella struttura ma sostituisce la stima MLR con
campionamento dal posterior. Si implementa con `blavaan`:

```r
library(blavaan)
fit4_bayes <- bsem(mod4, data = dati_full, cluster = "Field",
                   target = "stan", n.chains = 4, burnin = 1000, sample = 2000,
                   dp = dpriors(...))
```

### 10.2 Specificazione dei prior

Si possono usare le stime di M4-MLR come centri dei prior (**empirical Bayes
prior specification**):

```r
dp <- dpriors(
  lambda = "normal(0.82, 0.3)",   # centrato sulla stima ML del loading logN
  beta   = "normal(0, 1)",        # gamma: prior regolarizzante centrato su 0
  theta  = "exponential(2)",      # varianze specifiche: prior weakly informative
  psi    = "exponential(1)"       # varianza dei fattori
)
```

**È legittimo?** Sì, con alcune avvertenze:
- Se M4-MLR converge bene, il posterior è dominato dalla likelihood → i prior fanno
  poca differenza → l'empirical Bayes è sostanzialmente innocuo
- Il vantaggio è pratico: MCMC converge più velocemente, evita regioni patologiche
- Per le varianze usa sempre prior positivi (`exponential`, `half-normal`) indipendentemente
  dalle stime ML — le varianze ML al bordo zero sono instabili come punto di partenza

### 10.3 Vantaggi rispetto a M4-MLR

1. **Incertezza completa sui factor scores** — `lavPredict` dà solo la media del
   posterior condizionato; blavaan dà l'intera distribuzione
2. **Regularizzazione automatica** — prior debolmente informativi sui γ riducono
   la varianza di stima con 40 Field (equivalente a ridge regression)
3. **Integrazione naturale del GP spaziale** — il GP si aggiunge come componente
   nel modello latente senza riscrivere il framework
4. **Test d'ipotesi più ricchi** — P(λ > 0.3 | dati) invece di solo p-value
5. **Posterior predictive check** — `pp_check(fit4_bayes)` per ogni indicatore

### 10.4 Svantaggi

1. **Tempo di calcolo** — 30–90 minuti per 4 catene × 2000 iter
2. **Specificazione dei prior** — i prior contano; serve sensitivity analysis
3. **Diagnosi MCMC** — R-hat, ESS, trace plots invece di soli indici SEM

### 10.5 M4-Bayes con GP spaziale

Il modello completo integra struttura fattoriale e GP:

$$\boldsymbol{\eta}^B_j = \boldsymbol{\gamma}^B \mathbf{x}^B_j + \mathbf{w}_j + \boldsymbol{\zeta}^B_j$$

dove $\mathbf{w}_j \sim \mathcal{GP}(0, K_{\text{Matérn}}(\mathbf{s}_j, \mathbf{s}_{j'}; \ell, \sigma^2_w))$
è il processo gaussiano spaziale e $\boldsymbol{\zeta}^B_j \sim \mathcal{N}(0, \boldsymbol{\Psi}^B)$
è la variabilità between-Field non spaziale (nugget).

**Percorso implementativo consigliato:**

```
Step 1 → Moran's I sui BLUP di M4 (già eseguito: significativo per il fosforo)
Step 2 → brms multivariato con gp(Lat, Long) sul between-Field
          stima rho: a che distanza i Field diventano indipendenti?
Step 3 → (opzionale) blavaan per posterior completo su loadings e gamma
```
