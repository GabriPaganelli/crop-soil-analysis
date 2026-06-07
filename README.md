# Analisi gerarchica bayesiana dei profili verticali di nutrienti del suolo
### SOC, N e P in 40 campi agricoli in Tanzania

**Autori**: Bortoletto Davide · Faedo Piero · Paganelli Gabriele · Zannini Pietro  
**Corso**: Statistica per l'Iterazione — Università degli Studi di Padova, a.a. 2025/26  
**Report**: [`report.pdf`](report.pdf)

---

## Abstract

La distribuzione verticale dei nutrienti del suolo è un indicatore chiave della qualità agronomica, ma i meccanismi che determinano le differenze tra campi nella *forma* del profilo di profondità restano poco quantificati. In questo lavoro analizziamo i profili verticali di carbonio organico (SOC), azoto totale (N) e fosforo totale (P) in 40 campi agricoli in Tanzania, con misurazioni fino a sei profondità per campo (20–80 cm), per un totale di 220 osservazioni.

Il modello proposto — un modello bayesiano gerarchico con pendenza proporzionale e effetti casuali multivariati (M-SP-RIRS-MVRE) — stima esplicitamente come la variabilità inter-campo del profilo verticale cambi con la profondità. I risultati mostrano che **i campi con maggiore contenuto medio di carbonio organico presentano profili significativamente più uniformi**: la fertilità carboniosa si conserva meglio in profondità. Per l'azoto questo pattern è assente; per il fosforo è presente con evidenza moderata.

La selezione predittiva (*projpred*) identifica la **tessitura fine** (contrasto frazioni fini/sabbia, ILR) come principale predittore *within-field*. Le pratiche di gestione risultano non significative dopo aver controllato per la struttura fisica del suolo.

---

## Dataset

| Variabile | Dettaglio |
|---|---|
| Osservazioni | 220 (40 campi × max 6 profondità) |
| Profondità | 20, 30, 40, 50, 60, 80 cm |
| Variabili risposta | SOC (%), N totale (%), P totale (%) — scala log |
| Predittori *within-field* | log(profondità), Texture1, Texture2 (ILR), BulkDensity, PH |
| Predittori *between-field* | OnFarm, Irrigate, Fertilised, N_Natural |
| Fonte dati | `data/crop.csv` |

---

## Struttura del repository

```
crop-soil-analysis/
│
├── report.pdf              # Report finale
├── README.md
├── LICENSE                 # MIT
│
├── data/
│   ├── crop.csv            # Dati grezzi (fonte canonica)
│   ├── crop.rds            # Dataset pre-elaborato (da 01_preprocessing.R)
│   └── dati.rds            # Dataset con ILR, log-risposte, scaling
│
├── stan/                   # Modelli Stan — vedi stan/README.md
│   ├── model_mvre.stan     # Modello finale (M-SP-RIRS-MVRE)
│   └── ...
│
├── scripts/                # Pipeline R — vedi scripts/README.md
│   ├── 00_utilities.R      # Funzioni condivise
│   ├── 01_preprocessing.R
│   ├── 02_eda_texture.R    # → 22_robustezza_latlong.R
│   └── run_all.R           # Master script (esegue 07–22)
│
├── output/
│   ├── figures/            # Tutte le figure (EDA + modello + mappe)
│   └── tables/             # Tabelle CSV con risultati
│
├── report/
│   ├── report.tex          # Sorgente LaTeX
│   └── refs.bib
│
└── docs/
    └── crop.pdf            # Documentazione dati originali
```

---

## Dipendenze

### R (≥ 4.3)

```r
install.packages(c(
  "tidyverse",    # manipolazione dati e visualizzazione
  "here",         # path relativi
  "cmdstanr",     # interfaccia R → CmdStan
  "posterior",    # manipolazione draws MCMC
  "loo",          # LOO-CV e PSIS
  "projpred",     # variable selection predittiva
  "compositions", # ILR per dati composizionali (tessitura)
  "soiltexture",  # triangolo tessitura USDA (TT.plot)
  "patchwork",    # composizione figure ggplot
  "lme4",         # modelli LMM frequentisti (EDA)
  "lmerTest",     # p-value per lme4
  "nlme",         # modelli con struttura di correlazione (robustezza)
  "spdep",        # Moran's I (autocorrelazione spaziale)
  "ggbeeswarm",   # beeswarm plot (EDA)
  "MCMCglmm",     # LMM multivariato MCMC (esplorativo)
  "sommer"        # LMM multivariato frequentista (esplorativo)
))
```

### CmdStan (≥ 2.33)

```r
# Installazione CmdStan via cmdstanr
cmdstanr::install_cmdstan()

# Verificare versione
cmdstanr::cmdstan_version()
```

### Windows: Rtools

Su Windows è necessario [Rtools](https://cran.r-project.org/bin/windows/Rtools/) per compilare i modelli Stan.

---

## Quick start

```r
library(here)

# 1. Pre-elaborazione dati (una volta sola — produce data/dati.rds)
source(here("scripts", "01_preprocessing.R"))

# 2. EDA (opzionale, non richiesta per il modello)
source(here("scripts", "02_eda_texture.R"))

# 3. Pipeline completa: fit modelli 07–22 + figure
source(here("scripts", "run_all.R"))
```

> **Nota**: I fit Stan (file `.rds`, ~100–700 MB ciascuno) non sono inclusi nel repository (`.gitignore`). La prima esecuzione compila e campiona tutti i modelli — può richiedere diverse ore. I fit vengono salvati in `stan/` e riutilizzati nelle esecuzioni successive.

---

## Modello finale: M-SP-RIRS-MVRE

Il modello finale è **M-SP-RIRS-MVRE** (script `10_model_final_mvre.R`, Stan file `stan/model_mvre.stan`).

**Struttura**: per ogni campo *j* e profondità *i*:

```
logY_{ij} = α_r + V[2r-1, j] + V[2r, j] · logBottom_i + γ_r · X_W_i + β_r · X_B_j + ε_{ij}
```

dove `V[1:6, j]` è un vettore 6D di effetti casuali multivariati (intercette e slope per SOC, N, P), con struttura di covarianza LKJ che stima le correlazioni cross-risposta a livello di campo.

**Risultato chiave**: `ρ_int(SOC, N) = +0.386` — i campi ricchi di SOC tendono ad avere profili più uniformi anche per N.

Per la gerarchia completa dei modelli confrontati, vedi [`stan/README.md`](stan/README.md).

---

## Citazione

Se utilizzi questo codice o i risultati, cita il report:

> Bortoletto D., Faedo P., Paganelli G., Zannini P. (2026). *Analisi gerarchica bayesiana dei profili verticali di nutrienti del suolo: SOC, N e P in 40 campi agricoli in Tanzania*. Università degli Studi di Padova.

---

## Licenza

Codice rilasciato sotto licenza [MIT](LICENSE).
