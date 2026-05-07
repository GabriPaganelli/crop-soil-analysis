# Descrizione Dataset Crop

## Data Contracts

- Le variabili binarie (`OnFarm`, `Irrigate`, `Fertilised`, `N_Natural`) costituiscono una **combinazione lineare perfetta** per l'identificazione di ogni livello di `OwnId`. Non si tratta solo di una correlazione, ma di una mappatura univoca (7 combinazioni uniche per 7 classi).
- Il dataset non è un cubo perfetto. La profondità **60-80 cm** non presenta righe (record) per le classi `Offfarm IF`, `Offfarm INF`, `Offfarm NIF` e per uno solo dei due siti `Offfarm NINF`.
- **PercClay + PercSilt + PercSand ≈ 100%:** vincolo compositivo. Le tre frazioni non sono indipendenti — includerne tutte e tre nel modello crea collinearità perfetta. Verificare se la somma è esattamente 100% o se ci sono scostamenti (presenza di scheletro/materia organica?).
- **Class è deterministica dato (PercClay, PercSilt, PercSand):** la classe tessiturale è derivata dal triangolo tessiturale, non è una misura indipendente. Alcune righe (Field 21-40) hanno Class vuota — da verificare se ricostruibile.
- **BulkDensity ha un range sospetto (0.62 – 2.06):** valori > 1.7 indicano compattazione severa, valori < 0.7 sono rari per suoli minerali. Possibili outlier da verificare.
- **PercTotPhos ha distribuzione fortemente asimmetrica:** media (0.105) >> mediana (0.08), con un picco a 0.92% che è un valore "astronomico" — potenziale hotspot di fertilizzante concentrato o errore.
- **Bottom non è uniforme tra i tipi:** 6 profondità (20-80) per Farm e Natural, 5 profondità (20-60) per gli Offfarm e per il secondo blocco di NINF. Manca la profondità 80 cm.
- Field=1 e Plot=1 dentro Farm NI sono lo **stesso identico sito fisico**. Field=6 e Plot=1 dentro Farm I sono un altro sito fisico, cioè lo stesso pezzetto di terra ma numerato in due modi diversi. Plot si resetta, Field no. Nel modello usiamo Field come ID univoco.

## Implicazioni

- **Problema della multicollinearità:** l'uso simultaneo di `OwnId` e delle variabili dummy come esplicative causerebbe una singolarità nella matrice del disegno ($X'X$ non invertibile). Statisticamente, questo genera coefficienti `NA` poiché l'informazione portata dalle dummy è già interamente contenuta nel fattore `OwnId`.
- **Equivalenza delle previsioni:** un modello saturo con tutte le interazioni delle dummy ($N\_Natural + OnFarm \times Irrigate \times Fertilised$) produrrebbe gli stessi valori predetti (*fitted values*) del modello con `OwnId`. La differenza risiede solo nella scomposizione algebrica dei coefficienti, non nel potere esplicativo.
- **Assenza di record a 80 cm per alcuni siti Offfarm:** analisi comparative esclusivamente a parità di profondità oppure tenendone conto in modo opportuno, per evitare bias sistematici nel calcolo delle medie globali (i siti campionati più in profondità risulterebbero artificialmente impoveriti).
- **Vincolo compositivo clay + silt + sand:** impone di non usare tutte e tre come predittori. Si può usarne due (tipicamente clay e silt, omettendo sand) oppure trasformare con log-ratio (composizionale).
- **BulkDensity necessaria per convertire % in stock:** senza di essa non si può dire quale gestione sequestra più carbonio in termini assoluti (t/ha). Va usata come esplicativa o per calcolare gli stock come variabile risposta derivata.
- **Autocorrelazione spaziale:** potrebbe emergere nei residui. I due blocchi di NINF sono in zone geografiche distinte, e la posizione geografica potrebbe spiegare variabilità residua non catturata dal modello.

## Decisioni da prendere

- Scegliere se usare `OwnId` oppure le dummy.
- Come gestire la profondità 80 presente solo in alcuni tipi di gestione.
- Scelta delle covariate tessiturali: usare solo 2 su 3 tra `PercClay`/`PercSilt`/`PercSand` per il vincolo compositivo, oppure trasformare.
- Come trattare `Class`: usarla come fattore di stratificazione/blocco, come covariata alternativa alle percentuali, o eliminarla perché derivata? E come gestire i valori mancanti?
- Modellare le 3 risposte separatamente o congiuntamente? `PercSOC`, `PercTotNitro` e `PercTotPhos` sono fortemente correlate (C:N ≈ 12.4). Modelli univariati separati vs modello multivariato.
- `Bottom` come variabile continua o categorica? Continua permette interpolazione e meno parametri, categorica non impone forma funzionale.
- Verificare e gestire gli outlier di `PercTotPhos` (0.92%) e `BulkDensity` (estremi).
- Verificare l'autocorrelazione spaziale nei residui dei modelli — eventualmente includere `Lat`/`Long` o una struttura spaziale.

## Relazioni potenziali da teoria della fattoria

### Relazioni con le risposte (Y)

- `PercSOC` vs `Bottom` — calo atteso con profondità, probabilmente esponenziale.
- `PercTotNitro` vs `Bottom` — stesso pattern di SOC.
- `PercTotPhos` vs `Bottom` — fosforo poco mobile, concentrato in superficie.
- `PercSOC` vs `OwnId` — Natural atteso più alto, differenze per gestione.
- `PercTotNitro` vs `OwnId` — segue SOC.
- `PercTotPhos` vs `OwnId` — Farm fertilizzate attese più alte.
- `PercSOC` vs `PercClay` — correlazione positiva attesa (argilla protegge carbonio).
- `PercSOC` vs `PercSand` — correlazione negativa attesa (sabbia = meno protezione).
- `PercSOC` vs `pH` — correlazione inversa attesa (pH basso rallenta decomposizione).
- `PercSOC` vs `BulkDensity` — correlazione inversa forte (più organico = meno denso).
- **Confondimento gestione-tessitura:** se i siti Natural sono tutti argillosi e i Farm tutti sabbiosi, le differenze in SOC potrebbero essere spiegate dalla tessitura anziché dalla gestione. `PercClay` va usata come covariata per "pulire" l'effetto.

### Relazioni tra risposte (dipendenza)

- `PercTotNitro` vs `PercSOC` — correlazione altissima attesa (C:N ≈ 12), quasi deterministica.
- `PercTotPhos` vs `PercSOC` — correlazione più debole, fosforo ha anche componente minerale.
- C:N ratio per tipo di gestione — se Farm I ha C:N < 10, la fertilizzazione sta "bruciando" il carbonio.

### Relazioni tra covariate (da verificare per collinearità)

- `PercClay` + `PercSilt` + `PercSand` — vincolo compositivo, verificare somma ≈ 100%.
- `PercClay` vs `Bottom` — argilla attesa in aumento con profondità (illuviazione).
- `pH` vs `OwnId` — Farm atteso più neutro, Natural più acido.
- `BulkDensity` vs `Bottom` — aumento atteso con profondità.
- `pH` vs `PercClay` — argilla dà potere tampone, pH più stabile.

### Relazioni spaziali

- Mappa `Lat`/`Long` colorata per `OwnId` — verificare se Natural è isolato geograficamente.
- Mappa `Lat`/`Long` con `PercSOC` — verificare gradienti spaziali non spiegati dalla gestione.
- Confronto dei due blocchi NINF (Field 11-15 vs 36-40) — stesse proprietà nonostante la distanza?

## Business Layer

### Come si legano OwnId, Field e Plot

Il dataset ha una struttura gerarchica a 3 livelli. `OwnId` identifica il **tipo di gestione** del terreno (7 classi: Farm NI, Farm I, Offfarm NINF, Natural, ecc.). Ogni tipo contiene più **siti fisici** (pezzi di terra delle "fattorie"), ciascuno identificato da un numero univoco globale nella colonna `Field` (da 1 a 40). `Plot` è semplicemente la stessa unità di `Field`, ma numerata come posizione relativa dentro la propria classe (da 1 a 5, si resetta ad ogni cambio di `OwnId`). Infine, dentro ogni sito, sono state prelevate misurazioni a diverse profondità (`Bottom`: 20, 30, 40, 50, 60, 80 cm).

### Esempio concreto con 3 siti

```
OwnId: Farm NI                    OwnId: Farm I
(irrigata, non fertilizzata)       (irrigata, fertilizzata)
│                                  │
├── Field=1  (Plot=1)              ├── Field=6  (Plot=1)
│   ├── Bottom=20 → riga 1         │   ├── Bottom=20 → riga 31
│   ├── Bottom=30 → riga 2         │   ├── Bottom=30 → riga 32
│   ├── Bottom=40 → riga 3         │   └── ... (6 profondità)
│   ├── Bottom=50 → riga 4         │
│   ├── Bottom=60 → riga 5         ├── Field=7  (Plot=2)
│   └── Bottom=80 → riga 6         │   └── ... (6 profondità)
│                                  │
├── Field=2  (Plot=2)              └── ... (5 siti totali)
│   └── ... (6 profondità)
│
└── ... (5 siti totali)


OwnId: Natural
(foresta, riferimento ecologico)
│
├── Field=16 (Plot=1)
│   └── ... (6 profondità)
│
├── Field=17 (Plot=2)
│   └── ... (6 profondità)
│
└── ... (5 siti totali)
```

**Punto chiave:** `Field=1` e `Plot=1` dentro Farm NI sono lo **stesso identico sito fisico**. `Field=6` e `Plot=1` dentro Farm I sono un altro sito fisico, cioè lo stesso pezzetto di terra ma numerato in due modi diversi. `Plot` si resetta, `Field` no. Nel modello usiamo `Field` come ID univoco.

**Eccezione:** Offfarm NINF ha 10 siti (Field 11-15 e 36-40) anziché 5, in due zone geografiche distinte (verificato da Lat/Long). Questo causa lo sbilanciamento del dataset (55 righe vs 25-30 degli altri tipi).
