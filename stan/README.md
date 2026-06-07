# stan/

Modelli Stan per l'analisi bayesiana gerarchica. I file `.stan` sono il
sorgente compilabile; i file `.rds` con i fit MCMC sono esclusi dal repository
(`.gitignore`) per via delle dimensioni (100–700 MB ciascuno).

---

## Gerarchia dei modelli

I modelli sono stati sviluppati in progressione di complessità crescente:

```
M-RI  →  M-SP  →  M-SP-RIRS  →  M-SP-RIRS-MVRE  ← modello finale
                                        ↓
                              varianti A / B / full / gp
```

---

## File Stan

### Modelli principali (confronto)

| File | Sigla | Descrizione |
|---|---|---|
| `model_ri.stan` | **M-RI** | Random intercept puro (baseline). Un effetto casuale scalare per campo. |
| `model_ri_slope.stan` | **M-SP** | Aggiunge pendenza proporzionale: `slope_j = ρ_r · intercept_j`. La slope dipende linearmente dall'intercetta. |
| `model_rirs.stan` | **M-SP-RIRS** | Pendenza risposta-specifica: `(intercept_j, slope_j)` bivariati per ciascuna risposta (SOC, N, P) separatamente. |
| `model_mvre.stan` | **M-SP-RIRS-MVRE** ⭐ | **Modello finale.** Effetti casuali 6D multivariati `V[1:6, j]` = (int_SOC, slope_SOC, int_N, slope_N, int_P, slope_P). Stima le correlazioni cross-risposta a livello di campo tramite LKJ. |

### Modelli di robustezza e varianti

| File | Descrizione |
|---|---|
| `model_mvre_A.stan` | MVRE ridotto: senza predittori management (K_B = 0), Texture1 esclusa. Confronto in script 14. |
| `model_mvre_B.stan` | MVRE parsimonioso: N trattato come solo random intercept (5D invece di 6D). Confronto in script 14. |
| `model_mvre_full.stan` | MVRE + correlazioni residue tra risposte (MVNormal sui residui). Verifica indipendenza condizionale (script 16). |
| `model_mvre_gp.stan` | MVRE + Gaussian Process sugli intercetti di campo (kernel squared-exponential). Robustezza spaziale (script 21). |
| `model_rirs_mv.stan` | RIRS + residui multivariati. Step intermedio tra RIRS e MVRE. |
| `model_gp_fixed_slope.stan` | GP + pendenza fissa globale. Alternativa spaziale non gerarchica (script 21, confronto). |
| `model_gp_full.stan` | GP completo sulle 3 risposte. Confronto computazionalmente pesante. |

---

## Parametri chiave del modello finale (`model_mvre.stan`)

```
N, J, K_W, K_B      dimensioni dati (osservazioni, campi, predittori within/between)
field_id[N]         indice campo per ogni osservazione
logSOC, logN, logP  variabili risposta (log-scala)
X_W[N, K_W]         predittori within-field: logBottom, Texture1, Texture2, BulkDensity, PH
X_B[J, K_B]         predittori between-field: OnFarm, Irrigate, Fertilised, N_Natural

α_r                 intercetta globale per risposta r ∈ {SOC, N, P}
γ_r[K_W]            effetti fissi within-field per risposta r
β_r[K_B]            effetti fissi between-field per risposta r
V[6, J]             effetti casuali di campo (6D multivariati)
L_Omega[6×6]        Cholesky della matrice di correlazione LKJ
τ_alpha_r, τ_beta_r deviazioni standard delle componenti V
σ_r                 deviazione standard residua per risposta r
```

**Risultati principali**:
- `ρ_int(SOC, N) = +0.386` — correlazione tra intercette di campo SOC e N
- `ρ_slope(SOC, N) > 0` — i campi con profilo SOC più uniforme hanno anche N più uniforme
- ΔELPD MVRE vs RIRS = +0.4 (SE 0.9) — miglioramento marginale ma struttura interpretabile

---

## Compilazione e sampling

I modelli vengono compilati e campionati automaticamente dagli script R
corrispondenti (07–22) tramite `cmdstanr`. Tutti i fit usano `seed = 2024`
per garantire la riproducibilità.

```r
# Esempio: compilazione manuale
library(cmdstanr)
library(here)
mod <- cmdstan_model(here("stan", "model_mvre.stan"))
```

I file compilati (`.exe` su Windows) sono esclusi dal repository.
