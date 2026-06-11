# stan/
# *Bayesian Hierarchical Stan Models*

---

## Italiano

Modelli Stan per l'analisi bayesiana gerarchica dei profili verticali di nutrienti del suolo.
I file `.stan` sono la sorgente compilabile; i file `.rds` con i fit MCMC sono esclusi dal
repository (`.gitignore`) per via delle dimensioni (100–700 MB ciascuno).

### Gerarchia dei modelli

I modelli sono stati sviluppati in progressione di complessità crescente:

```
random intercept → proportional slope → response-specific random slopes → multivariate random effects (final)
                                                                                ↓
                                                                    varianti A / B / full / gp
```

### File Stan

**Modelli principali (confronto)**

| File | Descrizione |
|---|---|
| `model_ri.stan` | Random intercept puro (baseline). Un effetto casuale scalare per campo. |
| `model_ri_slope.stan` | Pendenza proporzionale: `slope_j = ρ_r · intercept_j`. La slope dipende linearmente dall'intercetta. |
| `model_rirs.stan` | Pendenza risposta-specifica: `(intercept_j, slope_j)` bivariati per ciascuna risposta (SOC, N, P) separatamente. |
| `model_mvre.stan` ⭐ | **Modello finale.** Effetti casuali 6D multivariati `V[1:6, j]`. Stima le correlazioni cross-risposta a livello di campo tramite LKJ. |

**Modelli di robustezza e varianti**

| File | Descrizione |
|---|---|
| `model_mvre_A.stan` | Modello finale ridotto: senza predittori management (K_B = 0), Texture1 esclusa. Confronto in script 14. |
| `model_mvre_B.stan` | Modello finale parsimonioso: N trattato come solo random intercept (5D invece di 6D). Confronto in script 14. |
| `model_mvre_full.stan` | Modello finale + correlazioni residue tra risposte (MVNormal sui residui). Verifica indipendenza condizionale (script 16). |
| `model_mvre_gp.stan` | Modello finale + Gaussian Process sugli intercetti di campo (kernel squared-exponential). Robustezza spaziale (script 21). |
| `model_rirs_mv.stan` | Pendenza risposta-specifica + residui multivariati. Step intermedio verso il modello finale. |
| `model_gp_fixed_slope.stan` | GP + pendenza fissa globale. Alternativa spaziale non gerarchica (script 21). |
| `model_gp_full.stan` | GP completo sulle 3 risposte. Confronto computazionalmente pesante. |

### Parametri chiave del modello finale

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
- ΔELPD modello finale vs pendenza risposta-specifica = +0.4 (SE 0.9) — miglioramento marginale ma struttura interpretabile

### Compilazione e sampling

I modelli vengono compilati e campionati automaticamente dagli script R
corrispondenti (07–22) tramite `cmdstanr`. Tutti i fit usano `seed = 2024`
per garantire la riproducibilità.

```r
# Compilazione manuale
library(cmdstanr)
library(here)
mod <- cmdstan_model(here("stan", "model_mvre.stan"))
```

I file compilati (`.exe` su Windows) sono esclusi dal repository.

---

## English

Stan models for the Bayesian hierarchical analysis of soil nutrient vertical profiles.
The `.stan` files are compilable source; the `.rds` MCMC fit files are excluded from the
repository (`.gitignore`) due to size (100–700 MB each).

### Model hierarchy

Models were developed in order of increasing complexity:

```
random intercept → proportional slope → response-specific random slopes → multivariate random effects (final)
                                                                                ↓
                                                                    variants A / B / full / gp
```

### Stan files

**Main models (comparison)**

| File | Description |
|---|---|
| `model_ri.stan` | Pure random intercept (baseline). One scalar random effect per field. |
| `model_ri_slope.stan` | Proportional slope: `slope_j = ρ_r · intercept_j`. Slope depends linearly on the intercept. |
| `model_rirs.stan` | Response-specific slope: bivariate `(intercept_j, slope_j)` for each response (SOC, N, P) separately. |
| `model_mvre.stan` ⭐ | **Final model.** 6D multivariate random effects `V[1:6, j]`. Estimates cross-response field-level correlations via LKJ. |

**Robustness models and variants**

| File | Description |
|---|---|
| `model_mvre_A.stan` | Reduced final model: no management predictors (K_B = 0), Texture1 excluded. Comparison in script 14. |
| `model_mvre_B.stan` | Parsimonious final model: N treated as random intercept only (5D instead of 6D). Comparison in script 14. |
| `model_mvre_full.stan` | Final model + residual cross-response correlations (MVNormal on residuals). Tests conditional independence (script 16). |
| `model_mvre_gp.stan` | Final model + Gaussian Process on field intercepts (squared-exponential kernel). Spatial robustness (script 21). |
| `model_rirs_mv.stan` | Response-specific random slopes + multivariate residuals. Intermediate step toward the final model. |
| `model_gp_fixed_slope.stan` | GP + fixed global slope. Non-hierarchical spatial alternative (script 21). |
| `model_gp_full.stan` | Full GP over all 3 responses. Computationally intensive comparison. |

### Key parameters of the final model

```
N, J, K_W, K_B      data dimensions (observations, fields, within/between predictors)
field_id[N]         field index for each observation
logSOC, logN, logP  response variables (log scale)
X_W[N, K_W]         within-field predictors: logBottom, Texture1, Texture2, BulkDensity, PH
X_B[J, K_B]         between-field predictors: OnFarm, Irrigate, Fertilised, N_Natural

α_r                 global intercept for response r ∈ {SOC, N, P}
γ_r[K_W]            within-field fixed effects for response r
β_r[K_B]            between-field fixed effects for response r
V[6, J]             field random effects (6D multivariate)
L_Omega[6×6]        Cholesky factor of the LKJ correlation matrix
τ_alpha_r, τ_beta_r standard deviations of the V components
σ_r                 residual standard deviation for response r
```

**Key results**:
- `ρ_int(SOC, N) = +0.386` — correlation between SOC and N field intercepts
- `ρ_slope(SOC, N) > 0` — fields with a flatter SOC profile also show a flatter N profile
- ΔELPD final model vs response-specific random slopes = +0.4 (SE 0.9) — marginal improvement but interpretable structure

### Compilation and sampling

Models are compiled and sampled automatically by the corresponding R scripts (07–22)
via `cmdstanr`. All fits use `seed = 2024` for reproducibility.

```r
# Manual compilation
library(cmdstanr)
library(here)
mod <- cmdstan_model(here("stan", "model_mvre.stan"))
```

Compiled files (`.exe` on Windows) are excluded from the repository.
