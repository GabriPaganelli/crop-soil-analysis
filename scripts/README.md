# scripts/
# *R Analysis Pipeline*

---

## Italiano

Pipeline di analisi in R, organizzata in 23 script numerati (00–22) più script EDA
supplementari. Tutti gli script usano `here::here()` per i path e caricano le funzioni
condivise da `00_utilities.R`.

### Esecuzione completa

```r
source(here::here("scripts", "run_all.R"))
```

`run_all.R` esegue in sequenza gli script **07–22** con gestione degli errori e caching
dei fit Stan: se il file `.rds` esiste già, il fit viene caricato invece di essere
ricalcolato.

### Pipeline numerata

| Script | Descrizione |
|---|---|
| `00_utilities.R` | Funzioni condivise: `setup_rtools()`, `carica_dati()`, utility EDA e visualizzazione |
| `01_preprocessing.R` | Carica `crop.csv`, calcola ILR tessitura (Texture1/2), produce `data/dati.rds` e `data/crop.rds` |
| `02_eda_texture.R` | EDA composizionale tessitura: CLR+PCA, ILR, triangolo USDA. Produce `fig_texture_*.pdf` |
| `03_eda_functional_form.R` | Forma funzionale (AIC), ICC, scatter intercetta-slope. Produce `fig_01`, `fig_02`, `fig_03` |
| `04_eda_depth_profiles.R` | Profili N/P per terzile, identificazione campi influenti, SNR |
| `05_eda_spatial.R` | Moran's I, semivariogramma, varianza between/within |
| `06_lmm_exploratory.R` | LMM multivariato frequentista (MCMCglmm, sommer) — esplorativo, non nel pipeline finale |
| `07_model_ri.R` | Fit **M-RI**: random intercept baseline (`model_ri.stan`) |
| `08_model_msp.R` | Fit **M-SP**: pendenza proporzionale (`model_ri_slope.stan`) |
| `09_model_msp_rirs.R` | Fit **M-SP-RIRS**: pendenza risposta-specifica (`model_rirs.stan`) |
| `10_model_final_mvre.R` | Fit **M-SP-RIRS-MVRE**: modello finale, effetti casuali 6D (`model_mvre.stan`) + LOO completo |
| `11_model_msp_landuse.R` | Variante RIRS con covariate Landuse — esplorativa |
| `12_selezione_variabili.R` | Selezione predittiva (*projpred*) su MVRE — produce `fig_08_projpred_*.pdf` |
| `13_validazione_projpred.R` | Confronto posterori proiettati vs MVRE. Produce `fig_13_proj_vs_msp.pdf` |
| `14_confronto_modelli.R` | Fit varianti A (ridotta) e B (parsimoniosa), LOO-CV — produce `fig_15_loo_ab_comparison.pdf` |
| `15_sensitivity_pareto.R` | Diagnostica Pareto-k, refit MVRE senza osservazioni influenti — produce `fig_16_sensitivity_*.pdf` |
| `16_robustezza_msp_mv.R` | MVRE-FULL: aggiunge correlazioni residue (stoichiometria residua). Verifica indipendenza condizionale |
| `17_frequentista_corexp.R` | Confronto frequentista: `nlme` con struttura di correlazione exponential (robustezza) |
| `18_figure_principali.R` | Genera figure principali per il report (`fig_04` → `fig_18`) in `output/figures/` |
| `19_figure_report.R` | Genera figure sensibilità/comparazione + copia in `report/images/` per compilazione LaTeX |
| `20_spatial_confounding.R` | Diagnostica confounding spaziale: Moran's I sulle covariate |
| `21_robustezza_gp.R` | Robustezza con Gaussian Process sugli intercetti (`model_mvre_gp.stan`) + LOO GP vs MVRE |
| `22_robustezza_latlong.R` | Robustezza con Lat/Long come covariate between-field aggiuntive |

### Script EDA supplementari

Generano figure salvate in `output/figures/`. Non vengono eseguiti da `run_all.R`.

| Script | Output |
|---|---|
| `eda_figures_profiles.R` | `eda_01_distribuzioni.pdf`, `eda_02_gestione.pdf`, `eda_03_traiettorie.pdf`, `eda_04_traiettorie_orig.pdf` |
| `eda_figures_trajectories.R` | `eda_traj_onfarm.pdf`, `eda_traj_irrigate.pdf`, `eda_traj_fertilised.pdf`, `eda_traj_gestione.pdf`, `eda_traj_gestione3.pdf` |
| `eda_figures_model_motivation.R` | `eda_mod_04_mappa.pdf` (mappa spaziale campi), `eda_mod_05_blup_cross.pdf` (correlazione BLUP intercette) |

### Dipendenze tra script

```
01 → {02, 03, 04, 05, 06}   (EDA usa dati.rds prodotto da 01)
01 → 07 → 08 → 09 → 10      (pipeline modelli, ogni fit si appoggia al precedente per confronto)
10 → {12, 13, 14, 15, 16}   (validazione e robustezza usano il fit del modello finale)
{07–16} → {18, 19}           (figure usano tutti i fit)
```

---

## English

R analysis pipeline with 23 numbered scripts (00–22) plus supplementary EDA scripts.
All scripts use `here::here()` for paths and load shared functions from `00_utilities.R`.

### Full execution

```r
source(here::here("scripts", "run_all.R"))
```

`run_all.R` runs scripts **07–22** sequentially with error handling and Stan fit caching:
if an `.rds` fit file already exists, it is loaded rather than resampled.

### Numbered pipeline

| Script | Description |
|---|---|
| `00_utilities.R` | Shared functions: `setup_rtools()`, `carica_dati()`, EDA and plotting utilities |
| `01_preprocessing.R` | Loads `crop.csv`, computes texture ILR (Texture1/2), produces `data/dati.rds` and `data/crop.rds` |
| `02_eda_texture.R` | Compositional texture EDA: CLR+PCA, ILR, USDA triangle. Produces `fig_texture_*.pdf` |
| `03_eda_functional_form.R` | Functional form (AIC), ICC, intercept-slope scatter. Produces `fig_01`, `fig_02`, `fig_03` |
| `04_eda_depth_profiles.R` | N/P profiles by tertile, influential field identification, SNR |
| `05_eda_spatial.R` | Moran's I, semivariogram, between/within variance |
| `06_lmm_exploratory.R` | Frequentist multivariate LMM (MCMCglmm, sommer) — exploratory, not in final pipeline |
| `07_model_ri.R` | Fit **M-RI**: random intercept baseline (`model_ri.stan`) |
| `08_model_msp.R` | Fit **M-SP**: proportional slope (`model_ri_slope.stan`) |
| `09_model_msp_rirs.R` | Fit **M-SP-RIRS**: response-specific slope (`model_rirs.stan`) |
| `10_model_final_mvre.R` | Fit **M-SP-RIRS-MVRE**: final model, 6D random effects (`model_mvre.stan`) + full LOO |
| `11_model_msp_landuse.R` | RIRS variant with Landuse covariates — exploratory |
| `12_selezione_variabili.R` | Predictive variable selection (*projpred*) on MVRE — produces `fig_08_projpred_*.pdf` |
| `13_validazione_projpred.R` | Projected vs MVRE posterior comparison. Produces `fig_13_proj_vs_msp.pdf` |
| `14_confronto_modelli.R` | Fit variants A (reduced) and B (parsimonious), LOO-CV — produces `fig_15_loo_ab_comparison.pdf` |
| `15_sensitivity_pareto.R` | Pareto-k diagnostics, MVRE refit without influential observations — produces `fig_16_sensitivity_*.pdf` |
| `16_robustezza_msp_mv.R` | MVRE-FULL: adds residual cross-response correlations. Tests conditional independence |
| `17_frequentista_corexp.R` | Frequentist comparison: `nlme` with exponential correlation structure (robustness) |
| `18_figure_principali.R` | Generates main report figures (`fig_04` → `fig_18`) in `output/figures/` |
| `19_figure_report.R` | Generates sensitivity/comparison figures + copies to `report/images/` for LaTeX |
| `20_spatial_confounding.R` | Spatial confounding diagnostics: Moran's I on covariates |
| `21_robustezza_gp.R` | Robustness with Gaussian Process on field intercepts (`model_mvre_gp.stan`) + LOO GP vs MVRE |
| `22_robustezza_latlong.R` | Robustness with Lat/Long as additional between-field covariates |

### Supplementary EDA scripts

Generate figures saved in `output/figures/`. Not executed by `run_all.R`.

| Script | Output |
|---|---|
| `eda_figures_profiles.R` | `eda_01_distribuzioni.pdf`, `eda_02_gestione.pdf`, `eda_03_traiettorie.pdf`, `eda_04_traiettorie_orig.pdf` |
| `eda_figures_trajectories.R` | `eda_traj_onfarm.pdf`, `eda_traj_irrigate.pdf`, `eda_traj_fertilised.pdf`, `eda_traj_gestione.pdf`, `eda_traj_gestione3.pdf` |
| `eda_figures_model_motivation.R` | `eda_mod_04_mappa.pdf` (spatial field map), `eda_mod_05_blup_cross.pdf` (BLUP intercept correlation) |

### Script dependencies

```
01 → {02, 03, 04, 05, 06}   (EDA uses dati.rds produced by 01)
01 → 07 → 08 → 09 → 10      (model pipeline, each fit used by the next for comparison)
10 → {12, 13, 14, 15, 16}   (validation and robustness use the final model fit)
{07–16} → {18, 19}           (figure scripts use all fits)
```
