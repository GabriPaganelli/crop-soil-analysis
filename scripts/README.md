# scripts/
# *R Analysis Pipeline*

---

## Italiano

Pipeline di analisi in R, organizzata in script numerati (00–22, più 04b e 05b) più
`run_all.R`. Tutti gli script usano `here::here()` per i path e caricano le funzioni
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
| `00_utilities.R` | Funzioni condivise: `setup_rtools()`, `load_data()`, utility EDA e visualizzazione |
| `01_preprocessing.R` | Carica `crop.csv`, calcola ILR tessitura (Texture1/2), produce `data/crop_analytic.rds` e `data/crop_full.rds` |
| `02_eda_texture.R` | EDA composizionale tessitura: CLR+PCA, ILR, triangolo USDA. Produce `fig_texture_*.pdf` |
| `03_eda_functional_form.R` | Forma funzionale (AIC), ICC, scatter intercetta-slope. Produce `fig_01`, `fig_02`, `fig_03` |
| `04_eda_depth_profiles.R` | Profili N/P per terzile, identificazione campi influenti, SNR |
| `04b_eda_profiles_summary.R` | EDA proprietà chimico-fisiche: distribuzioni, gestione, traiettorie. Produce `eda_01–04_*.pdf` |
| `05_eda_spatial.R` | Moran's I, semivariogramma, varianza between/within |
| `05b_eda_spatial_maps.R` | Mappa spaziale campi (SOC medio) e correlazioni BLUP cross-risposta. Produce `eda_mod_04–05_*.pdf` |
| `06_lmm_exploratory.R` | LMM multivariato frequentista (MCMCglmm, sommer) — esplorativo, non nel pipeline finale |
| `07_model_ri.R` | Fit **M-RI**: random intercept baseline (`model_ri.stan`) |
| `08_model_msp.R` | Fit **M-SP**: pendenza proporzionale (`model_ri_slope.stan`) |
| `09_model_msp_rirs.R` | Fit **M-SP-RIRS**: pendenza risposta-specifica (`model_rirs.stan`) |
| `10_model_final_mvre.R` | Fit **M-SP-RIRS-MVRE**: modello finale, effetti casuali 6D (`model_mvre.stan`) + LOO completo |
| `11_model_msp_landuse.R` | Variante RIRS con covariate Landuse — esplorativa |
| `12_variable_selection.R` | Selezione predittiva (*projpred*) su MVRE — produce `fig_08_projpred_*.pdf` |
| `13_projpred_validation.R` | Confronto posterori proiettati vs MVRE. Produce `fig_13_proj_vs_msp.pdf` |
| `14_model_comparison.R` | Fit varianti A (ridotta) e B (parsimoniosa), LOO-CV — produce `fig_15_loo_ab_comparison.pdf` |
| `15_sensitivity_pareto.R` | Diagnostica Pareto-k, refit MVRE senza osservazioni influenti — produce `fig_16_sensitivity_*.pdf` |
| `16_figures_main.R` | Genera figure principali per il report (`fig_04` → `fig_18`) in `output/figures/` |
| `17_figures_report.R` | Genera figure sensibilità/comparazione (`fig_15`, `fig_16`, `fig_17`) in `output/figures/` |
| `18_robustness_mvre_full.R` | MVRE-FULL: aggiunge correlazioni residue (stoichiometria residua). Verifica indipendenza condizionale |
| `19_frequentist_nlme.R` | Confronto frequentista: `nlme` con struttura di correlazione esponenziale cross-risposta (robustezza) |
| `20_spatial_confounding.R` | Diagnostica confounding spaziale: Moran's I sulle covariate |
| `21_robustness_gp.R` | Robustezza con Gaussian Process sugli intercetti (`model_mvre_gp.stan`) + LOO GP vs MVRE |
| `22_robustness_latlong.R` | Robustezza con Lat/Long come covariate between-field aggiuntive |

### Dipendenze tra script

```
01 → {02, 03, 04, 04b, 05, 05b, 06}   (EDA usa crop_analytic.rds prodotto da 01)
01 → 07 → 08 → 09 → 10                (pipeline modelli)
10 → {12, 13, 14, 15, 18}             (validazione e robustezza usano il fit finale)
{07–15} → {16, 17}                    (figure usano tutti i fit principali)
10 → {18, 19, 20, 21, 22}             (check di robustezza)
```

---

## English

R analysis pipeline with numbered scripts (00–22, plus 04b and 05b) and `run_all.R`.
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
| `00_utilities.R` | Shared functions: `setup_rtools()`, `load_data()`, EDA and plotting utilities |
| `01_preprocessing.R` | Loads `crop.csv`, computes texture ILR (Texture1/2), produces `data/crop_analytic.rds` and `data/crop_full.rds` |
| `02_eda_texture.R` | Compositional texture EDA: CLR+PCA, ILR, USDA triangle. Produces `fig_texture_*.pdf` |
| `03_eda_functional_form.R` | Functional form (AIC), ICC, intercept-slope scatter. Produces `fig_01`, `fig_02`, `fig_03` |
| `04_eda_depth_profiles.R` | N/P profiles by tertile, influential field identification, SNR |
| `04b_eda_profiles_summary.R` | EDA soil chemical-physical properties: distributions, management, trajectories. Produces `eda_01–04_*.pdf` |
| `05_eda_spatial.R` | Moran's I, semivariogram, between/within variance |
| `05b_eda_spatial_maps.R` | Spatial field map (mean SOC) and BLUP cross-response correlations. Produces `eda_mod_04–05_*.pdf` |
| `06_lmm_exploratory.R` | Frequentist multivariate LMM (MCMCglmm, sommer) — exploratory, not in final pipeline |
| `07_model_ri.R` | Fit **M-RI**: random intercept baseline (`model_ri.stan`) |
| `08_model_msp.R` | Fit **M-SP**: proportional slope (`model_ri_slope.stan`) |
| `09_model_msp_rirs.R` | Fit **M-SP-RIRS**: response-specific slope (`model_rirs.stan`) |
| `10_model_final_mvre.R` | Fit **M-SP-RIRS-MVRE**: final model, 6D random effects (`model_mvre.stan`) + full LOO |
| `11_model_msp_landuse.R` | RIRS variant with Landuse covariates — exploratory |
| `12_variable_selection.R` | Predictive variable selection (*projpred*) on MVRE — produces `fig_08_projpred_*.pdf` |
| `13_projpred_validation.R` | Projected vs MVRE posterior comparison. Produces `fig_13_proj_vs_msp.pdf` |
| `14_model_comparison.R` | Fit variants A (reduced) and B (parsimonious), LOO-CV — produces `fig_15_loo_ab_comparison.pdf` |
| `15_sensitivity_pareto.R` | Pareto-k diagnostics, MVRE refit without influential observations — produces `fig_16_sensitivity_*.pdf` |
| `16_figures_main.R` | Generates main report figures (`fig_04` → `fig_18`) in `output/figures/` |
| `17_figures_report.R` | Generates sensitivity/comparison figures (`fig_15`, `fig_16`, `fig_17`) in `output/figures/` |
| `18_robustness_mvre_full.R` | MVRE-FULL: adds residual cross-response correlations. Tests conditional independence |
| `19_frequentist_nlme.R` | Frequentist comparison: `nlme` with exponential cross-response correlation structure (robustness) |
| `20_spatial_confounding.R` | Spatial confounding diagnostics: Moran's I on covariates |
| `21_robustness_gp.R` | Robustness with Gaussian Process on field intercepts (`model_mvre_gp.stan`) + LOO GP vs MVRE |
| `22_robustness_latlong.R` | Robustness with Lat/Long as additional between-field covariates |

### Script dependencies

```
01 → {02, 03, 04, 04b, 05, 05b, 06}   (EDA uses crop_analytic.rds produced by 01)
01 → 07 → 08 → 09 → 10                (model pipeline)
10 → {12, 13, 14, 15, 18}             (validation and robustness use the final model fit)
{07–15} → {16, 17}                    (figure scripts use all main fits)
10 → {18, 19, 20, 21, 22}             (robustness checks)
```
