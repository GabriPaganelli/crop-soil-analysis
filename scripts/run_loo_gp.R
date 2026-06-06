# =============================================================================
# run_loo_gp.R  —  LOO-CV per M-SP-RIRS-MVRE-GP vs M-SP-RIRS-MVRE
#
# Computa e confronta ELPD per il modello con Gaussian Process (script 21)
# e il modello finale MVRE (script 10).
#
# Dipende da:
#   stan/fit_msp_rirs_mvre_gp.rds   (da script 21)
#   stan/fit_msp_rirs_mvre.rds      (da script 10)
# Produce:
#   output/loo_gp.rds
# =============================================================================

library(loo); library(here)
fit_gp   <- readRDS(here("stan","fit_msp_rirs_mvre_gp.rds"))
loo_gp   <- loo(fit_gp$draws("log_lik", format="matrix"), cores=2)
cat("ELPD_GP =", round(loo_gp$estimates["elpd_loo","Estimate"],1),
    "SE =", round(loo_gp$estimates["elpd_loo","SE"],1), "\n")
pk <- loo_gp$diagnostics$pareto_k
cat("k>0.7:", sum(pk>0.7), "k>1:", sum(pk>1), "\n")
saveRDS(loo_gp, here("output","loo_gp.rds"))
rm(fit_gp, loo_gp); gc()
fit_mvre  <- readRDS(here("stan","fit_msp_rirs_mvre.rds"))
loo_mvre  <- loo(fit_mvre$draws("log_lik", format="matrix"), cores=2)
cat("ELPD_MVRE =", round(loo_mvre$estimates["elpd_loo","Estimate"],1),
    "SE =", round(loo_mvre$estimates["elpd_loo","SE"],1), "\n")
rm(fit_mvre); gc()
loo_gp <- readRDS(here("output","loo_gp.rds"))
cat("DELTA =", round(loo_gp$estimates["elpd_loo","Estimate"] -
                     loo_mvre$estimates["elpd_loo","Estimate"], 1), "\n")
print(loo_compare(list(GP=loo_gp, MVRE=loo_mvre)))
