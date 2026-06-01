# =============================================================================
# 12_figure_principali.R  —  Export figure principali dal fit M-SP
#
# Produce in output/figures/:
#   fig_04_posterior_eta.pdf       — distribuzioni a posteriori di eta_r
#   fig_05_forest_gamma.pdf        — forest plot effetti within (gamma_r)
#   fig_06_forest_beta.pdf         — forest plot effetti between (beta_r)
#   fig_07_ppc.pdf                 — posterior predictive check (SOC, N, P)
#   fig_09_loo_comparison.pdf      — confronto LOO-CV 4 modelli
#   fig_10_trace_key.pdf           — trace plots parametri chiave
#   fig_11_posterior_psi_sigma.pdf — psi_r e sigma_r a posteriori (diagnostica)
#
# Tutto caricato da .rds (nessun MCMC). LOO e projpred summary usano cache.
# Dipende da: stan/fit_msp.rds (+ fit_m*.rds per LOO), data/dati.rds
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(posterior)
  library(bayesplot)
  library(loo)
  library(patchwork)
})
source(here("scripts", "00_utilities.R"))

dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "cache"),   recursive = TRUE, showWarnings = FALSE)
fig_dir   <- here("output", "figures")
cache_dir <- here("output", "cache")

save_fig <- function(fname, p, w = 16, h = 9, u = "cm") {
  ggplot2::ggsave(file.path(fig_dir, fname), plot = p,
                  width = w, height = h, units = u, device = "pdf")
  cat(sprintf("  [fig] Salvato: %s\n", fname))
}

# Etichette covariate
cov_W_labels <- c("logBottom", "Texture1 (ILR1)", "Texture2 (ILR2)", "BulkDensity")
cov_B_labels <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")
resp_levels  <- c("SOC", "N", "P")
resp_colors  <- c("SOC" = "#2166AC", "N" = "#1A9850", "P" = "#D73027")


# ── 1. DATI ───────────────────────────────────────────────────────────────────

cat("Carico dati...\n")
dati <- carica_dati()
field_levels <- sort(unique(as.integer(as.character(dati$Field))))
J <- length(field_levels)
N <- nrow(dati)

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))

X_W_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity")
X_B_cols <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")
X_W      <- as.matrix(dati_int[, X_W_cols])
X_B      <- dati_int |>
  distinct(field_int, across(all_of(X_B_cols))) |>
  arrange(field_int) |>
  select(all_of(X_B_cols)) |>
  as.matrix()
field_id <- dati_int$field_int


# ── 2. CARICA FIT M-SP ────────────────────────────────────────────────────────

cat("Carico fit_msp.rds...\n")
fit <- readRDS(here("stan", "fit_msp.rds"))

# Summary parametri chiave
smry <- fit$summary(c(
  paste0(c("alpha_", "psi_", "eta_", "b_", "sigma_"), rep(resp_levels, each = 5)),
  paste0("gamma_SOC[", 1:4, "]"),
  paste0("gamma_N[",   1:4, "]"),
  paste0("gamma_P[",   1:4, "]"),
  paste0("beta_SOC[",  1:4, "]"),
  paste0("beta_N[",    1:4, "]"),
  paste0("beta_P[",    1:4, "]")
)) |>
  select(variable, median, mean, sd, q5, q95, rhat, ess_bulk)

cat(sprintf("Fit caricato. Parametri nel summary: %d\n", nrow(smry)))


# ── 3. FIG_04: POSTERIOR di eta_r ─────────────────────────────────────────────

cat("\n[fig_04] Posterior eta_r...\n")

draws_eta <- fit$draws(c("eta_SOC", "eta_N", "eta_P"), format = "df") |>
  pivot_longer(cols = starts_with("eta_"),
               names_to  = "par",
               values_to = "valore") |>
  mutate(risposta = factor(sub("eta_", "", par), levels = resp_levels))

ci_eta <- draws_eta |>
  group_by(risposta) |>
  summarise(
    med = median(valore),
    lo  = quantile(valore, 0.05),
    hi  = quantile(valore, 0.95),
    .groups = "drop"
  )

p_eta <- ggplot(draws_eta, aes(x = valore, fill = risposta, color = risposta)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_segment(data = ci_eta,
               aes(x = lo, xend = hi, y = -Inf, yend = -Inf,
                   color = risposta),
               linewidth = 3, alpha = 0.8,
               inherit.aes = FALSE) +
  geom_point(data = ci_eta, aes(x = med, y = 0, color = risposta),
             size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values  = resp_colors, guide = "none") +
  scale_color_manual(values = resp_colors, guide = "none") +
  facet_wrap(~ risposta, scales = "free", nrow = 1) +
  labs(
    title    = expression("Distribuzione a posteriori di " * eta[r]),
    subtitle = "Punto = mediana | Segmento = IC 90% | Linea tratteggiata = 0",
    x        = expression(eta[r] ~ "(effetto slope proporzionale)"),
    y        = "Densità"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 12))

print(p_eta)
save_fig("fig_04_posterior_eta.pdf", p_eta, w = 18, h = 8)


# ── 4. FIG_05: FOREST PLOT gamma_r (effetti within-field) ─────────────────────

cat("\n[fig_05] Forest plot gamma_r...\n")

gamma_data <- smry |>
  filter(grepl("^gamma_", variable)) |>
  mutate(
    risposta = factor(
      case_when(
        grepl("SOC", variable) ~ "SOC",
        grepl("_N\\[", variable) ~ "N",
        grepl("_P\\[", variable) ~ "P"
      ),
      levels = resp_levels
    ),
    k         = as.integer(stringr::str_extract(variable, "\\d+")),
    covariata = factor(cov_W_labels[k], levels = rev(cov_W_labels))
  )

p_gamma <- ggplot(gamma_data,
                  aes(x = median, y = covariata, color = risposta)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_linerange(aes(xmin = q5, xmax = q95), linewidth = 1.5, alpha = 0.7,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = resp_colors, name = "Risposta") +
  labs(
    title    = expression("Effetti fissi within-field " * (gamma[r])),
    subtitle = "Mediana e IC 90% posteriori | Scala standardizzata",
    x        = "Effetto stimato (scala standardizzata)",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

print(p_gamma)
save_fig("fig_05_forest_gamma.pdf", p_gamma, w = 14, h = 9)


# ── 5. FIG_06: FOREST PLOT beta_r (effetti between-field / management) ────────

cat("\n[fig_06] Forest plot beta_r...\n")

beta_data <- smry |>
  filter(grepl("^beta_", variable)) |>
  mutate(
    risposta = factor(
      case_when(
        grepl("SOC", variable) ~ "SOC",
        grepl("_N\\[", variable) ~ "N",
        grepl("_P\\[", variable) ~ "P"
      ),
      levels = resp_levels
    ),
    k         = as.integer(stringr::str_extract(variable, "\\d+")),
    covariata = factor(cov_B_labels[k], levels = rev(cov_B_labels))
  )

p_beta <- ggplot(beta_data,
                 aes(x = median, y = covariata, color = risposta)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_linerange(aes(xmin = q5, xmax = q95), linewidth = 1.5, alpha = 0.7,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = resp_colors, name = "Risposta") +
  labs(
    title    = expression("Effetti fissi between-field — management " * (beta[r])),
    subtitle = "Mediana e IC 90% posteriori",
    x        = "Effetto stimato",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

print(p_beta)
save_fig("fig_06_forest_beta.pdf", p_beta, w = 14, h = 9)

# Combinata gamma + beta in un'unica figura
p_fixed_effects <- (p_gamma + labs(title = expression(gamma[r] ~ "(within-field)"))) /
                   (p_beta  + labs(title = expression(beta[r]  ~ "(between-field)"))) +
  plot_annotation(
    title    = "Effetti fissi del modello M-SP",
    subtitle = "Mediana e IC 90% posteriori"
  )
save_fig("fig_05b_forest_gamma_beta_combined.pdf", p_fixed_effects, w = 14, h = 16)


# ── 6. FIG_07: POSTERIOR PREDICTIVE CHECK ─────────────────────────────────────
# Calcolo y_rep = mu[d,n] + sigma[d] * N(0,1) per D_ppc draws thinned

cat("\n[fig_07] Posterior predictive check...\n")

draws_raw  <- fit$draws(format = "matrix")
draws_full <- matrix(as.double(draws_raw),
                     nrow = nrow(draws_raw), ncol = ncol(draws_raw),
                     dimnames = dimnames(draws_raw))
D_total   <- nrow(draws_full)
D_ppc     <- 200L
thin_idx  <- round(seq(1, D_total, length.out = D_ppc))
draws_ppc <- draws_full[thin_idx, , drop = FALSE]
rm(draws_raw, draws_full); gc()

extract_mat_ppc <- function(draws, prefix, K) {
  cols <- paste0(prefix, "[", seq_len(K), "]")
  draws[, cols, drop = FALSE]
}

compute_yrep <- function(draws_d, r) {
  alpha <- draws_d[, paste0("alpha_", r)]
  psi   <- draws_d[, paste0("psi_", r)]
  eta   <- draws_d[, paste0("eta_", r)]
  sigma <- draws_d[, paste0("sigma_", r)]
  z_nu  <- extract_mat_ppc(draws_d, paste0("z_nu_", r), J)
  gamma <- extract_mat_ppc(draws_d, paste0("gamma_", r), length(X_W_cols))
  beta  <- extract_mat_ppc(draws_d, paste0("beta_", r),  length(X_B_cols))

  X_B_n      <- X_B[field_id, ]
  z_nu_n     <- z_nu[, field_id]
  ri_factor  <- outer(psi, rep(1, N)) + outer(eta, X_W[, 1])
  rand_part  <- z_nu_n * ri_factor
  fix_within <- tcrossprod(gamma, X_W)
  fix_betw   <- tcrossprod(beta,  X_B_n)
  intercept  <- outer(alpha, rep(1, N))
  mu         <- intercept + rand_part + fix_within + fix_betw  # D_ppc × N

  # Simula y_rep
  eps  <- matrix(rnorm(D_ppc * N), nrow = D_ppc, ncol = N)
  yrep <- mu + outer(sigma, rep(1, N)) * eps
  yrep
}

set.seed(2024)
yrep_SOC <- compute_yrep(draws_ppc, "SOC")
yrep_N   <- compute_yrep(draws_ppc, "N")
yrep_P   <- compute_yrep(draws_ppc, "P")
rm(draws_ppc); gc()

# Usa bayesplot per PPC
color_scheme_set("mix-teal-pink")

ppc_list <- list(
  SOC = ppc_dens_overlay(y = dati_int$logSOC, yrep = yrep_SOC) +
          ggtitle("PPC — logSOC") + theme_minimal(base_size = 10),
  N   = ppc_dens_overlay(y = dati_int$logN,   yrep = yrep_N)   +
          ggtitle("PPC — logN")   + theme_minimal(base_size = 10),
  P   = ppc_dens_overlay(y = dati_int$logP,   yrep = yrep_P)   +
          ggtitle("PPC — logP")   + theme_minimal(base_size = 10)
)

p_ppc <- wrap_plots(ppc_list, nrow = 1) +
  plot_annotation(
    title    = "Posterior predictive check — M-SP",
    subtitle = paste0("Linea scura = dati osservati | ",
                      D_ppc, " curve chiare = repliche dalla posteriore")
  )
print(p_ppc)
save_fig("fig_07_ppc.pdf", p_ppc, w = 20, h = 8)
rm(yrep_SOC, yrep_N, yrep_P); gc()


# ── 7. FIG_09: LOO-CV COMPARISON ──────────────────────────────────────────────

cat("\n[fig_09] LOO comparison...\n")

loo_cache_path <- file.path(cache_dir, "loo_results.rds")

if (file.exists(loo_cache_path)) {
  cat("  Carico LOO cache...\n")
  loo_list <- readRDS(loo_cache_path)
} else {
  cat("  Calcolo LOO per tutti i modelli (può richiedere 2-5 min)...\n")
  loo_list <- list()

  # M-SP (già in memoria)
  cat("  LOO M-SP...\n")
  ll_msp         <- fit$draws("log_lik", format = "matrix")
  loo_list[["M-SP"]]  <- loo(ll_msp, cores = 4)
  rm(ll_msp); gc()

  # M-RI
  fit_mri_path <- here("stan", "fit_mri.rds")
  if (file.exists(fit_mri_path)) {
    cat("  LOO M-RI...\n")
    f <- readRDS(fit_mri_path)
    ll <- f$draws("log_lik", format = "matrix")
    loo_list[["M-RI"]] <- loo(ll, cores = 4)
    rm(f, ll); gc()
  }

  # M-GPS
  fit_mgps_path <- here("stan", "fit_mgps.rds")
  if (file.exists(fit_mgps_path)) {
    cat("  LOO M-GPS...\n")
    f <- readRDS(fit_mgps_path)
    ll <- f$draws("log_lik", format = "matrix")
    loo_list[["M-GPS"]] <- loo(ll, cores = 4)
    rm(f, ll); gc()
  }

  # M-GP
  fit_mgp_path <- here("stan", "fit_mgp.rds")
  if (file.exists(fit_mgp_path)) {
    cat("  LOO M-GP...\n")
    f <- readRDS(fit_mgp_path)
    ll <- f$draws("log_lik", format = "matrix")
    loo_list[["M-GP"]] <- loo(ll, cores = 4)
    rm(f, ll); gc()
  }

  saveRDS(loo_list, loo_cache_path)
  cat("  LOO cache salvato:", loo_cache_path, "\n")
}

# LOO comparison plot
if (length(loo_list) > 1) {
  cmp <- loo_compare(loo_list)
  cat("\nLOO comparison:\n")
  print(cmp)

  cmp_df <- as.data.frame(cmp) |>
    tibble::rownames_to_column("modello") |>
    mutate(
      modello = factor(modello,
                       levels = rev(c("M-SP", "M-RI", "M-GPS", "M-GP"))),
      label   = sprintf("elpd_diff = %.1f +/- %.1f", elpd_diff, se_diff)
    )

  p_loo <- ggplot(cmp_df,
                  aes(x = elpd_diff, y = modello)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.8) +
    geom_linerange(
      aes(xmin = elpd_diff - 2 * se_diff,
          xmax = elpd_diff + 2 * se_diff),
      linewidth = 2, color = "steelblue", alpha = 0.65
    ) +
    geom_point(size = 4, color = "steelblue") +
    geom_text(aes(label = label), hjust = -0.15, size = 3.2, color = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
    labs(
      title    = "Confronto LOO-CV tra i quattro modelli",
      subtitle = "ELPD diff rispetto a M-SP (riferimento) | Barre = +/- 2 SE",
      x        = "ΔELPD",
      y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank())

  print(p_loo)
  save_fig("fig_09_loo_comparison.pdf", p_loo, w = 15, h = 8)
} else {
  cat("  Meno di 2 modelli disponibili per il confronto.\n")
}


# ── 8. FIG_10: TRACE PLOTS parametri chiave ───────────────────────────────────

cat("\n[fig_10] Trace plots...\n")

params_trace <- c("alpha_SOC", "psi_SOC", "eta_SOC",
                  "alpha_N",   "psi_N",   "eta_N",
                  "alpha_P",   "psi_P",   "eta_P")

draws_trace <- fit$draws(variables = params_trace)

p_trace <- mcmc_trace(draws_trace, facet_args = list(ncol = 3)) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold")) +
  labs(title    = "Trace plots — parametri strutturali di M-SP",
       subtitle = "4 catene × 5000 iterazioni")

print(p_trace)
save_fig("fig_10_trace_key.pdf", p_trace, w = 20, h = 14)


# ── 9. FIG_11: POSTERIOR psi_r e sigma_r ──────────────────────────────────────

cat("\n[fig_11] Posterior psi_r e sigma_r...\n")

draws_psi_sigma <- fit$draws(
  c("psi_SOC", "psi_N", "psi_P", "sigma_SOC", "sigma_N", "sigma_P"),
  format = "df"
) |>
  pivot_longer(
    cols      = c(starts_with("psi_"), starts_with("sigma_")),
    names_to  = "par",
    values_to = "valore"
  ) |>
  mutate(
    tipo     = if_else(str_starts(par, "psi_"), "psi[r]", "sigma[r]"),
    risposta = factor(str_replace(par, "^(psi|sigma)_", ""), levels = resp_levels)
  )

p_psi_sigma <- ggplot(draws_psi_sigma,
                      aes(x = valore, fill = risposta, color = risposta)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  scale_fill_manual(values  = resp_colors, name = "Risposta") +
  scale_color_manual(values = resp_colors, name = "Risposta") +
  facet_grid(risposta ~ tipo, scales = "free",
             labeller = label_parsed) +
  labs(
    title    = expression("Distribuzioni a posteriori di " * psi[r] * " e " * sigma[r]),
    subtitle = expression(psi[r] ~ "= SD effetti casuali alla prof. media  |  " *
                          sigma[r] ~ "= SD residua"),
    x        = "Valore",
    y        = "Densità"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

print(p_psi_sigma)
save_fig("fig_11_posterior_psi_sigma.pdf", p_psi_sigma, w = 14, h = 12)


# ── 10. FIG_12: RIASSUNTO PARAMETRI PRINCIPALI (tabella visiva) ───────────────

cat("\n[fig_12] Panel riassuntivo parametri strutturali...\n")

struct_params <- smry |>
  filter(variable %in% c(
    "alpha_SOC", "psi_SOC", "eta_SOC", "b_SOC", "sigma_SOC",
    "alpha_N",   "psi_N",   "eta_N",   "b_N",   "sigma_N",
    "alpha_P",   "psi_P",   "eta_P",   "b_P",   "sigma_P"
  )) |>
  mutate(
    par      = str_extract(variable, "^[^_]+"),
    risposta = factor(str_extract(variable, "[A-Z]+$"), levels = resp_levels),
    par_label = factor(par, levels = c("alpha", "psi", "eta", "b", "sigma"),
                       labels = c("alpha[r]", "psi[r]", "eta[r]",
                                  "b[r]~'='~eta[r]/psi[r]", "sigma[r]"))
  )

p_struct <- ggplot(struct_params, aes(x = median, y = risposta, color = risposta)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_linerange(aes(xmin = q5, xmax = q95), linewidth = 2, alpha = 0.7) +
  geom_point(size = 3.5) +
  scale_color_manual(values = resp_colors, guide = "none") +
  facet_wrap(~ par_label, scales = "free_x", nrow = 1,
             labeller = label_parsed) +
  labs(
    title    = "Parametri strutturali del modello M-SP",
    subtitle = "Mediana e IC 90% posteriori per SOC, N, P",
    x        = "Valore stimato",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.major.y = element_blank())

print(p_struct)
save_fig("fig_12_struct_params_panel.pdf", p_struct, w = 22, h = 7)


# ── 11. EXPORT TABELLE CSV per LaTeX ──────────────────────────────────────────

cat("\n[tabelle] Export CSV per il report...\n")
dir.create(here("output", "tables"), recursive = TRUE, showWarnings = FALSE)
tab_dir <- here("output", "tables")

# Tabella 1: summary completo parametri
smry_export <- smry |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(smry_export, file.path(tab_dir, "tab_posterior_summary.csv"),
          row.names = FALSE)
cat("  Salvato: tab_posterior_summary.csv\n")

# Tabella 2: parametri strutturali (alpha, psi, eta, b, sigma)
tab_struct <- smry |>
  filter(grepl("^(alpha_|psi_|eta_|b_|sigma_)", variable)) |>
  mutate(
    par      = stringr::str_extract(variable, "^[^_]+"),
    risposta = stringr::str_extract(variable, "[A-Z]+$"),
    par      = factor(par, levels = c("alpha","psi","eta","b","sigma")),
    risposta = factor(risposta, levels = c("SOC","N","P"))
  ) |>
  arrange(risposta, par) |>
  select(risposta, par, median, sd, q5, q95, rhat, ess_bulk) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write.csv(tab_struct, file.path(tab_dir, "tab_struct_params.csv"),
          row.names = FALSE)
cat("  Salvato: tab_struct_params.csv\n")
cat("\n  Parametri strutturali:\n")
print(tab_struct, n = Inf)

# Tabella 3: gamma_r (effetti within)
tab_gamma <- smry |>
  filter(grepl("^gamma_", variable)) |>
  mutate(
    risposta  = factor(
      case_when(grepl("SOC", variable) ~ "SOC",
                grepl("_N\\[", variable) ~ "N",
                grepl("_P\\[", variable) ~ "P"),
      levels = resp_levels),
    k         = as.integer(stringr::str_extract(variable, "\\d+")),
    covariata = cov_W_labels[k]
  ) |>
  select(risposta, covariata, median, sd, q5, q95) |>
  arrange(risposta, covariata) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write.csv(tab_gamma, file.path(tab_dir, "tab_gamma_within.csv"),
          row.names = FALSE)
cat("  Salvato: tab_gamma_within.csv\n")

# Tabella 4: beta_r (effetti between / management)
tab_beta <- smry |>
  filter(grepl("^beta_", variable)) |>
  mutate(
    risposta  = factor(
      case_when(grepl("SOC", variable) ~ "SOC",
                grepl("_N\\[", variable) ~ "N",
                grepl("_P\\[", variable) ~ "P"),
      levels = resp_levels),
    k         = as.integer(stringr::str_extract(variable, "\\d+")),
    covariata = cov_B_labels[k]
  ) |>
  select(risposta, covariata, median, sd, q5, q95) |>
  arrange(risposta, covariata) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write.csv(tab_beta, file.path(tab_dir, "tab_beta_between.csv"),
          row.names = FALSE)
cat("  Salvato: tab_beta_between.csv\n")

# Tabella 5: LOO comparison (se cache disponibile)
if (file.exists(loo_cache_path)) {
  loo_list_tab <- readRDS(loo_cache_path)
  cmp_tab      <- as.data.frame(loo::loo_compare(loo_list_tab)) |>
    tibble::rownames_to_column("modello") |>
    mutate(across(where(is.numeric), ~ round(.x, 2)))
  write.csv(cmp_tab, file.path(tab_dir, "tab_loo_comparison.csv"),
            row.names = FALSE)
  cat("  Salvato: tab_loo_comparison.csv\n")
  cat("\n  LOO comparison:\n")
  print(cmp_tab[, 1:4])
}

cat(sprintf("\nTabelle salvate in: %s\n", tab_dir))

# ── RIEPILOGO ─────────────────────────────────────────────────────────────────

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("  FIGURE PRODOTTE in:", fig_dir, "\n")
cat("═══════════════════════════════════════════════════════════════════\n")
figs <- list.files(fig_dir, pattern = "\\.pdf$", full.names = FALSE)
for (f in sort(figs)) cat(sprintf("  %s\n", f))
cat(sprintf("\nTotale: %d figure\n", length(figs)))
cat("\n── Fine script 08 ───────────────────────────────────────────────\n")
