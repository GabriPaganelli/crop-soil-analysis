# =============================================================================
# 18_figure_principali.R  —  Export figure principali dal fit M-SP-RIRS-MVRE
#
# Produce in output/figures/:
#   fig_04_posterior_rho.pdf       — distribuzioni a posteriori di rho_r e tau_beta_r
#   fig_05_forest_gamma.pdf        — forest plot effetti within (gamma_r)
#   fig_06_forest_beta.pdf         — forest plot effetti between (beta_r)
#   fig_07_ppc.pdf                 — posterior predictive check (SOC, N, P)
#   fig_09_loo_comparison.pdf      — confronto LOO-CV modelli principali
#   fig_10_trace_key.pdf           — trace plots parametri chiave
#   fig_11_posterior_tau_sigma.pdf — tau_alpha_r, tau_beta_r, sigma_r a posteriori
#   fig_12_struct_params_panel.pdf — riassunto parametri strutturali
#   fig_18_cross_corr.pdf          — correlazioni cross-risposta (MVRE)
#
# Dipende da: stan/fit_msp_rirs_mvre.rds, data/dati.rds
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
dir.create(here("output", "tables"),  recursive = TRUE, showWarnings = FALSE)
fig_dir   <- here("output", "figures")
cache_dir <- here("output", "cache")
tab_dir   <- here("output", "tables")

save_fig <- function(fname, p, w = 16, h = 9, u = "cm") {
  ggplot2::ggsave(file.path(fig_dir, fname), plot = p,
                  width = w, height = h, units = u, device = "pdf")
  cat(sprintf("  [fig] Salvato: %s\n", fname))
}

cov_W_labels <- c("logBottom", "Texture1 (ILR1)", "Texture2 (ILR2)", "BulkDensity", "pH")
cov_B_labels <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")
resp_levels  <- c("SOC", "N", "P")
resp_colors  <- c("SOC" = "#2166AC", "N" = "#1A9850", "P" = "#D73027")


# ── 1. DATI ───────────────────────────────────────────────────────────────────

cat("Carico dati...\n")
dati <- readRDS(here("data", "dati.rds")) |>
  mutate(across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                ~ as.integer(as.character(.x)))) |>
  mutate(
    logSOC    = log(PercSOC),
    logN      = log(PercTotNitro),
    logP      = log(PercTotPhos),
    logBottom = log(Bottom)
  ) |>
  mutate(across(c(logBottom, Texture1, Texture2, BulkDensity, PH),
                ~ c(scale(.x)))) |>
  mutate(Field = factor(Field))

field_levels <- sort(unique(as.integer(as.character(dati$Field))))
J <- length(field_levels)
N <- nrow(dati)

dati_int <- dati |>
  mutate(field_int = as.integer(factor(as.integer(as.character(Field)),
                                       levels = field_levels)))

X_W_cols <- c("logBottom", "Texture1", "Texture2", "BulkDensity", "PH")
X_B_cols <- c("OnFarm", "Irrigate", "Fertilised", "N_Natural")
X_W      <- as.matrix(dati_int[, X_W_cols])
X_B      <- dati_int |>
  distinct(field_int, across(all_of(X_B_cols))) |>
  arrange(field_int) |>
  select(all_of(X_B_cols)) |>
  as.matrix()
field_id <- dati_int$field_int


# ── 2. CARICA FIT M-SP-RIRS ───────────────────────────────────────────────────

cat("Carico fit_msp_rirs_mvre.rds...\n")
fit <- readRDS(here("stan", "fit_msp_rirs_mvre.rds"))

smry <- fit$summary(c(
  paste0(c("alpha_","tau_alpha_","tau_beta_","rho_","sigma_"), rep(resp_levels, each = 5)),
  paste0("gamma_SOC[", 1:5, "]"),
  paste0("gamma_N[",   1:5, "]"),
  paste0("gamma_P[",   1:5, "]"),
  paste0("beta_SOC[",  1:4, "]"),
  paste0("beta_N[",    1:4, "]"),
  paste0("beta_P[",    1:4, "]"),
  "rho_int_SOC_N", "rho_int_SOC_P", "rho_int_N_P",
  "rho_slope_SOC_N", "rho_slope_SOC_P", "rho_slope_N_P"
)) |>
  select(variable, median, mean, sd, q5, q95, rhat, ess_bulk)

cat(sprintf("Fit caricato. Parametri nel summary: %d\n", nrow(smry)))


# ── 3. FIG_04: POSTERIOR rho_r e tau_beta_r ───────────────────────────────────

cat("\n[fig_04] Posterior rho_r e tau_beta_r...\n")

# Pannello sinistro: rho_r (correlazione RI-RS)
draws_rho <- fit$draws(c("rho_SOC", "rho_N", "rho_P"), format = "df") |>
  pivot_longer(cols = starts_with("rho_"),
               names_to  = "par",
               values_to = "valore") |>
  mutate(risposta = factor(sub("rho_", "", par), levels = resp_levels))

ci_rho <- draws_rho |>
  group_by(risposta) |>
  summarise(med = median(valore), lo = quantile(valore, .05), hi = quantile(valore, .95))

p_rho <- ggplot(draws_rho, aes(x = valore, fill = risposta, color = risposta)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_segment(data = ci_rho, aes(x = lo, xend = hi, y = -Inf, yend = -Inf, color = risposta),
               linewidth = 3, alpha = 0.8, inherit.aes = FALSE) +
  geom_point(data = ci_rho, aes(x = med, y = 0, color = risposta), size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = resp_colors, guide = "none") +
  scale_color_manual(values = resp_colors, guide = "none") +
  facet_wrap(~ risposta, scales = "free", nrow = 1) +
  labs(
    title    = expression("Correlazione RI–RS " * rho[r]),
    subtitle = expression(rho[r] %~~% "+1  " %->%  " M-SP confermato   |   " * rho[r] %~~% " 0  " %->% "  slope indip. da intercetta"),
    x        = expression(rho[r]),
    y        = "Densità"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 12))

# Pannello destro: tau_beta_r (SD delle slope tra campi)
draws_tb <- fit$draws(c("tau_beta_SOC","tau_beta_N","tau_beta_P"), format = "df") |>
  pivot_longer(cols = starts_with("tau_beta_"),
               names_to  = "par",
               values_to = "valore") |>
  mutate(risposta = factor(sub("tau_beta_", "", par), levels = resp_levels))

ci_tb <- draws_tb |>
  group_by(risposta) |>
  summarise(med = median(valore), lo = quantile(valore, .05), hi = quantile(valore, .95))

p_tb <- ggplot(draws_tb, aes(x = valore, fill = risposta, color = risposta)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  geom_segment(data = ci_tb, aes(x = lo, xend = hi, y = -Inf, yend = -Inf, color = risposta),
               linewidth = 3, alpha = 0.8, inherit.aes = FALSE) +
  geom_point(data = ci_tb, aes(x = med, y = 0, color = risposta), size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = resp_colors, guide = "none") +
  scale_color_manual(values = resp_colors, guide = "none") +
  facet_wrap(~ risposta, scales = "free", nrow = 1) +
  labs(
    title    = expression("SD slope tra campi " * tau[beta[r]]),
    subtitle = expression(tau[beta[r]] %~~% " 0  " %->% "  M-RI: solo intercetta casuale"),
    x        = expression(tau[beta[r]]),
    y        = "Densità"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 12))

p_rho_combined <- p_rho / p_tb +
  plot_annotation(title = "Struttura bivariata M-SP-RIRS-MVRE: correlazione RI–RS e variabilità slope")
print(p_rho_combined)
save_fig("fig_04_posterior_rho.pdf", p_rho_combined, w = 22, h = 14)


# ── 4. FIG_05: FOREST PLOT gamma_r (effetti within-field) ─────────────────────

cat("\n[fig_05] Forest plot gamma_r...\n")

gamma_data <- smry |>
  filter(grepl("^gamma_", variable)) |>
  mutate(
    risposta = factor(
      case_when(grepl("SOC", variable) ~ "SOC", grepl("_N\\[", variable) ~ "N",
                grepl("_P\\[", variable) ~ "P"),
      levels = resp_levels
    ),
    k         = as.integer(stringr::str_extract(variable, "\\d+")),
    covariata = factor(cov_W_labels[k], levels = rev(cov_W_labels))
  )

p_gamma <- ggplot(gamma_data, aes(x = median, y = covariata, color = risposta)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_linerange(aes(xmin = q5, xmax = q95), linewidth = 1.5, alpha = 0.7,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = resp_colors, name = "Risposta") +
  labs(
    title    = expression("Effetti fissi within-field " * (gamma[r])),
    subtitle = "Mediana e IC 90% posteriori | Scala standardizzata (M-SP-RIRS-MVRE)",
    x        = "Effetto stimato (scala standardizzata)",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

print(p_gamma)
save_fig("fig_05_forest_gamma.pdf", p_gamma, w = 14, h = 10)


# ── 5. FIG_06: FOREST PLOT beta_r ────────────────────────────────────────────

cat("\n[fig_06] Forest plot beta_r...\n")

beta_data <- smry |>
  filter(grepl("^beta_", variable)) |>
  mutate(
    risposta = factor(
      case_when(grepl("SOC", variable) ~ "SOC", grepl("_N\\[", variable) ~ "N",
                grepl("_P\\[", variable) ~ "P"),
      levels = resp_levels
    ),
    k         = as.integer(stringr::str_extract(variable, "\\d+")),
    covariata = factor(cov_B_labels[k], levels = rev(cov_B_labels))
  )

p_beta <- ggplot(beta_data, aes(x = median, y = covariata, color = risposta)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_linerange(aes(xmin = q5, xmax = q95), linewidth = 1.5, alpha = 0.7,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = resp_colors, name = "Risposta") +
  labs(
    title    = expression("Effetti fissi between-field — management " * (beta[r])),
    subtitle = "Mediana e IC 90% posteriori (M-SP-RIRS-MVRE)",
    x        = "Effetto stimato", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

print(p_beta)
save_fig("fig_06_forest_beta.pdf", p_beta, w = 14, h = 9)

p_fixed_effects <- (p_gamma + labs(title = expression(gamma[r] ~ "(within-field)"))) /
                   (p_beta  + labs(title = expression(beta[r]  ~ "(between-field)"))) +
  plot_annotation(title = "Effetti fissi del modello M-SP-RIRS-MVRE",
                  subtitle = "Mediana e IC 90% posteriori")
save_fig("fig_05b_forest_gamma_beta_combined.pdf", p_fixed_effects, w = 14, h = 16)


# ── 6. FIG_07: POSTERIOR PREDICTIVE CHECK ─────────────────────────────────────

cat("\n[fig_07] Posterior predictive check...\n")

draws_raw  <- fit$draws(format = "matrix")
draws_full <- matrix(as.double(draws_raw), nrow = nrow(draws_raw), ncol = ncol(draws_raw),
                     dimnames = dimnames(draws_raw))
D_total   <- nrow(draws_full)
D_ppc     <- 200L
thin_idx  <- round(seq(1, D_total, length.out = D_ppc))
draws_ppc <- draws_full[thin_idx, , drop = FALSE]
rm(draws_raw, draws_full); gc()

extract_V_row_ppc <- function(draws_d, row) {
  cols <- paste0("V[", row, ",", seq_len(J), "]")
  draws_d[, cols, drop = FALSE]
}
extract_mat_ppc <- function(draws_d, prefix, K) {
  cols <- paste0(prefix, "[", seq_len(K), "]")
  draws_d[, cols, drop = FALSE]
}

compute_yrep <- function(draws_d, r) {
  row_int <- c(SOC = 1L, N = 3L, P = 5L)[r]
  row_slo <- c(SOC = 2L, N = 4L, P = 6L)[r]
  alpha   <- draws_d[, paste0("alpha_", r)]
  sigma   <- draws_d[, paste0("sigma_", r)]
  u_int   <- extract_V_row_ppc(draws_d, row_int)  # D × J
  u_slope <- extract_V_row_ppc(draws_d, row_slo)  # D × J
  gamma   <- extract_mat_ppc(draws_d, paste0("gamma_", r), length(X_W_cols))
  beta    <- extract_mat_ppc(draws_d, paste0("beta_",  r), length(X_B_cols))

  X_B_n      <- X_B[field_id, ]
  rand_int   <- u_int[,   field_id]
  rand_slope <- u_slope[, field_id] * outer(rep(1, D_ppc), X_W[, 1])
  fix_within <- tcrossprod(gamma, X_W)
  fix_betw   <- tcrossprod(beta,  X_B_n)
  intercept  <- outer(alpha, rep(1, N))
  mu         <- intercept + rand_int + rand_slope + fix_within + fix_betw

  eps  <- matrix(rnorm(D_ppc * N), nrow = D_ppc, ncol = N)
  yrep <- mu + outer(sigma, rep(1, N)) * eps
  yrep
}

set.seed(2024)
yrep_SOC <- compute_yrep(draws_ppc, "SOC")
yrep_N   <- compute_yrep(draws_ppc, "N")
yrep_P   <- compute_yrep(draws_ppc, "P")
rm(draws_ppc); gc()

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
    title    = "Posterior predictive check — M-SP-RIRS-MVRE",
    subtitle = paste0("Linea scura = dati | ", D_ppc, " curve = repliche dalla posteriore")
  )
print(p_ppc)
save_fig("fig_07_ppc.pdf", p_ppc, w = 20, h = 8)
rm(yrep_SOC, yrep_N, yrep_P); gc()


# ── 7. FIG_09: LOO-CV COMPARISON ──────────────────────────────────────────────

cat("\n[fig_09] LOO comparison...\n")

loo_cache_path <- file.path(cache_dir, "loo_results_mvre.rds")

if (file.exists(loo_cache_path)) {
  loo_list <- readRDS(loo_cache_path)
} else {
  cat("  Calcolo LOO modelli principali...\n")
  loo_list <- list()
  ll_mvre <- fit$draws("log_lik", format = "matrix")
  loo_list[["M-SP-RIRS-MVRE"]] <- loo(ll_mvre, cores = 4)
  rm(ll_mvre); gc()

  model_map <- c("M-SP" = "fit_msp", "M-RI" = "fit_mri",
                  "M-GPS" = "fit_mgps", "M-GP" = "fit_mgp")
  for (label in names(model_map)) {
    p <- here("stan", paste0(model_map[[label]], ".rds"))
    if (file.exists(p)) {
      cat(sprintf("  LOO %s...\n", label))
      f <- readRDS(p)
      loo_list[[label]] <- loo(f$draws("log_lik", format = "matrix"), cores = 4)
      rm(f); gc()
    }
  }
  saveRDS(loo_list, loo_cache_path)
}

if (length(loo_list) > 1) {
  cmp    <- loo_compare(loo_list)
  cmp_df <- as.data.frame(cmp) |>
    tibble::rownames_to_column("modello") |>
    mutate(modello = factor(modello, levels = rev(rownames(cmp))))

  p_loo <- ggplot(cmp_df, aes(x = elpd_diff, y = modello)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    geom_linerange(aes(xmin = elpd_diff - 2*se_diff, xmax = elpd_diff + 2*se_diff),
                   linewidth = 2, color = "steelblue", alpha = 0.65) +
    geom_point(size = 4, color = "steelblue") +
    geom_text(aes(label = sprintf("%.1f ± %.1f", elpd_diff, se_diff)),
              hjust = -0.15, size = 3.2, color = "grey30") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
    labs(title    = "Confronto LOO-CV tra i modelli",
         subtitle = "ΔELPD rispetto al migliore | Barre = ±2 SE",
         x = "ΔELPD", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank())

  print(p_loo)
  save_fig("fig_09_loo_comparison.pdf", p_loo, w = 15, h = 8)
}


# ── 8. FIG_10: TRACE PLOTS ────────────────────────────────────────────────────

cat("\n[fig_10] Trace plots...\n")

params_trace <- c("alpha_SOC", "tau_alpha_SOC", "tau_beta_SOC", "rho_SOC",
                  "alpha_N",   "tau_alpha_N",   "tau_beta_N",   "rho_N",
                  "alpha_P",   "tau_alpha_P",   "tau_beta_P",   "rho_P")

p_trace <- mcmc_trace(fit$draws(variables = params_trace), facet_args = list(ncol = 4)) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold")) +
  labs(title    = "Trace plots — parametri strutturali di M-SP-RIRS-MVRE",
       subtitle = "4 catene × 5000 iterazioni")
print(p_trace)
save_fig("fig_10_trace_key.pdf", p_trace, w = 24, h = 14)


# ── 9. FIG_11: POSTERIOR tau_r e sigma_r ──────────────────────────────────────

cat("\n[fig_11] Posterior tau_alpha_r, tau_beta_r, sigma_r...\n")

draws_tau_sigma <- fit$draws(
  c("tau_alpha_SOC","tau_alpha_N","tau_alpha_P",
    "tau_beta_SOC", "tau_beta_N", "tau_beta_P",
    "sigma_SOC","sigma_N","sigma_P"),
  format = "df"
) |>
  pivot_longer(cols = -c(.chain, .iteration, .draw),
               names_to  = "par", values_to = "valore") |>
  mutate(
    tipo = case_when(
      str_starts(par, "tau_alpha_") ~ "tau_alpha[r]",
      str_starts(par, "tau_beta_")  ~ "tau_beta[r]",
      TRUE                          ~ "sigma[r]"
    ),
    risposta = factor(str_replace(par, "^(tau_alpha_|tau_beta_|sigma_)", ""),
                      levels = resp_levels)
  )

p_tau_sigma <- ggplot(draws_tau_sigma, aes(x = valore, fill = risposta, color = risposta)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  scale_fill_manual(values  = resp_colors, name = "Risposta") +
  scale_color_manual(values = resp_colors, name = "Risposta") +
  facet_grid(risposta ~ tipo, scales = "free", labeller = label_parsed) +
  labs(
    title    = expression("Distribuzioni a posteriori di " * tau[alpha[r]] *
                          ", " * tau[beta[r]] * " e " * sigma[r]),
    x = "Valore", y = "Densità"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
print(p_tau_sigma)
save_fig("fig_11_posterior_tau_sigma.pdf", p_tau_sigma, w = 16, h = 12)


# ── 10. FIG_12: PANEL PARAMETRI STRUTTURALI ────────────────────────────────────

cat("\n[fig_12] Panel strutturale...\n")

struct_params <- smry |>
  filter(variable %in% c(
    paste0(c("alpha_","tau_alpha_","tau_beta_","rho_","sigma_"),
           rep(resp_levels, each = 5))
  )) |>
  mutate(
    par      = str_extract(variable, "^[^_]+(?:_[^_]+)?"),
    risposta = factor(str_extract(variable, "[A-Z]+$"), levels = resp_levels),
    par_label = factor(par,
      levels = c("alpha","tau_alpha","tau_beta","rho","sigma"),
      labels = c("alpha[r]","tau[alpha[r]]","tau[beta[r]]","rho[r]","sigma[r]"))
  ) |>
  filter(!is.na(par_label))

p_struct <- ggplot(struct_params, aes(x = median, y = risposta, color = risposta)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_linerange(aes(xmin = q5, xmax = q95), linewidth = 2, alpha = 0.7) +
  geom_point(size = 3.5) +
  scale_color_manual(values = resp_colors, guide = "none") +
  facet_wrap(~ par_label, scales = "free_x", nrow = 1, labeller = label_parsed) +
  labs(title    = "Parametri strutturali del modello M-SP-RIRS-MVRE",
       subtitle = "Mediana e IC 90% posteriori per SOC, N, P",
       x = "Valore stimato", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"), panel.grid.major.y = element_blank())
print(p_struct)
save_fig("fig_12_struct_params_panel.pdf", p_struct, w = 26, h = 7)


# ── 10b. FIG_18: CROSS-RESPONSE CORRELATIONS (esclusivo di MVRE) ─────────────

cat("\n[fig_18] Cross-response correlations (MVRE)...\n")

cross_params <- c(
  "rho_int_SOC_N", "rho_int_SOC_P", "rho_int_N_P",
  "rho_slope_SOC_N", "rho_slope_SOC_P", "rho_slope_N_P"
)
cross_labels <- c(
  "rho[int]~SOC-N", "rho[int]~SOC-P", "rho[int]~N-P",
  "rho[slope]~SOC-N", "rho[slope]~SOC-P", "rho[slope]~N-P"
)
cross_tipo <- c(rep("Intercepts (rho[int])", 3), rep("Slopes (rho[slope])", 3))

cross_smry <- smry |>
  filter(variable %in% cross_params) |>
  left_join(
    data.frame(variable = cross_params, label = cross_labels, tipo = cross_tipo,
               stringsAsFactors = FALSE),
    by = "variable"
  ) |>
  mutate(label = factor(label, levels = rev(cross_labels)),
         tipo  = factor(tipo, levels = c("Intercepts (rho[int])", "Slopes (rho[slope])")))

if (nrow(cross_smry) > 0) {
  p_cross <- ggplot(cross_smry, aes(x = median, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    geom_linerange(aes(xmin = q5, xmax = q95), linewidth = 2, color = "steelblue", alpha = 0.7) +
    geom_point(size = 3.5, color = "steelblue") +
    geom_text(aes(label = sprintf("%.3f [%.3f,%.3f]", median, q5, q95)),
              hjust = -0.1, size = 2.9, color = "grey30") +
    facet_wrap(~ tipo, scales = "free_y", labeller = label_parsed) +
    scale_x_continuous(limits = c(-0.8, 1.2), expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_discrete(labels = function(x) parse(text = x)) +
    labs(
      title    = expression("Correlazioni cross-risposta — " * Omega[6] * " (M-SP-RIRS-MVRE)"),
      subtitle = expression("Mediana e IC 90% | " * rho[int]^"SOC-N" * " = +0.39 interpretabile biologicamente"),
      x        = "Correlazione posteriore", y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title  = element_text(face = "bold"),
      strip.text  = element_text(face = "bold"),
      panel.grid.major.y = element_blank()
    )
  print(p_cross)
  save_fig("fig_18_cross_corr.pdf", p_cross, w = 20, h = 10)
}


# ── 11. EXPORT TABELLE CSV ────────────────────────────────────────────────────

cat("\n[tabelle] Export CSV...\n")

smry_export <- smry |> mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(smry_export, file.path(tab_dir, "tab_posterior_summary.csv"), row.names = FALSE)

tab_struct <- smry |>
  filter(grepl("^(alpha_|tau_alpha_|tau_beta_|rho_|sigma_)", variable)) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(tab_struct, file.path(tab_dir, "tab_struct_params.csv"), row.names = FALSE)
cat("  Strutturali: tab_struct_params.csv\n")
print(tab_struct |> select(variable, median, sd, q5, q95), n = Inf)

tab_gamma <- smry |>
  filter(grepl("^gamma_", variable)) |>
  mutate(
    risposta  = factor(case_when(grepl("SOC",variable)~"SOC",grepl("_N\\[",variable)~"N",
                                 grepl("_P\\[",variable)~"P"), levels = resp_levels),
    k         = as.integer(stringr::str_extract(variable, "\\d+")),
    covariata = cov_W_labels[k]
  ) |>
  select(risposta, covariata, median, sd, q5, q95) |>
  arrange(risposta, covariata) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(tab_gamma, file.path(tab_dir, "tab_gamma_within.csv"), row.names = FALSE)
cat("  Within-field: tab_gamma_within.csv\n")

tab_beta <- smry |>
  filter(grepl("^beta_", variable)) |>
  mutate(
    risposta  = factor(case_when(grepl("SOC",variable)~"SOC",grepl("_N\\[",variable)~"N",
                                 grepl("_P\\[",variable)~"P"), levels = resp_levels),
    k         = as.integer(stringr::str_extract(variable, "\\d+")),
    covariata = cov_B_labels[k]
  ) |>
  select(risposta, covariata, median, sd, q5, q95) |>
  arrange(risposta, covariata) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(tab_beta, file.path(tab_dir, "tab_beta_between.csv"), row.names = FALSE)
cat("  Between-field: tab_beta_between.csv\n")

cat("\n── Fine script 18 ──────────────────────────────────────────────\n")
