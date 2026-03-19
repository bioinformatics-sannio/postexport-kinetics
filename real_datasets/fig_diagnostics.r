# =============================================================================
# Title: Diagnostic Plots and Bias Audits for Real-Dataset Results
# Description:
#   This script produces a set of diagnostic summaries for the real-dataset
#   nested-test results, including:
#     - distribution of the fit-improvement statistic T.obs,
#     - DeltaAIC versus estimated Sigma,
#     - QQ plots of p-values by dataset,
#     - bias audits against abundance and zero inflation,
#     - top-ranked gene tables per dataset.
#
#   The goal is to assess calibration, effect-size structure, and possible
#   dependencies of the model-selection signal on coverage or sparsity.
#
# Intended use:
#   Diagnostic visualization and supplementary table generation for the
#   real-data analysis.
#
# Copyright (c) 2026 Luigi Cerulo
#
# Permission is hereby granted to use, copy, modify, and distribute this
# software for academic and research purposes, provided that this notice is
# retained in all copies.
#
# This software is provided "as is", without warranty of any kind, express
# or implied, including but not limited to the warranties of merchantability,
# fitness for a particular purpose, and noninfringement. In no event shall
# the authors be liable for any claim, damages, or other liability arising
# from, out of, or in connection with the software or its use.
# =============================================================================


setwd("~/postexport-kinetics/real_datasets")
require(data.table)
require(ggplot2)

load("results_realdatasets.rdata", verbose = TRUE)


# =============================================================================
# Distribution of T.obs
# =============================================================================

# -----------------------------------------------------------------------------
# Prepare a log-transformed version of the observed fit-improvement statistic.
#
# T.obs is expected to be non-negative, so a small epsilon is added before
# taking log10.
# -----------------------------------------------------------------------------
eps <- 1e-8
results[!is.na(T.obs) & T.obs >= 0, logT := log10(T.obs + eps)]


# -----------------------------------------------------------------------------
# Dataset-level summary statistics for T.obs:
#   - sample size,
#   - probability that T.obs > 0,
#   - median,
#   - 95th percentile,
#   - 99th percentile.
# -----------------------------------------------------------------------------
T_stats <- results[!is.na(T.obs) & T.obs >= 0,
  .(
    n = .N,
    pr_pos = mean(T.obs > 0),
    q50 = quantile(T.obs, 0.50),
    q95 = quantile(T.obs, 0.95),
    q99 = quantile(T.obs, 0.99)
  ),
  by = dataset
]

# Convert key summary quantiles to log scale for plotting.
T_stats[, `:=`(
  v_med = log10(q50 + eps),
  v_95  = log10(q95 + eps),
  v_99  = log10(q99 + eps),
  ann_txt = sprintf(
    "n=%d\nPr(T>0)=%.3f\nmedian=%.3g\n95%%=%.3g\n99%%=%.3g",
    n, pr_pos, q50, q95, q99
  )
)]

# Determine per-panel label positions.
lab_pos_T <- results[!is.na(logT),
  .(
    x_lab = quantile(logT, 0.70),
    y_lab = Inf
  ),
  by = dataset
]

T_stats <- lab_pos_T[T_stats, on = "dataset"]


# -----------------------------------------------------------------------------
# Faceted histogram of log10(T.obs + eps) by dataset.
#
# Vertical lines show:
#   - median,
#   - 95th percentile,
#   - 99th percentile.
# -----------------------------------------------------------------------------
p_T_facet <- ggplot(results[!is.na(logT)], aes(x = logT)) +
  geom_histogram(
    bins = 60,
    linewidth = 0.2,
    fill = "grey80",
    color = "grey30"
  ) +
  geom_vline(data = T_stats, aes(xintercept = v_med), linetype = 2, linewidth = 0.5) +
  geom_vline(data = T_stats, aes(xintercept = v_95),  linetype = 3, linewidth = 0.5) +
  geom_vline(data = T_stats, aes(xintercept = v_99),  linetype = 1, linewidth = 0.5) +
  geom_label(
    data = T_stats,
    aes(x = x_lab, y = y_lab, label = ann_txt),
    hjust = 0,
    vjust = 1.1,
    size = 3.1,
    linewidth = 0.25
  ) +
  facet_wrap(~ dataset, scales = "free_y") +
  labs(
    title = "Distribution of fit improvement statistic by dataset",
    subtitle = expression(T[obs] == RSS[0] - RSS[1] ~ " (full vs null)"),
    x = expression(log[10](T[obs] + epsilon)),
    y = "Count"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_T_facet


# =============================================================================
# DeltaAIC versus Sigma
# =============================================================================

# -----------------------------------------------------------------------------
# Scatter plot of DeltaAIC against the estimated Sigma value.
#
# Horizontal reference lines at:
#   - 0
#   - 2
#   - 10
# are commonly used as rough model-comparison benchmarks.
# -----------------------------------------------------------------------------
p_simple_aic <- ggplot(results, aes(x = Sigma, y = deltaAIC)) +
  geom_hline(yintercept = 0,  linetype = 2, linewidth = 0.4) +
  geom_hline(yintercept = 2,  linetype = 3, linewidth = 0.4) +
  geom_hline(yintercept = 10, linetype = 3, linewidth = 0.4) +
  geom_point(size = 0.7, alpha = 0.6) +
  facet_wrap(~ dataset, ncol = 2, scales = "free_x") +
  labs(
    x = expression(hat(sigma)[c]),
    y = expression(Delta * AIC)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

p_simple_aic

ggsave(
  p_simple_aic,
  filename = "aic_by_dataset.pdf",
  width = 10,
  height = 8
)


# =============================================================================
# QQ plot of p-values by dataset
# =============================================================================

# -----------------------------------------------------------------------------
# Build dataset-specific QQ-plot data:
#   - observed p-values,
#   - expected p-values under Uniform(0,1),
#   - log10 transforms,
#   - Beta-based 95% envelope for order statistics.
# -----------------------------------------------------------------------------
qq_dt <- results[
  !is.na(p.value) & p.value > 0 & p.value <= 1,
  .(p = p.value),
  by = dataset
]

setorder(qq_dt, dataset, p)

# Expected p-values for each order statistic within each dataset.
qq_dt[, exp_p := (seq_len(.N) - 0.5) / .N, by = dataset]
qq_dt[, obs_logp := -log10(p)]
qq_dt[, exp_logp := -log10(exp_p)]

# Beta confidence band for ordered Uniform(0,1) p-values.
qq_dt[, i := seq_len(.N), by = dataset]
qq_dt[, m := .N, by = dataset]

alpha <- 0.05
qq_dt[, lo := qbeta(alpha / 2, i, m - i + 1)]
qq_dt[, hi := qbeta(1 - alpha / 2, i, m - i + 1)]
qq_dt[, lo_logp := -log10(lo)]
qq_dt[, hi_logp := -log10(hi)]


# -----------------------------------------------------------------------------
# Per-dataset facet labels:
#   - total number of tests,
#   - number of p-values < 0.05,
#   - number of p-values < 0.01.
# -----------------------------------------------------------------------------
lab_dt <- qq_dt[, .(
  n = .N,
  n_p05 = sum(p < 0.05),
  n_p01 = sum(p < 0.01),
  x_lab = quantile(exp_logp, 0.05),
  y_lab = quantile(obs_logp, 0.95)
), by = dataset]

lab_dt[, label := sprintf("n=%d\np<0.05: %d\np<0.01: %d", n, n_p05, n_p01)]

# Optional extension:
# If deltaAIC is available and desired, counts above a threshold could also
# be appended to the label.


# -----------------------------------------------------------------------------
# Faceted QQ plot:
#   - ribbon = 95% null envelope,
#   - points = observed ordered p-values,
#   - dashed line = exact null expectation.
# -----------------------------------------------------------------------------
p_qq_facet <- ggplot(qq_dt, aes(x = exp_logp, y = obs_logp)) +
  geom_ribbon(aes(ymin = lo_logp, ymax = hi_logp), alpha = 0.2) +
  geom_point(size = 0.6, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
  geom_label(
    data = lab_dt,
    aes(x = x_lab, y = y_lab, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -2,
    size = 3.1,
    linewidth = 0.25
  ) +
  facet_wrap(~ dataset, ncol = 2, scales = "fixed") +
  labs(
    x = "Expected -log10(p) under Uniform(0,1)",
    y = "Observed -log10(p)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

p_qq_facet

ggsave(
  p_qq_facet,
  filename = "qq_by_dataset.pdf",
  width = 10,
  height = 8
)

