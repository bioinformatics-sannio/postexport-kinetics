# =============================================================================
# Title: Calibration, Power, and Decision-Threshold Analysis on Synthetic Data
# Description:
#   This script evaluates inferential performance of the benchmarked methods
#   on the synthetic dataset, focusing on:
#     - calibration under the null hypothesis (type I error control),
#     - power under the alternative hypothesis,
#     - decision-level precision/recall trade-offs at different q-value cutoffs.
#
#   Main outputs:
#     - calibration curves,
#     - supplementary calibration summary tables,
#     - power curves as a function of true effect size,
#     - combined calibration/power figures,
#     - decision plots across significance thresholds.
#
# Intended use:
#   Benchmark diagnostics and figure generation for method validation.
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


# -----------------------------------------------------------------------------
# Set working directory and load required packages.
# -----------------------------------------------------------------------------
setwd("~/postexport-kinetics/synthetic_dataset")

require(data.table)
library(ggplot2)
library(scales)

# Load benchmark results.
load("test_sigma_2k_new.rdata")



# =============================================================================
# FIGURE 1A
# Calibration plot of the test under H0 (type I error control)
# =============================================================================

# -----------------------------------------------------------------------------
# Filter null genes and select the benchmark scenario used for calibration.
#
# Current choice:
#   - Nested-wild_whitened
#   - true negative genes only (Positive == 0)
#   - Lambda = 0.3
#   - Tsteps = 10
# -----------------------------------------------------------------------------
dt0 <- dt_test_results[
  Test_method %in% c("Nested-wild_whitened") &
  Positive == 0 &
  Lambda == 0.3 &
  Tsteps == 10
]

# Optional relabeling for plotting.
# Note: the factor levels here refer to older naming conventions and may be
# retained for backward compatibility with earlier plots/scripts.
dt0[, Test_method := factor(
  Test_method,
  levels = c("Nested-wild-rademacher", "Nested-wild-mammen"),
  labels = c("Nested test (r)", "Nested test (m)")
)]

# Remove non-finite p-values before calibration summaries.
dt0 <- dt0[is.finite(p.value)]

# Quick sanity checks.
dt0[, .N, by = .(Test_method)]
dt0[, .(pmin = min(p.value), pmax = max(p.value)), by = Test_method]

# Nominal alpha levels used for summary size evaluation.
alpha_grid <- c(0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.001)


# -----------------------------------------------------------------------------
# Compute empirical type I error ("size") at selected alpha levels.
#
# For each design stratum:
#   size(alpha) = P(p <= alpha | H0)
# -----------------------------------------------------------------------------
size_by_group <- function(DT, alphas = alpha_grid) {
  DT[, {
    pv <- p.value
    out <- data.table(alpha = alphas)
    out[, size := sapply(alpha, function(a) mean(pv <= a, na.rm = TRUE))]
    out
  }, by = .(
    Test_method,
    N_replicates, N_tsamples, Tsteps,
    Perturbation, Exprs_noise, Platform
  )]
}

dt_size <- size_by_group(dt0, alpha_grid)

# Deviation from perfect calibration, where empirical size should equal alpha.
dt_size[, diff := size - alpha]

# Identify worst calibration deviations at alpha = 0.05.
dt_size[alpha == 0.05][order(-abs(diff))][1:30]


# -----------------------------------------------------------------------------
# Compute the full calibration curve over a dense alpha grid.
# -----------------------------------------------------------------------------
calibration_curve <- function(DT, alphas = seq(0, 1, by = 0.01)) {
  DT[, {
    pv <- p.value
    data.table(
      alpha = alphas,
      size = sapply(alphas, function(a) mean(pv <= a))
    )
  }, by = .(
    Test_method,
    N_replicates, N_tsamples, Tsteps,
    Perturbation, Exprs_noise, Platform
  )]
}

dt_cal <- calibration_curve(dt0)


# -----------------------------------------------------------------------------
# Select the main scenario to display in the calibration figure.
#
# Current filter:
#   - N_replicates = 5
#   - Tsteps = 10
#   - Exprs_noise = Medium
# -----------------------------------------------------------------------------
dt_cal_sub <- dt_cal[
  N_replicates == 5 &
  Tsteps == 10 &
  Exprs_noise == "Medium"
]

# If needed, additional stratification by N_replicates can be enabled here.


# -----------------------------------------------------------------------------
# Summarize calibration across remaining scenarios using:
#   - median curve,
#   - interquartile range,
#   - 5% to 95% envelope.
# -----------------------------------------------------------------------------
dt_band <- dt_cal_sub[, .(
  N = .N,
  med = median(size, na.rm = TRUE),
  q25 = quantile(size, 0.25, na.rm = TRUE),
  q75 = quantile(size, 0.75, na.rm = TRUE),
  q05 = quantile(size, 0.05, na.rm = TRUE),
  q95 = quantile(size, 0.95, na.rm = TRUE)
), by = .(
  Test_method,
  N_tsamples, Tsteps,
  Perturbation, Exprs_noise, Platform, alpha
)]

dt_band[, N_tsamples := factor(N_tsamples, levels = c(3, 5, 10, 20))]

require(latex2exp)
require(scales)

# -----------------------------------------------------------------------------
# Calibration plot:
#   - dotted diagonal = ideal calibration,
#   - ribbons = variability across scenarios,
#   - lines = median empirical size.
# -----------------------------------------------------------------------------
p <- ggplot(dt_band, aes(x = alpha, group = N_tsamples)) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dotted", alpha = 0.3, linewidth = 0.8
  ) +
  geom_ribbon(
    aes(ymin = q05, ymax = q95, fill = N_tsamples),
    alpha = 0.08, colour = NA
  ) +
  geom_ribbon(
    aes(ymin = q25, ymax = q75, fill = N_tsamples),
    alpha = 0.12, colour = NA
  ) +
  geom_line(aes(y = med, color = N_tsamples), linewidth = 0.9) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    limits = c(0, 1),
    breaks = c(0.001, 0.05, 0.1, 0.2, 0.4),
    labels = c("0.001", "0.05", "0.1", "0.2", "0.4")
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    limits = c(0, 1),
    breaks = c(0.001, 0.05, 0.1, 0.2, 0.4),
    labels = c("0.001", "0.05", "0.1", "0.2", "0.4")
  ) +
  facet_grid(Platform ~ Perturbation) +
  labs(
    x = TeX("Nominal $\\; \\alpha$"),
    y = TeX("Empirical $\\; P(p \\leq \\alpha \\;|\\; H_0)$"),
    color = "# time points",
    fill  = "# time points"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.width = unit(1.6, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

p

ggsave("calibration-curve.pdf", p, width = 8, height = 12)


# =============================================================================
# Supplementary table:
# Calibration summaries at alpha = 0.05
# =============================================================================

library(data.table)

alpha0 <- 0.05

# -----------------------------------------------------------------------------
# Summarize empirical size at nominal alpha = 0.05 across simulation settings.
#
# Reported quantities:
#   - mean empirical size,
#   - standard deviation across design configurations,
#   - 5% and 95% quantiles,
#   - worst-case size,
#   - deviation from nominal alpha,
#   - inflation factor = mean size / alpha0.
# -----------------------------------------------------------------------------
dt_suppl <- dt_size[alpha == alpha0 & Tsteps == 10,
  .(
    n = .N,
    Mean = round(mean(size, na.rm = TRUE), 3),
    SD = round(sd(size, na.rm = TRUE), 3),
    q05 = round(quantile(size, 0.05, na.rm = TRUE), 3),
    q95 = round(quantile(size, 0.95, na.rm = TRUE), 3),
    Worst = round(max(size, na.rm = TRUE), 3),
    Deviation = round(mean(size, na.rm = TRUE) - alpha0, 3),
    Inflation = round(mean(size, na.rm = TRUE) / alpha0, 2)
  ),
  by = .(Test_method, Exprs_noise, Platform, Perturbation)
][order(Test_method, Exprs_noise, Platform, Perturbation)]

dt_suppl

# Minimal version for manuscript/supplementary reporting.
dt_suppl_min <- dt_suppl[, .(
  Platform,
  Exprs_noise,
  Perturbation,
  Mean,
  SD,
  Inflation
)]

# LaTeX table output.
kbl <- knitr::kable(
  dt_suppl_min,
  format = "latex",
  booktabs = TRUE,
  longtable = FALSE,
  escape = TRUE,
  digits = 3,
  caption = "Empirical type~I error rates at nominal level $\\alpha = 0.05$ across simulation settings. 
Mean size, standard deviation across design configurations, and inflation factor (Mean/$\\alpha$) are reported.",
  label = "tab:type1_supp",
  align = "llllrrr"
)

kbl



# =============================================================================
# FIGURE 1B
# Power curve as a function of alpha and true sigma_c
# =============================================================================

# -----------------------------------------------------------------------------
# Wilson confidence interval for a binomial proportion.
#
# Used to summarize uncertainty in empirical power estimates.
# -----------------------------------------------------------------------------
wilson_ci <- function(x, n, z = 1.96) {
  if (n == 0) return(c(lo = NA_real_, hi = NA_real_))

  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2*n)) / denom
  half <- (z * sqrt((p*(1-p) + z^2/(4*n)) / n)) / denom

  c(lo = max(0, center - half), hi = min(1, center + half))
}


# -----------------------------------------------------------------------------
# Filter positive genes and the scenario used for power analysis.
#
# Current choice:
#   - Tsteps = 10
#   - Lambda = 0.7
#   - Nested-wild_whitened
#   - SHUTOFF perturbation
#   - Positive genes only
# -----------------------------------------------------------------------------
w_dt <- dt_test_results[
  Tsteps == 10 &
  Lambda == 0.7 &
  Test_method == "Nested-wild_whitened" &
  Perturbation == "SHUTOFF" &
  Positive == 1L &
  is.finite(p.value) &
  is.finite(sigma_true_median)
]

# -----------------------------------------------------------------------------
# Bin genes by the true median sigma_c using quantile bins.
#
# This allows power to be visualized as a smooth function of effect size.
# -----------------------------------------------------------------------------
K <- 10
brks <- unique(quantile(
  w_dt$sigma_true_median,
  probs = seq(0, 1, length.out = K + 1),
  na.rm = TRUE
))

if (length(brks) < 4) {
  stop("Not enough unique sigma_true_median quantiles to bin. Reduce K or check sigma_true_median.")
}

w_dt[, sigma_bin := cut(sigma_true_median, breaks = brks, include.lowest = TRUE)]

# Representative x-coordinate per bin: median true sigma within the bin.
w_dt[, sigma_x := median(sigma_true_median),
     by = .(N_tsamples, N_replicates, Platform, Exprs_noise, sigma_bin)]


# -----------------------------------------------------------------------------
# Compute empirical power and Wilson confidence intervals at a given alpha.
# -----------------------------------------------------------------------------
make_pow <- function(alpha0) {
  dt <- w_dt[, .(
    rej = sum(p.value <= alpha0),
    n = .N,
    sigma = unique(sigma_x)
  ), by = .(N_tsamples, N_replicates, Platform, Exprs_noise, sigma_bin)]

  dt[, power := rej / n]

  dt[, c("lo", "hi") := {
    ci <- wilson_ci(rej, n)
    .(ci["lo"], ci["hi"])
  }, by = .(N_tsamples, N_replicates, Platform, Exprs_noise, sigma_bin)]

  dt[, alpha := factor(sprintf("%.2f", alpha0), levels = c("0.01", "0.05"))]
  dt[]
}

dt_pow <- rbindlist(list(make_pow(0.01), make_pow(0.05)), use.names = TRUE)

# Order rows by sigma for plotting consistency.
setorder(dt_pow, Platform, Exprs_noise, alpha, sigma)


# -----------------------------------------------------------------------------
# Main power plot for one representative design slice:
#   - N_tsamples = 10
#   - N_replicates = 5
# -----------------------------------------------------------------------------
p_power <- ggplot(
  dt_pow[N_tsamples == 10 & N_replicates == 5],
  aes(x = sigma, y = power, color = alpha, group = alpha)
) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.04, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.6) +
  facet_grid(Platform ~ Exprs_noise) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = TeX("True effect size ($\\sigma_c$)"),
    y = TeX("Power $\\; P(p \\leq \\alpha \\;|\\; H_1)$"),
    color = expression(alpha)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_power

ggsave("power-curve.pdf", p_power, width = 16, height = 10)


# -----------------------------------------------------------------------------
# Combine calibration and power figures into one joint panel.
# -----------------------------------------------------------------------------
library(patchwork)

p_combined <- p + p_power +
  plot_layout(ncol = 2, widths = c(1, 2)) &
  theme(legend.position = "top")

p_combined

ggsave("calibration_power_combined.pdf", p_combined, width = 20, height = 10)


# =============================================================================
# Supplementary power analysis
# =============================================================================

# -----------------------------------------------------------------------------
# Plot power curves at alpha = 0.05, fixing Exprs_noise = Medium and comparing
# platforms across numbers of time points and replicates.
# -----------------------------------------------------------------------------
dt_plot = dt_pow[
  alpha == 0.05 &
  Exprs_noise == "Medium"
]

require(latex2exp)

p_power <- ggplot(
  dt_plot,
  aes(x = sigma, y = power, colour = Platform)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = 0.005,
    linewidth = 0.1
  ) +
  facet_grid(
    N_tsamples ~ N_replicates,
    labeller = labeller(
      N_tsamples = function(x) paste0("T=", x),
      N_replicates = function(x) paste0("R=", x)
    )
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    y = TeX("Power $\\; P(p \\leq \\alpha \\;|\\; H_1)$"),
    x = TeX("True effect size ($\\sigma_c$)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_power

ggsave("power_curve_suppl.pdf", p_power, width = 8, height = 10)



# =============================================================================
# Decision plot at different q-value thresholds
# =============================================================================

# -----------------------------------------------------------------------------
# Filter one benchmark scenario for threshold-based decision analysis.
#
# Current choice:
#   - N_tsamples = 10
#   - N_replicates = 5
#   - Tsteps = 10
#   - Lambda = 0.5
#   - Nested-wild_whitened
# -----------------------------------------------------------------------------
w_dt <- dt_test_results[
  N_tsamples == 10 &
  N_replicates == 5 &
  Tsteps == 10 &
  Lambda == 0.5 &
  Test_method == "Nested-wild_whitened"
]

# Optional relabeling for older naming schemes.
w_dt[, Test_method := factor(
  Test_method,
  levels = c("Nested-wild-rademacher", "PSIBaseline"),
  labels = c("Nested test", "PSI Baseline")
)]

# Compute BH-adjusted q-values within each design stratum.
w_dt[, qvalue := p.value]
w_dt[, qvalue := p.adjust(p.value, "BH"),
     by = .(Platform, Exprs_noise, Perturbation,
            Test_method, N_tsamples, Tsteps, N_replicates)]

# Threshold grid used for decision analysis.
qs <- c(0.01, 0.02, 0.05, 0.1, 0.2)

w_dt = w_dt[!is.na(qvalue)]


# -----------------------------------------------------------------------------
# Build a table of binary decisions at each q-value threshold.
# -----------------------------------------------------------------------------
dt_qs_pred <- rbindlist(lapply(seq_along(qs), function(i) {

  qN <- qs[i]

  w_dt[
    , .(
      alpha = qN,
      pred = qvalue <= qN
    ),
    by = .(Platform, Exprs_noise, Perturbation,
           Test_method, N_tsamples, Tsteps, N_replicates, Gene, Positive)
  ]
}))


# -----------------------------------------------------------------------------
# Compute precision and recall from binary predictions.
# -----------------------------------------------------------------------------
pr_metrics <- function(pred, pos) {
  TP <- sum(pred == 1 & pos == 1L)
  FP <- sum(pred == 1 & pos == 0L)
  FN <- sum(pred == 0 & pos == 1L)

  P <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
  R <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)

  list(precision = P, recall = R)
}


# -----------------------------------------------------------------------------
# Bootstrap confidence intervals for decision-level precision and recall.
# -----------------------------------------------------------------------------
B <- 1000L

pr_ci_dt <- dt_qs_pred[, {

  pr0 = pr_metrics(pred, Positive)
  n <- .N

  boot_mat <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    pr_b = pr_metrics(pred[idx], Positive[idx])
    c(
      recall = pr_b$recall,
      precision = pr_b$precision
    )
  })

  recall_ci <- quantile(boot_mat["recall", ], c(0.025, 0.975), na.rm = TRUE)
  precision_ci <- quantile(boot_mat["precision", ], c(0.025, 0.975), na.rm = TRUE)

  .(
    Recall = pr0$recall,
    Recall_low = recall_ci[1],
    Recall_high = recall_ci[2],
    Precision = pr0$precision,
    Precision_low = precision_ci[1],
    Precision_high = precision_ci[2]
  )

}, by = .(alpha, Platform, Exprs_noise, Perturbation,
          Test_method, N_tsamples, Tsteps, N_replicates)]


library(grid)
library(scales)
library(ggrepel)

# Order q-value thresholds for plotting aesthetics.
pr_ci_dt$alpha_level <- factor(
  pr_ci_dt$alpha,
  levels = qs,
  ordered = TRUE
)

# Transparency scale across thresholds.
alpha_vals <- seq(1, 0.4, -0.6 / (length(qs) - 1))
names(alpha_vals) <- as.character(qs)


# -----------------------------------------------------------------------------
# Decision plot:
#   - x = recall
#   - y = precision
#   - fill = perturbation type
#   - alpha = q-value threshold
#
# Error bars represent bootstrap confidence intervals.
# -----------------------------------------------------------------------------
p_decision <- ggplot(
  pr_ci_dt,
  aes(x = Recall, y = Precision, fill = Perturbation, alpha = alpha_level)
) +
  geom_errorbar(
    aes(xmin = Recall_low, xmax = Recall_high),
    width = 0, linewidth = 0.30, alpha = 0.25, colour = "grey25"
  ) +
  geom_errorbar(
    aes(ymin = Precision_low, ymax = Precision_high),
    width = 0, linewidth = 0.30, alpha = 0.25, colour = "grey25"
  ) +
  geom_point(shape = 21, size = 2.8, stroke = 0.4, colour = "grey25") +
  scale_alpha_manual(
    name = expression(alpha),
    values = alpha_vals
  ) +
  guides(
    fill = guide_legend(
      title = NULL,
      order = 1
    ),
    alpha = guide_legend(
      order = 2,
      override.aes = list(
        fill = gray(seq(0.2, 0.85, length.out = length(alpha_vals))),
        colour = "grey25",
        shape = 21,
        size = 3
      )
    )
  ) +
  facet_grid(Platform ~ Exprs_noise) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

p_decision

ggsave("decision-plot.pdf", p_decision, width = 16, height = 10)