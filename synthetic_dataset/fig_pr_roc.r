# =============================================================================
# Title: Benchmark Visualization and Ranking Analysis on Synthetic Data
# Description:
#   This script loads the benchmark results produced on the synthetic dataset
#   and generates summary visualizations for prioritization performance.
#
#   Main outputs:
#     - PR-AUC and ROC-AUC summaries with bootstrap confidence intervals,
#     - precision-recall curves with confidence ribbons,
#     - ROC-style curves,
#     - top-K precision/recall operating points,
#     - publication-style faceted figures across platforms and noise levels.
#
#   The analysis is organized around different score definitions derived from
#   the benchmark output and compares methods across perturbation regimes.
#
# Intended use:
#   Downstream visualization and figure generation for benchmark comparisons.
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

# Load benchmark result table.
load("test_sigma_2k_new.rdata")


# -----------------------------------------------------------------------------
# Reshape selected ranking scores into long format.
#
# Each row will correspond to one score definition:
#   - score_siglog_IR_logq
#   - score_sig_IR_logq
#   - score_siglogq
#
# This allows easy filtering and comparison of prioritization strategies.
# -----------------------------------------------------------------------------
dt_test_results <- melt(
  dt_test_results,
  measure.vars = c("score_siglog_IR_logq", "score_sig_IR_logq", "score_siglogq"),
  variable.name = "Score_method",
  value.name    = "Score"
)


# =============================================================================
# FIGURE BLOCK 1
# Prioritization analysis:
#   - PR-AUC
#   - ROC-AUC
#   - Top-K summaries
#
# Goal:
#   Show the benefit of SHUTOFF and the dependence on platform and noise level.
#   Supplementary materials may include full ROC-AUC tables and additional
#   summaries.
# =============================================================================


# -----------------------------------------------------------------------------
# Helper: precision-recall AUC.
#
# Inputs:
#   - score: ranking score
#   - y01: binary truth labels (1 = positive, 0 = negative)
#
# Output:
#   PR-AUC computed using PRROC.
# -----------------------------------------------------------------------------
pr_auc <- function(score, y01) {
  s_pos <- score[y01 == 1L]
  s_neg <- score[y01 == 0L]
  if (length(s_pos) == 0L || length(s_neg) == 0L) return(NA_real_)

  PRROC::pr.curve(
    scores.class0 = s_pos,
    scores.class1 = s_neg,
    curve = FALSE
  )$auc.integral
}


# -----------------------------------------------------------------------------
# Helper: ROC-AUC.
#
# Uses pROC and returns NA if only one class is present.
# -----------------------------------------------------------------------------
roc_auc <- function(score, y01) {
  if (length(unique(y01)) < 2L) return(NA_real_)

  as.numeric(
    pROC::auc(
      pROC::roc(
        response = y01,
        predictor = score,
        quiet = TRUE
      )
    )
  )
}


# -----------------------------------------------------------------------------
# Helper: top-K precision and recall.
#
# The top K rows according to decreasing score are selected and basic
# classification summaries are computed.
# -----------------------------------------------------------------------------
top_pr <- function(score, y01, K = 100L) {
  topk_idx = order(score, decreasing = TRUE)[1:K]

  TP = sum(y01[topk_idx] == 1L)
  FP = sum(y01[topk_idx] == 0L)
  TN = sum(y01[-topk_idx] == 0L)
  FN = sum(y01[-topk_idx] == 1L)

  P = TP / (TP + FP)
  R = TP / (TP + FN)

  list(precision = P, recall = R)
}


# -----------------------------------------------------------------------------
# Filter one benchmark scenario for inspection.
#
# Current filter:
#   - N_replicates = 5
#   - score method = score_siglogq
#   - Tsteps = 10
#   - Lambda = 0.5
#   - Test_method = Nested-wild_whitened
#
# This slice is mainly used to inspect ROC behavior across platform, noise,
# and perturbation conditions.
# -----------------------------------------------------------------------------
w_dt <- dt_test_results[
  N_replicates == 5 &
  Score_method == "score_siglogq" &
  Tsteps == 10 &
  Lambda == 0.5 &
  Test_method %in% c("Nested-wild_whitened")
]


# -----------------------------------------------------------------------------
# Rename test methods for cleaner plotting labels.
# -----------------------------------------------------------------------------
w_dt[, Test_method := factor(
  Test_method,
  levels = c("Nested-wild_whitened", "PSIBaseline"),
  labels = c("Nested test", "PSI Baseline")
)]


# Check missing scores if needed.
w_dt[is.na(Score)]


# -----------------------------------------------------------------------------
# Bootstrap confidence intervals for:
#   - PR-AUC
#   - ROC-AUC
#   - Recall@K
#   - Precision@K
#
# Here K is set to 100 and bootstrap resampling is performed within each
# benchmark stratum.
# -----------------------------------------------------------------------------
B <- 1000L
K = 100L

auc_ci_dt <- w_dt[, {

  # Point estimates on the original sample.
  pr0  <- pr_auc(Score, Positive)
  roc0 <- roc_auc(Score, Positive)
  top_pr0 = top_pr(Score, Positive, K = K)

  n <- .N

  # Bootstrap resampling with replacement.
  boot_mat <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    top_pr_b = top_pr(Score[idx], Positive[idx], K = K)

    c(
      pr        = pr_auc(Score[idx], Positive[idx]),
      roc       = roc_auc(Score[idx], Positive[idx]),
      recall    = top_pr_b$recall,
      precision = top_pr_b$precision
    )
  })

  # Percentile confidence intervals.
  pr_ci  <- quantile(boot_mat["pr", ],        c(0.025, 0.975), na.rm = TRUE)
  roc_ci <- quantile(boot_mat["roc", ],       c(0.025, 0.975), na.rm = TRUE)
  recall_ci <- quantile(boot_mat["recall", ], c(0.025, 0.975), na.rm = TRUE)
  precision_ci <- quantile(boot_mat["precision", ], c(0.025, 0.975), na.rm = TRUE)

  .(
    PR_AUC = pr0,
    PR_AUC_low = pr_ci[1],
    PR_AUC_high = pr_ci[2],
    ROC_AUC = roc0,
    ROC_AUC_low = roc_ci[1],
    ROC_AUC_high = roc_ci[2],
    Recall = top_pr0$recall,
    Recall_low = recall_ci[1],
    Recall_high = recall_ci[2],
    Precision = top_pr0$precision,
    Precision_low = precision_ci[1],
    Precision_high = precision_ci[2]
  )

}, by = .(Score_method, N_replicates, Tsteps, Test_method, N_tsamples, Platform, Exprs_noise, Perturbation)]


# Save summary table with bootstrap intervals.
save(auc_ci_dt, file = "auc_ci_dt.rdata")


# -----------------------------------------------------------------------------
# Plot PR-AUC as a function of the number of sampled time points.
#
# One panel per platform/noise combination.
# Color indicates perturbation regime.
# -----------------------------------------------------------------------------
library(ggplot2)

auc_ci_dt[, N_tsamples_f := factor(N_tsamples, levels = c(3, 5, 10, 20))]

p_pr = ggplot(
  auc_ci_dt,
  aes(x = N_tsamples, y = PR_AUC, colour = Perturbation)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  geom_errorbar(
    aes(ymin = PR_AUC_low, ymax = PR_AUC_high),
    width = 0.5,
    linewidth = 0.45
  ) +
  facet_grid(Platform ~ Exprs_noise) +
  scale_x_continuous(breaks = c(3, 5, 10, 20)) +
  labs(
    x = "# time points",
    y = "PR-AUC"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_pr

ggsave("pr-auc_plot.pdf", p_pr, width = 16, height = 10, useDingbats = FALSE)


# -----------------------------------------------------------------------------
# Plot ROC-AUC as a function of the number of sampled time points.
# -----------------------------------------------------------------------------
p_roc = ggplot(
  auc_ci_dt,
  aes(x = N_tsamples, y = ROC_AUC, colour = Perturbation)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  geom_errorbar(
    aes(ymin = ROC_AUC_low, ymax = ROC_AUC_high),
    width = 0.5,
    linewidth = 0.45
  ) +
  facet_grid(Platform ~ Exprs_noise) +
  scale_x_continuous(breaks = c(3, 5, 10, 20)) +
  labs(
    x = "# time points",
    y = "ROC-AUC"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_roc

ggsave("roc-auc_plot.pdf", p_roc, width = 16, height = 10, useDingbats = FALSE)


# =============================================================================
# FIGURE BLOCK 2
# PR curves with confidence ribbons and top-K markers.
#
# Goal:
#   Compare the Nested test against the PSI baseline under a fixed SHUTOFF
#   scenario and controlled replicate/time-point settings.
# =============================================================================


# -----------------------------------------------------------------------------
# Helper: return full PR and ROC traversal coordinates.
#
# For a ranking score:
#   - recall and precision are computed along the sorted list,
#   - fpr is also computed for ROC-style plotting.
# -----------------------------------------------------------------------------
pr_xy <- function(score, y01) {
  n_y01 = 1 * (!y01)

  R = cumsum(y01[order(score, decreasing = TRUE)]) / sum(y01)
  P = cumsum(y01[order(score, decreasing = TRUE)]) / 1:length(score)
  FPR = cumsum(n_y01[order(score, decreasing = TRUE)]) / sum(n_y01)

  list(recall = R, precision = P, fpr = FPR)
}


# -----------------------------------------------------------------------------
# Helper: safer top-K precision/recall calculation.
#
# This version protects against empty inputs and K > n.
# -----------------------------------------------------------------------------
top_pr <- function(score, y01, K = 100L) {
  n <- length(score)
  if (n == 0L) return(list(precision = NA_real_, recall = NA_real_))

  K <- min(K, n)
  ord <- order(score, decreasing = TRUE)
  topk_idx <- ord[seq_len(K)]

  TP <- sum(y01[topk_idx] == 1L)
  FP <- sum(y01[topk_idx] == 0L)
  FN <- sum(y01[-topk_idx] == 1L)

  P <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  R <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_

  list(precision = P, recall = R)
}


# -----------------------------------------------------------------------------
# Filter the main comparison scenario.
#
# Current choice:
#   - N_tsamples = 10
#   - N_replicates = 5
#   - score_siglogq
#   - Tsteps = 10
#   - Perturbation = SHUTOFF
#   - Nested test with Lambda = 0.5
#   - PSI baseline with Lambda = 0
# -----------------------------------------------------------------------------
w_dt <- dt_test_results[
  N_tsamples == 10 &
  N_replicates == 5 &
  Score_method == "score_siglogq" &
  Tsteps == 10 &
  Perturbation == "SHUTOFF" &
  ((Lambda == 0.5 & Test_method == "Nested-wild_whitened") |
   (Lambda == 0   & Test_method == "PSIBaseline"))
]


# -----------------------------------------------------------------------------
# Flip PSI baseline score sign so that larger values consistently indicate
# stronger evidence for positives.
# -----------------------------------------------------------------------------
w_dt[Test_method == "PSIBaseline", Score := -Score]


# -----------------------------------------------------------------------------
# Set plotting labels for the compared methods.
# -----------------------------------------------------------------------------
w_dt[, Test_method := factor(
  Test_method,
  levels = c("Nested-parametric_whitened", "Nested-wild_whitened", "PSIBaseline"),
  labels = c("Nested test param", "Nested test wild", "PSI Baseline")
)]


# -----------------------------------------------------------------------------
# Bootstrap PR curves with pointwise confidence intervals.
#
# For each benchmark stratum:
#   - compute the empirical PR curve,
#   - bootstrap the precision profile B times,
#   - derive pointwise 2.5% and 97.5% quantiles.
# -----------------------------------------------------------------------------
B = 1000L

pr_curve_dt <- w_dt[, {
  pr0 <- pr_xy(Score, Positive)

  n = length(Score)
  boot_pr = matrix(NA_real_, nrow = length(pr0$recall), ncol = B)

  lapply(seq_len(B), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    pr_b <- pr_xy(Score[idx], Positive[idx])
    boot_pr[, b] <<- pr_b$precision
  })

  pr_hi = apply(boot_pr, 1, quantile, probs = 0.975, na.rm = TRUE)
  pr_lo = apply(boot_pr, 1, quantile, probs = 0.025, na.rm = TRUE)

  data.table(
    recall = pr0$recall,
    fpr = pr0$fpr,
    precision = pr0$precision,
    precision_low = pr_lo,
    precision_high = pr_hi
  )
}, by = .(N_replicates, Tsteps, Test_method, N_tsamples, Platform, Exprs_noise, Perturbation)]


# -----------------------------------------------------------------------------
# Compute top-K operating points for K = 50, 100, 200.
#
# These are later overlaid as markers on the PR curves.
# -----------------------------------------------------------------------------
topK_dt <- rbindlist(lapply(c(50L, 100L, 200L), function(K_top) {
  w_dt[, {
    tp <- top_pr(Score, Positive, K = K_top)
    data.table(recall = tp$recall, precision = tp$precision, K = K_top)
  }, by = .(N_replicates, Tsteps, Test_method, N_tsamples, Platform, Exprs_noise, Perturbation)]
}))

topK_dt <- topK_dt[is.finite(recall) & is.finite(precision)]

topK_dt[, K := factor(
  K,
  levels = c(50, 100, 200),
  labels = c("Top-50", "Top-100", "Top-200")
)]


# -----------------------------------------------------------------------------
# Color palette for method labels.
# -----------------------------------------------------------------------------
meth_cols <- c(
  "Nested test wild"  = "#E76F61",
  "Nested test param" = "#2a19bd",
  "PSI Baseline"      = "#2CA02C"
)


# -----------------------------------------------------------------------------
# Main PR-curve figure with confidence ribbons and Top-K markers.
# -----------------------------------------------------------------------------
p_pr_curves_main <- ggplot(
  pr_curve_dt,
  aes(x = recall, y = precision, color = Test_method)
) +
  geom_ribbon(
    aes(ymin = precision_low, ymax = precision_high, fill = Test_method),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 0.5) +
  geom_point(
    data = topK_dt,
    aes(x = recall, y = precision, fill = Test_method, shape = K),
    colour = "grey55",
    size = 2.6,
    stroke = 0.8,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  scale_color_manual(values = meth_cols) +
  scale_fill_manual(values = meth_cols) +
  scale_shape_manual(values = c("Top-50" = 21, "Top-100" = 24, "Top-200" = 22)) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linetype = 1, linewidth = 1.2)),
    fill  = "none",
    shape = guide_legend(order = 2, override.aes = list(fill = "white", colour = "grey40"))
  ) +
  facet_grid(Platform ~ Exprs_noise) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Recall", y = "Precision", color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

p_pr_curves_main

ggsave("pr-curves_with_topK.pdf", p_pr_curves_main, width = 8.5, height = 6, useDingbats = FALSE)


# -----------------------------------------------------------------------------
# ROC-style curve figure.
#
# Here:
#   x = false positive rate
#   y = recall (true positive rate)
# -----------------------------------------------------------------------------
p_roc_curves_main <- ggplot(
  pr_curve_dt,
  aes(x = fpr, y = recall, color = Test_method)
) +
  geom_line(linewidth = 0.5) +
  facet_grid(Platform ~ Exprs_noise) +
  labs(x = "TPR", y = "FPR", color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

p_roc_curves_main

ggsave("roc-curves.pdf", p_roc_curves_main, width = 8.5, height = 6, useDingbats = FALSE)