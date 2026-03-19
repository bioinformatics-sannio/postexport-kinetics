# =============================================================================
# Title: Representative Dynamics Plots for Real-Dataset Events
# Description:
#   This script loads the real-data nested-test results, derives a few
#   additional ranking scores, and generates representative model-fit plots
#   for selected example events.
#
#   Main workflow:
#     - load combined event-level data and nested-test results,
#     - compute significance-based and effect-based ranking scores,
#     - summarize dataset-level discovery statistics,
#     - select example positive events,
#     - refit the nested model for those events,
#     - visualize observed dynamics together with full/null model fits,
#     - export a combined multi-panel figure.
#
# Intended use:
#   Figure generation and qualitative illustration of representative
#   post-export dynamics in real datasets.
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
# Set working directory and load the nested-test and plotting utilities.
# -----------------------------------------------------------------------------
setwd("~/postexport-kinetics/real_datasets")

source("../commons/nested_test.r")
source("../commons/plot.r")

require(data.table)
require(ggplot2)


# -----------------------------------------------------------------------------
# Load:
#   - rmats_all: combined real-dataset event-level time-course table
#   - results_realdatasets: nested-test output table
# -----------------------------------------------------------------------------
load("rmats_all.rdata")
load("results_realdatasets.rdata", verbose = TRUE)


# -----------------------------------------------------------------------------
# Derive additional ranking scores from p-values, Sigma, and deltaAIC.
#
# These scores are exploratory summaries used for ranking and manual inspection.
# -----------------------------------------------------------------------------
eps <- 1e-10
results[, neglogp := -log10(pmax(p.value, eps))]

results[, score_siglogp := Sigma * pmin(neglogp, 6)]
results[, score_aiclogp := pmax(deltaAIC, 0) * pmin(neglogp, 6)]
results[, score_aic := deltaAIC]
results[, score_neglogp := -log10(pmax(p.value, 1e-12))]
results[, score_sigma := Sigma]


# =============================================================================
# Representative dynamics: preliminary dataset summaries
# =============================================================================

# -----------------------------------------------------------------------------
# Compute BH-adjusted q-values within each dataset.
# -----------------------------------------------------------------------------
results[, qvalue := p.adjust(p.value, method = "BH"), by = dataset]

# -----------------------------------------------------------------------------
# Dataset-level summary:
#   - number of tested events,
#   - smallest p-value,
#   - smallest q-value,
#   - number of very small p-values (< 1e-4).
# -----------------------------------------------------------------------------
results[, .(
  m = sum(!is.na(p.value)),
  min_p = min(p.value, na.rm = TRUE),
  min_q = min(qvalue),
  p_1e4 = sum(p.value < 1e-4, na.rm = TRUE)
), by = dataset]

library(data.table)
setDT(results)


# -----------------------------------------------------------------------------
# Quick manual inspection of candidate events in specific datasets.
#
# These lines are used interactively to identify interesting examples.
# -----------------------------------------------------------------------------
results[grepl("ESC", dataset) & p.value < 0.01 & Sigma > 0.01][order(-score_siglogq)]
results[grepl("Kc167", dataset) & p.value < 0.05][order(-IR)][1:5]

results[IR > 0.1 & p.value < 0.01]


# =============================================================================
# Example 1: representative positive event
# =============================================================================

# -----------------------------------------------------------------------------
# Reload plotting helpers explicitly if needed.
# -----------------------------------------------------------------------------
source("../commons/plot.r")

# Manually selected example event.
mygene = "FBgn0266173:chr3R:-:21356811:21356968:21357115:21357758"

# Extract the event-specific time-course data.
ts_gene <- rmats_all[event == mygene]

# Refit the nested model for this event.
res_nested = test_sigma_nested(ts_gene, t_star = 0, B_n = 5000, lambda_var = 0.7)

# Build readable title/subtitle from the event string and test output.
evento = unlist(strsplit(ts_gene$event[1], ":"))
evento = paste0(evento[2], ":", evento[3], ":", evento[5], "-", evento[6])

titolo = paste0("(", ts_gene$gene_symbol[1], ") ", evento)
sottotitolo = paste0(
  "(nested test p=",
  signif(res_nested$p.value, 2),
  ", Relative RSS=",
  signif(res_nested$IR, 2),
  ")"
)

# Plot observed data with full and null fitted trajectories.
p1 = plot_fit_two_models_paper(
  ts_gene,
  res_nested,
  t_star = 0,
  show_shutoff = FALSE,
  step = 1
) +
  ggtitle(titolo, subtitle = sottotitolo) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )


# =============================================================================
# Example 2: representative positive event
# =============================================================================

# -----------------------------------------------------------------------------
# Additional candidate examples are left here as commented alternatives.
# Only the final assignment is retained.
# -----------------------------------------------------------------------------
# Scaf4
# mygene="ENSMUSG00000022983.17:chr16:-:90039415:90039606:90039760:90039846"

# Tia1
mygene = "ENSMUSG00000071337.12:chr6:+:86395860:86395914:86397306:86397393"

# Wdr73
# mygene="ENSMUSG00000025722:chr7:-:80900015:80900082:80900556:80900603"

# Hp1bp3
# mygene="ENSMUSG00000028759.14:chr4:+:137953200:137953359:137955955:137956098"

# ZMYM5
# mygene="ENSG00000132950:chr13:-:19851371:19851464:19851705:19852206"

# mygene="ENSMUSG00000031939:chr9:+:15306846:15306963:15307701:15307794"

# Extract event-specific time-course data.
ts_gene <- rmats_all[event == mygene]

# Refit the nested model for this event.
res_nested = test_sigma_nested(ts_gene, t_star = 0, B_n = 5000, lambda_var = 0.7)

# Build readable title/subtitle.
evento = unlist(strsplit(ts_gene$event[1], ":"))
evento = paste0(evento[2], ":", evento[3], ":", evento[5], "-", evento[6])

titolo = paste0("(", ts_gene$gene_symbol[1], ") ", evento)
sottotitolo = paste0(
  "(nested test p=",
  signif(res_nested$p.value, 2),
  ", Relative RSS=",
  signif(res_nested$IR, 2),
  ")"
)

# Plot observed data with full and null fitted trajectories.
p2 = plot_fit_two_models_paper(
  ts_gene,
  res_nested,
  t_star = 0,
  show_shutoff = FALSE,
  step = 1
) +
  ggtitle(titolo, subtitle = sottotitolo) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )


# =============================================================================
# Combine and export representative plots
# =============================================================================

library(patchwork)

# Two-panel layout for the selected examples.
p_combined <- p1 + p2 +
  plot_layout(ncol = 2, widths = c(1, 1))

p_combined

# Export the combined figure.
ggsave("representative_dynamics.pdf", p_combined, width = 10, height = 5)