# =============================================================================
# Title: Nested-Test Analysis on Real Datasets
# Description:
#   This script applies the nested weighted NNLS test to multiple real
#   datasets, after merging event-level rMATS-derived measurements and
#   filtering low-information events.
#
#   Main workflow:
#     - load and merge multiple processed datasets,
#     - compute event-level QC summaries,
#     - define dataset-specific filtering thresholds,
#     - run the nested test on retained events,
#     - adjust p-values within dataset,
#     - derive prioritization scores from effect size, fit improvement,
#       and statistical evidence,
#     - export ranked result tables.
#
# Intended use:
#   Real-data application of the post-export kinetics framework.
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
# Set working directory and load the nested-test implementation.
# -----------------------------------------------------------------------------
setwd("~/postexport-kinetics/real_datasets")
source("../commons/nested_test.r")

require(data.table)


# -----------------------------------------------------------------------------
# Load each processed real dataset and append a human-readable dataset label.
# -----------------------------------------------------------------------------
load("GSE83620/rmats_dt_Kc167.rdata")
dt_kc167 <- copy(rmats_dt)
dt_kc167[, dataset := "Kc167 | GSE83620"]

load("GSE207924/rmats_dt_K562.rdata")
dt_k562 <- copy(rmats_dt)
dt_k562[, dataset := "K562 | GSE207924"]

load("GSE207924/rmats_dt_3T3.rdata")
dt_3t3 <- copy(rmats_dt)
dt_3t3[, dataset := "3T3 | GSE207924"]

load("GSE256335/rmats_dt_ECS.rdata")
dt_esc <- copy(rmats_dt)
dt_esc[, dataset := "mESC | GSE256335"]


# -----------------------------------------------------------------------------
# Merge all datasets into a single long table.
#
# fill = TRUE is used to tolerate small differences in column structure.
# -----------------------------------------------------------------------------
rmats_all <- rbindlist(
  list(dt_kc167, dt_k562, dt_3t3, dt_esc),
  use.names = TRUE,
  fill = TRUE
)

# Save combined object for reuse.
save(rmats_all, file = "rmats_all.rdata")


# -----------------------------------------------------------------------------
# Compute event-level QC metrics before statistical testing.
#
# Metrics include:
#   - total coverage across all compartments,
#   - cytoplasmic coverage,
#   - dynamic range of C and C_s,
#   - number and fraction of non-zero observations for C and C_s.
# -----------------------------------------------------------------------------
evt_qc <- rmats_all[, .(
  cov_total = sum(N + N_s + C + C_s, na.rm = TRUE),
  cov_cyt   = sum(C + C_s, na.rm = TRUE),
  rng_C     = max(C,   na.rm = TRUE) - min(C,   na.rm = TRUE),
  rng_Cs    = max(C_s, na.rm = TRUE) - min(C_s, na.rm = TRUE),
  nz_C      = sum(C   > 0, na.rm = TRUE),
  nz_Cs     = sum(C_s > 0, na.rm = TRUE),
  frac_nzC  = mean(C   > 0, na.rm = TRUE),
  frac_nzCs = mean(C_s > 0, na.rm = TRUE)
), by = .(dataset, event)]


# -----------------------------------------------------------------------------
# Define dataset-specific QC thresholds using low empirical quantiles.
#
# These thresholds are used as minimum information requirements before
# attempting the nested test.
# -----------------------------------------------------------------------------
thr <- evt_qc[, .(
  cov_total_min = quantile(cov_total, 0.01, na.rm = TRUE),
  cov_cyt_min   = quantile(cov_cyt,   0.01, na.rm = TRUE),
  rng_C_min     = quantile(rng_C,     0.01, na.rm = TRUE),
  rng_Cs_min    = quantile(rng_Cs,    0.01, na.rm = TRUE)
), by = dataset]


# -----------------------------------------------------------------------------
# Map number of replicates to variance-shrinkage strength.
#
# Heuristic choice:
#   R = 2  -> lambda = 0.7
#   R = 3  -> lambda = 0.5
#   R >= 4 -> lambda = 0.3
# -----------------------------------------------------------------------------
lambda = c("2" = 0.7, "3" = 0.5, "4" = 0.3, "5" = 0.3)


# -----------------------------------------------------------------------------
# Run the nested test event-by-event within each dataset.
#
# For each event:
#   1. apply QC and zero-inflation filters,
#   2. if the event is informative enough, run test_sigma_nested(),
#   3. return p-value, effect-size estimate, fit diagnostics, and annotations.
#
# Grouping variables:
#   - dataset
#   - event
#   - ensembl
# -----------------------------------------------------------------------------
results <- rmats_all[, {

  d <- .SD

  # ---------------------------------------------------------------------------
  # Hard QC thresholds.
  # These criteria enforce minimum information in both cytoplasmic and nuclear
  # compartments and exclude highly zero-inflated events.
  # ---------------------------------------------------------------------------
  min_pos      <- 3
  min_frac     <- 0.3
  min_pos_nuc  <- 3
  max_zero_frac_C  <- 0.7
  max_zero_frac_Cs <- 0.8

  # Dataset-specific quantile-based thresholds.
  th <- thr[dataset == dataset[1]]

  # Number of replicates available for this event.
  n_rep = length(unique(d$replicate))

  # ---------------------------------------------------------------------------
  # Determine whether the event should be discarded before testing.
  #
  # Reasons for exclusion include:
  #   - one compartment entirely zero,
  #   - too few positive observations,
  #   - too little total or cytoplasmic coverage,
  #   - too little variation in C or C_s,
  #   - excessive zero inflation.
  # ---------------------------------------------------------------------------
  bad <- (all(d$N == 0) || all(d$N_s == 0) || all(d$C == 0) || all(d$C_s == 0)) ||

    # Cytoplasmic informativeness
    (sum(d$C   > 0, na.rm = TRUE) < min_pos) ||
    (sum(d$C_s > 0, na.rm = TRUE) < min_pos) ||
    (mean(d$C  > 0, na.rm = TRUE) < min_frac) ||

    # Nuclear informativeness
    (sum(d$N   > 0, na.rm = TRUE) < min_pos_nuc) ||
    (sum(d$N_s > 0, na.rm = TRUE) < min_pos_nuc) ||

    # Coverage-based QC
    (sum(d$N + d$N_s + d$C + d$C_s, na.rm = TRUE) < th$cov_total_min) ||
    (sum(d$C + d$C_s, na.rm = TRUE) < th$cov_cyt_min) ||

    # Dynamic-range QC
    ((max(d$C,   na.rm = TRUE) - min(d$C,   na.rm = TRUE)) < th$rng_C_min) ||
    ((max(d$C_s, na.rm = TRUE) - min(d$C_s, na.rm = TRUE)) < th$rng_Cs_min) ||

    # Explicit zero-inflation filters
    (mean(d$C   == 0, na.rm = TRUE) > max_zero_frac_C) ||
    (mean(d$C_s == 0, na.rm = TRUE) > max_zero_frac_Cs)

  if (bad) {
    # Return no row for filtered-out events.
  } else {

    # -------------------------------------------------------------------------
    # Run the nested weighted NNLS test on the retained event.
    #
    # Notes:
    #   - scaling_A = TRUE enables column scaling,
    #   - parametric_whitened bootstrap is used here,
    #   - lambda_var depends on the number of replicates,
    #   - t_star = 0 effectively defines the split consistently for this setup.
    # -------------------------------------------------------------------------
    WWW <- test_sigma_nested(
      tsampled_data = d,
      scaling_A     = TRUE,
      boot_mode     = "parametric_whitened",
      B_n           = 5000,
      lambda_var    = lambda[n_rep],
      t_star        = 0
    )

    .(
      p.value = WWW$p.value,
      Sigma   = WWW$Sigma,
      Alpha   = WWW$Alpha,
      T.obs   = WWW$T.obs,
      RSS0    = WWW$RSS0,
      RSS1    = WWW$RSS1,
      IR      = WWW$IR,
      gene_symbol  = gene_symbol[1],
      description  = description[1],
      gene_biotype = gene_biotype[1]
    )
  }
}, by = .(dataset, event, ensembl)]


# -----------------------------------------------------------------------------
# Quick summary:
#   - number of tested events,
#   - number of unique genes,
#   - smallest p-value per dataset.
# -----------------------------------------------------------------------------
results[, .(
  m = .N,
  ng = length(unique(ensembl)),
  pmin = min(p.value, na.rm = TRUE)
), by = dataset]


# -----------------------------------------------------------------------------
# Estimate the fraction of events discarded by QC filtering.
# -----------------------------------------------------------------------------
a = results[, .(N = .N), by = dataset]
b = rmats_all[, .(N = length(unique(event))), by = dataset]
1 - a$N / b$N


# -----------------------------------------------------------------------------
# Within each dataset:
#   - adjust p-values by FDR,
#   - derive significance/evidence measures,
#   - compute fit-improvement diagnostics,
#   - define prioritization scores.
# -----------------------------------------------------------------------------
eps = 1e-10

results[, q.value := p.adjust(p.value, method = "fdr"), by = .(dataset)]
results[, neglogq := -log10(pmax(q.value, eps))]
results[, neglogp := -log10(pmax(p.value, eps))]

# Cap extreme evidence values.
cap <- 6
results[, neglogq_cap := pmin(neglogq, cap)]
results[, neglogp_cap := pmin(neglogp, cap)]


# -----------------------------------------------------------------------------
# Model-improvement metrics derived from RSS values.
# -----------------------------------------------------------------------------
results[, DeltaRSS := pmax(RSS0 - RSS1, 0)]
results[, IR := fifelse(RSS0 > 0, DeltaRSS / RSS0, 0)]
results[, logRSSratio := fifelse(RSS1 > 0, log(RSS0 / RSS1), 0)]


# -----------------------------------------------------------------------------
# Robust effect-size scaling.
#
# s0 is the median absolute Sigma within each dataset.
# This prevents very large Sigma values from dominating the score.
# -----------------------------------------------------------------------------
results[, s0 := median(abs(Sigma), na.rm = TRUE), by = .(dataset)]
results[s0 == 0, s0 := 1]

results[, sig_log := log1p(abs(Sigma + eps) / s0)]


# -----------------------------------------------------------------------------
# Ranking scores.
#
# score_siglogq:
#   effect size × evidence
#
# score_sig_IR_logq:
#   effect size × model improvement × evidence
#
# score_siglog_IR_logq:
#   robust effect size × model improvement × evidence
# -----------------------------------------------------------------------------
results[, score_siglogq := Sigma * neglogp_cap]
results[, score_sig_IR_logq := Sigma * IR * neglogp_cap]
results[, score_siglog_IR_logq := sig_log * IR * neglogp_cap]


# -----------------------------------------------------------------------------
# Quick summary of the final result table.
# -----------------------------------------------------------------------------
summary(results)

# Save analysis output.
save(results, file = "results_realdatasets.rdata")

# Reload if needed.
load("results_realdatasets.rdata")

# Inspect top-ranked events by the recommended score.
results[order(-score_sig_IR_logq)][1:10]


library(openxlsx)

# -----------------------------------------------------------------------------
# Create a complete supplementary-style table of tested events.
#
# Only events with non-missing p.value are included.
# -----------------------------------------------------------------------------
supp_all <- results[!is.na(p.value),
  .(
    dataset, event, ensembl,
    gene_symbol, gene_biotype,
    T.obs, RSS0, RSS1, score_sig_IR_logq,
    p.value, q.value
  )
]

# Sort within dataset by decreasing priority score.
setorder(supp_all, dataset, -score_sig_IR_logq)


# -----------------------------------------------------------------------------
# Quick manual inspection of selected high-confidence events.
# -----------------------------------------------------------------------------
results[IR > 0.1 & p.value < 0.01 & Sigma > 0.01]
results[order(-score_sig_IR_logq)][1:10]
results[order(-score_siglogq)][1:10]


# -----------------------------------------------------------------------------
# Build a dataset-wise ranked export table.
#
# For each dataset:
#   - sort by score_sig_IR_logq,
#   - retain selected columns,
#   - keep events with IR > 0.1 and p.value < 0.05.
# -----------------------------------------------------------------------------
supp_all = results[, {
  .SD[order(-score_sig_IR_logq),
      .(gene_symbol, gene_biotype, event, description, IR, p.value, Sigma)][IR > 0.1 & p.value < 0.05]
}, by = dataset]


# -----------------------------------------------------------------------------
# Export the ranked event table to Excel.
# -----------------------------------------------------------------------------
require(openxlsx)
write.xlsx(supp_all, file = "nested_test_results.xlsx", overwrite = TRUE)