# =============================================================================
# Title: Baseline PSI Comparison Using Logit-Transformed Pre/Post Testing
# Description:
#   This function implements a simple replicate-level comparison of PSI-like
#   proportions before and after a reference time point.
#
#   The procedure:
#     - aggregates C and C_s within pre- and post-windows,
#     - computes PSI = C / (C + C_s),
#     - applies a logit transform for stability,
#     - tests whether the pre/post difference is centered at zero,
#     - reports a simple effect-size summary on the PSI scale.
#
# Intended use:
#   Baseline reference method for comparison with the dynamical model.
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
# Compare pre- vs post-window PSI values using a t-test on logit-transformed
# proportions.
#
# Definition:
#   PSI = C / (C + C_s)
#
# Workflow:
#   1. identify pre and post time windows,
#   2. aggregate C and C_s within each window for each replicate,
#   3. compute replicate-level PSI in both windows,
#   4. compare logit(PSI_post) - logit(PSI_pre) to zero using a t-test,
#   5. return replicate-level summaries and a simple effect-size proxy.
#
# Notes:
#   - If t_star is provided:
#       pre  = times <= t_star
#       post = times >  t_star
#   - If t_star is not provided:
#       pre  = first time point only
#       post = all remaining time points
#   - Small eps offsets prevent division by zero and infinite logits.
#   - Replicates with too little total signal are filtered out.
#
# Output:
#   A list containing:
#     - Tab_coeff: replicate-level aggregated data and PSI summaries
#     - Sigma / Alpha: median decrease score on the PSI scale
#     - p.value: p-value from the t-test
#     - T.obs: observed t statistic
#     - additional placeholder fields for interface compatibility
# -----------------------------------------------------------------------------
baseline_psi_logit_ttest_pre_vs_post <- function(
  tsampled_data, ss_data,
  scaling_A = FALSE,
  usa_W = TRUE,
  add_SS = FALSE,
  B_n = -1,
  t_star = NULL,
  eps = 1e-6,
  min_total = 1L
){
  # Sort data by replicate and time for reproducibility/readability.
  tmain <- tsampled_data[order(tsampled_data$replicate, tsampled_data$time), ]

  # Check that required columns are present.
  req <- c("time","replicate","C","C_s")
  miss <- setdiff(req, names(tmain))
  if (length(miss) > 0) stop("tsampled_data is missing columns: ", paste(miss, collapse=", "))

  # Extract distinct observed time points.
  times_all <- sort(unique(tmain$time))
  if (length(times_all) < 2) stop("Need at least 2 distinct time points.")

  # Define pre/post windows.
  if (!is.null(t_star)) {
    # Pre window includes times up to and including t_star.
    pre_times  <- times_all[times_all <= t_star]

    # Post window includes times strictly after t_star.
    post_times <- times_all[times_all >  t_star]

    # Both windows must contain at least one time point.
    if (length(pre_times) < 1 || length(post_times) < 1) {
      stop("t_star yields empty pre or post window.")
    }
  } else {
    # Default split:
    #   pre  = first observed time point
    #   post = all subsequent time points
    pre_times  <- times_all[1]
    post_times <- times_all[-1]
  }

  # Aggregate C and C_s over the pre window for each replicate.
  dt_pre <- tmain[
    time %in% pre_times,
    .(C_pre = sum(C, na.rm=TRUE), Cs_pre = sum(C_s, na.rm=TRUE)),
    by = replicate
  ]

  # Aggregate C and C_s over the post window for each replicate.
  dt_post <- tmain[
    time %in% post_times,
    .(C_post = sum(C, na.rm=TRUE), Cs_post = sum(C_s, na.rm=TRUE)),
    by = replicate
  ]

  # Keep only replicates represented in both windows.
  dt <- merge(dt_pre, dt_post, by="replicate", all=FALSE)

  # Filter out replicates with insufficient total signal in either window.
  dt <- dt[((C_pre + Cs_pre) > min_total) & ((C_post + Cs_post) > min_total)]

  # If too few replicates remain, return a minimal NA result.
  if (nrow(dt) < 2) {
    return(list(
      Tab_coeff = as.data.frame(dt),
      Sigma = NA_real_,
      p.value = NA_real_,
      T.obs = NA_real_
    ))
  }

  # Compute PSI in the pre and post windows.
  psi_pre  <- dt$C_pre  / (dt$C_pre  + dt$Cs_pre  + eps)
  psi_post <- dt$C_post / (dt$C_post + dt$Cs_post + eps)

  # Stabilized logit transform.
  logit <- function(x) log((x + eps) / (1 - x + eps))

  # Replicate-level change on the logit scale.
  d <- logit(psi_post) - logit(psi_pre)

  # Test whether the average logit-scale change differs from zero.
  # Interpretation:
  #   d < 0 corresponds to PSI decreasing from pre to post.
  tt <- t.test(d, mu = 0, alternative = "two.sided")

  # Add replicate-level summaries on the original PSI scale.
  dt[, `:=`(
    PSI_pre = psi_pre,
    PSI_post = psi_post,
    dPSI = psi_post - psi_pre,

    # Positive score corresponds to a decrease in PSI after the split.
    score = -(psi_post - psi_pre)
  )]

  # Effect-size summary: median replicate-level decrease score.
  Sigma <- median(dt$score, na.rm=TRUE)

  # Return results using a structure aligned with the main pipeline.
  list(
    Tab_coeff = as.data.frame(dt),
    Sigma = Sigma,
    Alpha = Sigma,
    p.value = tt$p.value,

    # Placeholder values retained for compatibility with downstream code.
    RSS0 = 1,
    RSS1 = 1,
    IR = 1,

    # Observed t statistic.
    T.obs = as.numeric(tt$statistic),

    # Placeholder diagnostic retained for compatibility.
    deltaAIC = as.numeric(tt$statistic)
  )
}