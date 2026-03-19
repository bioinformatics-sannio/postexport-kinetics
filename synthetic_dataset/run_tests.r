# =============================================================================
# Title: Parallel Benchmark Evaluation of Statistical Tests on Synthetic Data
# Description:
#   This script evaluates multiple statistical testing strategies on the
#   synthetic ODE-based dataset across a wide range of experimental settings.
#
#   For each gene and each scenario, the script:
#     - selects the corresponding synthetic time-course data,
#     - simulates platform-specific technical noise,
#     - runs one or more statistical tests,
#     - stores p-values, effect-size estimates, and fit diagnostics,
#     - computes multiple-testing corrected q-values,
#     - derives ranking scores combining effect size, fit improvement,
#       and statistical significance.
#
#   The benchmark spans:
#     - multiple platforms (RT-qPCR, Gaussian noise, RNA-seq),
#     - baseline and shutoff perturbations,
#     - several technical noise regimes,
#     - multiple shrinkage settings,
#     - different numbers of replicates and sampled time points.
#
# Intended use:
#   Large-scale method comparison on simulated datasets for performance
#   assessment, calibration analysis, and prioritization benchmarking.
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
# Set working directory and load required analysis modules.
#
# nested_test.r:
#   weighted NNLS-based nested bootstrap test.
#
# psi_test.r:
#   simple PSI-based baseline comparator.
#
# platforms.r:
#   platform-specific technical noise simulation utilities.
# -----------------------------------------------------------------------------
setwd("~/postexport-kinetics/synthetic_dataset")
source("../commons/nested_test.r")
source("../commons/psi_test.r")
source("../commons/platforms.r")


require(data.table)
library(parallel)

# -----------------------------------------------------------------------------
# Load the synthetic benchmark dataset generated previously.
# -----------------------------------------------------------------------------
load("ode_states_2k_20p.rdata")


# -----------------------------------------------------------------------------
# Define the benchmark dimensions.
#
# range_platform:
#   measurement platform used to perturb the latent simulated data.
#
# range_perturbation:
#   scenario type:
#     - NONE    : baseline time course
#     - SHUTOFF : intervention time course
#
# range_technical_noise:
#   qualitative levels of measurement noise, later mapped to platform-specific
#   parameters.
#
# range_lambda:
#   shrinkage strengths used in the variance regularization step of the nested
#   test.
#
# scaleA:
#   whether to scale the columns of A before fitting.
#
# test_method:
#   list of testing procedures to evaluate.
# -----------------------------------------------------------------------------
range_platform <- c("RT-qPCR", "GAUSS", "RNA-seq")
range_perturbation <- c("NONE", "SHUTOFF")
range_technical_noise <- c("Very low", "Low", "Medium", "High")
range_lambda <-  c(0, 0.3, 0.5, 0.7)
scaleA <- T
test_method <- c("Nested-parametric_whitened", "Nested-wild_whitened", "PSIBaseline")


# -----------------------------------------------------------------------------
# Platform-specific Gaussian noise levels.
# Used only when Platform == "GAUSS".
# -----------------------------------------------------------------------------
range_gauss_noise <- c(0.02, 0.05, 0.1, 0.2)
names(range_gauss_noise) <- range_technical_noise


# -----------------------------------------------------------------------------
# Extract the design settings present in the synthetic dataset.
# -----------------------------------------------------------------------------
range_n_time_samples <- unique(ode_states$N_time_samples)
range_n_replicates <- unique(ode_states$N_replicates)

# Here T_step is fixed to 10 for the benchmark run.
# A more general option would be:
#   unique(ode_states$T_step)
range_tsteps <- 10

# All genes present in the dataset.
genes <- unique(ode_states$Gene)


# -----------------------------------------------------------------------------
# Run the full benchmark for one gene.
#
# For each gene, the function loops over all combinations of:
#   - number of sampled time points,
#   - perturbation type,
#   - shutoff step spacing,
#   - number of replicates,
#   - measurement platform,
#   - technical noise level,
#   - A-scaling setting,
#   - shrinkage parameter lambda,
#   - statistical test.
#
# At each combination:
#   1. the matching synthetic data are extracted,
#   2. platform-specific measurement noise is simulated,
#   3. the selected test is run,
#   4. the result is stored as one row of output.
#
# Arguments:
#   - gene_i: gene identifier
#   - N_boot: number of bootstrap replicates for nested tests
#
# Output:
#   A list of result rows for the selected gene, or NULL if no valid
#   combinations are found.
# -----------------------------------------------------------------------------
run_for_gene <- function(gene_i, N_boot=5000) {
  list_results_gene <- list()
  ki <- 1

  for (time_samples_i in range_n_time_samples) {
    for (perturbation_i in range_perturbation) {
      for (tstep_i in range_tsteps) {
        for (n_replicates_i in range_n_replicates) {
          for (platform_i in range_platform) {
            for (noise_i in range_technical_noise) {
              for (scaling_A_i in scaleA) {
                for (lambda_i in range_lambda) {
                  for (test_method_i in test_method) {

                    # ---------------------------------------------------------
                    # Parse test configuration.
                    #
                    # Examples:
                    #   "Nested-parametric_whitened"
                    #   "Nested-wild_whitened"
                    #   "PSIBaseline"
                    # ---------------------------------------------------------
                    s_test_i <- unlist(strsplit(test_method_i, "-"))
                    test_i <- s_test_i[1]
                    b_type_i <- if (length(s_test_i) > 1) s_test_i[2]

                    # ---------------------------------------------------------
                    # Select the subset of synthetic data for this gene and
                    # design combination.
                    # ---------------------------------------------------------
                    dt_results_sub <- ode_states[
                      Gene == gene_i &
                        N_replicates == n_replicates_i &
                        N_time_samples == time_samples_i &
                        T_step == tstep_i
                    ]

                    # Skip combinations not present in the dataset.
                    if (nrow(dt_results_sub) == 0L) next

                    # Separate time-course data from steady-state summaries.
                    dt_results_sub_timecourse <- dt_results_sub[Perturbation == perturbation_i]
                    dt_results_sub_steadystate <- dt_results_sub[Perturbation == "SteadyState"]

                    if (nrow(dt_results_sub_timecourse) == 0L ||
                        nrow(dt_results_sub_steadystate) == 0L) {
                      next
                    }

                    # State variables to which technical noise will be applied.
                    cols_to_noise <- c("N", "C", "C_s", "N_s")

                    # ---------------------------------------------------------
                    # Platform-specific technical noise simulation.
                    # ---------------------------------------------------------
                    if (platform_i == "GAUSS") {

                      # Add Gaussian noise proportional to signal magnitude.
                      noise_sd_i <- range_gauss_noise[noise_i]

                      dt_results_sub_timecourse <- add_gaussian_noise(
                        dt_results_sub_timecourse,
                        cols_to_noise,
                        noise_sd_i
                      )

                      dt_results_sub_steadystate <- add_gaussian_noise(
                        dt_results_sub_steadystate,
                        cols_to_noise,
                        noise_sd_i
                      )

                    } else if (platform_i == "RT-qPCR") {

                      # RT-qPCR-like measurement simulation:
                      #   - ct_sd controls technical Ct noise,
                      #   - scale_copies controls effective molecular scale.
                      cd_sd_i <- switch(noise_i,
                        "Very low" = 0.01,
                        "Low"      = 0.05,
                        "Medium"   = 0.1,
                        "High"     = 0.25
                      )

                      scale_copies_i <- switch(noise_i,
                        "Very low" = 20,
                        "Low"      = 15,
                        "Medium"   = 10,
                        "High"     = 5
                      )

                      dt_results_sub_timecourse <- simulate_rt_qpcr(
                        dt_results_sub_timecourse,
                        targets      = cols_to_noise,
                        ct_sd        = cd_sd_i,
                        scale_copies = scale_copies_i
                      )

                      dt_results_sub_steadystate <- simulate_rt_qpcr(
                        dt_results_sub_steadystate,
                        targets      = cols_to_noise,
                        ct_sd        = cd_sd_i,
                        scale_copies = scale_copies_i
                      )

                    } else if (platform_i == "RNA-seq") {

                      # RNA-seq-like simulation:
                      #   - scale_counts controls expected sequencing depth,
                      #   - mean_disp and cv_disp control overdispersion.
                      scale_counts_i <- switch(noise_i,
                        "Very low" = 10000,
                        "Low"      = 5000,
                        "Medium"   = 1000,
                        "High"     = 200
                      )

                      mean_disp_i <- switch(noise_i,
                        "Very low" = 0.01,
                        "Low"      = 0.05,
                        "Medium"   = 0.10,
                        "High"     = 0.25
                      ) * 0.25

                      cv_disp_i <- switch(noise_i,
                        "Very low" = 0.5,
                        "Low"      = 0.7,
                        "Medium"   = 0.8,
                        "High"     = 1.0
                      )

                      dt_results_sub_timecourse <- simulate_rnaseq(
                        dt_results_sub_timecourse,
                        targets = cols_to_noise,
                        scale_counts = scale_counts_i,
                        mean_disp = mean_disp_i,
                        cv_disp = cv_disp_i
                      )

                      dt_results_sub_steadystate <- simulate_rnaseq(
                        dt_results_sub_steadystate,
                        targets = cols_to_noise,
                        scale_counts = scale_counts_i,
                        mean_disp = mean_disp_i,
                        cv_disp = cv_disp_i
                      )
                    }

                    # ---------------------------------------------------------
                    # Run the selected statistical test.
                    # ---------------------------------------------------------
                    if (test_i == "Nested") {

                      # Nested weighted NNLS test with bootstrap.
                      # t_star is used only for shutoff trajectories.
                      WW <- test_sigma_nested(
                        tsampled_data = dt_results_sub_timecourse,
                        scaling_A = scaling_A_i,
                        boot_mode = b_type_i,
                        lambda_var = lambda_i,
                        B_n = N_boot,
                        t_star = if (perturbation_i == "NONE") NULL else unique(dt_results_sub_timecourse$T_star)
                      )

                    } else if (test_i == "PSIBaseline") {

                      # Baseline PSI method is only run for lambda = 0,
                      # since lambda is irrelevant for this approach.
                      if (lambda_i > 0) next

                      # Small pseudocount-like offset for numerical stability.
                      dt_results_sub_timecourse[, C := C + 1]
                      dt_results_sub_timecourse[, C_s := C_s + 1]

                      WW <- baseline_psi_logit_ttest_pre_vs_post(
                        tsampled_data = dt_results_sub_timecourse,
                        t_star = NULL
                        # Alternative:
                        # t_star = if (perturbation_i == "NONE") NULL else unique(dt_results_sub_timecourse$T_star)
                      )

                    } else {
                      next
                    }

                    # ---------------------------------------------------------
                    # Store benchmark result for the current combination.
                    # ---------------------------------------------------------
                    r_i <- data.frame(
                      Gene         = gene_i,
                      N_replicates = n_replicates_i,
                      N_tsamples   = time_samples_i,
                      Lambda       = lambda_i,
                      Tsteps       = tstep_i,
                      Perturbation = perturbation_i,
                      Exprs_noise  = noise_i,
                      Platform     = platform_i,
                      Scaling_A    = scaling_A_i,
                      Test_method  = test_method_i,
                      Positive     = unique(dt_results_sub$truth_pos),
                      p.value      = WW$p.value,
                      T.obs        = WW$T.obs,
                      Sigma        = WW$Sigma,
                      Alpha        = WW$Alpha,
                      RSS0         = WW$RSS0,
                      RSS1         = WW$RSS1,
                      IR           = WW$IR,
                      sigma_true_mean   = mean(dt_results_sub_steadystate$Base_sigma_c),
                      sigma_true_sd     = sd(dt_results_sub_steadystate$Base_sigma_c),
                      sigma_true_median = median(dt_results_sub_steadystate$Base_sigma_c),
                      sigma_true_min    = min(dt_results_sub_steadystate$Base_sigma_c),
                      sigma_true_max    = max(dt_results_sub_steadystate$Base_sigma_c)
                    )

                    list_results_gene[[ki]] <- r_i
                    ki <- ki + 1
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  # Return NULL if no valid rows were produced for this gene.
  if (length(list_results_gene) == 0L) {
    return(NULL)
  } else {
    return(list_results_gene)
  }
}


# -----------------------------------------------------------------------------
# Parallel execution across genes.
#
# The number of used cores is conservatively reduced by 90 relative to the
# machine total; this keeps the script lightweight on large shared systems,
# while still guaranteeing at least one core.
# -----------------------------------------------------------------------------
n_cores <- max(1L, detectCores() - 90)


# -----------------------------------------------------------------------------
# Quick timing check on one gene to estimate total runtime.
# -----------------------------------------------------------------------------
currentTs <- Sys.time()
run_for_gene(1)
elapsed <- Sys.time() - currentTs
elapsed * length(genes) / n_cores


# -----------------------------------------------------------------------------
# Launch the full benchmark in parallel.
# -----------------------------------------------------------------------------
res_par <- mclapply(genes, run_for_gene, mc.cores = n_cores)


# -----------------------------------------------------------------------------
# Combine all per-gene outputs into one flat results table.
# -----------------------------------------------------------------------------
res_par <- Filter(Negate(is.null), res_par)
res_flat <- unlist(res_par, recursive = FALSE)

# Total number of tests executed.
length(res_flat)

dt_test_results <- data.table::rbindlist(res_flat, use.names = TRUE, fill = TRUE)


# -----------------------------------------------------------------------------
# Standardize factor levels for downstream summaries and plotting.
# -----------------------------------------------------------------------------
dt_test_results$Exprs_noise <- factor(
  dt_test_results$Exprs_noise,
  levels = c("Very low", "Low", "Medium", "High"),
  ordered = TRUE
)

dt_test_results$Perturbation <- factor(
  dt_test_results$Perturbation,
  levels = c("NONE", "SHUTOFF")
)

dt_test_results$Platform <- factor(
  dt_test_results$Platform,
  levels = c("RT-qPCR", "GAUSS", "RNA-seq")
)

dt_test_results$Test_method <- factor(dt_test_results$Test_method)


# -----------------------------------------------------------------------------
# Multiple-testing correction and evidence summaries.
#
# q-values are computed separately within each benchmark stratum defined by:
#   Lambda, N_replicates, Tsteps, Test_method, N_tsamples,
#   Platform, Exprs_noise, Perturbation
#
# Additional transformed evidence measures are then derived.
# -----------------------------------------------------------------------------
eps = 1e-10

dt_test_results[, q.value := p.adjust(p.value, method="fdr"),
  by = .(Lambda, N_replicates, Tsteps, Test_method, N_tsamples, Platform, Exprs_noise, Perturbation)
]

dt_test_results[, neglogq := -log10(pmax(q.value, eps))]
dt_test_results[, neglogp := -log10(pmax(p.value, eps))]

# Cap extreme significance values to avoid a few points dominating ranking.
cap <- 6
dt_test_results[, neglogq_cap := pmin(neglogq, cap)]


# -----------------------------------------------------------------------------
# Model-improvement summaries derived from RSS values.
#
# DeltaRSS:
#   absolute reduction in residual sum of squares.
#
# IR:
#   improvement ratio in [0, 1], normalized by RSS0.
#
# logRSSratio:
#   log-scale fit improvement diagnostic.
# -----------------------------------------------------------------------------
dt_test_results[, DeltaRSS := pmax(RSS0 - RSS1, 0)]
dt_test_results[, IR := fifelse(RSS0 > 0, DeltaRSS / RSS0, 0)]
dt_test_results[, logRSSratio := fifelse(RSS1 > 0, log(RSS0 / RSS1), 0)]


# -----------------------------------------------------------------------------
# Robust scaling for Sigma.
#
# Within each benchmark stratum, s0 is defined as the median absolute Sigma.
# This allows the construction of a log-scaled version of Sigma that is
# less sensitive to outliers.
# -----------------------------------------------------------------------------
dt_test_results[, s0 := median(abs(Sigma), na.rm = TRUE),
  by = .(Lambda, N_replicates, Tsteps, Test_method, N_tsamples, Platform, Exprs_noise, Perturbation)
]

dt_test_results[s0 == 0, s0 := 1]

dt_test_results[, sig_log := log1p(abs(Sigma + eps) / s0)]


# -----------------------------------------------------------------------------
# Composite prioritization scores.
#
# score_siglogq:
#   combines raw Sigma with significance evidence.
#
# score_sig_IR_logq:
#   recommended score combining:
#     - estimated effect size (Sigma),
#     - fit improvement (IR),
#     - statistical evidence (-log10 q-value).
#
# score_siglog_IR_logq:
#   more robust variant replacing Sigma with a log-scaled version.
# -----------------------------------------------------------------------------
dt_test_results[, score_siglogq := Sigma * neglogq_cap]
dt_test_results[, score_sig_IR_logq := Sigma * IR * neglogq_cap]
dt_test_results[, score_siglog_IR_logq := sig_log * IR * neglogq_cap]


# -----------------------------------------------------------------------------
# Quick global summary of the final benchmark table.
# -----------------------------------------------------------------------------
summary(dt_test_results)


# -----------------------------------------------------------------------------
# Save the full benchmark results.
# -----------------------------------------------------------------------------
save(dt_test_results, file = "test_sigma_2k_new.rdata")