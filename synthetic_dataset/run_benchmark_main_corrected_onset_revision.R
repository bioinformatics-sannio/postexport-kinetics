# =============================================================================
# Title:
# Main Synthetic Benchmark for the Bootstrap-Calibrated Weighted NNLS Test
#
# Description:
#   This script performs the principal synthetic-data benchmark used to
#   evaluate the bootstrap-calibrated weighted NNLS model-comparison framework
#   for post-export RNA conversion.
#
#   The benchmark evaluates two principal experimental regimes:
#
#     1. NONE
#        Continuous transcription with no intervention.
#
#     2. SHUTOFF
#        Complete transcriptional shutoff at T_star.
#
#   The benchmark spans:
#     - true-null genes (sigma_c = 0);
#     - true-positive genes (sigma_c > 0);
#     - RT-qPCR, Gaussian, and RNA-seq observation models;
#     - multiple measurement-noise regimes;
#     - multiple numbers of destructive biological replicates;
#     - multiple numbers of sampled time points;
#     - multiple temporal spacings after shutoff.
#
# Statistical method:
#   The inference procedure:
#     - treats biological replicates as destructive independent samples at
#       each time point;
#     - estimates within-time multivariate covariance;
#     - propagates the full covariance of interval differences;
#     - accounts for uncertainty in both A and b through replicate-level
#       parametric bootstrap datasets generated under H0;
#     - reconstructs A*, b*, and Sigma_b* at every bootstrap iteration;
#     - retains non-negative least-squares estimation;
#     - uses the deterministic add-one bootstrap p-value
#
#           p = (1 + sum(T* >= Tobs)) / (B + 1).
#
# Primary objectives:
#   1. Empirical Type-I error under the null.
#   2. Statistical power under sigma_c > 0.
#   3. Parameter recovery and practical identifiability.
#   4. Calibration under realistic sparse experimental designs.
#   5. Quantification of the benefits and limitations of transcriptional
#      shutoff.
#
# Restart and checkpoint strategy:
#   This implementation is explicitly designed for long HPC runs.
#
#   Existing legacy batch checkpoints:
#
#       benchmark_main_batch_001.rdata
#       benchmark_main_batch_002.rdata
#       ...
#
#   are automatically inspected and their completed genes are recovered.
#
#   Newly computed genes are checkpointed INDIVIDUALLY:
#
#       benchmark_main_checkpoints/gene_000401.rds
#       benchmark_main_checkpoints/gene_000402.rds
#       ...
#
#   A gene is considered complete only if its checkpoint contains the full
#   expected number of benchmark scenarios.
#
#   Importantly, a worker writes its gene checkpoint BEFORE returning its
#   small completion message to the PSOCK master. Therefore, if communication
#   between workers and the master becomes blocked, completed calculations are
#   still recoverable on restart.
#
#   PSOCK workers return only lightweight status records; the large per-gene
#   benchmark tables are never transmitted back through the socket connection.
#   This substantially reduces the risk of the socket-buffer deadlock observed
#   with the previous batch-level implementation.
#
# Scientific scope:
#   sigma_c represents a generic post-export conversion rate in the kinetic
#   model. Positive evidence should be interpreted as a kinetic pattern
#   consistent with an additional post-export conversion component, rather
#   than experimental confirmation of a particular molecular mechanism.
#
# Input:
#   ode_states_2k_20p_corrected_onset.rdata
#
# Expected object:
#   ode_states
#
# Perturbation labels used by this benchmark:
#   NONE
#   SHUTOFF
#   SteadyState
#
# Main outputs:
#   benchmark_main_revision_raw.rdata
#   benchmark_main_revision_raw.tsv
#   benchmark_main_revision_summary.rdata
#   benchmark_main_revision_summary.tsv
#   benchmark_main_sessionInfo.txt
#
# Checkpoint output:
#   benchmark_main_checkpoints/
#
# Progress log:
#   benchmark_main_progress.tsv
#
# Author:
#   Luigi Cerulo
#
# Copyright:
#   Copyright (c) 2026 Luigi Cerulo
#
# License:
#   Permission is hereby granted to use, copy, modify, and distribute this
#   software for academic and research purposes, provided that this copyright
#   notice and permission notice are retained in all copies or substantial
#   portions of the software.
#
# Disclaimer:
#   This software is provided "as is", without warranty of any kind, express
#   or implied, including but not limited to the warranties of merchantability,
#   fitness for a particular purpose, and noninfringement. In no event shall
#   the author be liable for any claim, damages, or other liability arising
#   from, out of, or in connection with the software or its use.
#
# Version:
#   Major-revision benchmark, restart-safe implementation, 2026.
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/synthetic_dataset")

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

library(data.table)
library(parallel)

data.table::setDTthreads(1)

source("../commons/nested_test2.r")
source("../commons/psi_test.r")
source("../commons/platforms.r")


# =============================================================================
# 1. User-configurable settings
# =============================================================================

N_BOOT <- 1999L

# Version tag prevents accidental reuse of checkpoints/results generated from
# the historical observation-shift simulator.
BENCHMARK_VERSION <- "corrected_onset_v1"

INPUT_FILE <- "ode_states_2k_20p_corrected_onset.rdata"

# Set TRUE only if you deliberately want to inspect/reuse historical batch
# files from exactly the same simulator/version. For this revision rerun it
# MUST remain FALSE.
ALLOW_LEGACY_RECOVERY <- FALSE

# Number of simultaneously active PSOCK workers.
N_WORKERS <- 100L

# Number of genes submitted to one parLapplyLB call.
#
# This is NOT the checkpoint size: checkpoints are gene-specific.
#
# Keeping this approximately equal to N_WORKERS gives one initial gene to
# each worker while preventing one single PSOCK call from covering the entire
# 2000-gene experiment.
TASK_BATCH_SIZE <- 100L

# Socket inactivity timeout.
#
# The previous cluster used the default 30-day timeout. A blocked worker could
# therefore stall the benchmark almost indefinitely.
#
# Twelve hours is deliberately conservative relative to the normal runtime of
# one gene while still allowing a genuinely dead connection to terminate.
PSOCK_TIMEOUT_SECONDS <- 12L * 60L * 60L

CHECKPOINT_DIR <- "benchmark_main_corrected_onset_checkpoints"

PROGRESS_FILE <- "benchmark_main_corrected_onset_progress.tsv"

RAW_RDATA_FILE <- "benchmark_main_corrected_onset_raw.rdata"
RAW_TSV_FILE <- "benchmark_main_corrected_onset_raw.tsv"
SUMMARY_RDATA_FILE <- "benchmark_main_corrected_onset_summary.rdata"
SUMMARY_TSV_FILE <- "benchmark_main_corrected_onset_summary.tsv"
SESSION_FILE <- "benchmark_main_corrected_onset_sessionInfo.txt"
RUN_METADATA_FILE <- "benchmark_main_corrected_onset_run_metadata.tsv"


# =============================================================================
# 2. Load synthetic data
# =============================================================================

cat("\nLoading synthetic dataset...\n")

load(INPUT_FILE)

stopifnot(
  exists("ode_states"),
  data.table::is.data.table(ode_states)
)


# -----------------------------------------------------------------------------
# PRE-FLIGHT VALIDATION OF THE CORRECTED SYNTHETIC DATASET
# -----------------------------------------------------------------------------

required_input_columns <- c(
  "Gene",
  "Perturbation",
  "truth_pos",
  "Base_R",
  "Base_tau",
  "Base_tau_s",
  "Base_sigma_c",
  "Base_sigma_n",
  "Base_alpha",
  "Base_alpha_s",
  "N_time_samples",
  "N_replicates",
  "T_step",
  "T_star",
  "Post_R_fraction",
  "Onset_shift",
  "Onset_time",
  "time",
  "replicate",
  "N",
  "N_s",
  "C",
  "C_s"
)

missing_input_columns <- setdiff(
  required_input_columns,
  names(ode_states)
)

if (
  length(missing_input_columns) > 0L
) {
  stop(
    paste0(
      "Corrected synthetic dataset is missing required columns:\n  ",
      paste(
        missing_input_columns,
        collapse = "\n  "
      )
    )
  )
}

# One transcriptional onset per gene across ALL designs and perturbations.
timing_by_gene <- ode_states[
  ,
  .(
    n_onset_shift = uniqueN(Onset_shift),
    n_onset_time = uniqueN(Onset_time),
    n_T_star = uniqueN(T_star)
  ),
  by = Gene
]

if (
  any(
    timing_by_gene$n_onset_shift != 1L |
      timing_by_gene$n_onset_time != 1L |
      timing_by_gene$n_T_star != 1L
  )
) {
  stop(
    "Input timing validation failed: every gene must have exactly one onset and one common T_star."
  )
}

# Complete SHUTOFF must really encode rho = 0.
bad_shutoff_rho <- ode_states[
  Perturbation == "SHUTOFF" &
    (
      !is.finite(Post_R_fraction) |
        abs(Post_R_fraction) > 1e-12
    )
]

if (
  nrow(bad_shutoff_rho) > 0L
) {
  stop(
    "Input validation failed: SHUTOFF rows with Post_R_fraction != 0 were found."
  )
}

# NONE must represent uninterrupted transcription.
bad_none_rho <- ode_states[
  Perturbation == "NONE" &
    (
      !is.finite(Post_R_fraction) |
        abs(Post_R_fraction - 1) > 1e-12
    )
]

if (
  nrow(bad_none_rho) > 0L
) {
  stop(
    "Input validation failed: NONE rows with Post_R_fraction != 1 were found."
  )
}

# Ground truth must be invariant within a gene.
truth_by_gene <- ode_states[
  ,
  .(
    n_truth = uniqueN(truth_pos),
    n_sigma = uniqueN(Base_sigma_c)
  ),
  by = Gene
]

if (
  any(
    truth_by_gene$n_truth != 1L |
      truth_by_gene$n_sigma != 1L
  )
) {
  stop(
    "Input validation failed: truth_pos or Base_sigma_c varies within at least one gene."
  )
}

# Null/alternative labels must agree with sigma_c.
truth_mismatch <- unique(
  ode_states[
    ,
    .(
      Gene,
      truth_pos,
      Base_sigma_c
    )
  ]
)[
  (
    truth_pos == 0 &
      abs(Base_sigma_c) > 1e-12
  ) |
    (
      truth_pos == 1 &
        Base_sigma_c <= 0
    )
]

if (
  nrow(truth_mismatch) > 0L
) {
  stop(
    "Input validation failed: truth_pos is inconsistent with Base_sigma_c."
  )
}

cat(
  "\nCorrected-input preflight checks PASSED.\n",
  "Genes: ",
  uniqueN(ode_states$Gene),
  "\n",
  "Common T_star: ",
  paste(
    sort(unique(ode_states$T_star)),
    collapse = ", "
  ),
  "\n",
  "Onset shift range: ",
  min(ode_states$Onset_shift),
  " to ",
  max(ode_states$Onset_shift),
  "\n",
  sep = ""
)


# =============================================================================
# 3. Benchmark design
# =============================================================================

RANGE_PLATFORM <- c(
  "RT-qPCR",
  "GAUSS",
  "RNA-seq"
)

RANGE_PERTURBATION <- c(
  "NONE",
  "SHUTOFF"
)

RANGE_NOISE <- c(
  "Very low",
  "Low",
  "Medium",
  "High"
)

RANGE_N_TIME <- sort(
  unique(
    ode_states$N_time_samples
  )
)

RANGE_N_REP <- sort(
  unique(
    ode_states$N_replicates
  )
)

RANGE_TSTEP <- sort(
  unique(
    ode_states$T_step
  )
)

SCALING_A <- TRUE

LAMBDA_TIME <- 0.5
LAMBDA_DIAG <- 0.1
REL_FLOOR <- 1e-8

TRUNCATE_NONNEGATIVE_BOOT <- FALSE
MAX_FAILURE_RATE <- 0.05


# Expected number of benchmark scenarios for one gene.

EXPECTED_TESTS_PER_GENE <- (
  length(RANGE_PLATFORM) *
    length(RANGE_PERTURBATION) *
    length(RANGE_NOISE) *
    length(RANGE_N_TIME) *
    length(RANGE_N_REP) *
    length(RANGE_TSTEP)
)

cat(
  "Expected tests per gene:",
  EXPECTED_TESTS_PER_GENE,
  "\n"
)


# =============================================================================
# 4. Gaussian noise settings
# =============================================================================

RANGE_GAUSS_NOISE <- c(
  "Very low" = 0.02,
  "Low"      = 0.05,
  "Medium"   = 0.10,
  "High"     = 0.20
)


# =============================================================================
# 5. Restrict synthetic dataset to required perturbations
# =============================================================================

ode_states_main <- ode_states[
  Perturbation %in% c(
    "NONE",
    "SHUTOFF",
    "SteadyState"
  )
]

setkey(
  ode_states_main,
  Gene,
  N_time_samples,
  N_replicates,
  T_step,
  Perturbation
)

genes <- sort(
  unique(
    ode_states_main$Gene
  )
)

N_GENES_TOTAL <- length(genes)

cat(
  "Total genes:",
  N_GENES_TOTAL,
  "\n"
)

rm(ode_states)
gc()


# =============================================================================
# 6. Checkpoint directory
# =============================================================================

if (!dir.exists(CHECKPOINT_DIR)) {

  dir.create(
    CHECKPOINT_DIR,
    recursive = TRUE
  )
}

CHECKPOINT_DIR <- normalizePath(
  CHECKPOINT_DIR,
  mustWork = TRUE
)


# =============================================================================
# 7. Platform simulation helper
# =============================================================================

add_platform_noise_main <- function(
  dt,
  platform,
  noise
) {

  dt <- copy(
    as.data.table(dt)
  )

  targets <- c(
    "N",
    "C",
    "C_s",
    "N_s"
  )


  if (platform == "GAUSS") {

    noise_sd <- unname(
      RANGE_GAUSS_NOISE[noise]
    )

    return(
      add_gaussian_noise(
        dt,
        cols = targets,
        noise_sd = noise_sd
      )
    )
  }


  if (platform == "RT-qPCR") {

    ct_sd <- switch(
      noise,
      "Very low" = 0.01,
      "Low"      = 0.05,
      "Medium"   = 0.10,
      "High"     = 0.25
    )

    scale_copies <- switch(
      noise,
      "Very low" = 20,
      "Low"      = 15,
      "Medium"   = 10,
      "High"     = 5
    )

    return(
      simulate_rt_qpcr(
        dt,
        targets = targets,
        ct_sd = ct_sd,
        scale_copies = scale_copies
      )
    )
  }


  if (platform == "RNA-seq") {

    scale_counts <- switch(
      noise,
      "Very low" = 10000,
      "Low"      = 5000,
      "Medium"   = 1000,
      "High"     = 200
    )

    mean_disp <- switch(
      noise,
      "Very low" = 0.01,
      "Low"      = 0.05,
      "Medium"   = 0.10,
      "High"     = 0.25
    ) * 0.25

    cv_disp <- switch(
      noise,
      "Very low" = 0.5,
      "Low"      = 0.7,
      "Medium"   = 0.8,
      "High"     = 1.0
    )

    return(
      simulate_rnaseq(
        dt,
        targets = targets,
        scale_counts = scale_counts,
        mean_disp = mean_disp,
        cv_disp = cv_disp
      )
    )
  }


  stop(
    paste(
      "Unknown platform:",
      platform
    )
  )
}


# =============================================================================
# 8. Reproducible seed helper
# =============================================================================

# Stable deterministic hash.
#
# Unlike the previous arithmetic seed construction, this avoids accidental
# collisions between different combinations of platform, noise, perturbation,
# and experimental design.

stable_seed_from_string <- function(x) {

  ints <- utf8ToInt(
    enc2utf8(x)
  )

  h <- 104729

  for (ii in ints) {

    h <- (
      h * 1000003 +
        ii
    ) %% 2000000000
  }

  seed <- as.integer(
    h
  )

  if (
    !is.finite(seed) ||
      seed <= 0L
  ) {
    seed <- 1L
  }

  seed
}


make_main_seed <- function(
  gene,
  ntime,
  nrep,
  tstep,
  perturbation,
  platform,
  noise,
  stream = c(
    "measurement",
    "bootstrap"
  )
) {

  stream <- match.arg(
    stream
  )

  key <- paste(
    gene,
    ntime,
    nrep,
    tstep,
    perturbation,
    platform,
    noise,
    stream,
    sep = "|"
  )

  stable_seed_from_string(
    key
  )
}


# =============================================================================
# 9. Safe extraction helpers
# =============================================================================

safe_value <- function(
  x,
  field,
  default = NA_real_
) {

  if (
    is.null(x) ||
      is.null(x[[field]]) ||
      length(x[[field]]) == 0L
  ) {

    return(
      default
    )
  }

  x[[field]][1]
}


safe_coef <- function(
  x,
  name
) {

  if (
    is.null(x) ||
      is.null(names(x)) ||
      !(name %in% names(x))
  ) {

    return(
      NA_real_
    )
  }

  unname(
    x[name]
  )
}


# =============================================================================
# 10. Worker computation for one gene
# =============================================================================

run_gene_main <- function(
  gene_i,
  N_boot
) {

  out <- list()

  kk <- 1L


  for (ntime_i in RANGE_N_TIME) {

    for (nrep_i in RANGE_N_REP) {

      for (tstep_i in RANGE_TSTEP) {


        dt_design <- ode_states_main[
          Gene == gene_i &
            N_time_samples == ntime_i &
            N_replicates == nrep_i &
            T_step == tstep_i
        ]


        if (nrow(dt_design) == 0L) {
          next
        }


        dt_ss <- dt_design[
          Perturbation == "SteadyState"
        ]


        if (nrow(dt_ss) == 0L) {
          next
        }


        for (perturb_i in RANGE_PERTURBATION) {


          dt_tc <- dt_design[
            Perturbation == perturb_i
          ]


          if (nrow(dt_tc) == 0L) {
            next
          }


          truth_i <- unique(
            dt_tc$truth_pos
          )[1]


          sigma_true <- unique(
            dt_tc$Base_sigma_c
          )[1]


          if (perturb_i == "SHUTOFF") {

            tstar <- unique(
              dt_tc$T_star
            )

            tstar <- tstar[
              is.finite(tstar)
            ]

            if (length(tstar) != 1L) {
              next
            }

            t_star_fit <- tstar[1]

          } else {

            t_star_fit <- NULL
          }


          for (platform_i in RANGE_PLATFORM) {

            for (noise_i in RANGE_NOISE) {


              measurement_seed <- make_main_seed(
                gene = gene_i,
                ntime = ntime_i,
                nrep = nrep_i,
                tstep = tstep_i,
                perturbation = perturb_i,
                platform = platform_i,
                noise = noise_i,
                stream = "measurement"
              )


              bootstrap_seed <- make_main_seed(
                gene = gene_i,
                ntime = ntime_i,
                nrep = nrep_i,
                tstep = tstep_i,
                perturbation = perturb_i,
                platform = platform_i,
                noise = noise_i,
                stream = "bootstrap"
              )


              set.seed(
                measurement_seed
              )


              dt_noisy <- tryCatch(

                add_platform_noise_main(
                  dt = dt_tc,
                  platform = platform_i,
                  noise = noise_i
                ),

                error = function(e) NULL
              )


              if (is.null(dt_noisy)) {
                next
              }


              test_start <- Sys.time()


              WW <- tryCatch(

                test_sigma_nested(

                  tsampled_data =
                    dt_noisy,

                  scaling_A =
                    SCALING_A,

                  t_star =
                    t_star_fit,

                  B_n =
                    N_boot,

                  seed =
                    bootstrap_seed,

                  lambda_time =
                    LAMBDA_TIME,

                  lambda_diag =
                    LAMBDA_DIAG,

                  rel_floor =
                    REL_FLOOR,

                  truncate_nonnegative_boot =
                    TRUNCATE_NONNEGATIVE_BOOT,

                  max_failure_rate =
                    MAX_FAILURE_RATE,

                  return_boot =
                    FALSE,

                  verbose =
                    FALSE
                ),

                error = function(e) {

                  list(
                    p.value = NA_real_,
                    status = "test_error",
                    error.message = conditionMessage(e)
                  )
                }
              )


              elapsed_seconds <- as.numeric(
                difftime(
                  Sys.time(),
                  test_start,
                  units = "secs"
                )
              )


              coef_full_i <- WW$coef_full
              coef_null_i <- WW$coef_null


              out[[kk]] <- data.frame(

                Gene =
                  gene_i,

                Positive =
                  truth_i,

                sigma_true =
                  sigma_true,

                Platform =
                  platform_i,

                Perturbation =
                  perturb_i,

                Exprs_noise =
                  noise_i,

                N_tsamples =
                  ntime_i,

                N_replicates =
                  nrep_i,

                Tsteps =
                  tstep_i,

                T_star =
                  ifelse(
                    is.null(t_star_fit),
                    NA_real_,
                    t_star_fit
                  ),

                N_boot =
                  N_boot,

                Benchmark_version =
                  BENCHMARK_VERSION,

                Onset_shift =
                  unique(
                    dt_tc$Onset_shift
                  )[1],

                Onset_time =
                  unique(
                    dt_tc$Onset_time
                  )[1],

                Post_R_fraction =
                  unique(
                    dt_tc$Post_R_fraction
                  )[1],


                # -------------------------------------------------------------
                # Inference
                # -------------------------------------------------------------

                p.value =
                  safe_value(
                    WW,
                    "p.value"
                  ),

                T.obs =
                  safe_value(
                    WW,
                    "T.obs"
                  ),

                Sigma =
                  safe_value(
                    WW,
                    "Sigma"
                  ),

                Alpha =
                  safe_value(
                    WW,
                    "Alpha"
                  ),

                RSS0 =
                  safe_value(
                    WW,
                    "RSS0"
                  ),

                RSS1 =
                  safe_value(
                    WW,
                    "RSS1"
                  ),

                IR =
                  safe_value(
                    WW,
                    "IR"
                  ),


                # -------------------------------------------------------------
                # Boundary diagnostics
                # -------------------------------------------------------------

                atom_zero =
                  safe_value(
                    WW,
                    "atom.zero"
                  ),


                # -------------------------------------------------------------
                # Practical-identifiability diagnostics
                # -------------------------------------------------------------

                condition_number =
                  safe_value(
                    WW,
                    "condition.number"
                  ),

                min_singular_value =
                  safe_value(
                    WW,
                    "min.singular.value"
                  ),

                rank_full =
                  safe_value(
                    WW,
                    "rank.full"
                  ),

                rank_null =
                  safe_value(
                    WW,
                    "rank.null"
                  ),


                # -------------------------------------------------------------
                # Bootstrap stability
                # -------------------------------------------------------------

                bootstrap_failure_rate =
                  safe_value(
                    WW,
                    "bootstrap.failure.rate"
                  ),

                n_bootstrap_valid =
                  safe_value(
                    WW,
                    "n.bootstrap.valid"
                  ),

                boot_condition_median =
                  safe_value(
                    WW,
                    "bootstrap.condition.median"
                  ),

                boot_condition_q95 =
                  safe_value(
                    WW,
                    "bootstrap.condition.q95"
                  ),

                boot_condition_max =
                  safe_value(
                    WW,
                    "bootstrap.condition.max"
                  ),

                boot_rank_deficient_fraction =
                  safe_value(
                    WW,
                    "bootstrap.rank.deficient.fraction"
                  ),


                # -------------------------------------------------------------
                # Full-model fitted parameters
                # -------------------------------------------------------------

                R_hat =
                  safe_coef(
                    coef_full_i,
                    "R"
                  ),

                tau_hat =
                  safe_coef(
                    coef_full_i,
                    "tau"
                  ),

                tau_s_hat =
                  safe_coef(
                    coef_full_i,
                    "tau_s"
                  ),

                sigma_c_hat =
                  safe_coef(
                    coef_full_i,
                    "sigma_c"
                  ),

                sigma_n_hat =
                  safe_coef(
                    coef_full_i,
                    "sigma_n"
                  ),

                alpha_hat =
                  safe_coef(
                    coef_full_i,
                    "alpha"
                  ),

                alpha_s_hat =
                  safe_coef(
                    coef_full_i,
                    "alpha_s"
                  ),


                # -------------------------------------------------------------
                # Null-model fitted parameters
                # -------------------------------------------------------------

                R_hat_null =
                  safe_coef(
                    coef_null_i,
                    "R"
                  ),

                tau_hat_null =
                  safe_coef(
                    coef_null_i,
                    "tau"
                  ),

                tau_s_hat_null =
                  safe_coef(
                    coef_null_i,
                    "tau_s"
                  ),

                sigma_n_hat_null =
                  safe_coef(
                    coef_null_i,
                    "sigma_n"
                  ),

                alpha_hat_null =
                  safe_coef(
                    coef_null_i,
                    "alpha"
                  ),

                alpha_s_hat_null =
                  safe_coef(
                    coef_null_i,
                    "alpha_s"
                  ),


                # -------------------------------------------------------------
                # Ground-truth parameters
                # -------------------------------------------------------------

                R_true =
                  unique(
                    dt_tc$Base_R
                  )[1],

                tau_true =
                  unique(
                    dt_tc$Base_tau
                  )[1],

                tau_s_true =
                  unique(
                    dt_tc$Base_tau_s
                  )[1],

                sigma_n_true =
                  unique(
                    dt_tc$Base_sigma_n
                  )[1],

                alpha_true =
                  unique(
                    dt_tc$Base_alpha
                  )[1],

                alpha_s_true =
                  unique(
                    dt_tc$Base_alpha_s
                  )[1],


                # -------------------------------------------------------------
                # Reproducibility
                # -------------------------------------------------------------

                measurement_seed =
                  measurement_seed,

                bootstrap_seed =
                  bootstrap_seed,

                test_seconds =
                  elapsed_seconds,

                status =
                  if (
                    is.null(WW$status)
                  ) {
                    NA_character_
                  } else {
                    as.character(
                      WW$status
                    )
                  },

                error_message =
                  if (
                    is.null(WW$error.message)
                  ) {
                    NA_character_
                  } else {
                    as.character(
                      WW$error.message
                    )
                  },

                stringsAsFactors =
                  FALSE
              )


              kk <- kk + 1L
            }
          }
        }
      }
    }
  }


  if (length(out) == 0L) {
    return(NULL)
  }


  data.table::rbindlist(
    out,
    use.names = TRUE,
    fill = TRUE
  )
}


# =============================================================================
# 11. Per-gene checkpoint helper
# =============================================================================

gene_checkpoint_file <- function(
  gene_i
) {

  file.path(
    CHECKPOINT_DIR,
    sprintf(
      "gene_%06d.rds",
      as.integer(gene_i)
    )
  )
}


valid_gene_result <- function(
  dt,
  gene_i,
  expected_tests = EXPECTED_TESTS_PER_GENE
) {

  if (
    is.null(dt) ||
      !data.table::is.data.table(dt)
  ) {
    return(FALSE)
  }


  required_checkpoint_columns <- c(
    "Gene",
    "Benchmark_version",
    "Onset_shift",
    "Onset_time",
    "Post_R_fraction"
  )

  if (
    any(
      !required_checkpoint_columns %in%
        names(dt)
    )
  ) {
    return(FALSE)
  }

  if (
    data.table::uniqueN(
      dt$Benchmark_version
    ) != 1L ||
      unique(
        dt$Benchmark_version
      )[1] != BENCHMARK_VERSION
  ) {
    return(FALSE)
  }

  if (
    data.table::uniqueN(
      dt$Onset_shift
    ) != 1L ||
      data.table::uniqueN(
        dt$Onset_time
      ) != 1L
  ) {
    return(FALSE)
  }


  if (
    nrow(dt) != expected_tests
  ) {
    return(FALSE)
  }


  if (
    length(unique(dt$Gene)) != 1L ||
      unique(dt$Gene)[1] != gene_i
  ) {
    return(FALSE)
  }


  # Verify uniqueness of benchmark design combinations.

  scenario_key <- paste(
    dt$Platform,
    dt$Perturbation,
    dt$Exprs_noise,
    dt$N_tsamples,
    dt$N_replicates,
    dt$Tsteps,
    sep = "|"
  )


  if (
    data.table::uniqueN(scenario_key) != expected_tests
  ) {
    return(FALSE)
  }


  TRUE
}


run_gene_and_checkpoint <- function(
  gene_i,
  N_boot
) {

  pid <- Sys.getpid()

  t0 <- Sys.time()

  checkpoint_file <- gene_checkpoint_file(
    gene_i
  )


  # ---------------------------------------------------------------------------
  # If a valid checkpoint already exists, do not recompute.
  # ---------------------------------------------------------------------------

  if (file.exists(checkpoint_file)) {

    existing <- tryCatch(
      readRDS(
        checkpoint_file
      ),
      error = function(e) NULL
    )


    if (
      valid_gene_result(
        existing,
        gene_i
      )
    ) {

      return(
        data.frame(
          Gene = gene_i,
          pid = pid,
          status = "already_checkpointed",
          n_tests = nrow(existing),
          elapsed_seconds = 0,
          checkpoint = checkpoint_file,
          stringsAsFactors = FALSE
        )
      )
    }
  }


  # ---------------------------------------------------------------------------
  # Compute gene.
  # ---------------------------------------------------------------------------

  result <- tryCatch(

    run_gene_main(
      gene_i = gene_i,
      N_boot = N_boot
    ),

    error = function(e) {

      structure(
        list(
          error_message =
            conditionMessage(e)
        ),
        class =
          "gene_benchmark_error"
      )
    }
  )


  if (
    inherits(
      result,
      "gene_benchmark_error"
    )
  ) {

    return(
      data.frame(
        Gene = gene_i,
        pid = pid,
        status = "gene_error",
        n_tests = 0L,
        elapsed_seconds =
          as.numeric(
            difftime(
              Sys.time(),
              t0,
              units = "secs"
            )
          ),
        checkpoint = NA_character_,
        error_message = result$error_message,
        stringsAsFactors = FALSE
      )
    )
  }


  if (
    !valid_gene_result(
      result,
      gene_i
    )
  ) {

    return(
      data.frame(
        Gene = gene_i,
        pid = pid,
        status = "incomplete_gene_result",
        n_tests =
          if (
            is.null(result)
          ) {
            0L
          } else {
            nrow(result)
          },
        elapsed_seconds =
          as.numeric(
            difftime(
              Sys.time(),
              t0,
              units = "secs"
            )
          ),
        checkpoint = NA_character_,
        error_message =
          sprintf(
            "Expected %d tests for gene %d.",
            EXPECTED_TESTS_PER_GENE,
            gene_i
          ),
        stringsAsFactors = FALSE
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Atomic checkpoint.
  #
  # Write first to a temporary file. Only a successfully completed file is
  # renamed to the canonical checkpoint name.
  # ---------------------------------------------------------------------------

  tmp_file <- paste0(
    checkpoint_file,
    ".tmp.",
    pid
  )


  saveRDS(
    result,
    file = tmp_file,
    compress = "gzip"
  )


  rename_ok <- file.rename(
    tmp_file,
    checkpoint_file
  )


  if (!rename_ok) {

    if (file.exists(tmp_file)) {
      unlink(tmp_file)
    }

    return(
      data.frame(
        Gene = gene_i,
        pid = pid,
        status = "checkpoint_write_failed",
        n_tests = nrow(result),
        elapsed_seconds =
          as.numeric(
            difftime(
              Sys.time(),
              t0,
              units = "secs"
            )
          ),
        checkpoint = NA_character_,
        error_message =
          "Unable to atomically rename checkpoint file.",
        stringsAsFactors = FALSE
      )
    )
  }


  elapsed <- as.numeric(
    difftime(
      Sys.time(),
      t0,
      units = "secs"
    )
  )


  message(
    sprintf(
      "[PID %d] CHECKPOINT SAVED gene %d with %d tests (%.1f min)",
      pid,
      gene_i,
      nrow(result),
      elapsed / 60
    )
  )


  # ---------------------------------------------------------------------------
  # IMPORTANT:
  # Return only a tiny status record to the PSOCK master.
  #
  # The large benchmark result remains on disk.
  # ---------------------------------------------------------------------------

  data.frame(
    Gene = gene_i,
    pid = pid,
    status = "checkpoint_saved",
    n_tests = nrow(result),
    elapsed_seconds = elapsed,
    checkpoint = checkpoint_file,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# 12. Recover completed genes from legacy batch checkpoints
# =============================================================================

legacy_batch_files <- if (
  isTRUE(
    ALLOW_LEGACY_RECOVERY
  )
) {
  list.files(
    path = ".",
    pattern =
      "^benchmark_main_corrected_onset_batch_[0-9]+\\.rdata$",
    full.names = TRUE
  )
} else {
  character(0)
}


completed_from_legacy <- integer(0)


if (length(legacy_batch_files) > 0L) {

  cat(
    "\nInspecting",
    length(legacy_batch_files),
    "legacy batch checkpoint(s)...\n"
  )


  for (ff in legacy_batch_files) {

    env_tmp <- new.env(
      parent = emptyenv()
    )


    load_ok <- tryCatch(
      {

        load(
          ff,
          envir = env_tmp
        )

        TRUE
      },
      error = function(e) FALSE
    )


    if (!load_ok) {

      warning(
        paste(
          "Could not read legacy checkpoint:",
          ff
        )
      )

      next
    }


    if (!exists(
      "dt_bb",
      envir = env_tmp,
      inherits = FALSE
    )) {

      warning(
        paste(
          "Legacy checkpoint does not contain dt_bb:",
          ff
        )
      )

      next
    }


    dt_tmp <- data.table::as.data.table(
      get(
        "dt_bb",
        envir = env_tmp
      )
    )


    if (!"Gene" %in% names(dt_tmp)) {

      warning(
        paste(
          "Legacy checkpoint has no Gene column:",
          ff
        )
      )

      next
    }


    gene_counts <- dt_tmp[
      ,
      .N,
      by = Gene
    ]


    good_genes <- gene_counts[
      N == EXPECTED_TESTS_PER_GENE,
      Gene
    ]


    completed_from_legacy <- c(
      completed_from_legacy,
      good_genes
    )


    cat(
      basename(ff),
      ":",
      length(good_genes),
      "complete genes\n"
    )


    rm(
      env_tmp,
      dt_tmp,
      gene_counts
    )

    gc()
  }
}


completed_from_legacy <- sort(
  unique(
    as.integer(
      completed_from_legacy
    )
  )
)


# =============================================================================
# 13. Recover completed genes from per-gene checkpoints
# =============================================================================

checkpoint_files <- list.files(
  path = CHECKPOINT_DIR,
  pattern = "^gene_[0-9]+\\.rds$",
  full.names = TRUE
)


completed_from_gene_files <- integer(0)


if (length(checkpoint_files) > 0L) {

  cat(
    "\nInspecting",
    length(checkpoint_files),
    "per-gene checkpoint(s)...\n"
  )


  for (ff in checkpoint_files) {

    gene_from_file <- suppressWarnings(
      as.integer(
        sub(
          "^gene_([0-9]+)\\.rds$",
          "\\1",
          basename(ff)
        )
      )
    )


    if (!is.finite(gene_from_file)) {
      next
    }


    dt_tmp <- tryCatch(
      readRDS(
        ff
      ),
      error = function(e) NULL
    )


    if (
      valid_gene_result(
        dt_tmp,
        gene_from_file
      )
    ) {

      completed_from_gene_files <- c(
        completed_from_gene_files,
        gene_from_file
      )
    } else {

      warning(
        paste(
          "Invalid or incomplete gene checkpoint:",
          ff,
          "- this gene will be recomputed."
        )
      )
    }


    rm(
      dt_tmp
    )
  }
}


completed_from_gene_files <- sort(
  unique(
    completed_from_gene_files
  )
)


# =============================================================================
# 14. Determine genes still requiring computation
# =============================================================================

completed_genes <- sort(
  unique(
    c(
      completed_from_legacy,
      completed_from_gene_files
    )
  )
)


genes_remaining <- setdiff(
  genes,
  completed_genes
)


cat(
  "\n============================================================\n"
)

cat(
  "RESTART STATUS\n"
)

cat(
  "============================================================\n"
)

cat(
  "Total genes:                ",
  length(genes),
  "\n"
)

cat(
  "Recovered from old batches: ",
  length(completed_from_legacy),
  "\n"
)

cat(
  "Recovered gene checkpoints: ",
  length(completed_from_gene_files),
  "\n"
)

cat(
  "Total completed genes:      ",
  length(completed_genes),
  "\n"
)

cat(
  "Genes still to run:         ",
  length(genes_remaining),
  "\n"
)


if (length(completed_genes) > 0L) {

  cat(
    "Completed range:            ",
    min(completed_genes),
    "-",
    max(completed_genes),
    "\n"
  )
}


if (length(genes_remaining) > 0L) {

  cat(
    "First remaining genes:      ",
    paste(
      head(
        genes_remaining,
        20
      ),
      collapse = ", "
    ),
    "\n"
  )
}


cat(
  "============================================================\n\n"
)


# =============================================================================
# 14b. Fast preflight smoke test
# =============================================================================
#
# This deliberately uses a small bootstrap. Its purpose is NOT inferential;
# it checks that one complete gene traverses every design/platform/noise branch
# before the long HPC run starts.
#
SMOKE_BOOT <- 19L

smoke_gene <- genes[1]

cat(
  "\n============================================================\n",
  "PRE-RUN SMOKE TEST\n",
  "============================================================\n",
  "Gene: ",
  smoke_gene,
  "\n",
  "Bootstrap iterations/test: ",
  SMOKE_BOOT,
  "\n",
  sep = ""
)

smoke_result <- run_gene_main(
  gene_i = smoke_gene,
  N_boot = SMOKE_BOOT
)

if (
  is.null(smoke_result) ||
    nrow(smoke_result) != EXPECTED_TESTS_PER_GENE
) {
  stop(
    paste0(
      "Smoke test failed: expected ",
      EXPECTED_TESTS_PER_GENE,
      " rows but obtained ",
      if (is.null(smoke_result)) 0L else nrow(smoke_result),
      "."
    )
  )
}

smoke_key <- paste(
  smoke_result$Platform,
  smoke_result$Perturbation,
  smoke_result$Exprs_noise,
  smoke_result$N_tsamples,
  smoke_result$N_replicates,
  smoke_result$Tsteps,
  sep = "|"
)

if (
  uniqueN(smoke_key) != EXPECTED_TESTS_PER_GENE
) {
  stop(
    "Smoke test failed: duplicate or missing benchmark scenarios."
  )
}

if (
  uniqueN(smoke_result$Onset_time) != 1L ||
    uniqueN(smoke_result$Onset_shift) != 1L
) {
  stop(
    "Smoke test failed: onset metadata changed across scenarios for one gene."
  )
}

cat(
  "\nSmoke-test status counts:\n"
)

print(
  as.data.table(smoke_result)[
    ,
    .N,
    by = status
  ]
)

cat(
  "\nSMOKE TEST PASSED.\n"
)

rm(smoke_result)
gc()


# =============================================================================
# 15. PSOCK cluster and restart-safe execution
# =============================================================================

overall_start <- Sys.time()


if (length(genes_remaining) > 0L) {


  # ---------------------------------------------------------------------------
  # Create PSOCK cluster.
  # ---------------------------------------------------------------------------

  cat(
    "Starting",
    N_WORKERS,
    "PSOCK workers...\n"
  )


  cl <- parallel::makePSOCKcluster(
    N_WORKERS,
    outfile = "",
    timeout = PSOCK_TIMEOUT_SECONDS,
    setup_timeout = 120
  )


  # ---------------------------------------------------------------------------
  # Entire parallel phase enclosed in try/finally so workers are terminated
  # when the R script itself encounters a normal error.
  # ---------------------------------------------------------------------------

  tryCatch(

    {


      parallel::clusterEvalQ(
        cl,
        {

          Sys.setenv(
            OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1",
            MKL_NUM_THREADS = "1",
            VECLIB_MAXIMUM_THREADS = "1",
            NUMEXPR_NUM_THREADS = "1"
          )

          library(data.table)
          library(nnls)
          library(MASS)

          data.table::setDTthreads(1)

          NULL
        }
      )


      objects_to_export <- c(

        "ode_states_main",

        "RANGE_PLATFORM",
        "RANGE_PERTURBATION",
        "RANGE_NOISE",
        "RANGE_N_TIME",
        "RANGE_N_REP",
        "RANGE_TSTEP",

        "SCALING_A",
        "LAMBDA_TIME",
        "LAMBDA_DIAG",
        "REL_FLOOR",
        "TRUNCATE_NONNEGATIVE_BOOT",
        "MAX_FAILURE_RATE",

        "RANGE_GAUSS_NOISE",

        "EXPECTED_TESTS_PER_GENE",
        "CHECKPOINT_DIR",
        "BENCHMARK_VERSION",

        "add_platform_noise_main",
        "stable_seed_from_string",
        "make_main_seed",
        "safe_value",
        "safe_coef",
        "run_gene_main",
        "gene_checkpoint_file",
        "valid_gene_result",
        "run_gene_and_checkpoint",

        "add_gaussian_noise",
        "simulate_rt_qpcr",
        "simulate_rnaseq",
        "sample_dispersion_gamma",

        "KINETIC_VARS",
        "PARAM_NAMES",

        "make_spd",
        "inverse_sqrt_matrix",
        "time_summary_cov_shrink",
        "build_sigma_means",
        "build_difference_matrix",
        "build_Ab_fullcov",
        "fit_nnls_nested_once",
        "kinetic_matrix",
        "cn_interval",
        "predict_null_cn",
        "simulate_destructive_null",
        "test_sigma_nested"
      )


      parallel::clusterExport(
        cl,
        varlist = objects_to_export,
        envir = .GlobalEnv
      )


      # -----------------------------------------------------------------------
      # Split remaining genes into manageable PSOCK calls.
      #
      # Gene-level checkpoints are written inside each task, so even if the
      # master process becomes blocked before parLapplyLB returns, completed
      # genes remain recoverable.
      # -----------------------------------------------------------------------

      task_batches <- split(
        genes_remaining,
        ceiling(
          seq_along(genes_remaining) /
            TASK_BATCH_SIZE
        )
      )


      for (bb in seq_along(task_batches)) {


        genes_bb <- task_batches[[bb]]


        cat(
          "\n============================================================\n"
        )

        cat(
          "Running task batch",
          bb,
          "/",
          length(task_batches),
          "| genes",
          min(genes_bb),
          "-",
          max(genes_bb),
          "| N =",
          length(genes_bb),
          "\n"
        )

        cat(
          "============================================================\n"
        )


        batch_start <- Sys.time()


        status_bb <- parallel::parLapplyLB(

          cl,

          genes_bb,

          fun =
            run_gene_and_checkpoint,

          N_boot =
            N_BOOT
        )


        status_bb <- data.table::rbindlist(
          status_bb,
          use.names = TRUE,
          fill = TRUE
        )


        elapsed_batch <- as.numeric(
          difftime(
            Sys.time(),
            batch_start,
            units = "mins"
          )
        )


        status_bb[
          ,
          `:=`(
            task_batch = bb,
            task_batch_elapsed_minutes = elapsed_batch,
            master_timestamp = as.character(
              Sys.time()
            )
          )
        ]


        # ---------------------------------------------------------------------
        # Progress log.
        # ---------------------------------------------------------------------

        if (!file.exists(PROGRESS_FILE)) {

          data.table::fwrite(
            status_bb,
            file = PROGRESS_FILE,
            sep = "\t"
          )

        } else {

          data.table::fwrite(
            status_bb,
            file = PROGRESS_FILE,
            sep = "\t",
            append = TRUE,
            col.names = FALSE
          )
        }


        cat(
          "\nTask batch returned to master.\n"
        )

        print(
          status_bb[
            ,
            .N,
            by = status
          ]
        )


        cat(
          "Elapsed:",
          round(
            elapsed_batch,
            2
          ),
          "minutes\n"
        )


        # ---------------------------------------------------------------------
        # Count all currently durable completed genes.
        # ---------------------------------------------------------------------

        current_checkpoints <- list.files(
          path = CHECKPOINT_DIR,
          pattern = "^gene_[0-9]+\\.rds$",
          full.names = FALSE
        )


        cat(
          "Per-gene checkpoint files currently present:",
          length(current_checkpoints),
          "\n"
        )


        rm(
          status_bb
        )

        gc()
      }
    },


    finally = {

      cat(
        "\nStopping PSOCK cluster...\n"
      )

      try(
        parallel::stopCluster(
          cl
        ),
        silent = TRUE
      )
    }
  )

} else {

  cat(
    "All genes already have valid checkpointed results.\n"
  )
}


# =============================================================================
# 16. Re-scan checkpoints after computation
# =============================================================================

cat(
  "\nRe-scanning completed results...\n"
)


checkpoint_files <- list.files(
  path = CHECKPOINT_DIR,
  pattern = "^gene_[0-9]+\\.rds$",
  full.names = TRUE
)


valid_gene_checkpoint_map <- list()


for (ff in checkpoint_files) {

  gene_from_file <- suppressWarnings(
    as.integer(
      sub(
        "^gene_([0-9]+)\\.rds$",
        "\\1",
        basename(ff)
      )
    )
  )


  if (!is.finite(gene_from_file)) {
    next
  }


  dt_tmp <- tryCatch(
    readRDS(
      ff
    ),
    error = function(e) NULL
  )


  if (
    valid_gene_result(
      dt_tmp,
      gene_from_file
    )
  ) {

    valid_gene_checkpoint_map[[
      as.character(
        gene_from_file
      )
    ]] <- ff
  }


  rm(
    dt_tmp
  )
}


genes_checkpointed_final <- as.integer(
  names(
    valid_gene_checkpoint_map
  )
)


# =============================================================================
# 17. Load final results
# =============================================================================

cat(
  "\nCombining final benchmark table...\n"
)


result_list <- list()

rr <- 1L

genes_loaded <- integer(0)


# -----------------------------------------------------------------------------
# First load valid per-gene checkpoints.
#
# These take precedence over legacy batch results if both exist.
# -----------------------------------------------------------------------------

for (gene_i in sort(genes_checkpointed_final)) {

  ff <- valid_gene_checkpoint_map[[
    as.character(
      gene_i
    )
  ]]


  dt_gene <- readRDS(
    ff
  )


  result_list[[rr]] <- dt_gene

  rr <- rr + 1L

  genes_loaded <- c(
    genes_loaded,
    gene_i
  )
}


# -----------------------------------------------------------------------------
# Recover remaining genes from legacy batch checkpoints.
# -----------------------------------------------------------------------------

genes_needed_from_legacy <- setdiff(
  genes,
  genes_loaded
)


if (
  length(genes_needed_from_legacy) > 0L &&
    length(legacy_batch_files) > 0L
) {


  for (ff in legacy_batch_files) {


    env_tmp <- new.env(
      parent = emptyenv()
    )


    load(
      ff,
      envir = env_tmp
    )


    if (!exists(
      "dt_bb",
      envir = env_tmp,
      inherits = FALSE
    )) {
      next
    }


    dt_tmp <- data.table::as.data.table(
      get(
        "dt_bb",
        envir = env_tmp
      )
    )


    genes_here <- intersect(
      unique(
        dt_tmp$Gene
      ),
      genes_needed_from_legacy
    )


    if (length(genes_here) == 0L) {
      next
    }


    for (gene_i in genes_here) {


      dt_gene <- dt_tmp[
        Gene == gene_i
      ]


      if (
        valid_gene_result(
          dt_gene,
          gene_i
        )
      ) {


        result_list[[rr]] <- copy(
          dt_gene
        )

        rr <- rr + 1L


        genes_loaded <- c(
          genes_loaded,
          gene_i
        )
      }
    }


    rm(
      env_tmp,
      dt_tmp
    )

    gc()
  }
}


genes_loaded <- sort(
  unique(
    genes_loaded
  )
)


missing_final_genes <- setdiff(
  genes,
  genes_loaded
)


if (length(missing_final_genes) > 0L) {

  stop(
    paste0(
      "\nBenchmark incomplete.\n",
      "Missing ",
      length(missing_final_genes),
      " gene(s): ",
      paste(
        head(
          missing_final_genes,
          50
        ),
        collapse = ", "
      ),
      if (
        length(missing_final_genes) > 50
      ) {
        " ..."
      } else {
        ""
      },
      "\nRestart this script. Existing per-gene checkpoints will be reused."
    )
  )
}


dt_main <- data.table::rbindlist(
  result_list,
  use.names = TRUE,
  fill = TRUE
)


rm(
  result_list
)

gc()


# =============================================================================
# 18. Final integrity checks
# =============================================================================

gene_counts_final <- dt_main[
  ,
  .N,
  by = Gene
]


bad_gene_counts <- gene_counts_final[
  N != EXPECTED_TESTS_PER_GENE
]


if (nrow(bad_gene_counts) > 0L) {

  print(
    bad_gene_counts
  )

  stop(
    "Final benchmark contains genes with an unexpected number of tests."
  )
}


scenario_key_final <- paste(
  dt_main$Gene,
  dt_main$Platform,
  dt_main$Perturbation,
  dt_main$Exprs_noise,
  dt_main$N_tsamples,
  dt_main$N_replicates,
  dt_main$Tsteps,
  sep = "|"
)


if (
  data.table::uniqueN(
    scenario_key_final
  ) != nrow(dt_main)
) {

  stop(
    "Duplicate gene/scenario combinations found in final benchmark."
  )
}


EXPECTED_TOTAL_ROWS <- (
  length(genes) *
    EXPECTED_TESTS_PER_GENE
)


if (nrow(dt_main) != EXPECTED_TOTAL_ROWS) {

  stop(
    sprintf(
      "Final row count mismatch: found %d, expected %d.",
      nrow(dt_main),
      EXPECTED_TOTAL_ROWS
    )
  )
}


# Timing/provenance must remain invariant within gene after the benchmark.
timing_final <- dt_main[
  ,
  .(
    n_onset_shift = uniqueN(Onset_shift),
    n_onset_time = uniqueN(Onset_time),
    n_version = uniqueN(Benchmark_version)
  ),
  by = Gene
]

if (
  any(
    timing_final$n_onset_shift != 1L |
      timing_final$n_onset_time != 1L |
      timing_final$n_version != 1L
  )
) {
  stop(
    "Final benchmark timing/provenance integrity check failed."
  )
}

if (
  any(
    dt_main$Benchmark_version != BENCHMARK_VERSION
  )
) {
  stop(
    "Final benchmark contains results from a different benchmark version."
  )
}

cat(
  "\nFinal integrity checks passed.\n"
)

cat(
  "Genes:",
  data.table::uniqueN(dt_main$Gene),
  "\n"
)

cat(
  "Rows:",
  nrow(dt_main),
  "\n"
)


# =============================================================================
# 19. Valid-test indicator
# =============================================================================

dt_main[
  ,
  valid_test :=
    status == "ok" &
    is.finite(
      p.value
    )
]


# =============================================================================
# 20. Derived quantities
# =============================================================================

dt_main[
  ,
  q.value :=
    p.adjust(
      p.value,
      method = "BH"
    ),
  by = .(
    Platform,
    Perturbation,
    Exprs_noise,
    N_tsamples,
    N_replicates,
    Tsteps
  )
]


dt_main[
  ,
  DeltaRSS :=
    pmax(
      RSS0 - RSS1,
      0
    )
]


dt_main[
  ,
  logRSSratio :=
    fifelse(
      RSS1 > 0,
      log(
        RSS0 / RSS1
      ),
      NA_real_
    )
]


# -----------------------------------------------------------------------------
# Parameter errors
# -----------------------------------------------------------------------------

dt_main[
  ,
  sigma_c_error :=
    sigma_c_hat -
    sigma_true
]


dt_main[
  ,
  tau_error :=
    tau_hat -
    tau_true
]


dt_main[
  ,
  tau_s_error :=
    tau_s_hat -
    tau_s_true
]


dt_main[
  ,
  sigma_n_error :=
    sigma_n_hat -
    sigma_n_true
]


dt_main[
  ,
  alpha_error :=
    alpha_hat -
    alpha_true
]


dt_main[
  ,
  alpha_s_error :=
    alpha_s_hat -
    alpha_s_true
]


eps_param <- 1e-12


dt_main[
  ,
  sigma_c_rel_error :=
    fifelse(
      abs(sigma_true) > eps_param,
      (
        sigma_c_hat -
          sigma_true
      ) /
        sigma_true,
      NA_real_
    )
]


# =============================================================================
# 21. Summary helper
# =============================================================================

safe_median_finite <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  median(
    x
  )
}


safe_quantile_finite <- function(
  x,
  prob
) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  as.numeric(
    quantile(
      x,
      prob,
      na.rm = TRUE,
      names = FALSE
    )
  )
}


wilson_interval <- function(
  successes,
  n,
  conf = 0.95
) {

  if (
    n <= 0L ||
      successes < 0L ||
      successes > n
  ) {
    return(
      c(
        low = NA_real_,
        high = NA_real_
      )
    )
  }

  z <- qnorm(
    1 -
      (1 - conf) / 2
  )

  phat <- successes / n

  denom <- 1 + z^2 / n

  centre <- (
    phat +
      z^2 / (2 * n)
  ) / denom

  half <- (
    z *
      sqrt(
        phat *
          (1 - phat) /
          n +
          z^2 /
            (4 * n^2)
      )
  ) / denom

  c(
    low =
      max(
        0,
        centre - half
      ),
    high =
      min(
        1,
        centre + half
      )
  )
}


# =============================================================================
# 22. Summary for tables and figures
# =============================================================================

benchmark_summary <- dt_main[
  ,
  {

    valid <- valid_test


    p0 <- p.value[
      valid &
        Positive == 0
    ]


    p1 <- p.value[
      valid &
        Positive == 1
    ]


    n_null_valid <- length(
      p0
    )

    typeI_successes_005 <- sum(
      p0 <= 0.05,
      na.rm = TRUE
    )

    typeI_ci_005 <- wilson_interval(
      successes = typeI_successes_005,
      n = n_null_valid
    )


    cond_valid <- condition_number[
      valid
    ]


    list(

      N_total =
        .N,

      N_valid =
        sum(
          valid
        ),

      Valid_fraction =
        mean(
          valid
        ),

      N_null =
        sum(
          valid &
            Positive == 0
        ),

      N_alt =
        sum(
          valid &
            Positive == 1
        ),


      # ---------------------------------------------------------------
      # Type-I error
      # ---------------------------------------------------------------

      TypeI_001 =
        mean(
          p0 <= 0.01,
          na.rm = TRUE
        ),

      TypeI_005 =
        mean(
          p0 <= 0.05,
          na.rm = TRUE
        ),

      TypeI_010 =
        mean(
          p0 <= 0.10,
          na.rm = TRUE
        ),

      Inflation_005 =
        mean(
          p0 <= 0.05,
          na.rm = TRUE
        ) /
        0.05,

      TypeI_005_Wilson_low =
        unname(
          typeI_ci_005["low"]
        ),

      TypeI_005_Wilson_high =
        unname(
          typeI_ci_005["high"]
        ),

      TypeI_005_CI_contains_005 =
        is.finite(
          typeI_ci_005["low"]
        ) &&
        typeI_ci_005["low"] <= 0.05 &&
        typeI_ci_005["high"] >= 0.05,


      # ---------------------------------------------------------------
      # Power
      # ---------------------------------------------------------------

      Power_001 =
        mean(
          p1 <= 0.01,
          na.rm = TRUE
        ),

      Power_005 =
        mean(
          p1 <= 0.05,
          na.rm = TRUE
        ),

      Power_010 =
        mean(
          p1 <= 0.10,
          na.rm = TRUE
        ),


      # ---------------------------------------------------------------
      # BH discovery rates
      # ---------------------------------------------------------------

      Null_q005 =
        mean(
          q.value[
            valid &
              Positive == 0
          ] <= 0.05,
          na.rm = TRUE
        ),

      Alt_q005 =
        mean(
          q.value[
            valid &
              Positive == 1
          ] <= 0.05,
          na.rm = TRUE
        ),


      # ---------------------------------------------------------------
      # Parameter recovery
      # ---------------------------------------------------------------

      sigma_c_bias_alt =
        mean(
          sigma_c_error[
            valid &
              Positive == 1
          ],
          na.rm = TRUE
        ),

      sigma_c_MAE_alt =
        mean(
          abs(
            sigma_c_error[
              valid &
                Positive == 1
            ]
          ),
          na.rm = TRUE
        ),

      sigma_c_RMSE_alt =
        sqrt(
          mean(
            sigma_c_error[
              valid &
                Positive == 1
            ]^2,
            na.rm = TRUE
          )
        ),

      tau_MAE =
        mean(
          abs(
            tau_error[
              valid
            ]
          ),
          na.rm = TRUE
        ),

      alpha_MAE =
        mean(
          abs(
            alpha_error[
              valid
            ]
          ),
          na.rm = TRUE
        ),


      # ---------------------------------------------------------------
      # Practical identifiability
      # ---------------------------------------------------------------

      Finite_condition_fraction =
        mean(
          is.finite(
            cond_valid
          )
        ),

      Median_condition =
        safe_median_finite(
          cond_valid
        ),

      Condition_q95 =
        safe_quantile_finite(
          cond_valid,
          0.95
        ),

      Median_atom_zero_null =
        median(
          atom_zero[
            valid &
              Positive == 0
          ],
          na.rm = TRUE
        ),

      Mean_boot_failure =
        mean(
          bootstrap_failure_rate[
            valid
          ],
          na.rm = TRUE
        ),

      Mean_boot_rank_deficient =
        mean(
          boot_rank_deficient_fraction[
            valid
          ],
          na.rm = TRUE
        ),

      Median_test_seconds =
        median(
          test_seconds[
            valid
          ],
          na.rm = TRUE
        )
    )
  },

  by = .(
    Platform,
    Perturbation,
    Exprs_noise,
    N_tsamples,
    N_replicates,
    Tsteps
  )
]


# =============================================================================
# 23. Save final benchmark
# =============================================================================

cat(
  "\nSaving complete raw benchmark...\n"
)


save(
  dt_main,
  file =
    RAW_RDATA_FILE
)


data.table::fwrite(
  dt_main,
  file =
    RAW_TSV_FILE,
  sep = "\t"
)


cat(
  "Saving benchmark summary...\n"
)


save(
  benchmark_summary,
  file =
    SUMMARY_RDATA_FILE
)


data.table::fwrite(
  benchmark_summary,
  file =
    SUMMARY_TSV_FILE,
  sep = "\t"
)


# =============================================================================
# 23b. Run metadata
# =============================================================================

run_metadata <- data.table(
  benchmark_version =
    BENCHMARK_VERSION,
  input_file =
    INPUT_FILE,
  N_BOOT =
    N_BOOT,
  N_WORKERS =
    N_WORKERS,
  N_genes =
    uniqueN(
      dt_main$Gene
    ),
  N_tests =
    nrow(
      dt_main
    ),
  onset_shift_min =
    min(
      dt_main$Onset_shift
    ),
  onset_shift_max =
    max(
      dt_main$Onset_shift
    ),
  common_T_star =
    paste(
      sort(
        unique(
          ode_states_main$T_star
        )
      ),
      collapse = ","
    ),
  completed_at =
    as.character(
      Sys.time()
    )
)

fwrite(
  run_metadata,
  RUN_METADATA_FILE,
  sep = "\t"
)


# =============================================================================
# 24. Session information
# =============================================================================

sink(
  SESSION_FILE
)


cat(
  "Benchmark completed:\n"
)

print(
  Sys.time()
)


cat(
  "\nN_BOOT:\n"
)

print(
  N_BOOT
)


cat(
  "\nN_WORKERS:\n"
)

print(
  N_WORKERS
)


cat(
  "\nExpected tests per gene:\n"
)

print(
  EXPECTED_TESTS_PER_GENE
)


cat(
  "\nTotal benchmark rows:\n"
)

print(
  nrow(dt_main)
)


cat(
  "\nSession information:\n"
)

print(
  sessionInfo()
)


sink()


# =============================================================================
# 25. Completion summary
# =============================================================================

overall_elapsed_hours <- as.numeric(
  difftime(
    Sys.time(),
    overall_start,
    units = "hours"
  )
)


cat(
  "\n============================================================\n"
)

cat(
  "MAIN BENCHMARK COMPLETE\n"
)

cat(
  "============================================================\n"
)

cat(
  "Genes:",
  data.table::uniqueN(
    dt_main$Gene
  ),
  "\n"
)

cat(
  "Tests:",
  nrow(dt_main),
  "\n"
)

cat(
  "Valid tests:",
  sum(
    dt_main$valid_test
  ),
  "\n"
)

cat(
  "Total elapsed:",
  round(
    overall_elapsed_hours,
    2
  ),
  "hours\n"
)

cat(
  "\nSaved:\n"
)

cat(
  "  ", RAW_RDATA_FILE, "\n"
)

cat(
  "  ", RAW_TSV_FILE, "\n"
)

cat(
  "  ", SUMMARY_RDATA_FILE, "\n"
)

cat(
  "  ", SUMMARY_TSV_FILE, "\n"
)

cat(
  "  ", SESSION_FILE, "\n"
)

cat(
  "  ",
  CHECKPOINT_DIR,
  "/gene_*.rds\n",
  sep = ""
)

cat(
  "============================================================\n"
)