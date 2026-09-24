# =============================================================================
# Title:
# Targeted Pseudo-Shutoff Misspecification Benchmark
#
# Purpose:
#   Quantify robustness of the constrained nested-model test when inference
#   assumes COMPLETE transcriptional shutoff but the generating process retains
#   residual transcription after the common intervention time:
#
#       R_post = rho * R_pre
#
#   with
#
#       rho = 0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 1.00.
#
#   rho = 0 is the correctly specified complete-SHUTOFF reference.
#   rho = 1 is the limiting case of no effective transcriptional suppression,
#   sampled on the intervention-centered time grid.
#
# Scientific question:
#   How rapidly do null calibration, power, parameter recovery, and boundary
#   behavior deteriorate when the complete-shutoff assumption is violated?
#
# Design:
#   This is a TARGETED robustness analysis, not a repetition of the complete
#   1,152-cell factorial benchmark.
#
#   Two representative SHUTOFF designs are used:
#
#     Sparse:
#       5 sampled time points x 3 biological replicates x T_step = 10
#
#     Replicated:
#       5 sampled time points x 10 biological replicates x T_step = 10
#
#   These designs keep the temporal structure fixed while increasing
#   replication, allowing residual-transcription sensitivity to be separated
#   from a major change in temporal design.
#
#   All three observation platforms and all four noise regimes are retained.
#
#   The same synthetic genes and the same deterministic measurement-noise
#   random streams are reused across rho values, producing a paired sensitivity
#   analysis.
#
# Input:
#   ode_states_2k_20p_corrected_onset.rdata
#
# Required object:
#   ode_states
#
# Primary outputs:
#   pseudoshutoff_benchmark/
#       pseudoshutoff_raw.tsv
#       pseudoshutoff_summary.tsv
#       pseudoshutoff_design_summary.tsv
#       pseudoshutoff_progress.tsv
#       Fig_pseudoshutoff_typeI.pdf
#       Fig_pseudoshutoff_power.pdf
#       Fig_pseudoshutoff_sigma_bias.pdf
#       pseudoshutoff_benchmark.rdata
#       sessionInfo.txt
#
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
library(ggplot2)

data.table::setDTthreads(1)

source("../commons/nested_test2.r")
source("../commons/platforms.r")


# =============================================================================
# 1. User settings
# =============================================================================

INPUT_FILE <- "ode_states_2k_20p_corrected_onset.rdata"

OUTPUT_DIR <- "pseudoshutoff_benchmark"

CHECKPOINT_DIR <- file.path(
  OUTPUT_DIR,
  "checkpoints"
)

PROGRESS_FILE <- file.path(
  OUTPUT_DIR,
  "pseudoshutoff_progress.tsv"
)

RAW_FILE <- file.path(
  OUTPUT_DIR,
  "pseudoshutoff_raw.tsv"
)

SUMMARY_FILE <- file.path(
  OUTPUT_DIR,
  "pseudoshutoff_summary.tsv"
)

DESIGN_SUMMARY_FILE <- file.path(
  OUTPUT_DIR,
  "pseudoshutoff_design_summary.tsv"
)

RDATA_FILE <- file.path(
  OUTPUT_DIR,
  "pseudoshutoff_benchmark.rdata"
)

SESSION_FILE <- file.path(
  OUTPUT_DIR,
  "sessionInfo.txt"
)


# -----------------------------------------------------------------------------
# Inferential settings: identical to the principal factorial benchmark.
# -----------------------------------------------------------------------------

N_BOOT <- 1999L

SCALING_A <- TRUE

LAMBDA_TIME <- 0.5
LAMBDA_DIAG <- 0.1
REL_FLOOR <- 1e-8

TRUNCATE_NONNEGATIVE_BOOT <- FALSE
MAX_FAILURE_RATE <- 0.05


# -----------------------------------------------------------------------------
# HPC settings.
# -----------------------------------------------------------------------------

N_WORKERS <- 100L

TASK_BATCH_SIZE <- 100L

PSOCK_TIMEOUT_SECONDS <- 12L * 60L * 60L


# -----------------------------------------------------------------------------
# Number of genes.
#
# The primary misspecification endpoint is Type-I error. Five hundred null
# genes give an approximate binomial SE of 0.01 at p = 0.05.
#
# Alternative genes are retained for the secondary power/recovery analysis.
# -----------------------------------------------------------------------------

N_NULL <- 500L
N_ALT <- 250L


# -----------------------------------------------------------------------------
# Residual-transcription fractions.
# -----------------------------------------------------------------------------

RHO_LEVELS <- c(
  0.00,
  0.05,
  0.10,
  0.25,
  0.50,
  0.75,
  0.90,
  1.00
)


# -----------------------------------------------------------------------------
# Representative designs.
# -----------------------------------------------------------------------------

DESIGNS <- data.table(
  Design = c(
    "5tp_3rep_dt10",
    "5tp_10rep_dt10"
  ),
  N_time_samples = c(
    5L,
    5L
  ),
  N_replicates = c(
    3L,
    10L
  ),
  T_step = c(
    10L,
    10L
  )
)


# -----------------------------------------------------------------------------
# Platforms and noise.
# -----------------------------------------------------------------------------

RANGE_PLATFORM <- c(
  "RT-qPCR",
  "GAUSS",
  "RNA-seq"
)

RANGE_NOISE <- c(
  "Very low",
  "Low",
  "Medium",
  "High"
)

RANGE_GAUSS_NOISE <- c(
  "Very low" = 0.02,
  "Low"      = 0.05,
  "Medium"   = 0.10,
  "High"     = 0.20
)

BENCHMARK_VERSION <- "pseudoshutoff_corrected_onset_v1"


# =============================================================================
# 2. Prepare directories
# =============================================================================

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE
  )
}

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
# 3. Load corrected synthetic dataset
# =============================================================================

cat(
  "\nLoading corrected synthetic dataset...\n"
)

load(
  INPUT_FILE
)

if (
  !exists("ode_states") ||
    !data.table::is.data.table(ode_states)
) {
  stop(
    "Object 'ode_states' was not found as a data.table."
  )
}


# =============================================================================
# 4. Input validation
# =============================================================================

required_columns <- c(
  "Gene",
  "Perturbation",
  "truth_pos",
  "Base_sigma_c",
  "N_time_samples",
  "N_replicates",
  "T_step",
  "T_star",
  "Post_R_fraction",
  "Residual_transcription_pct",
  "Onset_shift",
  "Onset_time",
  "time",
  "replicate",
  "N",
  "N_s",
  "C",
  "C_s"
)

missing_columns <- setdiff(
  required_columns,
  names(
    ode_states
  )
)

if (
  length(
    missing_columns
  ) > 0L
) {
  stop(
    paste0(
      "Input dataset is missing required columns:\n  ",
      paste(
        missing_columns,
        collapse = "\n  "
      )
    )
  )
}


# -----------------------------------------------------------------------------
# Check pseudo-shutoff levels.
# -----------------------------------------------------------------------------

observed_rho <- sort(
  unique(
    ode_states[
      Perturbation %in%
        c(
          "SHUTOFF",
          "PSEUDO_SHUTOFF"
        ),
      Post_R_fraction
    ]
  )
)

missing_rho <- setdiff(
  RHO_LEVELS,
  observed_rho
)

if (
  length(
    missing_rho
  ) > 0L
) {
  stop(
    paste0(
      "The synthetic dataset does not contain all requested rho values. ",
      "Missing: ",
      paste(
        missing_rho,
        collapse = ", "
      )
    )
  )
}


# -----------------------------------------------------------------------------
# Check corrected timing semantics.
# -----------------------------------------------------------------------------

timing_check <- ode_states[
  ,
  .(
    n_onset =
      uniqueN(
        Onset_time
      ),
    n_shift =
      uniqueN(
        Onset_shift
      ),
    n_tstar =
      uniqueN(
        T_star
      )
  ),
  by =
    Gene
]

if (
  any(
    timing_check$n_onset != 1L |
      timing_check$n_shift != 1L |
      timing_check$n_tstar != 1L
  )
) {
  stop(
    paste0(
      "Timing integrity check failed: each gene must retain one onset, ",
      "one onset shift, and one common T_star."
    )
  )
}


# -----------------------------------------------------------------------------
# Validate representative designs.
# -----------------------------------------------------------------------------

for (
  ii in seq_len(
    nrow(
      DESIGNS
    )
  )
) {

  dd <- DESIGNS[ii]

  n_rows_i <- ode_states[
    N_time_samples ==
      dd$N_time_samples &
      N_replicates ==
        dd$N_replicates &
      T_step ==
        dd$T_step &
      Perturbation %in%
        c(
          "SHUTOFF",
          "PSEUDO_SHUTOFF"
        ),
    .N
  ]

  if (
    n_rows_i == 0L
  ) {
    stop(
      paste(
        "Requested design is absent:",
        dd$Design
      )
    )
  }
}


cat(
  "\nInput preflight checks PASSED.\n"
)


# =============================================================================
# 5. Restrict to intervention conditions and selected designs
# =============================================================================

ode_pseudo <- ode_states[
  Perturbation %in%
    c(
      "SHUTOFF",
      "PSEUDO_SHUTOFF"
    ) &
    Post_R_fraction %in%
      RHO_LEVELS
]

design_keys <- DESIGNS[
  ,
  .(
    N_time_samples,
    N_replicates,
    T_step
  )
]

ode_pseudo <- merge(
  ode_pseudo,
  design_keys,
  by = c(
    "N_time_samples",
    "N_replicates",
    "T_step"
  ),
  all = FALSE
)

setkey(
  ode_pseudo,
  Gene,
  N_time_samples,
  N_replicates,
  T_step,
  Post_R_fraction
)


# =============================================================================
# 6. Select paired genes
# =============================================================================

gene_truth <- unique(
  ode_pseudo[
    ,
    .(
      Gene,
      truth_pos
    )
  ]
)

null_genes_all <- sort(
  gene_truth[
    truth_pos == 0,
    Gene
  ]
)

alt_genes_all <- sort(
  gene_truth[
    truth_pos == 1,
    Gene
  ]
)

if (
  length(
    null_genes_all
  ) <
    N_NULL
) {
  stop(
    paste(
      "Requested",
      N_NULL,
      "null genes but only",
      length(
        null_genes_all
      ),
      "are available."
    )
  )
}

if (
  length(
    alt_genes_all
  ) <
    N_ALT
) {
  stop(
    paste(
      "Requested",
      N_ALT,
      "alternative genes but only",
      length(
        alt_genes_all
      ),
      "are available."
    )
  )
}


# Deterministic random selection.
set.seed(
  20260922L
)

NULL_GENES <- sort(
  sample(
    null_genes_all,
    N_NULL,
    replace = FALSE
  )
)

ALT_GENES <- sort(
  sample(
    alt_genes_all,
    N_ALT,
    replace = FALSE
  )
)

GENES_KEEP <- sort(
  c(
    NULL_GENES,
    ALT_GENES
  )
)

ode_pseudo <- ode_pseudo[
  Gene %in%
    GENES_KEEP
]


cat(
  "\nSelected genes:\n",
  "  Null:        ",
  length(
    NULL_GENES
  ),
  "\n",
  "  Alternative: ",
  length(
    ALT_GENES
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 7. Platform simulation helper
#
# Identical observation-model settings to the principal factorial benchmark.
# =============================================================================

add_platform_noise <- function(
  dt,
  platform,
  noise
) {

  dt <- copy(
    as.data.table(
      dt
    )
  )

  targets <- c(
    "N",
    "C",
    "C_s",
    "N_s"
  )


  if (
    platform ==
      "GAUSS"
  ) {

    noise_sd <- unname(
      RANGE_GAUSS_NOISE[
        noise
      ]
    )

    return(
      add_gaussian_noise(
        dt,
        cols =
          targets,
        noise_sd =
          noise_sd
      )
    )
  }


  if (
    platform ==
      "RT-qPCR"
  ) {

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
        targets =
          targets,
        ct_sd =
          ct_sd,
        scale_copies =
          scale_copies
      )
    )
  }


  if (
    platform ==
      "RNA-seq"
  ) {

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
    ) *
      0.25

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
        targets =
          targets,
        scale_counts =
          scale_counts,
        mean_disp =
          mean_disp,
        cv_disp =
          cv_disp
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
# 8. Reproducible seed helpers
#
# rho is deliberately EXCLUDED from the measurement/bootstrap seed.
#
# Therefore the same gene/design/platform/noise combination uses paired random
# streams across residual-transcription levels. This reduces Monte Carlo noise
# when assessing the effect of rho.
# =============================================================================

stable_seed_from_string <- function(
  x
) {

  ints <- utf8ToInt(
    enc2utf8(
      x
    )
  )

  h <- 104729

  for (ii in ints) {

    h <- (
      h *
        1000003 +
        ii
    ) %%
      2000000000
  }

  out <- as.integer(
    h
  )

  if (
    !is.finite(
      out
    ) ||
      out <= 0L
  ) {
    out <- 1L
  }

  out
}


make_seed <- function(
  gene,
  design,
  platform,
  noise,
  stream
) {

  stable_seed_from_string(
    paste(
      BENCHMARK_VERSION,
      gene,
      design,
      platform,
      noise,
      stream,
      sep = "|"
    )
  )
}


# =============================================================================
# 9. Safe extraction helpers
# =============================================================================

safe_value <- function(
  x,
  field,
  default =
    NA_real_
) {

  if (
    is.null(
      x
    ) ||
      is.null(
        x[[field]]
      ) ||
      length(
        x[[field]]
      ) == 0L
  ) {
    return(
      default
    )
  }

  suppressWarnings(
    as.numeric(
      x[[field]][1]
    )
  )
}


safe_text <- function(
  x,
  field,
  default =
    NA_character_
) {

  if (
    is.null(
      x
    ) ||
      is.null(
        x[[field]]
      ) ||
      length(
        x[[field]]
      ) == 0L
  ) {
    return(
      default
    )
  }

  as.character(
    x[[field]][1]
  )
}


safe_coef <- function(
  x,
  name
) {

  if (
    is.null(
      x
    ) ||
      is.null(
        names(
          x
        )
      ) ||
      !name %in%
        names(
          x
        )
  ) {
    return(
      NA_real_
    )
  }

  unname(
    x[
      name
    ]
  )
}


# =============================================================================
# 10. Run one gene
# =============================================================================

run_gene_pseudoshutoff <- function(
  gene_i,
  N_boot
) {

  out <- list()

  kk <- 1L


  for (
    di in seq_len(
      nrow(
        DESIGNS
      )
    )
  ) {

    design_i <- DESIGNS[di]


    for (
      rho_i in RHO_LEVELS
    )
    {

      dt_tc <- ode_pseudo[
        Gene ==
          gene_i &
          N_time_samples ==
            design_i$N_time_samples &
          N_replicates ==
            design_i$N_replicates &
          T_step ==
            design_i$T_step &
          abs(
            Post_R_fraction -
              rho_i
          ) <
            1e-12
      ]


      if (
        nrow(
          dt_tc
        ) == 0L
      ) {
        stop(
          paste(
            "Missing latent trajectory for gene",
            gene_i,
            "design",
            design_i$Design,
            "rho",
            rho_i
          )
        )
      }


      truth_values <- unique(
        dt_tc$truth_pos
      )

      if (
        length(
          truth_values
        ) != 1L
      ) {
        stop(
          "truth_pos is not unique within gene/design/rho."
        )
      }

      truth_i <- truth_values[1]


      sigma_values <- unique(
        dt_tc$Base_sigma_c
      )

      if (
        length(
          sigma_values
        ) != 1L
      ) {
        stop(
          "Base_sigma_c is not unique within gene/design/rho."
        )
      }

      sigma_true <- sigma_values[1]


      tstar_values <- unique(
        dt_tc$T_star
      )

      tstar_values <- tstar_values[
        is.finite(
          tstar_values
        )
      ]

      if (
        length(
          tstar_values
        ) != 1L
      ) {
        stop(
          "T_star is not unique within gene/design/rho."
        )
      }

      # CRITICAL MISSPECIFICATION:
      # inference always assumes a complete shutoff at this common T_star,
      # even when rho > 0 in the generating process.
      t_star_fit <- tstar_values[1]


      for (
        platform_i in RANGE_PLATFORM
      ) {

        for (
          noise_i in RANGE_NOISE
        ) {


          measurement_seed <- make_seed(
            gene =
              gene_i,
            design =
              design_i$Design,
            platform =
              platform_i,
            noise =
              noise_i,
            stream =
              "measurement"
          )


          bootstrap_seed <- make_seed(
            gene =
              gene_i,
            design =
              design_i$Design,
            platform =
              platform_i,
            noise =
              noise_i,
            stream =
              "bootstrap"
          )


          set.seed(
            measurement_seed
          )


          dt_noisy <- tryCatch(

            add_platform_noise(
              dt =
                dt_tc,
              platform =
                platform_i,
              noise =
                noise_i
            ),

            error = function(e) {

              structure(
                list(
                  error_message =
                    conditionMessage(
                      e
                    )
                ),
                class =
                  "measurement_error"
              )
            }
          )


          if (
            inherits(
              dt_noisy,
              "measurement_error"
            )
          ) {

            out[[kk]] <- data.table(

              Benchmark_version =
                BENCHMARK_VERSION,

              Gene =
                gene_i,

              Positive =
                truth_i,

              sigma_true =
                sigma_true,

              Design =
                design_i$Design,

              N_tsamples =
                design_i$N_time_samples,

              N_replicates =
                design_i$N_replicates,

              Tsteps =
                design_i$T_step,

              Platform =
                platform_i,

              Exprs_noise =
                noise_i,

              Post_R_fraction =
                rho_i,

              Residual_transcription_pct =
                100 *
                  rho_i,

              T_star =
                t_star_fit,

              N_boot =
                N_boot,

              status =
                "measurement_error",

              error_message =
                dt_noisy$error_message
            )

            kk <- kk + 1L

            next
          }


          test_start <- Sys.time()


          WW <- tryCatch(

            test_sigma_nested(

              tsampled_data =
                dt_noisy,

              scaling_A =
                SCALING_A,

              # The fitted model assumes complete shutoff at T_star.
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
                p.value =
                  NA_real_,

                status =
                  "test_error",

                error.message =
                  conditionMessage(
                    e
                  )
              )
            }
          )


          elapsed_seconds <- as.numeric(
            difftime(
              Sys.time(),
              test_start,
              units =
                "secs"
            )
          )


          coef_full_i <- WW$coef_full
          coef_null_i <- WW$coef_null


          out[[kk]] <- data.table(

            Benchmark_version =
              BENCHMARK_VERSION,

            Gene =
              gene_i,

            Positive =
              truth_i,

            sigma_true =
              sigma_true,

            Design =
              design_i$Design,

            N_tsamples =
              design_i$N_time_samples,

            N_replicates =
              design_i$N_replicates,

            Tsteps =
              design_i$T_step,

            Platform =
              platform_i,

            Exprs_noise =
              noise_i,

            Perturbation =
              ifelse(
                rho_i ==
                  0,
                "SHUTOFF",
                "PSEUDO_SHUTOFF"
              ),

            Post_R_fraction =
              rho_i,

            Residual_transcription_pct =
              100 *
                rho_i,

            T_star =
              t_star_fit,

            N_boot =
              N_boot,

            Onset_shift =
              unique(
                dt_tc$Onset_shift
              )[1],

            Onset_time =
              unique(
                dt_tc$Onset_time
              )[1],

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

            IR =
              safe_value(
                WW,
                "IR"
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

            atom_zero =
              safe_value(
                WW,
                "atom.zero"
              ),

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

            boot_rank_deficient_fraction =
              safe_value(
                WW,
                "bootstrap.rank.deficient.fraction"
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

            status =
              safe_text(
                WW,
                "status",
                default =
                  "ok"
              ),

            error_message =
              safe_text(
                WW,
                "error.message"
              ),

            runtime_seconds =
              elapsed_seconds
          )


          kk <- kk + 1L
        }
      }
    }
  }


  rbindlist(
    out,
    use.names =
      TRUE,
    fill =
      TRUE
  )
}


# =============================================================================
# 11. Checkpoint helpers
# =============================================================================

EXPECTED_TESTS_PER_GENE <- (
  nrow(
    DESIGNS
  ) *
    length(
      RHO_LEVELS
    ) *
    length(
      RANGE_PLATFORM
    ) *
    length(
      RANGE_NOISE
    )
)


gene_checkpoint_file <- function(
  gene_i
) {

  file.path(
    CHECKPOINT_DIR,
    sprintf(
      "gene_%06d.rds",
      gene_i
    )
  )
}


valid_gene_result <- function(
  x,
  gene_i
) {

  if (
    is.null(
      x
    ) ||
      !is.data.table(
        x
      )
  ) {
    return(
      FALSE
    )
  }


  required <- c(
    "Benchmark_version",
    "Gene",
    "Design",
    "Platform",
    "Exprs_noise",
    "Post_R_fraction",
    "Positive",
    "p.value",
    "status"
  )

  if (
    any(
      !required %in%
        names(
          x
        )
    )
  ) {
    return(
      FALSE
    )
  }


  if (
    nrow(
      x
    ) !=
      EXPECTED_TESTS_PER_GENE
  ) {
    return(
      FALSE
    )
  }


  if (
    uniqueN(
      x$Benchmark_version
    ) !=
      1L ||
      x$Benchmark_version[1] !=
        BENCHMARK_VERSION
  ) {
    return(
      FALSE
    )
  }


  if (
    uniqueN(
      x$Gene
    ) !=
      1L ||
      x$Gene[1] !=
        gene_i
  ) {
    return(
      FALSE
    )
  }


  scenario_key <- paste(
    x$Design,
    x$Platform,
    x$Exprs_noise,
    sprintf(
      "%.2f",
      x$Post_R_fraction
    ),
    sep =
      "|"
  )


  if (
    uniqueN(
      scenario_key
    ) !=
      EXPECTED_TESTS_PER_GENE
  ) {
    return(
      FALSE
    )
  }


  TRUE
}


run_gene_and_checkpoint <- function(
  gene_i,
  N_boot
) {

  checkpoint <- gene_checkpoint_file(
    gene_i
  )


  if (
    file.exists(
      checkpoint
    )
  ) {

    old <- tryCatch(
      readRDS(
        checkpoint
      ),
      error =
        function(e) NULL
    )

    if (
      valid_gene_result(
        old,
        gene_i
      )
    ) {

      return(
        data.table(
          Gene =
            gene_i,
          status =
            "already_checkpointed",
          n_tests =
            nrow(
              old
            ),
          error_message =
            NA_character_
        )
      )
    }
  }


  t0 <- Sys.time()


  result <- tryCatch(

    run_gene_pseudoshutoff(
      gene_i =
        gene_i,
      N_boot =
        N_boot
    ),

    error = function(e) {

      structure(
        list(
          error_message =
            conditionMessage(
              e
            )
        ),
        class =
          "gene_error"
      )
    }
  )


  if (
    inherits(
      result,
      "gene_error"
    )
  ) {

    return(
      data.table(
        Gene =
          gene_i,
        status =
          "gene_error",
        n_tests =
          0L,
        elapsed_seconds =
          as.numeric(
            difftime(
              Sys.time(),
              t0,
              units =
                "secs"
            )
          ),
        error_message =
          result$error_message
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
      data.table(
        Gene =
          gene_i,
        status =
          "invalid_result",
        n_tests =
          nrow(
            result
          ),
        elapsed_seconds =
          as.numeric(
            difftime(
              Sys.time(),
              t0,
              units =
                "secs"
            )
          ),
        error_message =
          "Gene result failed integrity validation."
      )
    )
  }


  tmp <- paste0(
    checkpoint,
    ".tmp.",
    Sys.getpid()
  )


  saveRDS(
    result,
    tmp,
    compress =
      "gzip"
  )


  if (
    !file.rename(
      tmp,
      checkpoint
    )
  ) {

    unlink(
      tmp
    )

    return(
      data.table(
        Gene =
          gene_i,
        status =
          "checkpoint_write_failed",
        n_tests =
          nrow(
            result
          ),
        error_message =
          "Atomic checkpoint rename failed."
      )
    )
  }


  data.table(
    Gene =
      gene_i,
    status =
      "completed",
    n_tests =
      nrow(
        result
      ),
    elapsed_seconds =
      as.numeric(
        difftime(
          Sys.time(),
          t0,
          units =
            "secs"
        )
      ),
    error_message =
      NA_character_
  )
}


# =============================================================================
# 12. Fast serial smoke test
# =============================================================================

SMOKE_BOOT <- 19L

smoke_gene <- GENES_KEEP[1]


cat(
  "\n============================================================\n",
  "SERIAL PSEUDO-SHUTOFF SMOKE TEST\n",
  "============================================================\n",
  "Gene: ",
  smoke_gene,
  "\n",
  "Expected scenarios/gene: ",
  EXPECTED_TESTS_PER_GENE,
  "\n",
  "Bootstrap/test: ",
  SMOKE_BOOT,
  "\n",
  sep = ""
)


smoke <- run_gene_pseudoshutoff(
  gene_i =
    smoke_gene,
  N_boot =
    SMOKE_BOOT
)


if (
  nrow(
    smoke
  ) !=
    EXPECTED_TESTS_PER_GENE
) {
  stop(
    paste(
      "Smoke test returned",
      nrow(
        smoke
      ),
      "rows; expected",
      EXPECTED_TESTS_PER_GENE
    )
  )
}


smoke_key <- paste(
  smoke$Design,
  smoke$Platform,
  smoke$Exprs_noise,
  sprintf(
    "%.2f",
    smoke$Post_R_fraction
  ),
  sep =
    "|"
)


if (
  uniqueN(
    smoke_key
  ) !=
    EXPECTED_TESTS_PER_GENE
) {
  stop(
    "Smoke test contains duplicate/missing scenarios."
  )
}


if (
  uniqueN(
    smoke$Onset_time
  ) !=
    1L
) {
  stop(
    "Smoke test found multiple onset times for one gene."
  )
}


cat(
  "\nSmoke status counts:\n"
)

print(
  smoke[
    ,
    .N,
    by =
      status
  ]
)


cat(
  "\nSMOKE TEST PASSED.\n"
)


rm(
  smoke
)

gc()


# =============================================================================
# 13. Determine tasks remaining
# =============================================================================

done <- vapply(
  GENES_KEEP,
  function(gene_i) {

    ff <- gene_checkpoint_file(
      gene_i
    )

    if (
      !file.exists(
        ff
      )
    ) {
      return(
        FALSE
      )
    }

    xx <- tryCatch(
      readRDS(
        ff
      ),
      error =
        function(e) NULL
    )

    valid_gene_result(
      xx,
      gene_i
    )
  },
  logical(1)
)


genes_remaining <- GENES_KEEP[
  !done
]


cat(
  "\nAlready complete: ",
  sum(
    done
  ),
  " / ",
  length(
    GENES_KEEP
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 14. PSOCK cluster
# =============================================================================

overall_start <- Sys.time()


if (
  length(
    genes_remaining
  ) > 0L
) {

  workers_use <- min(
    N_WORKERS,
    length(
      genes_remaining
    )
  )


  cat(
    "\nStarting PSOCK cluster with ",
    workers_use,
    " workers...\n",
    sep = ""
  )


  cl <- makePSOCKcluster(
    workers_use,
    outfile = "",
    timeout =
      PSOCK_TIMEOUT_SECONDS,
    setup_timeout =
      120
  )


  tryCatch(
    {

      clusterEvalQ(
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

          data.table::setDTthreads(1)

          source("../commons/nested_test2.r")
          source("../commons/platforms.r")

          NULL
        }
      )


      export_names <- c(
        "ode_pseudo",
        "DESIGNS",
        "RHO_LEVELS",
        "RANGE_PLATFORM",
        "RANGE_NOISE",
        "RANGE_GAUSS_NOISE",
        "BENCHMARK_VERSION",
        "SCALING_A",
        "LAMBDA_TIME",
        "LAMBDA_DIAG",
        "REL_FLOOR",
        "TRUNCATE_NONNEGATIVE_BOOT",
        "MAX_FAILURE_RATE",
        "CHECKPOINT_DIR",
        "EXPECTED_TESTS_PER_GENE",
        "add_platform_noise",
        "stable_seed_from_string",
        "make_seed",
        "safe_value",
        "safe_text",
        "safe_coef",
        "run_gene_pseudoshutoff",
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


      clusterExport(
        cl,
        export_names,
        envir =
          .GlobalEnv
      )


      batch_ids <- split(
        seq_along(
          genes_remaining
        ),
        ceiling(
          seq_along(
            genes_remaining
          ) /
            TASK_BATCH_SIZE
        )
      )


      for (
        bb in seq_along(
          batch_ids
        )
      ) {

        genes_batch <- genes_remaining[
          batch_ids[[bb]]
        ]


        cat(
          "\n============================================================\n",
          "Running batch ",
          bb,
          " / ",
          length(
            batch_ids
          ),
          " | genes = ",
          length(
            genes_batch
          ),
          "\n",
          "============================================================\n",
          sep = ""
        )


        batch_start <- Sys.time()


        stat <- parLapplyLB(
          cl,
          genes_batch,
          run_gene_and_checkpoint,
          N_boot =
            N_BOOT
        )


        stat <- rbindlist(
          stat,
          use.names =
            TRUE,
          fill =
            TRUE
        )


        stat[
          ,
          `:=`(
            task_batch =
              bb,

            batch_elapsed_minutes =
              as.numeric(
                difftime(
                  Sys.time(),
                  batch_start,
                  units =
                    "mins"
                )
              ),

            master_timestamp =
              as.character(
                Sys.time()
              )
          )
        ]


        fwrite(
          stat,
          PROGRESS_FILE,
          sep =
            "\t",
          append =
            file.exists(
              PROGRESS_FILE
            ),
          col.names =
            !file.exists(
              PROGRESS_FILE
            )
        )


        print(
          stat[
            ,
            .N,
            by =
              status
          ]
        )


        errors <- stat[
          !is.na(
            error_message
          )
        ]


        if (
          nrow(
            errors
          ) > 0L
        ) {

          cat(
            "\nErrors in current batch:\n"
          )

          print(
            unique(
              errors[
                ,
                .(
                  status,
                  error_message
                )
              ]
            )
          )
        }
      }

    },

    finally = {

      cat(
        "\nStopping PSOCK cluster...\n"
      )

      try(
        stopCluster(
          cl
        ),
        silent =
          TRUE
      )
    }
  )
}


# =============================================================================
# 15. Collect checkpoints
# =============================================================================

cat(
  "\n============================================================\n",
  "COLLECTING CHECKPOINTS\n",
  "============================================================\n",
  sep = ""
)


result_list <- list()

rr <- 1L


for (
  gene_i in GENES_KEEP
) {

  ff <- gene_checkpoint_file(
    gene_i
  )

  if (
    !file.exists(
      ff
    )
  ) {
    next
  }


  xx <- tryCatch(
    readRDS(
      ff
    ),
    error =
      function(e) NULL
  )


  if (
    valid_gene_result(
      xx,
      gene_i
    )
  ) {

    result_list[[rr]] <- xx

    rr <- rr + 1L
  }
}


if (
  length(
    result_list
  ) == 0L
) {
  stop(
    "No valid checkpoints were recovered."
  )
}


results <- rbindlist(
  result_list,
  use.names =
    TRUE,
  fill =
    TRUE
)


setorder(
  results,
  Positive,
  Gene,
  Design,
  Platform,
  Exprs_noise,
  Post_R_fraction
)


# =============================================================================
# 16. Final raw-result validation
# =============================================================================

results[
  ,
  valid :=
    is.finite(
      p.value
    ) &
    p.value >=
      0 &
    p.value <=
      1 &
    !status %in%
      c(
        "test_error",
        "measurement_error"
      )
]


expected_rows <- (
  length(
    GENES_KEEP
  ) *
    EXPECTED_TESTS_PER_GENE
)


if (
  nrow(
    results
  ) !=
    expected_rows
) {
  stop(
    paste(
      "Final row-count mismatch. Expected",
      expected_rows,
      "observed",
      nrow(
        results
      )
    )
  )
}


# =============================================================================
# 17. Wilson interval helper
# =============================================================================

wilson_interval <- function(
  successes,
  n,
  conf =
    0.95
) {

  if (
    n <= 0L
  ) {

    return(
      c(
        low =
          NA_real_,
        high =
          NA_real_
      )
    )
  }


  z <- qnorm(
    1 -
      (
        1 -
          conf
      ) /
        2
  )


  phat <- successes /
    n


  denom <- 1 +
    z^2 /
      n


  centre <- (
    phat +
      z^2 /
        (
          2 *
            n
        )
  ) /
    denom


  half <- (
    z *
      sqrt(
        phat *
          (
            1 -
              phat
          ) /
          n +
          z^2 /
            (
              4 *
                n^2
            )
      )
  ) /
    denom


  c(
    low =
      max(
        0,
        centre -
          half
      ),

    high =
      min(
        1,
        centre +
          half
      )
  )
}


# =============================================================================
# 18. Scenario summary
# =============================================================================

summary_dt <- results[
  ,
  {

    null <- .SD[
      valid &
        Positive ==
          0
    ]

    alt <- .SD[
      valid &
        Positive ==
          1
    ]


    n0 <- nrow(
      null
    )

    n1 <- nrow(
      alt
    )


    k0 <- if (
      n0 >
        0L
    ) {
      sum(
        null$p.value <=
          0.05
      )
    } else {
      0L
    }


    ci <- wilson_interval(
      successes =
        k0,
      n =
        n0
    )


    typeI <- if (
      n0 >
        0L
    ) {
      k0 /
        n0
    } else {
      NA_real_
    }


    power <- if (
      n1 >
        0L
    ) {
      mean(
        alt$p.value <=
          0.05
      )
    } else {
      NA_real_
    }


    sigma_bias_null <- if (
      n0 >
        0L
    ) {
      mean(
        null$Sigma,
        na.rm =
          TRUE
      )
    } else {
      NA_real_
    }


    sigma_bias_alt <- if (
      n1 >
        0L
    ) {
      mean(
        alt$Sigma -
          alt$sigma_true,
        na.rm =
          TRUE
      )
    } else {
      NA_real_
    }


    sigma_rmse_alt <- if (
      n1 >
        0L
    ) {
      sqrt(
        mean(
          (
            alt$Sigma -
              alt$sigma_true
          )^2,
          na.rm =
            TRUE
        )
      )
    } else {
      NA_real_
    }


    .(
      N_null =
        n0,

      N_alt =
        n1,

      Valid_fraction_null =
        mean(
          .SD[
            Positive ==
              0,
            valid
          ]
        ),

      Valid_fraction_alt =
        mean(
          .SD[
            Positive ==
              1,
            valid
          ]
        ),

      TypeI_005 =
        typeI,

      TypeI_Wilson_low =
        unname(
          ci[
            "low"
          ]
        ),

      TypeI_Wilson_high =
        unname(
          ci[
            "high"
          ]
        ),

      TypeI_CI_contains_005 =
        is.finite(
          ci[
            "low"
          ]
        ) &&
        ci[
          "low"
        ] <=
          0.05 &&
        ci[
          "high"
        ] >=
          0.05,

      Inflation_005 =
        typeI /
          0.05,

      Power_005 =
        power,

      Null_sigma_hat_mean =
        sigma_bias_null,

      Alt_sigma_bias =
        sigma_bias_alt,

      Alt_sigma_RMSE =
        sigma_rmse_alt,

      Null_boundary_fraction =
        mean(
          null$Sigma <=
            1e-12,
          na.rm =
            TRUE
        ),

      Alt_boundary_fraction =
        mean(
          alt$Sigma <=
            1e-12,
          na.rm =
            TRUE
        ),

      Median_atom_zero_null =
        median(
          null$atom_zero,
          na.rm =
            TRUE
        ),

      Finite_condition_fraction =
        mean(
          is.finite(
            .SD$condition_number
          )
        ),

      Median_condition =
        median(
          .SD$condition_number[
            is.finite(
              .SD$condition_number
            )
          ],
          na.rm =
            TRUE
        ),

      Mean_bootstrap_failure_rate =
        mean(
          .SD$bootstrap_failure_rate,
          na.rm =
            TRUE
        ),

      Median_test_seconds =
        median(
          .SD$runtime_seconds,
          na.rm =
            TRUE
        )
    )
  },
  by = .(
    Design,
    N_tsamples,
    N_replicates,
    Tsteps,
    Platform,
    Exprs_noise,
    Post_R_fraction,
    Residual_transcription_pct
  )
]


setorder(
  summary_dt,
  Design,
  Platform,
  Exprs_noise,
  Post_R_fraction
)


# =============================================================================
# 19. Design-level summary across platform x noise
# =============================================================================

design_summary <- summary_dt[
  ,
  .(
    N_scenarios =
      .N,

    TypeI_median =
      median(
        TypeI_005,
        na.rm =
          TRUE
      ),

    TypeI_q25 =
      quantile(
        TypeI_005,
        0.25,
        na.rm =
          TRUE,
        names =
          FALSE
      ),

    TypeI_q75 =
      quantile(
        TypeI_005,
        0.75,
        na.rm =
          TRUE,
        names =
          FALSE
      ),

    TypeI_max =
      max(
        TypeI_005,
        na.rm =
          TRUE
      ),

    Fraction_CI_contains_005 =
      mean(
        TypeI_CI_contains_005
      ),

    Power_median =
      median(
        Power_005,
        na.rm =
          TRUE
      ),

    Null_sigma_hat_mean_median =
      median(
        Null_sigma_hat_mean,
        na.rm =
          TRUE
      ),

    Alt_sigma_bias_median =
      median(
        Alt_sigma_bias,
        na.rm =
          TRUE
      ),

    Alt_sigma_RMSE_median =
      median(
        Alt_sigma_RMSE,
        na.rm =
          TRUE
      ),

    Null_boundary_fraction_median =
      median(
        Null_boundary_fraction,
        na.rm =
          TRUE
      ),

    Mean_bootstrap_failure_rate =
      mean(
        Mean_bootstrap_failure_rate,
        na.rm =
          TRUE
      )
  ),
  by = .(
    Design,
    N_tsamples,
    N_replicates,
    Tsteps,
    Post_R_fraction,
    Residual_transcription_pct
  )
]


setorder(
  design_summary,
  Design,
  Post_R_fraction
)


# =============================================================================
# 20. Save tables
# =============================================================================

fwrite(
  results,
  RAW_FILE,
  sep =
    "\t"
)

fwrite(
  summary_dt,
  SUMMARY_FILE,
  sep =
    "\t"
)

fwrite(
  design_summary,
  DESIGN_SUMMARY_FILE,
  sep =
    "\t"
)


# =============================================================================
# 21. Figures
# =============================================================================

theme_pseudo <- theme_classic(
  base_size =
    11
) +
  theme(
    panel.grid.major =
      element_line(
        colour =
          "grey93",
        linewidth =
          0.25
      ),
    panel.grid.minor =
      element_blank(),
    strip.background =
      element_blank(),
    strip.text =
      element_text(
        face =
          "bold"
      ),
    legend.position =
      "bottom"
  )


# -----------------------------------------------------------------------------
# Type-I sensitivity.
# -----------------------------------------------------------------------------

p_typeI <- ggplot(
  summary_dt,
  aes(
    x =
      Residual_transcription_pct,
    y =
      TypeI_005,
    group =
      Exprs_noise,
    linetype =
      Exprs_noise,
    shape =
      Exprs_noise
  )
) +

  geom_hline(
    yintercept =
      0.05,
    linetype =
      "dashed",
    linewidth =
      0.55,
    colour =
      "grey35"
  ) +

  annotate(
    "rect",
    xmin =
      -Inf,
    xmax =
      Inf,
    ymin =
      0.025,
    ymax =
      0.075,
    alpha =
      0.05
  ) +

  geom_line(
    linewidth =
      0.8
  ) +

  geom_point(
    size =
      1.8
  ) +

  facet_grid(
    Platform ~ Design
  ) +

  scale_x_continuous(
    breaks =
      100 *
        RHO_LEVELS
  ) +

  labs(
    x =
      "Residual transcription after nominal shutoff (%)",
    y =
      "Empirical Type-I error",
    linetype =
      "Noise",
    shape =
      "Noise"
  ) +

  theme_pseudo


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Fig_pseudoshutoff_typeI.pdf"
    ),
  plot =
    p_typeI,
  width =
    10.5,
  height =
    8,
  units =
    "in",
  device =
    cairo_pdf
)


# -----------------------------------------------------------------------------
# Power sensitivity.
# -----------------------------------------------------------------------------

p_power <- ggplot(
  summary_dt,
  aes(
    x =
      Residual_transcription_pct,
    y =
      Power_005,
    group =
      Exprs_noise,
    linetype =
      Exprs_noise,
    shape =
      Exprs_noise
  )
) +

  geom_line(
    linewidth =
      0.8
  ) +

  geom_point(
    size =
      1.8
  ) +

  facet_grid(
    Platform ~ Design
  ) +

  scale_x_continuous(
    breaks =
      100 *
        RHO_LEVELS
  ) +

  labs(
    x =
      "Residual transcription after nominal shutoff (%)",
    y =
      "Empirical power",
    linetype =
      "Noise",
    shape =
      "Noise"
  ) +

  theme_pseudo


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Fig_pseudoshutoff_power.pdf"
    ),
  plot =
    p_power,
  width =
    10.5,
  height =
    8,
  units =
    "in",
  device =
    cairo_pdf
)


# -----------------------------------------------------------------------------
# Null sigma estimate / false-positive conversion signal.
# -----------------------------------------------------------------------------

p_sigma <- ggplot(
  summary_dt,
  aes(
    x =
      Residual_transcription_pct,
    y =
      Null_sigma_hat_mean,
    group =
      Exprs_noise,
    linetype =
      Exprs_noise,
    shape =
      Exprs_noise
  )
) +

  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.5,
    colour =
      "grey35"
  ) +

  geom_line(
    linewidth =
      0.8
  ) +

  geom_point(
    size =
      1.8
  ) +

  facet_grid(
    Platform ~ Design,
    scales =
      "free_y"
  ) +

  scale_x_continuous(
    breaks =
      100 *
        RHO_LEVELS
  ) +

  labs(
    x =
      "Residual transcription after nominal shutoff (%)",
    y =
      expression(
        "Mean " *
          hat(
            sigma
          )[c] *
          " under " *
          H[0]
      ),
    linetype =
      "Noise",
    shape =
      "Noise"
  ) +

  theme_pseudo


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Fig_pseudoshutoff_sigma_bias.pdf"
    ),
  plot =
    p_sigma,
  width =
    10.5,
  height =
    8,
  units =
    "in",
  device =
    cairo_pdf
)


# =============================================================================
# 22. Save analysis object
# =============================================================================

pseudoshutoff_benchmark <- list(

  settings =
    list(
      benchmark_version =
        BENCHMARK_VERSION,
      input_file =
        INPUT_FILE,
      N_boot =
        N_BOOT,
      N_null =
        N_NULL,
      N_alt =
        N_ALT,
      rho_levels =
        RHO_LEVELS,
      designs =
        DESIGNS,
      platforms =
        RANGE_PLATFORM,
      noise_levels =
        RANGE_NOISE
    ),

  selected_null_genes =
    NULL_GENES,

  selected_alt_genes =
    ALT_GENES,

  raw =
    results,

  summary =
    summary_dt,

  design_summary =
    design_summary
)


save(
  pseudoshutoff_benchmark,
  file =
    RDATA_FILE
)


# =============================================================================
# 23. Session info
# =============================================================================

sink(
  SESSION_FILE
)

print(
  sessionInfo()
)

sink()


# =============================================================================
# 24. Final report
# =============================================================================

cat(
  "\n============================================================\n",
  "PSEUDO-SHUTOFF BENCHMARK COMPLETE\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nDesign-level sensitivity to residual transcription:\n"
)

print(
  design_summary
)


cat(
  "\nScenario-level summary:\n"
)

print(
  summary_dt
)


cat(
  "\nSaved in:\n  ",
  OUTPUT_DIR,
  "\n",
  sep = ""
)


cat(
  "\nInterpretation:\n",
  "  rho = 0 is the correctly specified complete-shutoff reference.\n",
  "  rho > 0 is generated with residual transcription but fitted assuming\n",
  "  complete shutoff at the common T_star.\n",
  "  Primary endpoint: empirical Type-I error under sigma_c = 0.\n",
  "  Secondary endpoints: power, sigma_c bias/RMSE, boundary behavior,\n",
  "  numerical conditioning, and bootstrap stability.\n",
  sep = ""
)
