# =============================================================================
# Title:
# Fraction-Separation Robustness Benchmark
#
# Purpose:
#   Quantify robustness of the constrained nested-model test to two forms of
#   compartment-separation misspecification:
#
#   (1) Nuclear/cytoplasmic relative scaling mismatch
#
#         C_obs   = s_C * C
#         C_s_obs = s_C * C_s
#
#       while N and N_s are unchanged.
#
#   (2) Symmetric nuclear/cytoplasmic cross-contamination
#
#         N_obs   = (1-eps) * N   + eps * C
#         N_s_obs = (1-eps) * N_s + eps * C_s
#         C_obs   = eps * N       + (1-eps) * C
#         C_s_obs = eps * N_s     + (1-eps) * C_s
#
#   In both analyses, inference is performed with the ORIGINAL model, i.e. the
#   misspecification is deliberately NOT represented in the fitted model.
#
# Scientific question:
#   Can systematic fraction scaling or cross-fraction contamination generate
#   spurious evidence for sigma_c > 0, alter power, or destabilize inference?
#
# Design:
#   Uses the corrected synthetic dataset under TRUE complete SHUTOFF only.
#
#   Representative designs:
#       5 time points x 3 biological replicates x T_step = 10
#       5 time points x 10 biological replicates x T_step = 10
#
#   Platforms:
#       RT-qPCR, GAUSS, RNA-seq
#
#   Noise:
#       Very low, Low, Medium, High
#
#   Genes:
#       500 null + 250 alternative
#
#   Bootstrap:
#       1,999 replicates/test
#
# Outputs:
#   fraction_robustness_benchmark/
#       fraction_robustness_raw.tsv
#       fraction_robustness_summary.tsv
#       fraction_robustness_design_summary.tsv
#       fraction_robustness_progress.tsv
#       Fig_fraction_scaling_typeI.pdf
#       Fig_fraction_contamination_typeI.pdf
#       Fig_fraction_scaling_sigma_bias.pdf
#       Fig_fraction_contamination_sigma_bias.pdf
#       fraction_robustness_benchmark.rdata
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
# 1. Files and settings
# =============================================================================

INPUT_FILE <- "ode_states_2k_20p_corrected_onset.rdata"

OUTPUT_DIR <- "fraction_robustness_benchmark"

CHECKPOINT_DIR <- file.path(
  OUTPUT_DIR,
  "checkpoints"
)

PROGRESS_FILE <- file.path(
  OUTPUT_DIR,
  "fraction_robustness_progress.tsv"
)

RAW_FILE <- file.path(
  OUTPUT_DIR,
  "fraction_robustness_raw.tsv"
)

SUMMARY_FILE <- file.path(
  OUTPUT_DIR,
  "fraction_robustness_summary.tsv"
)

DESIGN_SUMMARY_FILE <- file.path(
  OUTPUT_DIR,
  "fraction_robustness_design_summary.tsv"
)

RDATA_FILE <- file.path(
  OUTPUT_DIR,
  "fraction_robustness_benchmark.rdata"
)

SESSION_FILE <- file.path(
  OUTPUT_DIR,
  "sessionInfo.txt"
)

BENCHMARK_VERSION <- "fraction_robustness_corrected_onset_v1"


# -----------------------------------------------------------------------------
# Inferential settings: identical to the principal benchmark.
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
# Gene counts.
# -----------------------------------------------------------------------------

N_NULL <- 500L
N_ALT <- 250L


# -----------------------------------------------------------------------------
# Misspecification levels.
#
# Scaling:
#   1.0 is the correctly specified reference.
#
# Contamination:
#   0.0 is the correctly specified reference.
#
# The two axes are studied SEPARATELY, not factorially combined.
# -----------------------------------------------------------------------------

SCALING_LEVELS <- c(
  0.50,
  0.75,
  1.00,
  1.25,
  1.50,
  2.00
)

CONTAMINATION_LEVELS <- c(
  0.000,
  0.010,
  0.025,
  0.050,
  0.100,
  0.200
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
# 4. Input checks
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
# Verify corrected gene-level timing semantics.
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
    "Timing integrity check failed."
  )
}


# -----------------------------------------------------------------------------
# TRUE complete-shutoff trajectories only.
# -----------------------------------------------------------------------------

ode_shutoff <- ode_states[
  Perturbation ==
    "SHUTOFF" &
    abs(
      Post_R_fraction
    ) <
      1e-12
]

if (
  nrow(
    ode_shutoff
  ) == 0L
) {
  stop(
    "No complete-SHUTOFF trajectories found."
  )
}


# -----------------------------------------------------------------------------
# Restrict to selected designs.
# -----------------------------------------------------------------------------

design_keys <- DESIGNS[
  ,
  .(
    N_time_samples,
    N_replicates,
    T_step
  )
]

ode_shutoff <- merge(
  ode_shutoff,
  design_keys,
  by = c(
    "N_time_samples",
    "N_replicates",
    "T_step"
  ),
  all = FALSE
)

if (
  nrow(
    ode_shutoff
  ) == 0L
) {
  stop(
    "Selected designs were not found in complete-SHUTOFF data."
  )
}

setkey(
  ode_shutoff,
  Gene,
  N_time_samples,
  N_replicates,
  T_step
)

cat(
  "\nInput preflight checks PASSED.\n"
)


# =============================================================================
# 5. Select paired genes
# =============================================================================

gene_truth <- unique(
  ode_shutoff[
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

ode_shutoff <- ode_shutoff[
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
# 6. Misspecification transforms
# =============================================================================

apply_fraction_scaling <- function(
  dt,
  scale_cyt
) {

  out <- copy(
    as.data.table(
      dt
    )
  )

  out[
    ,
    C :=
      C *
        scale_cyt
  ]

  out[
    ,
    C_s :=
      C_s *
        scale_cyt
  ]

  out
}


apply_cross_contamination <- function(
  dt,
  epsilon
) {

  out <- copy(
    as.data.table(
      dt
    )
  )

  N0 <- out$N
  Ns0 <- out$N_s
  C0 <- out$C
  Cs0 <- out$C_s

  out[
    ,
    N :=
      (
        1 -
          epsilon
      ) *
        N0 +
        epsilon *
          C0
  ]

  out[
    ,
    N_s :=
      (
        1 -
          epsilon
      ) *
        Ns0 +
        epsilon *
          Cs0
  ]

  out[
    ,
    C :=
      epsilon *
        N0 +
        (
          1 -
            epsilon
        ) *
          C0
  ]

  out[
    ,
    C_s :=
      epsilon *
        Ns0 +
        (
          1 -
            epsilon
        ) *
          Cs0
  ]

  out
}


# =============================================================================
# 7. Platform simulation helper
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
# 8. Reproducible seeds
#
# The misspecification level is deliberately EXCLUDED from measurement and
# bootstrap seeds. This creates paired Monte Carlo streams across levels.
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
  analysis,
  stream
) {

  stable_seed_from_string(
    paste(
      BENCHMARK_VERSION,
      gene,
      design,
      platform,
      noise,
      analysis,
      stream,
      sep = "|"
    )
  )
}


# =============================================================================
# 9. Safe helpers
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
# 10. Run one misspecified scenario
# =============================================================================

run_one_test <- function(
  dt_latent,
  gene_i,
  truth_i,
  sigma_true,
  design_name,
  n_time,
  n_rep,
  tstep,
  platform_i,
  noise_i,
  analysis_type,
  misspec_value,
  t_star_fit,
  N_boot
) {

  if (
    analysis_type ==
      "SCALING"
  ) {

    dt_misspec <- apply_fraction_scaling(
      dt_latent,
      scale_cyt =
        misspec_value
    )

    scale_cyt <- misspec_value
    contamination <- 0

  } else if (
    analysis_type ==
      "CONTAMINATION"
  ) {

    dt_misspec <- apply_cross_contamination(
      dt_latent,
      epsilon =
        misspec_value
    )

    scale_cyt <- 1
    contamination <- misspec_value

  } else {

    stop(
      paste(
        "Unknown analysis_type:",
        analysis_type
      )
    )
  }


  measurement_seed <- make_seed(
    gene =
      gene_i,
    design =
      design_name,
    platform =
      platform_i,
    noise =
      noise_i,
    analysis =
      analysis_type,
    stream =
      "measurement"
  )


  bootstrap_seed <- make_seed(
    gene =
      gene_i,
    design =
      design_name,
    platform =
      platform_i,
    noise =
      noise_i,
    analysis =
      analysis_type,
    stream =
      "bootstrap"
  )


  set.seed(
    measurement_seed
  )


  dt_noisy <- tryCatch(

    add_platform_noise(
      dt =
        dt_misspec,
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

    return(
      data.table(
        Benchmark_version =
          BENCHMARK_VERSION,
        Gene =
          gene_i,
        Positive =
          truth_i,
        sigma_true =
          sigma_true,
        Design =
          design_name,
        N_tsamples =
          n_time,
        N_replicates =
          n_rep,
        Tsteps =
          tstep,
        Platform =
          platform_i,
        Exprs_noise =
          noise_i,
        Analysis =
          analysis_type,
        Misspecification_value =
          misspec_value,
        Cytoplasmic_scale =
          scale_cyt,
        Contamination_fraction =
          contamination,
        T_star =
          t_star_fit,
        N_boot =
          N_boot,
        status =
          "measurement_error",
        error_message =
          dt_noisy$error_message
      )
    )
  }


  t0 <- Sys.time()


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
      t0,
      units =
        "secs"
    )
  )


  coef_full_i <- WW$coef_full


  data.table(
    Benchmark_version =
      BENCHMARK_VERSION,
    Gene =
      gene_i,
    Positive =
      truth_i,
    sigma_true =
      sigma_true,
    Design =
      design_name,
    N_tsamples =
      n_time,
    N_replicates =
      n_rep,
    Tsteps =
      tstep,
    Platform =
      platform_i,
    Exprs_noise =
      noise_i,
    Analysis =
      analysis_type,
    Misspecification_value =
      misspec_value,
    Cytoplasmic_scale =
      scale_cyt,
    Contamination_fraction =
      contamination,
    T_star =
      t_star_fit,
    N_boot =
      N_boot,
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
}


# =============================================================================
# 11. Run one gene
# =============================================================================

run_gene_fraction_robustness <- function(
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


    dt_latent <- ode_shutoff[
      Gene ==
        gene_i &
        N_time_samples ==
          design_i$N_time_samples &
        N_replicates ==
          design_i$N_replicates &
        T_step ==
          design_i$T_step
    ]


    if (
      nrow(
        dt_latent
      ) == 0L
    ) {
      stop(
        paste(
          "Missing latent data for gene",
          gene_i,
          "design",
          design_i$Design
        )
      )
    }


    truth_values <- unique(
      dt_latent$truth_pos
    )

    sigma_values <- unique(
      dt_latent$Base_sigma_c
    )

    tstar_values <- unique(
      dt_latent$T_star
    )

    tstar_values <- tstar_values[
      is.finite(
        tstar_values
      )
    ]


    if (
      length(
        truth_values
      ) !=
        1L ||
        length(
          sigma_values
        ) !=
          1L ||
        length(
          tstar_values
        ) !=
          1L
    ) {
      stop(
        paste(
          "Non-unique truth/sigma/T_star for gene",
          gene_i
        )
      )
    }


    truth_i <- truth_values[1]
    sigma_true <- sigma_values[1]
    t_star_fit <- tstar_values[1]


    # -------------------------------------------------------------------------
    # A. Relative fraction scaling
    # -------------------------------------------------------------------------

    for (
      scale_i in SCALING_LEVELS
    ) {

      for (
        platform_i in RANGE_PLATFORM
      ) {

        for (
          noise_i in RANGE_NOISE
        ) {

          out[[kk]] <- run_one_test(
            dt_latent =
              dt_latent,
            gene_i =
              gene_i,
            truth_i =
              truth_i,
            sigma_true =
              sigma_true,
            design_name =
              design_i$Design,
            n_time =
              design_i$N_time_samples,
            n_rep =
              design_i$N_replicates,
            tstep =
              design_i$T_step,
            platform_i =
              platform_i,
            noise_i =
              noise_i,
            analysis_type =
              "SCALING",
            misspec_value =
              scale_i,
            t_star_fit =
              t_star_fit,
            N_boot =
              N_boot
          )

          kk <- kk + 1L
        }
      }
    }


    # -------------------------------------------------------------------------
    # B. Cross-fraction contamination
    # -------------------------------------------------------------------------

    for (
      eps_i in CONTAMINATION_LEVELS
    ) {

      for (
        platform_i in RANGE_PLATFORM
      ) {

        for (
          noise_i in RANGE_NOISE
        ) {

          out[[kk]] <- run_one_test(
            dt_latent =
              dt_latent,
            gene_i =
              gene_i,
            truth_i =
              truth_i,
            sigma_true =
              sigma_true,
            design_name =
              design_i$Design,
            n_time =
              design_i$N_time_samples,
            n_rep =
              design_i$N_replicates,
            tstep =
              design_i$T_step,
            platform_i =
              platform_i,
            noise_i =
              noise_i,
            analysis_type =
              "CONTAMINATION",
            misspec_value =
              eps_i,
            t_star_fit =
              t_star_fit,
            N_boot =
              N_boot
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
# 12. Checkpoints
# =============================================================================

EXPECTED_TESTS_PER_GENE <- (
  nrow(
    DESIGNS
  ) *
    length(
      RANGE_PLATFORM
    ) *
    length(
      RANGE_NOISE
    ) *
    (
      length(
        SCALING_LEVELS
      ) +
        length(
          CONTAMINATION_LEVELS
        )
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
    "Analysis",
    "Misspecification_value",
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
    x$Analysis,
    sprintf(
      "%.3f",
      x$Misspecification_value
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

    run_gene_fraction_robustness(
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
# 13. Serial smoke test
# =============================================================================

SMOKE_BOOT <- 19L
smoke_gene <- GENES_KEEP[1]


cat(
  "\n============================================================\n",
  "SERIAL FRACTION-ROBUSTNESS SMOKE TEST\n",
  "============================================================\n",
  "Gene: ",
  smoke_gene,
  "\n",
  "Expected tests/gene: ",
  EXPECTED_TESTS_PER_GENE,
  "\n",
  sep = ""
)


smoke <- run_gene_fraction_robustness(
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
  smoke$Analysis,
  sprintf(
    "%.3f",
    smoke$Misspecification_value
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
# 14. Determine remaining genes
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
# 15. Parallel run
# =============================================================================

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
        "ode_shutoff",
        "DESIGNS",
        "SCALING_LEVELS",
        "CONTAMINATION_LEVELS",
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
        "apply_fraction_scaling",
        "apply_cross_contamination",
        "add_platform_noise",
        "stable_seed_from_string",
        "make_seed",
        "safe_value",
        "safe_text",
        "safe_coef",
        "run_one_test",
        "run_gene_fraction_robustness",
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
# 16. Collect checkpoints
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
  Analysis,
  Design,
  Platform,
  Exprs_noise,
  Misspecification_value
)


# =============================================================================
# 17. Raw validation
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
# 18. Wilson helper
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
# 19. Scenario summary
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
      k0,
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


    .(
      N_null =
        n0,

      N_alt =
        n1,

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
        if (
          n1 >
            0L
        ) {
          mean(
            alt$p.value <=
              0.05
          )
        } else {
          NA_real_
        },

      Null_sigma_hat_mean =
        if (
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
        },

      Alt_sigma_bias =
        if (
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
        },

      Alt_sigma_RMSE =
        if (
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
        },

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
    Analysis,
    Design,
    N_tsamples,
    N_replicates,
    Tsteps,
    Platform,
    Exprs_noise,
    Misspecification_value,
    Cytoplasmic_scale,
    Contamination_fraction
  )
]


setorder(
  summary_dt,
  Analysis,
  Design,
  Platform,
  Exprs_noise,
  Misspecification_value
)


# =============================================================================
# 20. Design-level summary
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
    Analysis,
    Design,
    N_tsamples,
    N_replicates,
    Tsteps,
    Misspecification_value,
    Cytoplasmic_scale,
    Contamination_fraction
  )
]


setorder(
  design_summary,
  Analysis,
  Design,
  Misspecification_value
)


# =============================================================================
# 21. Save tables
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
# 22. Figures
# =============================================================================

theme_fraction <- theme_classic(
  base_size =
    10.5
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
# A. Scaling: Type-I
# -----------------------------------------------------------------------------

scaling_dt <- summary_dt[
  Analysis ==
    "SCALING"
]


p_scaling_typeI <- ggplot(
  scaling_dt,
  aes(
    x =
      Cytoplasmic_scale,
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

  geom_hline(
    yintercept =
      0.05,
    linetype =
      "dashed",
    linewidth =
      0.5,
    colour =
      "grey35"
  ) +

  geom_vline(
    xintercept =
      1,
    linetype =
      "dotted",
    linewidth =
      0.45,
    colour =
      "grey50"
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
      SCALING_LEVELS
  ) +

  labs(
    x =
      "Cytoplasmic relative scaling factor",
    y =
      "Empirical Type-I error",
    linetype =
      "Noise",
    shape =
      "Noise"
  ) +

  theme_fraction


ggsave(
  file.path(
    OUTPUT_DIR,
    "Fig_fraction_scaling_typeI.pdf"
  ),
  p_scaling_typeI,
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
# B. Contamination: Type-I
# -----------------------------------------------------------------------------

contam_dt <- summary_dt[
  Analysis ==
    "CONTAMINATION"
]


p_contam_typeI <- ggplot(
  contam_dt,
  aes(
    x =
      100 *
        Contamination_fraction,
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

  geom_hline(
    yintercept =
      0.05,
    linetype =
      "dashed",
    linewidth =
      0.5,
    colour =
      "grey35"
  ) +

  geom_vline(
    xintercept =
      0,
    linetype =
      "dotted",
    linewidth =
      0.45,
    colour =
      "grey50"
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
        CONTAMINATION_LEVELS
  ) +

  labs(
    x =
      "Cross-fraction contamination (%)",
    y =
      "Empirical Type-I error",
    linetype =
      "Noise",
    shape =
      "Noise"
  ) +

  theme_fraction


ggsave(
  file.path(
    OUTPUT_DIR,
    "Fig_fraction_contamination_typeI.pdf"
  ),
  p_contam_typeI,
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
# C. Scaling: spurious sigma_c under H0
# -----------------------------------------------------------------------------

p_scaling_sigma <- ggplot(
  scaling_dt,
  aes(
    x =
      Cytoplasmic_scale,
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

  geom_vline(
    xintercept =
      1,
    linetype =
      "dotted",
    linewidth =
      0.45,
    colour =
      "grey50"
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
      SCALING_LEVELS
  ) +

  labs(
    x =
      "Cytoplasmic relative scaling factor",
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

  theme_fraction


ggsave(
  file.path(
    OUTPUT_DIR,
    "Fig_fraction_scaling_sigma_bias.pdf"
  ),
  p_scaling_sigma,
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
# D. Contamination: spurious sigma_c under H0
# -----------------------------------------------------------------------------

p_contam_sigma <- ggplot(
  contam_dt,
  aes(
    x =
      100 *
        Contamination_fraction,
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
        CONTAMINATION_LEVELS
  ) +

  labs(
    x =
      "Cross-fraction contamination (%)",
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

  theme_fraction


ggsave(
  file.path(
    OUTPUT_DIR,
    "Fig_fraction_contamination_sigma_bias.pdf"
  ),
  p_contam_sigma,
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
# 23. Save analysis object
# =============================================================================

fraction_robustness_benchmark <- list(

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
      scaling_levels =
        SCALING_LEVELS,
      contamination_levels =
        CONTAMINATION_LEVELS,
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
  fraction_robustness_benchmark,
  file =
    RDATA_FILE
)


# =============================================================================
# 24. Session information
# =============================================================================

sink(
  SESSION_FILE
)

print(
  sessionInfo()
)

sink()


# =============================================================================
# 25. Final report
# =============================================================================

cat(
  "\n============================================================\n",
  "FRACTION-SEPARATION ROBUSTNESS BENCHMARK COMPLETE\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nDesign-level summary:\n"
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
  "  SCALING studies a multiplicative cytoplasmic-vs-nuclear mismatch while\n",
  "  inference assumes the original compartment scale.\n",
  "  CONTAMINATION studies symmetric cross-fraction mixing while inference\n",
  "  assumes pure nuclear and cytoplasmic fractions.\n",
  "  Primary endpoint: empirical Type-I error under sigma_c = 0.\n",
  "  Secondary endpoints: power, sigma_c bias/RMSE, NNLS boundary behavior,\n",
  "  conditioning, and bootstrap stability.\n",
  sep = ""
)
