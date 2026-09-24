# =============================================================================
# Title:
# Revised Bootstrap-Calibrated NNLS Analysis of Real Compartment-Resolved
# RNA-seq Datasets
#
# Description:
#   This script applies the revised compartment-resolved kinetic model and
#   bootstrap-calibrated constrained NNLS test to the real RNA-seq datasets
#   analyzed in the manuscript.
#
#   The analysis was revised in response to peer-review concerns regarding:
#
#     - uncertainty in both A and b;
#     - correlation among adjacent interval differences;
#     - destructive biological sampling across time points;
#     - non-negativity and the sigma_c = 0 boundary;
#     - multiple-testing correction;
#     - distinction between true pharmacological shutoff and pseudo-shutoff;
#     - practical identifiability and numerical conditioning.
#
# Statistical inference:
#   The revised test:
#
#     1. treats biological replicates as independent destructive samples at
#        each time point and does not pair replicate identifiers longitudinally;
#
#     2. estimates within-time covariance matrices for
#        N, N_s, C, and C_s;
#
#     3. propagates the full covariance structure of interval differences,
#        including covariance between adjacent intervals sharing a time point;
#
#     4. fits null (sigma_c = 0) and full (sigma_c >= 0) models by generalized
#        weighted non-negative least squares;
#
#     5. generates replicate-level parametric bootstrap datasets under H0 and
#        reconstructs A*, b*, and Sigma_b* in every bootstrap replicate;
#
#     6. uses the deterministic add-one bootstrap p-value
#
#          p = (1 + #{T* >= Tobs}) / (B + 1).
#
# Biological interpretation:
#   sigma_c is interpreted as a phenomenological post-export conversion rate.
#   Statistical evidence for sigma_c > 0 identifies kinetic trajectories
#   consistent with an additional post-export conversion component and does
#   not by itself establish cytoplasmic splicing or another specific molecular
#   mechanism.
#
# Experimental regimes:
#
#   Kc167 | GSE83620
#       PSEUDO_SHUTOFF
#       Pre-existing RNA was experimentally isolated after 4sU labeling by
#       depletion of newly synthesized labeled RNA.
#
#   K562 | GSE207924
#       PSEUDO_SHUTOFF
#       Event-level abundances were combined with gene-level estimates of the
#       pre-existing RNA fraction.
#
#   3T3 | GSE207924
#       PSEUDO_SHUTOFF
#       Same conceptual reconstruction as above.
#
#   mESC | GSE256335
#       TRUE_SHUTOFF
#       Direct pharmacological transcriptional inhibition.
#
# Multiple testing:
#   Benjamini-Hochberg adjusted q-values are calculated separately within each
#   dataset. Raw p-values, q-values, sigma_c, RSS improvement, and practical
#   identifiability diagnostics are all retained.
#
# Candidate interpretation:
#   IR and sigma_c are retained as effect/model-improvement descriptors.
#   They are NOT used as substitutes for FDR control. Results not surviving
#   conventional FDR thresholds should be interpreted as exploratory rankings.
#
# Inputs:
#   GSE83620/rmats_dt_Kc167.rdata
#   GSE207924/rmats_dt_K562.rdata
#   GSE207924/rmats_dt_3T3.rdata
#   GSE256335/rmats_dt_ECS.rdata
#
# Outputs:
#   results_realdatasets_revision.rdata
#   results_realdatasets_revision.tsv
#   results_realdatasets_revision.xlsx
#   real_dataset_qc_summary.tsv
#   real_dataset_discovery_summary.tsv
#   real_dataset_sessionInfo.txt
#
# Author:
#   Luigi Cerulo
#
# Copyright:
#   Copyright (c) 2026 Luigi Cerulo
#
# License:
#   Permission is hereby granted to use, copy, modify, and distribute this
#   software for academic and research purposes, provided that this notice is
#   retained in all copies or substantial portions of the software.
#
# Disclaimer:
#   This software is provided "as is", without warranty of any kind, express
#   or implied, including but not limited to the warranties of merchantability,
#   fitness for a particular purpose, and noninfringement. In no event shall
#   the author be liable for any claim, damages, or other liability arising
#   from, out of, or in connection with the software or its use.
#
# Version:
#   Major-revision real-data analysis, 2026.
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/real_datasets")

# Avoid hidden nested BLAS/OpenMP parallelism.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

library(data.table)
library(parallel)
library(openxlsx)

data.table::setDTthreads(1)

# IMPORTANT: revised test.
source("../commons/nested_test2.r")


# =============================================================================
# 1. Analysis settings
# =============================================================================

N_BOOT <- 1999

SCALING_A <- TRUE

# Frozen covariance-shrinkage settings used in the revised benchmark.
LAMBDA_TIME <- 0.5
LAMBDA_DIAG <- 0.1

REL_FLOOR <- 1e-8

TRUNCATE_NONNEGATIVE_BOOT <- FALSE

MAX_FAILURE_RATE <- 0.05

# Use only CPUs that are currently available.
# Change this according to the machine load caused by the synthetic benchmark.
N_WORKERS <- 30L

# Intervention/reference time.
#
# IMPORTANT:
# This assumes that the time variables in all four processed rmats_dt objects
# are already expressed relative to the beginning of the true/pseudo shutoff.
#
# Verify this below before running the analysis.
T_STAR_DEFAULT <- 0


# =============================================================================
# 2. Load datasets
# =============================================================================

load("GSE83620/rmats_dt_Kc167.rdata")

dt_kc167 <- copy(rmats_dt)

dt_kc167[
  ,
  `:=`(
    dataset = "Kc167 | GSE83620",
    perturbation_type = "PSEUDO_SHUTOFF",
    pseudo_method = "4sU_unlabelled_preexisting"
  )
]


load("GSE207924/rmats_dt_K562.rdata")

dt_k562 <- copy(rmats_dt)

dt_k562[
  ,
  `:=`(
    dataset = "K562 | GSE207924",
    perturbation_type = "PSEUDO_SHUTOFF",
    pseudo_method = "gene_fraction_reconstruction"
  )
]


load("GSE207924/rmats_dt_3T3.rdata")

dt_3t3 <- copy(rmats_dt)

dt_3t3[
  ,
  `:=`(
    dataset = "3T3 | GSE207924",
    perturbation_type = "PSEUDO_SHUTOFF",
    pseudo_method = "gene_fraction_reconstruction"
  )
]


load("GSE256335/rmats_dt_ECS.rdata")

dt_esc <- copy(rmats_dt)

dt_esc[
  ,
  `:=`(
    dataset = "mESC | GSE256335",
    perturbation_type = "TRUE_SHUTOFF",
    pseudo_method = NA_character_
  )
]


rm(rmats_dt)


# =============================================================================
# 3. Merge
# =============================================================================

rmats_all <- rbindlist(
  list(
    dt_kc167,
    dt_k562,
    dt_3t3,
    dt_esc
  ),
  use.names = TRUE,
  fill = TRUE
)


# Required columns.
required_cols <- c(
  "dataset",
  "event",
  "ensembl",
  "time",
  "replicate",
  "N",
  "N_s",
  "C",
  "C_s"
)

missing_cols <- setdiff(
  required_cols,
  names(rmats_all)
)

if (length(missing_cols) > 0L) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}


# =============================================================================
# 4. Critical time-scale inspection
# =============================================================================

cat("\n============================================================\n")
cat("TIME VARIABLE CHECK\n")
cat("============================================================\n")

print(
  rmats_all[
    ,
    .(
      times = paste(
        sort(unique(time)),
        collapse = ", "
      ),
      min_time = min(time, na.rm = TRUE),
      max_time = max(time, na.rm = TRUE),
      n_timepoints = uniqueN(time),
      n_replicates = uniqueN(replicate)
    ),
    by = .(
      dataset,
      perturbation_type
    )
  ]
)

cat(
  "\nVerify that these times are on the intended common kinetic time scale\n",
  "and that time zero corresponds to the start/reference of the\n",
  "true/pseudo-shutoff analysis.\n\n",
  sep = ""
)


# =============================================================================
# 5. Save merged input for reproducibility
# =============================================================================

save(
  rmats_all,
  file = "rmats_all_revision.rdata"
)


# =============================================================================
# 6. Event-level QC
# =============================================================================

evt_qc <- rmats_all[
  ,
  {

    safe_range <- function(x) {

      z <- x[
        is.finite(x)
      ]

      if (length(z) == 0L) {
        return(NA_real_)
      }

      max(z) - min(z)
    }


    list(

      cov_total =
        sum(
          N + N_s + C + C_s,
          na.rm = TRUE
        ),

      cov_cyt =
        sum(
          C + C_s,
          na.rm = TRUE
        ),

      rng_C =
        safe_range(C),

      rng_Cs =
        safe_range(C_s),

      nz_N =
        sum(
          N > 0,
          na.rm = TRUE
        ),

      nz_Ns =
        sum(
          N_s > 0,
          na.rm = TRUE
        ),

      nz_C =
        sum(
          C > 0,
          na.rm = TRUE
        ),

      nz_Cs =
        sum(
          C_s > 0,
          na.rm = TRUE
        ),

      frac_nzN =
        mean(
          N > 0,
          na.rm = TRUE
        ),

      frac_nzNs =
        mean(
          N_s > 0,
          na.rm = TRUE
        ),

      frac_nzC =
        mean(
          C > 0,
          na.rm = TRUE
        ),

      frac_nzCs =
        mean(
          C_s > 0,
          na.rm = TRUE
        ),

      n_time =
        uniqueN(time),

      n_rep =
        uniqueN(replicate)
    )
  },

  by = .(
    dataset,
    event
  )
]


# =============================================================================
# 7. Dataset-specific empirical QC thresholds
# =============================================================================

thr <- evt_qc[
  ,
  .(

    cov_total_min =
      as.numeric(
        quantile(
          cov_total,
          0.01,
          na.rm = TRUE
        )
      ),

    cov_cyt_min =
      as.numeric(
        quantile(
          cov_cyt,
          0.01,
          na.rm = TRUE
        )
      ),

    rng_C_min =
      as.numeric(
        quantile(
          rng_C,
          0.01,
          na.rm = TRUE
        )
      ),

    rng_Cs_min =
      as.numeric(
        quantile(
          rng_Cs,
          0.01,
          na.rm = TRUE
        )
      )
  ),

  by = dataset
]


# =============================================================================
# 8. Explicit event-filter function
# =============================================================================

filter_event <- function(
  d,
  threshold_row
) {

  min_pos <- 3L
  min_frac <- 0.30
  min_pos_nuc <- 3L

  max_zero_frac_C <- 0.70
  max_zero_frac_Cs <- 0.80


  bad_all_zero <- (
    all(d$N == 0, na.rm = TRUE) ||
    all(d$N_s == 0, na.rm = TRUE) ||
    all(d$C == 0, na.rm = TRUE) ||
    all(d$C_s == 0, na.rm = TRUE)
  )


  bad_cyt <- (
    sum(d$C > 0, na.rm = TRUE) < min_pos ||
    sum(d$C_s > 0, na.rm = TRUE) < min_pos ||
    mean(d$C > 0, na.rm = TRUE) < min_frac
  )


  bad_nuc <- (
    sum(d$N > 0, na.rm = TRUE) < min_pos_nuc ||
    sum(d$N_s > 0, na.rm = TRUE) < min_pos_nuc
  )


  cov_total_i <- sum(
    d$N + d$N_s + d$C + d$C_s,
    na.rm = TRUE
  )


  cov_cyt_i <- sum(
    d$C + d$C_s,
    na.rm = TRUE
  )


  rng_C_i <- diff(
    range(
      d$C,
      na.rm = TRUE
    )
  )


  rng_Cs_i <- diff(
    range(
      d$C_s,
      na.rm = TRUE
    )
  )


  bad_coverage <- (
    cov_total_i < threshold_row$cov_total_min ||
    cov_cyt_i < threshold_row$cov_cyt_min
  )


  bad_range <- (
    rng_C_i < threshold_row$rng_C_min ||
    rng_Cs_i < threshold_row$rng_Cs_min
  )


  bad_zero <- (
    mean(
      d$C == 0,
      na.rm = TRUE
    ) > max_zero_frac_C ||
    mean(
      d$C_s == 0,
      na.rm = TRUE
    ) > max_zero_frac_Cs
  )


  bad_replication <- any(
    d[
      ,
      .N,
      by = time
    ]$N < 2L
  )


  bad <- (
    bad_all_zero ||
    bad_cyt ||
    bad_nuc ||
    bad_coverage ||
    bad_range ||
    bad_zero ||
    bad_replication
  )


  reasons <- character()


  if (bad_all_zero) {
    reasons <- c(
      reasons,
      "state_all_zero"
    )
  }

  if (bad_cyt) {
    reasons <- c(
      reasons,
      "cytoplasmic_information"
    )
  }

  if (bad_nuc) {
    reasons <- c(
      reasons,
      "nuclear_information"
    )
  }

  if (bad_coverage) {
    reasons <- c(
      reasons,
      "coverage"
    )
  }

  if (bad_range) {
    reasons <- c(
      reasons,
      "dynamic_range"
    )
  }

  if (bad_zero) {
    reasons <- c(
      reasons,
      "zero_inflation"
    )
  }

  if (bad_replication) {
    reasons <- c(
      reasons,
      "insufficient_within_time_replication"
    )
  }


  list(
    pass = !bad,
    reason = if (
      length(reasons) == 0L
    ) {
      "pass"
    } else {
      paste(
        reasons,
        collapse = ";"
      )
    }
  )
}


# =============================================================================
# 9. Construct event task table
# =============================================================================

tasks <- unique(
  rmats_all[
    ,
    .(
      dataset,
      perturbation_type,
      pseudo_method,
      event,
      ensembl
    )
  ]
)

tasks[
  ,
  task_id := .I
]


# =============================================================================
# 10. Reproducible seed
# =============================================================================

make_real_seed <- function(
  dataset,
  event
) {

  chars <- utf8ToInt(
    paste0(
      dataset,
      "::",
      event
    )
  )

  value <- sum(
    chars *
      seq_along(chars)
  )

  as.integer(
    (
      700000000 +
        value
    ) %%
      2000000000
  )
}


# =============================================================================
# 11. Single-event worker
# =============================================================================

run_real_event <- function(
  task_row,
  N_boot
) {

  dataset_i <- as.character(
    task_row$dataset
  )

  event_i <- as.character(
    task_row$event
  )

  ensembl_i <- as.character(
    task_row$ensembl
  )


  d <- rmats_all[
    dataset == dataset_i &
      event == event_i &
      ensembl == ensembl_i
  ]


  if (nrow(d) == 0L) {
    return(NULL)
  }


  # ---------------------------------------------------------------------------
  # Correct dataset-specific threshold selection.
  # ---------------------------------------------------------------------------

  threshold_row <- thr[
    dataset == dataset_i
  ]


  if (nrow(threshold_row) != 1L) {

    return(
      data.table(
        dataset = dataset_i,
        event = event_i,
        ensembl = ensembl_i,
        filter_pass = FALSE,
        filter_reason = "threshold_lookup_failed",
        p.value = NA_real_,
        status = "not_tested"
      )
    )
  }


  filter_result <- filter_event(
    d,
    threshold_row
  )


  perturbation_i <- d$perturbation_type[1]

  pseudo_method_i <- d$pseudo_method[1]


  if (!filter_result$pass) {

    return(
      data.table(

        dataset =
          dataset_i,

        perturbation_type =
          perturbation_i,

        pseudo_method =
          pseudo_method_i,

        event =
          event_i,

        ensembl =
          ensembl_i,

        gene_symbol =
          d$gene_symbol[1],

        description =
          d$description[1],

        gene_biotype =
          d$gene_biotype[1],

        n_timepoints =
          uniqueN(
            d$time
          ),

        n_replicates =
          uniqueN(
            d$replicate
          ),

        filter_pass =
          FALSE,

        filter_reason =
          filter_result$reason,

        p.value =
          NA_real_,

        status =
          "filtered"
      )
    )
  }


  seed_i <- make_real_seed(
    dataset_i,
    event_i
  )


  start_i <- Sys.time()


  WWW <- tryCatch(

    test_sigma_nested(

      tsampled_data =
        d,

      scaling_A =
        SCALING_A,

      t_star =
        T_STAR_DEFAULT,

      B_n =
        N_boot,

      seed =
        seed_i,

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


  elapsed_i <- as.numeric(
    difftime(
      Sys.time(),
      start_i,
      units = "secs"
    )
  )


  coef_full_i <- WWW$coef_full

  coef_null_i <- WWW$coef_null


  get_coef <- function(
    x,
    nm
  ) {

    if (
      is.null(x) ||
      is.null(names(x)) ||
      !nm %in% names(x)
    ) {
      return(NA_real_)
    }

    unname(
      x[nm]
    )
  }


  get_field <- function(
    x,
    nm,
    default = NA_real_
  ) {

    if (
      is.null(x[[nm]]) ||
      length(x[[nm]]) == 0L
    ) {
      return(default)
    }

    x[[nm]][1]
  }


  data.table(

    dataset =
      dataset_i,

    perturbation_type =
      perturbation_i,

    pseudo_method =
      pseudo_method_i,

    event =
      event_i,

    ensembl =
      ensembl_i,

    gene_symbol =
      d$gene_symbol[1],

    description =
      d$description[1],

    gene_biotype =
      d$gene_biotype[1],

    n_timepoints =
      uniqueN(
        d$time
      ),

    n_replicates =
      uniqueN(
        d$replicate
      ),

    filter_pass =
      TRUE,

    filter_reason =
      "pass",


    # -------------------------------------------------------------------------
    # Inference
    # -------------------------------------------------------------------------

    p.value =
      get_field(
        WWW,
        "p.value"
      ),

    Sigma =
      get_field(
        WWW,
        "Sigma"
      ),

    Alpha =
      get_field(
        WWW,
        "Alpha"
      ),

    T.obs =
      get_field(
        WWW,
        "T.obs"
      ),

    RSS0 =
      get_field(
        WWW,
        "RSS0"
      ),

    RSS1 =
      get_field(
        WWW,
        "RSS1"
      ),

    IR =
      get_field(
        WWW,
        "IR"
      ),


    # -------------------------------------------------------------------------
    # Boundary / bootstrap diagnostics
    # -------------------------------------------------------------------------

    atom_zero =
      get_field(
        WWW,
        "atom.zero"
      ),

    bootstrap_failure_rate =
      get_field(
        WWW,
        "bootstrap.failure.rate"
      ),

    n_bootstrap_valid =
      get_field(
        WWW,
        "n.bootstrap.valid"
      ),

    bootstrap_condition_median =
      get_field(
        WWW,
        "bootstrap.condition.median"
      ),

    bootstrap_condition_q95 =
      get_field(
        WWW,
        "bootstrap.condition.q95"
      ),

    bootstrap_condition_max =
      get_field(
        WWW,
        "bootstrap.condition.max"
      ),

    bootstrap_rank_deficient_fraction =
      get_field(
        WWW,
        "bootstrap.rank.deficient.fraction"
      ),


    # -------------------------------------------------------------------------
    # Observed-system identifiability diagnostics
    # -------------------------------------------------------------------------

    condition_number =
      get_field(
        WWW,
        "condition.number"
      ),

    min_singular_value =
      get_field(
        WWW,
        "min.singular.value"
      ),

    rank_full =
      get_field(
        WWW,
        "rank.full"
      ),

    rank_null =
      get_field(
        WWW,
        "rank.null"
      ),


    # -------------------------------------------------------------------------
    # Full-model kinetic estimates
    # -------------------------------------------------------------------------

    R_hat =
      get_coef(
        coef_full_i,
        "R"
      ),

    tau_hat =
      get_coef(
        coef_full_i,
        "tau"
      ),

    tau_s_hat =
      get_coef(
        coef_full_i,
        "tau_s"
      ),

    sigma_c_hat =
      get_coef(
        coef_full_i,
        "sigma_c"
      ),

    sigma_n_hat =
      get_coef(
        coef_full_i,
        "sigma_n"
      ),

    alpha_hat =
      get_coef(
        coef_full_i,
        "alpha"
      ),

    alpha_s_hat =
      get_coef(
        coef_full_i,
        "alpha_s"
      ),


    # -------------------------------------------------------------------------
    # Null-model estimates
    # -------------------------------------------------------------------------

    R_hat_null =
      get_coef(
        coef_null_i,
        "R"
      ),

    tau_hat_null =
      get_coef(
        coef_null_i,
        "tau"
      ),

    tau_s_hat_null =
      get_coef(
        coef_null_i,
        "tau_s"
      ),

    sigma_n_hat_null =
      get_coef(
        coef_null_i,
        "sigma_n"
      ),

    alpha_hat_null =
      get_coef(
        coef_null_i,
        "alpha"
      ),

    alpha_s_hat_null =
      get_coef(
        coef_null_i,
        "alpha_s"
      ),


    # -------------------------------------------------------------------------
    # Reproducibility
    # -------------------------------------------------------------------------

    seed =
      seed_i,

    test_seconds =
      elapsed_i,

    status =
      if (
        is.null(WWW$status)
      ) {
        NA_character_
      } else {
        as.character(
          WWW$status
        )
      },

    error_message =
      if (
        is.null(WWW$error.message)
      ) {
        NA_character_
      } else {
        as.character(
          WWW$error.message
        )
      }
  )
}


# =============================================================================
# 12. Parallel execution
# =============================================================================

cat(
  "\nTotal RI-event tasks:",
  nrow(tasks),
  "\n"
)

cat(
  "Workers:",
  N_WORKERS,
  "\n"
)


cl <- parallel::makeCluster(
  N_WORKERS,
  type = "PSOCK",
  outfile = ""
)


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


# Export all objects/functions required by nested_test2.r.
parallel::clusterExport(
  cl,
  varlist = c(

    "rmats_all",
    "thr",

    "N_BOOT",
    "SCALING_A",
    "LAMBDA_TIME",
    "LAMBDA_DIAG",
    "REL_FLOOR",
    "TRUNCATE_NONNEGATIVE_BOOT",
    "MAX_FAILURE_RATE",
    "T_STAR_DEFAULT",

    "filter_event",
    "make_real_seed",
    "run_real_event",

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
  ),
  envir = .GlobalEnv
)


# Split each row into an independent task.
task_list <- split(
  tasks,
  seq_len(
    nrow(tasks)
  )
)


start_all <- Sys.time()


res_list <- parallel::parLapplyLB(
  cl,
  task_list,
  fun = run_real_event,
  N_boot = N_BOOT
)


parallel::stopCluster(
  cl
)


cat(
  "\nElapsed:",
  difftime(
    Sys.time(),
    start_all,
    units = "hours"
  ),
  "hours\n"
)


# =============================================================================
# 13. Combine results
# =============================================================================

results_all <- rbindlist(
  res_list,
  use.names = TRUE,
  fill = TRUE
)


# Keep both tested and filtered events.
results_tested <- results_all[
  filter_pass == TRUE
]


# =============================================================================
# 14. Test status
# =============================================================================

results_tested[
  ,
  valid_test :=
    status == "ok" &
    is.finite(
      p.value
    )
]


cat(
  "\n============================================================\n"
)

cat(
  "TEST STATUS\n"
)

cat(
  "============================================================\n"
)


print(
  results_tested[
    ,
    .(
      N_tested = .N,
      N_valid = sum(valid_test),
      Valid_fraction = mean(valid_test),
      Mean_boot_failure = mean(
        bootstrap_failure_rate,
        na.rm = TRUE
      ),
      Median_condition = median(
        condition_number,
        na.rm = TRUE
      )
    ),
    by = .(
      dataset,
      perturbation_type
    )
  ]
)


# =============================================================================
# 15. Multiple testing
# =============================================================================

results_tested[
  valid_test == TRUE,
  q.value :=
    p.adjust(
      p.value,
      method = "BH"
    ),
  by = dataset
]


# =============================================================================
# 16. Derived quantities
# =============================================================================

eps <- 1e-10


results_tested[
  ,
  DeltaRSS :=
    pmax(
      RSS0 - RSS1,
      0
    )
]


results_tested[
  ,
  IR :=
    fifelse(
      is.finite(RSS0) &
        RSS0 > 0,
      DeltaRSS / RSS0,
      NA_real_
    )
]


results_tested[
  ,
  logRSSratio :=
    fifelse(
      is.finite(RSS1) &
        RSS1 > 0,
      log(
        RSS0 / RSS1
      ),
      NA_real_
    )
]


results_tested[
  ,
  neglogp :=
    -log10(
      pmax(
        p.value,
        eps
      )
    )
]


results_tested[
  ,
  neglogq :=
    -log10(
      pmax(
        q.value,
        eps
      )
    )
]


cap <- 6


results_tested[
  ,
  neglogp_cap :=
    pmin(
      neglogp,
      cap
    )
]


results_tested[
  ,
  neglogq_cap :=
    pmin(
      neglogq,
      cap
    )
]


# =============================================================================
# 17. Effect-size normalization / exploratory ranking
# =============================================================================

results_tested[
  ,
  s0 :=
    median(
      abs(Sigma),
      na.rm = TRUE
    ),
  by = dataset
]


results_tested[
  !is.finite(s0) |
    s0 <= 0,
  s0 := 1
]


results_tested[
  ,
  sig_log :=
    log1p(
      abs(Sigma) /
        s0
    )
]


# Keep the historical ranking definitions for comparison/reproducibility.
#
# IMPORTANT:
# q-value is the inferential multiple-testing quantity.
# Ranking scores are exploratory prioritization tools only.

results_tested[
  ,
  score_sigma_p :=
    Sigma *
    neglogp_cap
]


results_tested[
  ,
  score_sigma_IR_p :=
    Sigma *
    IR *
    neglogp_cap
]


results_tested[
  ,
  score_siglog_IR_p :=
    sig_log *
    IR *
    neglogp_cap
]


# FDR-aware versions, useful for exploratory ranking.
results_tested[
  ,
  score_sigma_q :=
    Sigma *
    neglogq_cap
]


results_tested[
  ,
  score_sigma_IR_q :=
    Sigma *
    IR *
    neglogq_cap
]


results_tested[
  ,
  score_siglog_IR_q :=
    sig_log *
    IR *
    neglogq_cap
]


# =============================================================================
# 18. Characteristic kinetic time scales
# =============================================================================

# Assumes the input time variable has been standardized to minutes.
#
# DO NOT interpret this column as minutes if the input datasets have not yet
# been converted to a common minute scale.

results_tested[
  ,
  sigma_c_half_time :=
    fifelse(
      is.finite(Sigma) &
        Sigma > 0,
      log(2) / Sigma,
      NA_real_
    )
]


# =============================================================================
# 19. Dataset-level discovery summaries
# =============================================================================

discovery_summary <- results_tested[
  valid_test == TRUE,
  .(

    N_tested =
      .N,

    N_unique_genes =
      uniqueN(
        ensembl
      ),

    N_p005 =
      sum(
        p.value < 0.05,
        na.rm = TRUE
      ),

    Fraction_p005 =
      mean(
        p.value < 0.05,
        na.rm = TRUE
      ),

    N_q020 =
      sum(
        q.value < 0.20,
        na.rm = TRUE
      ),

    N_q010 =
      sum(
        q.value < 0.10,
        na.rm = TRUE
      ),

    N_q005 =
      sum(
        q.value < 0.05,
        na.rm = TRUE
      ),

    Min_p =
      min(
        p.value,
        na.rm = TRUE
      ),

    Min_q =
      min(
        q.value,
        na.rm = TRUE
      ),

    Median_sigma =
      median(
        Sigma,
        na.rm = TRUE
      ),

    Median_IR =
      median(
        IR,
        na.rm = TRUE
      ),

    Median_condition =
      median(
        condition_number,
        na.rm = TRUE
      ),

    Median_atom_zero =
      median(
        atom_zero,
        na.rm = TRUE
      )
  ),

  by = .(
    dataset,
    perturbation_type,
    pseudo_method
  )
]


cat(
  "\n============================================================\n"
)

cat(
  "DISCOVERY SUMMARY\n"
)

cat(
  "============================================================\n"
)

print(
  discovery_summary
)


# =============================================================================
# 20. Filtering/QC summary
# =============================================================================

qc_summary <- results_all[
  ,
  .(
    N_events_input = .N,
    N_pass = sum(
      filter_pass,
      na.rm = TRUE
    ),
    N_filtered = sum(
      !filter_pass,
      na.rm = TRUE
    ),
    fraction_pass = mean(
      filter_pass,
      na.rm = TRUE
    )
  ),
  by = .(
    dataset,
    perturbation_type
  )
]


filter_reason_summary <- results_all[
  filter_pass == FALSE,
  .N,
  by = .(
    dataset,
    filter_reason
  )
]


# =============================================================================
# 21. Ranked tables
# =============================================================================

# Complete table: all valid tests.
supp_all <- results_tested[
  valid_test == TRUE
]


setorder(
  supp_all,
  dataset,
  q.value,
  -score_sigma_IR_q
)


# ---------------------------------------------------------------------------
# FDR-controlled tables
# ---------------------------------------------------------------------------

fdr05 <- supp_all[
  q.value < 0.05
]


fdr10 <- supp_all[
  q.value < 0.10
]


fdr20 <- supp_all[
  q.value < 0.20
]


# ---------------------------------------------------------------------------
# Exploratory ranking.
#
# No claim of "high-confidence" is attached to this table.
# ---------------------------------------------------------------------------

exploratory_rank <- supp_all[
  order(
    dataset,
    -score_sigma_IR_q
  )
]


# Historical nominal criterion retained only for comparison with the original
# manuscript, NOT as the revised discovery definition.

historical_nominal <- supp_all[
  p.value < 0.05 &
    IR > 0.10
]


# =============================================================================
# 22. Save R objects
# =============================================================================

save(
  results_all,
  results_tested,
  discovery_summary,
  qc_summary,
  filter_reason_summary,
  supp_all,
  fdr05,
  fdr10,
  fdr20,
  exploratory_rank,
  historical_nominal,
  file = "results_realdatasets_revision.rdata"
)


# =============================================================================
# 23. TSV outputs
# =============================================================================

fwrite(
  results_tested,
  "results_realdatasets_revision.tsv",
  sep = "\t"
)


fwrite(
  discovery_summary,
  "real_dataset_discovery_summary.tsv",
  sep = "\t"
)


fwrite(
  qc_summary,
  "real_dataset_qc_summary.tsv",
  sep = "\t"
)


fwrite(
  filter_reason_summary,
  "real_dataset_filter_reasons.tsv",
  sep = "\t"
)


# =============================================================================
# 24. Excel workbook
# =============================================================================

wb <- createWorkbook()


addWorksheet(
  wb,
  "All_valid_tests"
)

writeData(
  wb,
  "All_valid_tests",
  supp_all
)


addWorksheet(
  wb,
  "FDR_005"
)

writeData(
  wb,
  "FDR_005",
  fdr05
)


addWorksheet(
  wb,
  "FDR_010"
)

writeData(
  wb,
  "FDR_010",
  fdr10
)


addWorksheet(
  wb,
  "FDR_020"
)

writeData(
  wb,
  "FDR_020",
  fdr20
)


addWorksheet(
  wb,
  "Exploratory_rank"
)

writeData(
  wb,
  "Exploratory_rank",
  exploratory_rank
)


addWorksheet(
  wb,
  "Original_nominal_rule"
)

writeData(
  wb,
  "Original_nominal_rule",
  historical_nominal
)


addWorksheet(
  wb,
  "Discovery_summary"
)

writeData(
  wb,
  "Discovery_summary",
  discovery_summary
)


addWorksheet(
  wb,
  "QC_summary"
)

writeData(
  wb,
  "QC_summary",
  qc_summary
)


saveWorkbook(
  wb,
  "results_realdatasets_revision.xlsx",
  overwrite = TRUE
)


# =============================================================================
# 25. Session info
# =============================================================================

sink(
  "real_dataset_sessionInfo.txt"
)

print(
  sessionInfo()
)

sink()


# =============================================================================
# 26. Final console summary
# =============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "REAL-DATA ANALYSIS COMPLETE\n"
)

cat(
  "============================================================\n"
)


print(
  discovery_summary
)


cat(
  "\nFiles written:\n",
  "  results_realdatasets_revision.rdata\n",
  "  results_realdatasets_revision.tsv\n",
  "  results_realdatasets_revision.xlsx\n",
  "  real_dataset_discovery_summary.tsv\n",
  "  real_dataset_qc_summary.tsv\n",
  "  real_dataset_filter_reasons.tsv\n",
  "  real_dataset_sessionInfo.txt\n",
  sep = ""
)