# =============================================================================
# Title:
# High-Resolution Bootstrap Analysis of the mESC Pharmacological-Shutoff Dataset
#
# Description:
#   This script performs the definitive high-resolution real-data analysis of
#   retained-intron (RI) events in the mouse embryonic stem-cell dataset
#   GSE256335.
#
#   The dataset represents the true pharmacological transcriptional-shutoff
#   experiment used in the real-data component of the post-export RNA kinetics
#   study.
#
#   This analysis is restricted to mESC because the initial real-data analysis
#   showed FDR-controlled evidence only in this dataset. The pseudo-shutoff
#   datasets are retained as exploratory analyses and are not recomputed at
#   this higher bootstrap resolution.
#
# Statistical procedure:
#   For each RI event passing the same preprocessing and quality-control
#   criteria used in the complete real-data analysis, the script compares:
#
#       H0: sigma_c = 0
#
#   against the constrained full model:
#
#       H1: sigma_c >= 0
#
#   using the bootstrap-calibrated weighted NNLS model-comparison procedure.
#
#   The bootstrap:
#     - operates on destructive replicate-level observations;
#     - generates replicate-level datasets under H0;
#     - reconstructs A*, b*, and the associated covariance structure;
#     - uses the complete interval-difference covariance;
#     - accounts for uncertainty in both A and b;
#     - uses non-negative least squares;
#     - uses the add-one bootstrap p-value:
#
#           p = (1 + sum(T* >= Tobs)) / (B + 1).
#
# Bootstrap resolution:
#
#       B = 19,999
#
#   giving:
#
#       p_min = 1 / 20,000 = 5e-5.
#
# Restart strategy:
#   Each RI event is checkpointed independently immediately after completion:
#
#       mesc_20k_checkpoints/event_<hash>.rds
#
#   Workers return only lightweight status records to the PSOCK master.
#   Therefore:
#
#     - completed events survive a PSOCK communication failure;
#     - rerunning this script automatically skips valid checkpoints;
#     - only incomplete or corrupted events are recomputed.
#
# Main outputs:
#
#   mesc_20k_results.rdata
#   mesc_20k_results.tsv
#   mesc_20k_FDR05.tsv
#   mesc_20k_FDR10.tsv
#   mesc_20k_summary.tsv
#   mesc_20k_progress.tsv
#   mesc_20k_sessionInfo.txt
#
# Scientific interpretation:
#   sigma_c is treated as a phenomenological post-export conversion parameter.
#   FDR-supported events represent kinetic patterns consistent with an
#   additional post-export conversion component. They should not be
#   interpreted as direct experimental evidence of cytoplasmic splicing or of
#   a uniquely determined molecular mechanism.
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
#   Definitive mESC high-resolution bootstrap analysis, 2026.
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/real_datasets")

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


# =============================================================================
# 1. Main settings
# =============================================================================

DATASET_NAME <- "mESC | GSE256335"

N_BOOT <- 19999L

# Use only CPUs that are really free while the synthetic benchmark is running.
# Change this according to available resources.
N_WORKERS <- 80L

# Number of events passed to one PSOCK call.
# Checkpoints are nevertheless PER EVENT.
TASK_BATCH_SIZE <- N_WORKERS

PSOCK_TIMEOUT_SECONDS <- 12L * 60L * 60L

SCALING_A <- TRUE

LAMBDA_TIME <- 0.5
LAMBDA_DIAG <- 0.1

REL_FLOOR <- 1e-8

TRUNCATE_NONNEGATIVE_BOOT <- FALSE

MAX_FAILURE_RATE <- 0.05


CHECKPOINT_DIR <- "mesc_20k_checkpoints"

PROGRESS_FILE <- "mesc_20k_progress.tsv"


# =============================================================================
# 2. Load mESC data
# =============================================================================

cat("\nLoading mESC dataset...\n")

load("GSE256335/rmats_dt_ECS.rdata")

stopifnot(
  exists("rmats_dt")
)

dt_mesc <- copy(
  as.data.table(
    rmats_dt
  )
)

dt_mesc[
  ,
  dataset :=
    DATASET_NAME
]


cat(
  "Rows:",
  nrow(dt_mesc),
  "\n"
)

cat(
  "RI events:",
  uniqueN(dt_mesc$event),
  "\n"
)

cat(
  "Genes:",
  uniqueN(dt_mesc$ensembl),
  "\n"
)

cat(
  "Times:",
  paste(
    sort(unique(dt_mesc$time)),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Replicates:",
  uniqueN(dt_mesc$replicate),
  "\n"
)


# =============================================================================
# 3. Event-level QC metrics
# =============================================================================

cat("\nComputing event-level QC...\n")

evt_qc <- dt_mesc[
  ,
  .(

    cov_total =
      sum(
        N +
          N_s +
          C +
          C_s,
        na.rm = TRUE
      ),

    cov_cyt =
      sum(
        C +
          C_s,
        na.rm = TRUE
      ),

    rng_C =
      max(
        C,
        na.rm = TRUE
      ) -
      min(
        C,
        na.rm = TRUE
      ),

    rng_Cs =
      max(
        C_s,
        na.rm = TRUE
      ) -
      min(
        C_s,
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

    frac_nzC =
      mean(
        C > 0,
        na.rm = TRUE
      ),

    frac_nzCs =
      mean(
        C_s > 0,
        na.rm = TRUE
      )

  ),
  by = .(
    event
  )
]


# =============================================================================
# 4. Dataset-specific QC thresholds
# =============================================================================

thr <- evt_qc[
  ,
  .(

    cov_total_min =
      quantile(
        cov_total,
        0.01,
        na.rm = TRUE
      ),

    cov_cyt_min =
      quantile(
        cov_cyt,
        0.01,
        na.rm = TRUE
      ),

    rng_C_min =
      quantile(
        rng_C,
        0.01,
        na.rm = TRUE
      ),

    rng_Cs_min =
      quantile(
        rng_Cs,
        0.01,
        na.rm = TRUE
      )
  )
]


print(
  thr
)


# =============================================================================
# 5. Apply exactly the same event filters as the primary real-data analysis
# =============================================================================

MIN_POS <- 3L

MIN_FRAC <- 0.30

MIN_POS_NUC <- 3L

MAX_ZERO_FRAC_C <- 0.70

MAX_ZERO_FRAC_CS <- 0.80


event_filter <- dt_mesc[
  ,
  {

    d <- .SD


    bad <- (

      all(
        d$N == 0
      ) ||

      all(
        d$N_s == 0
      ) ||

      all(
        d$C == 0
      ) ||

      all(
        d$C_s == 0
      ) ||

      sum(
        d$C > 0,
        na.rm = TRUE
      ) < MIN_POS ||

      sum(
        d$C_s > 0,
        na.rm = TRUE
      ) < MIN_POS ||

      mean(
        d$C > 0,
        na.rm = TRUE
      ) < MIN_FRAC ||

      sum(
        d$N > 0,
        na.rm = TRUE
      ) < MIN_POS_NUC ||

      sum(
        d$N_s > 0,
        na.rm = TRUE
      ) < MIN_POS_NUC ||

      sum(
        d$N +
          d$N_s +
          d$C +
          d$C_s,
        na.rm = TRUE
      ) < thr$cov_total_min ||

      sum(
        d$C +
          d$C_s,
        na.rm = TRUE
      ) < thr$cov_cyt_min ||

      (
        max(
          d$C,
          na.rm = TRUE
        ) -
          min(
            d$C,
            na.rm = TRUE
          )
      ) < thr$rng_C_min ||

      (
        max(
          d$C_s,
          na.rm = TRUE
        ) -
          min(
            d$C_s,
            na.rm = TRUE
          )
      ) < thr$rng_Cs_min ||

      mean(
        d$C == 0,
        na.rm = TRUE
      ) > MAX_ZERO_FRAC_C ||

      mean(
        d$C_s == 0,
        na.rm = TRUE
      ) > MAX_ZERO_FRAC_CS
    )


    .(
      filter_pass =
        !bad,

      n_timepoints =
        uniqueN(
          d$time
        ),

      n_replicates =
        uniqueN(
          d$replicate
        ),

      gene_symbol =
        gene_symbol[1],

      gene_biotype =
        gene_biotype[1],

      description =
        description[1],

      ensembl =
        ensembl[1]
    )
  },

  by = event
]


cat(
  "\nQC summary:\n"
)

print(
  event_filter[
    ,
    .N,
    by = filter_pass
  ]
)


events_to_test <- event_filter[
  filter_pass == TRUE,
  event
]


cat(
  "\nEvents to test:",
  length(events_to_test),
  "\n"
)


# =============================================================================
# 6. Deterministic seed per event
# =============================================================================

stable_event_seed <- function(
  event,
  offset = 900000000L
) {

  ints <- utf8ToInt(
    enc2utf8(
      event
    )
  )


  h <- 104729


  for (ii in ints) {

    h <- (
      h * 1000003 +
        ii
    ) %% 1000000000
  }


  seed <- as.integer(
    (
      h +
        offset
    ) %%
      2000000000L
  )


  if (
    !is.finite(seed) ||
      seed <= 0
  ) {

    seed <- 1L
  }


  seed
}


# =============================================================================
# 7. Safe helper functions
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
# 8. Event checkpoint filename
# =============================================================================

event_hash <- function(
  event
) {

  sprintf(
    "%010d",
    stable_event_seed(
      event,
      offset = 0L
    )
  )
}


event_checkpoint_file <- function(
  event
) {

  file.path(
    CHECKPOINT_DIR,
    paste0(
      "event_",
      event_hash(event),
      ".rds"
    )
  )
}


# =============================================================================
# 9. Checkpoint directory
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
# 10. Run one mESC event
# =============================================================================

run_mesc_event <- function(
  event_i,
  N_boot
) {

  d <- copy(
    dt_mesc[
      event == event_i
    ]
  )


  if (nrow(d) == 0L) {

    stop(
      paste(
        "No rows found for event:",
        event_i
      )
    )
  }


  n_rep <- uniqueN(
    d$replicate
  )


  seed_i <- stable_event_seed(
    event_i
  )


  t0 <- Sys.time()


  WW <- test_sigma_nested(

    tsampled_data =
      d,

    scaling_A =
      SCALING_A,

    t_star =
      0,

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
  )


  elapsed_seconds <- as.numeric(
    difftime(
      Sys.time(),
      t0,
      units = "secs"
    )
  )


  coef_full <- WW$coef_full

  coef_null <- WW$coef_null


  data.table(

    dataset =
      DATASET_NAME,

    perturbation_type =
      "TRUE_SHUTOFF",

    event =
      event_i,

    ensembl =
      d$ensembl[1],

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
      n_rep,

    N_boot =
      N_boot,

    seed =
      seed_i,


    # -------------------------------------------------------------------------
    # Inference
    # -------------------------------------------------------------------------

    p.value =
      safe_value(
        WW,
        "p.value"
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

    T.obs =
      safe_value(
        WW,
        "T.obs"
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


    # -------------------------------------------------------------------------
    # Boundary / bootstrap diagnostics
    # -------------------------------------------------------------------------

    atom_zero =
      safe_value(
        WW,
        "atom.zero"
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

    bootstrap_condition_median =
      safe_value(
        WW,
        "bootstrap.condition.median"
      ),

    bootstrap_condition_q95 =
      safe_value(
        WW,
        "bootstrap.condition.q95"
      ),

    bootstrap_condition_max =
      safe_value(
        WW,
        "bootstrap.condition.max"
      ),

    bootstrap_rank_deficient_fraction =
      safe_value(
        WW,
        "bootstrap.rank.deficient.fraction"
      ),


    # -------------------------------------------------------------------------
    # Practical-identifiability diagnostics
    # -------------------------------------------------------------------------

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


    # -------------------------------------------------------------------------
    # Full model
    # -------------------------------------------------------------------------

    R_hat =
      safe_coef(
        coef_full,
        "R"
      ),

    tau_hat =
      safe_coef(
        coef_full,
        "tau"
      ),

    tau_s_hat =
      safe_coef(
        coef_full,
        "tau_s"
      ),

    sigma_c_hat =
      safe_coef(
        coef_full,
        "sigma_c"
      ),

    sigma_n_hat =
      safe_coef(
        coef_full,
        "sigma_n"
      ),

    alpha_hat =
      safe_coef(
        coef_full,
        "alpha"
      ),

    alpha_s_hat =
      safe_coef(
        coef_full,
        "alpha_s"
      ),


    # -------------------------------------------------------------------------
    # Null model
    # -------------------------------------------------------------------------

    R_hat_null =
      safe_coef(
        coef_null,
        "R"
      ),

    tau_hat_null =
      safe_coef(
        coef_null,
        "tau"
      ),

    tau_s_hat_null =
      safe_coef(
        coef_null,
        "tau_s"
      ),

    sigma_n_hat_null =
      safe_coef(
        coef_null,
        "sigma_n"
      ),

    alpha_hat_null =
      safe_coef(
        coef_null,
        "alpha"
      ),

    alpha_s_hat_null =
      safe_coef(
        coef_null,
        "alpha_s"
      ),


    test_seconds =
      elapsed_seconds,

    status =
      if (
        is.null(WW$status)
      ) {
        "ok"
      } else {
        as.character(
          WW$status
        )
      }
  )
}


# =============================================================================
# 11. Run event and checkpoint immediately
# =============================================================================

run_event_and_checkpoint <- function(
  event_i,
  N_boot
) {

  pid <- Sys.getpid()


  checkpoint_file <- event_checkpoint_file(
    event_i
  )


  # ---------------------------------------------------------------------------
  # Existing checkpoint
  # ---------------------------------------------------------------------------

  if (file.exists(checkpoint_file)) {

    old <- tryCatch(
      readRDS(
        checkpoint_file
      ),
      error = function(e) NULL
    )


    if (
      !is.null(old) &&
        data.table::is.data.table(old) &&
        nrow(old) == 1L &&
        identical(
          as.character(old$event),
          as.character(event_i)
        ) &&
        old$N_boot == N_boot
    ) {

      return(
        data.table(
          event =
            event_i,

          status =
            "already_checkpointed",

          pid =
            pid,

          elapsed_seconds =
            0,

          checkpoint =
            checkpoint_file
        )
      )
    }
  }


  # ---------------------------------------------------------------------------
  # Compute event
  # ---------------------------------------------------------------------------

  t0 <- Sys.time()


  result <- tryCatch(

    run_mesc_event(
      event_i =
        event_i,

      N_boot =
        N_boot
    ),

    error = function(e) {

      structure(
        list(
          error_message =
            conditionMessage(e)
        ),
        class =
          "mesc_event_error"
      )
    }
  )


  if (
    inherits(
      result,
      "mesc_event_error"
    )
  ) {

    return(
      data.table(

        event =
          event_i,

        status =
          "event_error",

        pid =
          pid,

        elapsed_seconds =
          as.numeric(
            difftime(
              Sys.time(),
              t0,
              units = "secs"
            )
          ),

        checkpoint =
          NA_character_,

        error_message =
          result$error_message
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Atomic checkpoint write
  # ---------------------------------------------------------------------------

  tmp_file <- paste0(
    checkpoint_file,
    ".tmp.",
    pid
  )


  saveRDS(
    result,
    file =
      tmp_file,

    compress =
      "gzip"
  )


  rename_ok <- file.rename(
    tmp_file,
    checkpoint_file
  )


  if (!rename_ok) {

    if (file.exists(tmp_file)) {

      unlink(
        tmp_file
      )
    }


    return(
      data.table(

        event =
          event_i,

        status =
          "checkpoint_error",

        pid =
          pid,

        elapsed_seconds =
          as.numeric(
            difftime(
              Sys.time(),
              t0,
              units = "secs"
            )
          ),

        checkpoint =
          NA_character_,

        error_message =
          "Checkpoint rename failed."
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
      "[PID %d] CHECKPOINT SAVED %s | p=%s | sigma_c=%s | %.1f min",
      pid,
      result$gene_symbol,
      format(
        result$p.value,
        digits = 4
      ),
      format(
        result$Sigma,
        digits = 4
      ),
      elapsed / 60
    )
  )


  data.table(

    event =
      event_i,

    status =
      "checkpoint_saved",

    pid =
      pid,

    elapsed_seconds =
      elapsed,

    checkpoint =
      checkpoint_file,

    error_message =
      NA_character_
  )
}


# =============================================================================
# 12. Detect valid existing checkpoints
# =============================================================================

valid_checkpoint_events <- character(0)


for (event_i in events_to_test) {


  ff <- event_checkpoint_file(
    event_i
  )


  if (!file.exists(ff)) {
    next
  }


  x <- tryCatch(
    readRDS(
      ff
    ),
    error = function(e) NULL
  )


  if (
    !is.null(x) &&
      data.table::is.data.table(x) &&
      nrow(x) == 1L &&
      identical(
        as.character(x$event),
        as.character(event_i)
      ) &&
      x$N_boot == N_BOOT
  ) {

    valid_checkpoint_events <- c(
      valid_checkpoint_events,
      event_i
    )
  }
}


valid_checkpoint_events <- unique(
  valid_checkpoint_events
)


events_remaining <- setdiff(
  events_to_test,
  valid_checkpoint_events
)


cat(
  "\n============================================================\n"
)

cat(
  "mESC 20K RESTART STATUS\n"
)

cat(
  "============================================================\n"
)

cat(
  "Events passing QC:      ",
  length(events_to_test),
  "\n"
)

cat(
  "Existing checkpoints:   ",
  length(valid_checkpoint_events),
  "\n"
)

cat(
  "Events still to test:   ",
  length(events_remaining),
  "\n"
)

cat(
  "Bootstrap replicates:   ",
  N_BOOT,
  "\n"
)

cat(
  "Minimum possible p:     ",
  1 / (
    N_BOOT +
      1
  ),
  "\n"
)

cat(
  "============================================================\n"
)


# =============================================================================
# 13. PSOCK execution
# =============================================================================

overall_start <- Sys.time()


if (length(events_remaining) > 0L) {


  workers_use <- min(
    N_WORKERS,
    length(events_remaining)
  )


  cat(
    "\nStarting",
    workers_use,
    "PSOCK workers...\n"
  )


  cl <- parallel::makePSOCKcluster(

    workers_use,

    outfile =
      "",

    timeout =
      PSOCK_TIMEOUT_SECONDS,

    setup_timeout =
      120
  )


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

        "DATASET_NAME",

        "dt_mesc",

        "N_BOOT",

        "SCALING_A",

        "LAMBDA_TIME",

        "LAMBDA_DIAG",

        "REL_FLOOR",

        "TRUNCATE_NONNEGATIVE_BOOT",

        "MAX_FAILURE_RATE",

        "CHECKPOINT_DIR",

        "stable_event_seed",

        "event_hash",

        "event_checkpoint_file",

        "safe_value",

        "safe_coef",

        "run_mesc_event",

        "run_event_and_checkpoint",

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
        varlist =
          objects_to_export,
        envir =
          .GlobalEnv
      )


      task_batches <- split(
        events_remaining,
        ceiling(
          seq_along(events_remaining) /
            TASK_BATCH_SIZE
        )
      )


      for (bb in seq_along(task_batches)) {


        events_bb <- task_batches[[bb]]


        cat(
          "\n============================================================\n"
        )

        cat(
          "Running event batch",
          bb,
          "/",
          length(task_batches),
          "|",
          length(events_bb),
          "events\n"
        )

        cat(
          "============================================================\n"
        )


        batch_start <- Sys.time()


        status_bb <- parallel::parLapplyLB(

          cl,

          events_bb,

          fun =
            run_event_and_checkpoint,

          N_boot =
            N_BOOT
        )


        status_bb <- rbindlist(
          status_bb,
          fill = TRUE
        )


        status_bb[
          ,
          `:=`(

            task_batch =
              bb,

            timestamp =
              as.character(
                Sys.time()
              )
          )
        ]


        if (!file.exists(PROGRESS_FILE)) {

          fwrite(
            status_bb,
            PROGRESS_FILE,
            sep = "\t"
          )

        } else {

          fwrite(
            status_bb,
            PROGRESS_FILE,
            sep = "\t",
            append = TRUE,
            col.names = FALSE
          )
        }


        cat(
          "\nBatch status:\n"
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
            as.numeric(
              difftime(
                Sys.time(),
                batch_start,
                units = "mins"
              )
            ),
            1
          ),
          "min\n"
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
}


# =============================================================================
# 14. Load all completed event checkpoints
# =============================================================================

cat(
  "\nLoading final mESC event checkpoints...\n"
)


res_list <- vector(
  "list",
  length(events_to_test)
)


missing_events <- character(0)


for (ii in seq_along(events_to_test)) {


  event_i <- events_to_test[ii]


  ff <- event_checkpoint_file(
    event_i
  )


  if (!file.exists(ff)) {

    missing_events <- c(
      missing_events,
      event_i
    )

    next
  }


  x <- tryCatch(
    readRDS(
      ff
    ),
    error = function(e) NULL
  )


  if (
    is.null(x) ||
      !data.table::is.data.table(x) ||
      nrow(x) != 1L ||
      x$N_boot != N_BOOT
  ) {

    missing_events <- c(
      missing_events,
      event_i
    )

    next
  }


  res_list[[ii]] <- x
}


if (length(missing_events) > 0L) {

  stop(
    paste0(
      "\nAnalysis incomplete. ",
      length(missing_events),
      " event(s) are missing.\n",
      "Restart the same script; valid event checkpoints will be reused."
    )
  )
}


results_mesc_20k <- rbindlist(
  res_list,
  use.names = TRUE,
  fill = TRUE
)


# =============================================================================
# 15. Multiple-testing correction
# =============================================================================

results_mesc_20k[
  ,
  q.value :=
    p.adjust(
      p.value,
      method = "BH"
    )
]


# =============================================================================
# 16. Derived model-comparison quantities
# =============================================================================

results_mesc_20k[
  ,
  DeltaRSS :=
    pmax(
      RSS0 -
        RSS1,
      0
    )
]


results_mesc_20k[
  ,
  IR :=
    fifelse(
      RSS0 > 0,
      DeltaRSS /
        RSS0,
      0
    )
]


results_mesc_20k[
  ,
  logRSSratio :=
    fifelse(
      RSS1 > 0,
      log(
        RSS0 /
          RSS1
      ),
      NA_real_
    )
]


# =============================================================================
# 17. Conversion timescale
# =============================================================================

results_mesc_20k[
  ,
  sigma_c_half_time :=
    fifelse(
      is.finite(
        Sigma
      ) &
        Sigma > 0,

      log(2) /
        Sigma,

      NA_real_
    )
]


# =============================================================================
# 18. Boundary diagnostics
# =============================================================================

results_mesc_20k[
  ,
  n_boundary_nuisance :=
    rowSums(
      cbind(

        tau_hat <= 1e-12,

        tau_s_hat <= 1e-12,

        sigma_n_hat <= 1e-12,

        alpha_hat <= 1e-12,

        alpha_s_hat <= 1e-12
      ),

      na.rm = TRUE
    )
]


# =============================================================================
# 19. Add experimental/QC summaries useful for future validation
# =============================================================================

event_validation_qc <- dt_mesc[
  ,
  .(

    mean_total =
      mean(
        N +
          N_s +
          C +
          C_s,
        na.rm = TRUE
      ),

    mean_cyt =
      mean(
        C +
          C_s,
        na.rm = TRUE
      ),

    range_N =
      diff(
        range(
          N,
          na.rm = TRUE
        )
      ),

    range_Ns =
      diff(
        range(
          N_s,
          na.rm = TRUE
        )
      ),

    range_C =
      diff(
        range(
          C,
          na.rm = TRUE
        )
      ),

    range_Cs =
      diff(
        range(
          C_s,
          na.rm = TRUE
        )
      ),

    frac_nonzero_N =
      mean(
        N > 0,
        na.rm = TRUE
      ),

    frac_nonzero_Ns =
      mean(
        N_s > 0,
        na.rm = TRUE
      ),

    frac_nonzero_C =
      mean(
        C > 0,
        na.rm = TRUE
      ),

    frac_nonzero_Cs =
      mean(
        C_s > 0,
        na.rm = TRUE
      )

  ),
  by = event
]


results_mesc_20k <- merge(
  results_mesc_20k,
  event_validation_qc,
  by = "event",
  all.x = TRUE
)


# =============================================================================
# 20. Replicate consistency
# =============================================================================

long_qc <- melt(

  dt_mesc,

  id.vars = c(
    "event",
    "time",
    "replicate"
  ),

  measure.vars = c(
    "N",
    "N_s",
    "C",
    "C_s"
  ),

  variable.name =
    "state",

  value.name =
    "value"
)


rep_qc <- long_qc[
  ,
  .(

    state_mean =
      mean(
        value,
        na.rm = TRUE
      ),

    state_sd =
      sd(
        value,
        na.rm = TRUE
      )

  ),
  by = .(
    event,
    time,
    state
  )
]


rep_qc[
  ,
  cv :=
    fifelse(
      state_mean > 0,
      state_sd /
        state_mean,
      NA_real_
    )
]


rep_event_qc <- rep_qc[
  ,
  .(

    median_replicate_cv =
      median(
        cv,
        na.rm = TRUE
      ),

    q90_replicate_cv =
      as.numeric(
        quantile(
          cv,
          0.90,
          na.rm = TRUE,
          names = FALSE
        )
      )

  ),
  by = event
]


results_mesc_20k <- merge(
  results_mesc_20k,
  rep_event_qc,
  by = "event",
  all.x = TRUE
)


# =============================================================================
# 21. Validation-oriented ranking
#
# IMPORTANT:
# This is an experimental-prioritization score, NOT a statistical test.
# Inferential conclusions remain based on BH q-values.
# =============================================================================

rank01 <- function(
  x,
  higher_better = TRUE
) {

  z <- frank(
    x,
    ties.method = "average",
    na.last = "keep"
  )


  if (
    max(
      z,
      na.rm = TRUE
    ) > 0
  ) {

    z <- z /
      max(
        z,
        na.rm = TRUE
      )
  }


  if (!higher_better) {

    z <- 1 -
      z
  }


  z
}


results_mesc_20k[
  ,
  validation_score := {

    score_q <-
      rank01(
        -log10(
          pmax(
            q.value,
            1e-12
          )
        ),
        TRUE
      )

    score_IR <-
      rank01(
        IR,
        TRUE
      )

    score_sigma <-
      rank01(
        log1p(
          pmax(
            Sigma,
            0
          )
        ),
        TRUE
      )

    score_boundary <-
      rank01(
        n_boundary_nuisance,
        FALSE
      )

    score_abundance <-
      rank01(
        log1p(
          mean_cyt
        ),
        TRUE
      )

    score_dynamic <-
      rank01(
        log1p(
          range_C +
            range_Cs
        ),
        TRUE
      )

    score_replicates <-
      rank01(
        median_replicate_cv,
        FALSE
      )


    0.25 *
      score_q +

      0.20 *
      score_IR +

      0.10 *
      score_sigma +

      0.15 *
      score_boundary +

      0.10 *
      score_abundance +

      0.10 *
      score_dynamic +

      0.10 *
      score_replicates
  }
]


# =============================================================================
# 22. Final ordering
# =============================================================================

setorder(
  results_mesc_20k,
  q.value,
  -validation_score,
  -IR
)


# =============================================================================
# 23. Summary
# =============================================================================

mesc_summary <- results_mesc_20k[
  ,
  .(

    N_tested =
      .N,

    p_lt_005 =
      sum(
        p.value < 0.05,
        na.rm = TRUE
      ),

    q_lt_020 =
      sum(
        q.value < 0.20,
        na.rm = TRUE
      ),

    q_lt_010 =
      sum(
        q.value < 0.10,
        na.rm = TRUE
      ),

    q_lt_005 =
      sum(
        q.value < 0.05,
        na.rm = TRUE
      ),

    sigma_zero =
      sum(
        Sigma <= 1e-12,
        na.rm = TRUE
      ),

    frac_sigma_zero =
      mean(
        Sigma <= 1e-12,
        na.rm = TRUE
      ),

    p_one =
      sum(
        p.value == 1,
        na.rm = TRUE
      ),

    frac_p_one =
      mean(
        p.value == 1,
        na.rm = TRUE
      ),

    median_atom_zero =
      median(
        atom_zero,
        na.rm = TRUE
      ),

    median_boot_failure =
      median(
        bootstrap_failure_rate,
        na.rm = TRUE
      ),

    median_test_seconds =
      median(
        test_seconds,
        na.rm = TRUE
      )
  )
]


# =============================================================================
# 24. FDR tables
# =============================================================================

FDR05 <- results_mesc_20k[
  q.value < 0.05
]


FDR10 <- results_mesc_20k[
  q.value < 0.10
]


# =============================================================================
# 25. Save
# =============================================================================

save(
  results_mesc_20k,
  FDR05,
  FDR10,
  mesc_summary,
  file =
    "mesc_20k_results.rdata"
)


fwrite(
  results_mesc_20k,
  "mesc_20k_results.tsv",
  sep = "\t"
)


fwrite(
  FDR05,
  "mesc_20k_FDR05.tsv",
  sep = "\t"
)


fwrite(
  FDR10,
  "mesc_20k_FDR10.tsv",
  sep = "\t"
)


fwrite(
  mesc_summary,
  "mesc_20k_summary.tsv",
  sep = "\t"
)


# =============================================================================
# 26. Session information
# =============================================================================

sink(
  "mesc_20k_sessionInfo.txt"
)


cat(
  "mESC high-resolution analysis\n"
)

cat(
  "Date:\n"
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
  "\nMinimum bootstrap p-value:\n"
)

print(
  1 /
    (
      N_BOOT +
        1
    )
)


cat(
  "\nSummary:\n"
)

print(
  mesc_summary
)


cat(
  "\nSession information:\n"
)

print(
  sessionInfo()
)


sink()


# =============================================================================
# 27. Final console output
# =============================================================================

elapsed_hours <- as.numeric(
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
  "mESC 20K ANALYSIS COMPLETE\n"
)

cat(
  "============================================================\n"
)


print(
  mesc_summary
)


cat(
  "\nFDR < 0.05 events:\n"
)


print(
  FDR05[
    ,
    .(
      gene_symbol,
      event,
      p.value,
      q.value,
      Sigma,
      sigma_c_half_time,
      IR,
      n_boundary_nuisance,
      validation_score
    )
  ]
)


cat(
  "\nTotal runtime:",
  round(
    elapsed_hours,
    2
  ),
  "hours\n"
)


cat(
  "\nSaved:\n",
  "  mesc_20k_results.rdata\n",
  "  mesc_20k_results.tsv\n",
  "  mesc_20k_FDR05.tsv\n",
  "  mesc_20k_FDR10.tsv\n",
  "  mesc_20k_summary.tsv\n",
  "  mesc_20k_progress.tsv\n",
  "  mesc_20k_sessionInfo.txt\n",
  sep = ""
)


cat(
  "============================================================\n"
)