# =============================================================================
# Title:
# Residual-Transcription and Pseudo-Shutoff Misspecification Benchmark for the
# Revised Bootstrap-Calibrated NNLS Test
#
# Description:
#   This script evaluates robustness of the revised weighted NNLS test to
#   violations of the complete transcriptional-shutoff assumption.
#
#   Synthetic trajectories are generated directly from the ODE model with:
#
#       R_post = rho * R_pre
#
#   where rho is the residual-transcription fraction stored in:
#
#       Post_R_fraction
#
#   The analysis includes complete shutoff (rho=0) and all simulated
#   pseudo-shutoff conditions.
#
#   Crucially, the statistical test is deliberately fitted assuming complete
#   shutoff at T_star for every pseudo-shutoff trajectory. The experiment
#   therefore quantifies model misspecification rather than simply fitting the
#   known simulated residual-transcription level.
#
# Primary objectives:
#   1. Determine empirical Type-I error as residual transcription increases.
#   2. Quantify statistical power under partial shutoff.
#   3. Identify regimes where incomplete shutoff induces false-positive
#      evidence for sigma_c.
#   4. Evaluate whether misspecification is absorbed by sigma_c, export, or
#      degradation-rate estimates.
#   5. Quantify changes in practical identifiability and matrix conditioning.
#
# Collected outputs include:
#   - p-values and q-values;
#   - sigma_c and all kinetic parameter estimates;
#   - true parameters;
#   - RSS improvement;
#   - condition numbers;
#   - NNLS boundary mass;
#   - bootstrap failure and rank diagnostics;
#   - platform, residual-transcription level, noise, design, and seeds.
#
# Intended use:
#   This analysis directly supports reviewer requests concerning:
#     - pseudo-shutoff assumptions;
#     - residual transcription;
#     - model misspecification;
#     - false-positive rates;
#     - practical identifiability.
#
# Input:
#   ode_states_2k_20p.rdata
#
# Expected perturbation labels:
#       SHUTOFF
#       PSEUDO_SHUTOFF
#
# Output:
#   benchmark_pseudoshutoff_revision_raw.rdata
#   benchmark_pseudoshutoff_revision_raw.tsv
#   benchmark_pseudoshutoff_revision_summary.rdata
#   benchmark_pseudoshutoff_revision_summary.tsv
#   benchmark_pseudoshutoff_sessionInfo.txt
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
#   Major-revision misspecification benchmark, 2026.
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
source("../commons/platforms.r")


# =============================================================================
# 1. Load data
# =============================================================================

load("ode_states_2k_20p.rdata")


# =============================================================================
# 2. Benchmark design
# =============================================================================

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

# For the misspecification experiment I would initially use the experimentally
# relevant sparse designs rather than every enormous crossing.
RANGE_N_TIME <- c(
  3,
  5
)

RANGE_N_REP <- c(
  3,
  5
)

RANGE_TSTEP <- c(
  5,
  10,
  20
)

SCALING_A <- TRUE

LAMBDA_TIME <- 0.5
LAMBDA_DIAG <- 0.1
REL_FLOOR <- 1e-8

TRUNCATE_NONNEGATIVE_BOOT <- FALSE
MAX_FAILURE_RATE <- 0.05

N_BOOT <- 1999L

N_WORKERS <- 100L

BATCH_SIZE <- 25L


# =============================================================================
# 3. Restrict pseudo-shutoff dataset
# =============================================================================

ode_states_pseudo <- ode_states[
  Perturbation %in% c(
    "SHUTOFF",
    "PSEUDO_SHUTOFF"
  ) &
    N_time_samples %in%
      RANGE_N_TIME &
    N_replicates %in%
      RANGE_N_REP &
    T_step %in%
      RANGE_TSTEP
]


setkey(
  ode_states_pseudo,
  Gene,
  N_time_samples,
  N_replicates,
  T_step,
  Perturbation,
  Post_R_fraction
)


genes <- sort(
  unique(
    ode_states_pseudo$Gene
  )
)


rho_values <- sort(
  unique(
    ode_states_pseudo$Post_R_fraction
  )
)


cat(
  "Residual transcription fractions:",
  paste(
    rho_values,
    collapse = ", "
  ),
  "\n"
)


rm(ode_states)
gc()


# =============================================================================
# 4. Noise settings
# =============================================================================

RANGE_GAUSS_NOISE <- c(
  "Very low" = 0.02,
  "Low"      = 0.05,
  "Medium"   = 0.10,
  "High"     = 0.20
)


add_platform_noise_pseudo <- function(
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

    return(
      add_gaussian_noise(
        dt,
        cols = targets,
        noise_sd =
          unname(
            RANGE_GAUSS_NOISE[noise]
          )
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


  stop("Unknown platform.")
}


# =============================================================================
# 5. Seed helper
# =============================================================================

make_pseudo_seed <- function(
  gene,
  rho,
  ntime,
  nrep,
  tstep,
  platform,
  noise,
  stream = c(
    "measurement",
    "bootstrap"
  )
) {

  stream <- match.arg(stream)


  rho_code <- as.integer(
    round(
      rho *
        1000
    )
  )


  stream_id <- if (
    stream == "measurement"
  ) {
    1L
  } else {
    2L
  }


  seed <- (
    500000L +
      gene * 100000L +
      rho_code * 1000L +
      ntime * 100L +
      nrep * 10L +
      tstep +
      match(
        platform,
        RANGE_PLATFORM
      ) * 10000000L +
      match(
        noise,
        RANGE_NOISE
      ) * 20000000L +
      stream_id
  )


  as.integer(
    seed %% 2000000000L
  )
}


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
    return(default)
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
    return(NA_real_)
  }

  unname(
    x[name]
  )
}


# =============================================================================
# 6. Gene worker
# =============================================================================

run_gene_pseudo <- function(
  gene_i,
  N_boot
) {

  out <- list()

  kk <- 1L


  dt_gene <- ode_states_pseudo[
    Gene == gene_i
  ]


  if (nrow(dt_gene) == 0L) {
    return(NULL)
  }


  for (ntime_i in RANGE_N_TIME) {

    for (nrep_i in RANGE_N_REP) {

      for (tstep_i in RANGE_TSTEP) {


        dt_design <- dt_gene[
          N_time_samples == ntime_i &
            N_replicates == nrep_i &
            T_step == tstep_i
        ]


        if (nrow(dt_design) == 0L) {
          next
        }


        for (rho_i in rho_values) {


          dt_rho <- dt_design[
            Post_R_fraction == rho_i
          ]


          if (nrow(dt_rho) == 0L) {
            next
          }


          truth_i <- unique(
            dt_rho$truth_pos
          )[1]


          tstar <- unique(
            dt_rho$T_star
          )

          tstar <- tstar[
            is.finite(tstar)
          ]


          if (length(tstar) != 1L) {
            next
          }


          t_star_fit <- tstar[1]


          perturb_i <- unique(
            dt_rho$Perturbation
          )[1]


          for (platform_i in RANGE_PLATFORM) {

            for (noise_i in RANGE_NOISE) {


              measurement_seed <- make_pseudo_seed(
                gene = gene_i,
                rho = rho_i,
                ntime = ntime_i,
                nrep = nrep_i,
                tstep = tstep_i,
                platform = platform_i,
                noise = noise_i,
                stream = "measurement"
              )


              bootstrap_seed <- make_pseudo_seed(
                gene = gene_i,
                rho = rho_i,
                ntime = ntime_i,
                nrep = nrep_i,
                tstep = tstep_i,
                platform = platform_i,
                noise = noise_i,
                stream = "bootstrap"
              )


              set.seed(
                measurement_seed
              )


              dt_noisy <- tryCatch(

                add_platform_noise_pseudo(
                  dt = dt_rho,
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


              cf <- WW$coef_full


              out[[kk]] <- data.frame(

                Gene =
                  gene_i,

                Positive =
                  truth_i,

                Perturbation =
                  perturb_i,

                Post_R_fraction =
                  rho_i,

                Residual_transcription_pct =
                  100 *
                    rho_i,

                Platform =
                  platform_i,

                Exprs_noise =
                  noise_i,

                N_tsamples =
                  ntime_i,

                N_replicates =
                  nrep_i,

                Tsteps =
                  tstep_i,

                N_boot =
                  N_boot,


                # -------------------------------------------------------------
                # Inference
                # -------------------------------------------------------------

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

                atom_zero =
                  safe_value(
                    WW,
                    "atom.zero"
                  ),


                # -------------------------------------------------------------
                # Identifiability
                # -------------------------------------------------------------

                condition_number =
                  safe_value(
                    WW,
                    "condition.number"
                  ),

                rank_full =
                  safe_value(
                    WW,
                    "rank.full"
                  ),

                bootstrap_failure_rate =
                  safe_value(
                    WW,
                    "bootstrap.failure.rate"
                  ),

                boot_condition_q95 =
                  safe_value(
                    WW,
                    "bootstrap.condition.q95"
                  ),

                boot_rank_deficient_fraction =
                  safe_value(
                    WW,
                    "bootstrap.rank.deficient.fraction"
                  ),


                # -------------------------------------------------------------
                # Parameter estimates
                # -------------------------------------------------------------

                R_hat =
                  safe_coef(
                    cf,
                    "R"
                  ),

                tau_hat =
                  safe_coef(
                    cf,
                    "tau"
                  ),

                tau_s_hat =
                  safe_coef(
                    cf,
                    "tau_s"
                  ),

                sigma_c_hat =
                  safe_coef(
                    cf,
                    "sigma_c"
                  ),

                sigma_n_hat =
                  safe_coef(
                    cf,
                    "sigma_n"
                  ),

                alpha_hat =
                  safe_coef(
                    cf,
                    "alpha"
                  ),

                alpha_s_hat =
                  safe_coef(
                    cf,
                    "alpha_s"
                  ),


                # -------------------------------------------------------------
                # Ground truth
                # -------------------------------------------------------------

                R_true =
                  unique(
                    dt_rho$Base_R
                  )[1],

                tau_true =
                  unique(
                    dt_rho$Base_tau
                  )[1],

                tau_s_true =
                  unique(
                    dt_rho$Base_tau_s
                  )[1],

                sigma_c_true =
                  unique(
                    dt_rho$Base_sigma_c
                  )[1],

                sigma_n_true =
                  unique(
                    dt_rho$Base_sigma_n
                  )[1],

                alpha_true =
                  unique(
                    dt_rho$Base_alpha
                  )[1],

                alpha_s_true =
                  unique(
                    dt_rho$Base_alpha_s
                  )[1],


                # -------------------------------------------------------------
                # Runtime / status
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

  out
}


# =============================================================================
# 7. Parallel execution
# =============================================================================

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


objects_to_export <- c(

  "ode_states_pseudo",
  "rho_values",

  "RANGE_PLATFORM",
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

  "add_platform_noise_pseudo",
  "make_pseudo_seed",
  "safe_value",
  "safe_coef",
  "run_gene_pseudo",

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


gene_batches <- split(
  genes,
  ceiling(
    seq_along(genes) /
      BATCH_SIZE
  )
)


all_batches <- vector(
  "list",
  length(gene_batches)
)


for (bb in seq_along(gene_batches)) {

  cat(
    "\nPseudo-shutoff batch",
    bb,
    "/",
    length(gene_batches),
    "\n"
  )


  res_bb <- parallel::parLapplyLB(
    cl,
    gene_batches[[bb]],
    fun = run_gene_pseudo,
    N_boot = N_BOOT
  )


  res_bb <- Filter(
    Negate(is.null),
    res_bb
  )


  flat_bb <- unlist(
    res_bb,
    recursive = FALSE
  )


  dt_bb <- data.table::rbindlist(
    flat_bb,
    use.names = TRUE,
    fill = TRUE
  )


  all_batches[[bb]] <- dt_bb


  save(
    dt_bb,
    file = sprintf(
      "benchmark_pseudoshutoff_batch_%03d.rdata",
      bb
    )
  )


  data.table::fwrite(
    dt_bb,
    file = sprintf(
      "benchmark_pseudoshutoff_batch_%03d.tsv",
      bb
    ),
    sep = "\t"
  )
}


parallel::stopCluster(
  cl
)


# =============================================================================
# 8. Combine and derive quantities
# =============================================================================

dt_pseudo <- data.table::rbindlist(
  all_batches,
  use.names = TRUE,
  fill = TRUE
)


dt_pseudo[
  ,
  valid_test :=
    status == "ok" &
    is.finite(
      p.value
    )
]


dt_pseudo[
  ,
  q.value :=
    p.adjust(
      p.value,
      method = "BH"
    ),
  by = .(
    Post_R_fraction,
    Platform,
    Exprs_noise,
    N_tsamples,
    N_replicates,
    Tsteps
  )
]


dt_pseudo[
  ,
  DeltaRSS :=
    pmax(
      RSS0 - RSS1,
      0
    )
]


dt_pseudo[
  ,
  sigma_c_error :=
    sigma_c_hat -
    sigma_c_true
]


dt_pseudo[
  ,
  tau_error :=
    tau_hat -
    tau_true
]


dt_pseudo[
  ,
  alpha_error :=
    alpha_hat -
    alpha_true
]


# =============================================================================
# 9. Summary
# =============================================================================

pseudo_summary <- dt_pseudo[
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


    list(

      N_total =
        .N,

      N_valid =
        sum(
          valid
        ),

      TypeI_005 =
        mean(
          p0 <= 0.05,
          na.rm = TRUE
        ),

      Inflation_005 =
        mean(
          p0 <= 0.05,
          na.rm = TRUE
        ) /
        0.05,

      Power_005 =
        mean(
          p1 <= 0.05,
          na.rm = TRUE
        ),

      Median_sigma_null =
        median(
          sigma_c_hat[
            valid &
              Positive == 0
          ],
          na.rm = TRUE
        ),

      sigma_c_bias_alt =
        mean(
          sigma_c_error[
            valid &
              Positive == 1
          ],
          na.rm = TRUE
        ),

      Median_tau_error =
        median(
          tau_error[
            valid
          ],
          na.rm = TRUE
        ),

      Median_alpha_error =
        median(
          alpha_error[
            valid
          ],
          na.rm = TRUE
        ),

      Median_condition =
        median(
          condition_number[
            valid
          ],
          na.rm = TRUE
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
        )
    )
  },

  by = .(
    Post_R_fraction,
    Residual_transcription_pct,
    Platform,
    Exprs_noise,
    N_tsamples,
    N_replicates,
    Tsteps
  )
]


# =============================================================================
# 10. Save
# =============================================================================

save(
  dt_pseudo,
  file =
    "benchmark_pseudoshutoff_revision_raw.rdata"
)


data.table::fwrite(
  dt_pseudo,
  file =
    "benchmark_pseudoshutoff_revision_raw.tsv",
  sep = "\t"
)


save(
  pseudo_summary,
  file =
    "benchmark_pseudoshutoff_revision_summary.rdata"
)


data.table::fwrite(
  pseudo_summary,
  file =
    "benchmark_pseudoshutoff_revision_summary.tsv",
  sep = "\t"
)


sink(
  "benchmark_pseudoshutoff_sessionInfo.txt"
)

print(
  sessionInfo()
)

sink()


cat(
  "\nPseudo-shutoff benchmark complete.\n"
)