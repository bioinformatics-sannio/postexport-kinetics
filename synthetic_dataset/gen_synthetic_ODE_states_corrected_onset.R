# =============================================================================
# Title:
# Synthetic Dataset Generation for ODE-Based RNA Kinetics with Complete,
# Partial, and Continuous Transcription Regimes
#
# Description:
#   This script generates the synthetic benchmark dataset used to evaluate
#   the statistical framework for detecting post-export RNA conversion from
#   compartment-resolved time-course data.
#
#   For each simulated gene, the script:
#     - samples a gene-specific kinetic parameter set;
#     - assigns the gene to the null or alternative class;
#     - sets sigma_c = 0 under the null and sigma_c > 0 under the alternative;
#     - simulates biological replicate variability;
#     - generates trajectories under multiple sampling designs;
#     - generates continuous-transcription, complete-shutoff, and
#       partial/pseudo-shutoff conditions;
#     - stores pre-intervention steady-state values;
#     - attaches the ground-truth kinetic parameters and experimental-design
#       metadata required for calibration, power, parameter-recovery, and
#       model-misspecification analyses.
#
# Experimental regimes:
#
#   1. NONE
#      Continuous transcription throughout the experiment.
#
#   2. SHUTOFF
#      Complete transcriptional shutoff at T_star:
#
#          R_post = 0
#
#   3. PSEUDO_SHUTOFF
#      Incomplete transcriptional shutoff at T_star:
#
#          R_post = rho * R_pre
#
#      where rho is the residual-transcription fraction. The default simulated
#      values are:
#
#          rho = 0.05, 0.10, 0.25, 0.50
#
#      These scenarios are included specifically to assess robustness to
#      violations of the complete-shutoff assumption and to support the
#      sensitivity analyses requested during peer review.
#
#   4. SteadyState
#      Pre-intervention steady-state quantities generated from the same
#      gene-specific kinetic parameters.
#
# Sampling model:
#   Biological replicates represent destructive biological samples. Replicate
#   identifiers are used for simulation bookkeeping but do not imply
#   longitudinal pairing across sampling times. This reflects experimental
#   designs in which a culture is harvested and destroyed at each time point.
#
# Ground truth:
#   Each simulated gene is assigned:
#
#       truth_pos = 0  -> sigma_c = 0
#       truth_pos = 1  -> sigma_c > 0
#
#   The original kinetic parameters are stored using the prefix "Base_".
#
# Output:
#   The complete synthetic benchmark is stored in the object:
#
#       ode_states
#
#   and saved to:
#
#       ode_states_2k_20p.rdata
#
#   The same file contains NONE, SHUTOFF, PSEUDO_SHUTOFF, and SteadyState
#   observations. This filename is intentionally retained for compatibility
#   with the existing downstream analysis scripts.
#
#   Additional columns introduced for the robustness analysis include:
#
#       Post_R_fraction
#       Residual_transcription_pct
#
# Reproducibility:
#   Gene-specific random seeds are used so that kinetic parameters and class
#   assignments are reproducible across runs. Simulation-design metadata and
#   R session information are saved separately.
#
# Intended use:
#   This code supports:
#     - empirical Type-I error assessment;
#     - statistical power analysis;
#     - parameter-recovery and practical-identifiability analyses;
#     - comparison across experimental platforms and sampling designs;
#     - complete versus partial transcriptional-shutoff sensitivity analysis;
#     - model-misspecification analyses requested during peer review.
#
# Scientific scope:
#   The simulated post-export conversion parameter sigma_c represents a
#   generic post-export conversion process within the kinetic model. Simulation
#   results should therefore be interpreted as evaluation of the statistical
#   and kinetic inference framework rather than as experimental confirmation
#   of a specific molecular mechanism.
#
# Dependencies:
#   - R
#   - data.table
#   - parallel
#   - deSolve
#   - local ODE utilities in ../ode_model/ode.r
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
#   The software may be redistributed as part of research repositories,
#   supplementary materials, and archival research releases, including
#   services such as GitHub, Zenodo, Figshare, or Software Heritage, subject
#   to the conditions above.
#
# Disclaimer:
#   This software is provided "as is", without warranty of any kind, express
#   or implied, including but not limited to the warranties of merchantability,
#   fitness for a particular purpose, and noninfringement. In no event shall
#   the author be liable for any claim, damages, or other liability arising
#   from, out of, or in connection with the software or its use.
#
# Version:
#   Revision developed for the 2026 major revision of the manuscript:
#
#   "A kinetic modeling framework to detect post-export RNA processing from
#    time-resolved transcriptomic data"
#
# Peer-review motivation:
#   The partial/pseudo-shutoff simulations were added to quantify the effect
#   of residual transcription and violations of the idealized complete-shutoff
#   assumption on statistical calibration and inference for sigma_c.
# =============================================================================

library(data.table)
library(parallel)


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/synthetic_dataset")

source("../ode_model/ode.r")


# Prevent nested numerical parallelism.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

data.table::setDTthreads(1)


# =============================================================================
# 1. Sampling-design grid
# =============================================================================

range_n_replicates <- c(
  3,
  5,
  10
)

range_time_samples <- c(
  3,
  5,
  10,
  20
)

range_tstep <- c(
  5,
  10,
  20,
  50
)


# =============================================================================
# 2. Residual-transcription levels
# =============================================================================

# Fraction of transcription remaining after T_star.
#
# 0.00 = complete shutoff
# 0.05 =  5% residual transcription
# 0.10 = 10%
# 0.25 = 25%
# 0.50 = 50%
#
# NONE already corresponds to continuous transcription, so rho=1
# is not generated again among intervention conditions.
range_post_R_fraction <- c(
  0,
  0.05,
  0.10,
  0.25,
  0.50,
  0.75,
  0.90,
  1.00
)

# =============================================================================
# 3. Global benchmark settings
# =============================================================================

ngenes <- 2000L

frac_pos <- 0.20


# =============================================================================
# 4. Global simulation time horizon
# =============================================================================

tmax <- 3 / min(
  r_sigma_n_min + r_tau_min,
  r_tau_s_min,
  r_sigma_c_min + r_alpha_min,
  r_alpha_s_min
)

times <- seq(
  0,
  tmax,
  by = 1
)


# Intervention time used throughout the benchmark.
T_star <- times[
  floor(
    length(times) / 3
  )
]


cat(
  "Simulation horizon:",
  min(times),
  "-",
  max(times),
  "\n"
)

cat(
  "T_star:",
  T_star,
  "\n"
)


# =============================================================================
# 4b. Gene-specific transcription-onset heterogeneity
# =============================================================================
#
# The intervention time T_star is common to every gene. A single gene-specific
# transcription onset is generated once and reused across every sampling design
# and every perturbation regime for that gene.
#
MAX_ONSET_SHIFT <- floor(
  0.1 * max(times)
)


# =============================================================================
# 5. Design grid
# =============================================================================

grid <- CJ(
  N_replicates = range_n_replicates,
  N_time_samples = range_time_samples,
  T_step = range_tstep
)


cat(
  "Number of design configurations:",
  nrow(grid),
  "\n"
)


# =============================================================================
# 6. Simulate one gene
# =============================================================================

simulate_one_gene <- function(
  gene_id,
  grid,
  frac_pos,
  times,
  T_star,
  range_post_R_fraction,
  max_onset_shift
) {

  # ---------------------------------------------------------------------------
  # Gene-specific reproducibility
  # ---------------------------------------------------------------------------

  set.seed(
    1000000L +
      as.integer(gene_id)
  )


  # ---------------------------------------------------------------------------
  # Sample kinetic parameters
  # ---------------------------------------------------------------------------

  r_param <- random_params()


  # ---------------------------------------------------------------------------
  # Assign true class
  # ---------------------------------------------------------------------------

  truth_pos <- as.integer(
    runif(1) <= frac_pos
  )


  if (truth_pos == 1L) {

    r_param$sigma_c <- runif(
      1,
      r_sigma_c_min,
      r_sigma_c_max
    )

  } else {

    r_param$sigma_c <- 0
  }


  # ---------------------------------------------------------------------------
  # Sample transcription rate
  # ---------------------------------------------------------------------------

  ALPHA_R <- 5
  BETA_R <- 20

  r_param$R <- rgamma(
    1,
    shape = ALPHA_R,
    scale = BETA_R
  )


  # ---------------------------------------------------------------------------
  # Constrain alpha_s relative to alpha
  # ---------------------------------------------------------------------------

  upper_alpha_s <- min(
    r_param$alpha,
    r_alpha_s_max
  )


  if (upper_alpha_s <= r_alpha_s_min) {

    r_param$alpha_s <- upper_alpha_s

  } else {

    r_param$alpha_s <- runif(
      1,
      r_alpha_s_min,
      upper_alpha_s
    )
  }


  # ---------------------------------------------------------------------------
  # Base parameter metadata
  # ---------------------------------------------------------------------------

  base_param_vec <- unlist(
    r_param
  )

  names(base_param_vec) <- paste0(
    "Base_",
    names(base_param_vec)
  )


  # ---------------------------------------------------------------------------
  # Gene-specific transcription onset
  # ---------------------------------------------------------------------------
  #
  # Use a dedicated seed and restore the RNG state afterwards. This preserves
  # the gene-specific kinetic/truth draws produced above.
  #
  rng_state_after_parameters <- .Random.seed

  set.seed(
    2000000L +
      as.integer(gene_id)
  )

  onset_shift_gene <- sample(
    -as.integer(max_onset_shift):
      as.integer(max_onset_shift),
    1L
  )

  .Random.seed <- rng_state_after_parameters

  # Nominal onset is zero.
  onset_time_gene <- onset_shift_gene


  # ---------------------------------------------------------------------------
  # Output container
  # ---------------------------------------------------------------------------

  out_list <- list()

  out_index <- 1L


  # ===========================================================================
  # Loop over designs
  # ===========================================================================

  for (gi in seq_len(nrow(grid))) {

    n_replicates_i <- grid[
      gi,
      N_replicates
    ]

    time_samples_i <- grid[
      gi,
      N_time_samples
    ]

    tstep_i <- grid[
      gi,
      T_step
    ]


    # =========================================================================
    # Baseline sampling times
    # =========================================================================

    u <- seq(
      0,
      1,
      length.out = time_samples_i + 2
    )

    p <- 3

    idx <- floor(
      max(times) /
        3 *
        (
          u ^ p
        )
    )


    sampled_times <- sort(
      unique(idx)
    )


    # Remove first and last points as in the original generator.
    if (length(sampled_times) >= 3L) {

      sampled_times <- sampled_times[-1]

      sampled_times <- sampled_times[
        -length(sampled_times)
      ]
    }


    # =========================================================================
    # Shutoff / pseudo-shutoff sampling times
    # =========================================================================

    if (time_samples_i < 2L) {
      stop(
        "N_time_samples must be >= 2."
      )
    }


    if (time_samples_i == 2L) {

      shutoff_sampled_times <- c(
        T_star - tstep_i,
        T_star
      )

    } else {

      step_adders <- cumsum(
        rep(
          tstep_i,
          time_samples_i - 2L
        )
      )


      shutoff_sampled_times <- sort(
        c(
          T_star - tstep_i,
          T_star,
          T_star + step_adders
        )
      )
    }


    # Keep only values within simulation horizon.
    shutoff_sampled_times <- shutoff_sampled_times[
      shutoff_sampled_times >= min(times) &
        shutoff_sampled_times <= max(times)
    ]


    # =========================================================================
    # Initial conditions
    # =========================================================================

    y0 <- c(
      N = 0,
      N_s = 0,
      C = 0,
      C_s = 0
    )


    # =========================================================================
    # A. NONE
    # =========================================================================

    dd_none <- generate_ODE_states(

      base_params = r_param,

      y0 = y0,

      times = times,

      n_replicates = n_replicates_i,

      model_kinetics = rna_kinetics,

      stimes = sampled_times,

      shutofftimes = shutoff_sampled_times,

      t_star = NULL,

      post_R_fraction = 1,

      # Fixed gene-specific onset; do not draw another shift inside ode.r.
      use_onset_shift = FALSE,
      max_shift = 0,
      nominal_onset_time = onset_time_gene
    )


    dt_none <- as.data.table(
      dd_none$tsampled_data
    )


    dt_none[
      ,
      `:=`(
        Perturbation = "NONE",
        Post_R_fraction = 1,
        Residual_transcription_pct = 100
      )
    ]


    dt_none[
      ,
      `:=`(
        Gene = gene_id,
        truth_pos = truth_pos,
        N_replicates = n_replicates_i,
        N_time_samples = time_samples_i,
        T_step = tstep_i,
        T_star = T_star,
        Onset_shift = onset_shift_gene,
        Onset_time = onset_time_gene
      )
    ]


    dt_none[
      ,
      (names(base_param_vec)) :=
        as.list(base_param_vec)
    ]


    out_list[[out_index]] <- dt_none

    out_index <- out_index + 1L


    # =========================================================================
    # B. SHUTOFF + PSEUDO_SHUTOFF
    # =========================================================================

    for (rho_i in range_post_R_fraction) {


      dd_rho <- generate_ODE_states(

        base_params = r_param,

        y0 = y0,

        times = times,

        n_replicates = n_replicates_i,

        model_kinetics = rna_kinetics,

        stimes = sampled_times,

        shutofftimes = shutoff_sampled_times,

        t_star = T_star,

        post_R_fraction = rho_i,

        # Same transcriptional history as NONE for this gene.
        # The pharmacological intervention remains at the common T_star.
        use_onset_shift = FALSE,
        max_shift = 0,
        nominal_onset_time = onset_time_gene
      )


      dt_rho <- as.data.table(
        dd_rho$intervention_tsampled_data
      )


      perturb_label <- if (
        rho_i == 0
      ) {

        "SHUTOFF"

      } else {

        "PSEUDO_SHUTOFF"
      }


      dt_rho[
        ,
        `:=`(
          Perturbation = perturb_label,
          Post_R_fraction = rho_i,
          Residual_transcription_pct = 100 * rho_i
        )
      ]


      dt_rho[
        ,
        `:=`(
          Gene = gene_id,
          truth_pos = truth_pos,
          N_replicates = n_replicates_i,
          N_time_samples = time_samples_i,
          T_step = tstep_i,
          T_star = T_star,
          Onset_shift = onset_shift_gene,
          Onset_time = onset_time_gene
        )
      ]


      dt_rho[
        ,
        (names(base_param_vec)) :=
          as.list(base_param_vec)
    ]


      out_list[[out_index]] <- dt_rho

      out_index <- out_index + 1L
    }


    # =========================================================================
    # C. SteadyState
    # =========================================================================

    dt_ss <- as.data.table(
      dd_none$ss_data
    )


    dt_ss[
      ,
      `:=`(
        Perturbation = "SteadyState",
        time = 0,
        Post_R_fraction = NA_real_,
        Residual_transcription_pct = NA_real_
      )
    ]


    dt_ss[
      ,
      `:=`(
        Gene = gene_id,
        truth_pos = truth_pos,
        N_replicates = n_replicates_i,
        N_time_samples = time_samples_i,
        T_step = tstep_i,
        T_star = T_star,
        Onset_shift = onset_shift_gene,
        Onset_time = onset_time_gene
      )
    ]


    dt_ss[
      ,
      (names(base_param_vec)) :=
        as.list(base_param_vec)
    ]


    out_list[[out_index]] <- dt_ss

    out_index <- out_index + 1L
  }


  # ---------------------------------------------------------------------------
  # Combine all scenarios for the gene
  # ---------------------------------------------------------------------------

  data.table::rbindlist(
    out_list,
    use.names = TRUE,
    fill = TRUE
  )
}


# =============================================================================
# 7. Syntax-safe small sanity test
# =============================================================================

cat(
  "\nRunning one-gene sanity test...\n"
)


test_grid <- grid[
  N_replicates == 3 &
    N_time_samples == 5 &
    T_step == 10
]


test_gene <- simulate_one_gene(

  gene_id = 1L,

  grid = test_grid,

  frac_pos = frac_pos,

  times = times,

  T_star = T_star,

  range_post_R_fraction =
    range_post_R_fraction,

  max_onset_shift =
    MAX_ONSET_SHIFT
)


cat(
  "\nSanity-test rows by condition:\n"
)


print(
  test_gene[
    ,
    .N,
    by = .(
      Perturbation,
      Post_R_fraction
    )
  ][
    order(
      Post_R_fraction
    )
  ]
)


# =============================================================================
# 8. Sanity checks
# =============================================================================

# ---------------------------------------------------------------------------
# SHUTOFF must exist
# ---------------------------------------------------------------------------

if (
  test_gene[
    Perturbation == "SHUTOFF",
    .N
  ] == 0L
) {

  stop(
    "Sanity check failed: SHUTOFF rows missing."
  )
}


# ---------------------------------------------------------------------------
# PSEUDO_SHUTOFF must exist at every requested nonzero rho
# ---------------------------------------------------------------------------

expected_pseudo <- range_post_R_fraction[
  range_post_R_fraction > 0
]


observed_pseudo <- sort(
  unique(
    test_gene[
      Perturbation == "PSEUDO_SHUTOFF",
      Post_R_fraction
    ]
  )
)


if (
  !identical(
    as.numeric(observed_pseudo),
    as.numeric(expected_pseudo)
  )
) {

  stop(
    "Sanity check failed: pseudo-shutoff levels are incomplete."
  )
}


# ---------------------------------------------------------------------------
# SHUTOFF must have rho=0
# ---------------------------------------------------------------------------

if (
  any(
    test_gene[
      Perturbation == "SHUTOFF",
      Post_R_fraction
    ] != 0
  )
) {

  stop(
    "Sanity check failed: SHUTOFF does not have Post_R_fraction=0."
  )
}


# ---------------------------------------------------------------------------
# NONE must have rho=1
# ---------------------------------------------------------------------------

if (
  any(
    test_gene[
      Perturbation == "NONE",
      Post_R_fraction
    ] != 1
  )
) {

  stop(
    "Sanity check failed: NONE does not have Post_R_fraction=1."
  )
}


# ---------------------------------------------------------------------------
# Timing consistency: one onset and one T_star per gene
# ---------------------------------------------------------------------------

if (
  uniqueN(
    test_gene$Onset_time
  ) != 1L
) {
  stop(
    "Sanity check failed: one gene has multiple transcription-onset times."
  )
}

if (
  uniqueN(
    test_gene$T_star
  ) != 1L
) {
  stop(
    "Sanity check failed: one gene has multiple T_star values."
  )
}

cat(
  "\nGene-specific timing check:\n"
)

print(
  unique(
    test_gene[
      ,
      .(
        Gene,
        Onset_shift,
        Onset_time,
        T_star
      )
    ]
  )
)

cat(
  "\nOne-gene sanity test PASSED.\n"
)


# =============================================================================
# 9. Parallel generation
# =============================================================================

detected_cores <- parallel::detectCores()


# You can set this explicitly if desired.
# Example:
#
# ncores <- 100L

ncores <- max(
  1L,
  detected_cores - 80L
)


cat(
  "\n============================================================\n"
)

cat(
  "STARTING SYNTHETIC DATA GENERATION\n"
)

cat(
  "Genes:",
  ngenes,
  "\n"
)

cat(
  "Workers:",
  ncores,
  "\n"
)

cat(
  "Pseudo-shutoff levels:",
  paste(
    range_post_R_fraction[
      range_post_R_fraction > 0
    ],
    collapse = ", "
  ),
  "\n"
)

cat(
  "============================================================\n"
)


generation_start <- Sys.time()


res_list <- parallel::mclapply(

  X = seq_len(ngenes),

  FUN = simulate_one_gene,

  grid = grid,

  frac_pos = frac_pos,

  times = times,

  T_star = T_star,

  range_post_R_fraction =
    range_post_R_fraction,

  max_onset_shift =
    MAX_ONSET_SHIFT,

  mc.cores = ncores,

  mc.set.seed = TRUE
)


generation_elapsed <- difftime(
  Sys.time(),
  generation_start,
  units = "hours"
)


# =============================================================================
# 10. Combine all genes
# =============================================================================

cat(
  "\nCombining synthetic data...\n"
)


ode_states <- data.table::rbindlist(
  res_list,
  use.names = TRUE,
  fill = TRUE
)


rm(res_list)

gc()


# =============================================================================
# 11. Integrity checks
# =============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "DATASET INTEGRITY CHECKS\n"
)

cat(
  "============================================================\n"
)


cat(
  "Rows:",
  format(
    nrow(ode_states),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Genes:",
  uniqueN(
    ode_states$Gene
  ),
  "\n"
)


# ---------------------------------------------------------------------------
# Timing consistency
# ---------------------------------------------------------------------------

timing_consistency <- ode_states[
  ,
  .(
    n_onset_shift = uniqueN(
      Onset_shift
    ),
    n_onset_time = uniqueN(
      Onset_time
    ),
    n_tstar = uniqueN(
      T_star
    )
  ),
  by = Gene
]

if (
  any(
    timing_consistency$n_onset_shift != 1L |
      timing_consistency$n_onset_time != 1L |
      timing_consistency$n_tstar != 1L
  )
) {
  stop(
    "Timing integrity check failed: onset or T_star is inconsistent within at least one gene."
  )
}


# ---------------------------------------------------------------------------
# Truth consistency
# ---------------------------------------------------------------------------

truth_consistency <- ode_states[
  ,
  .(
    n_truth = uniqueN(
      truth_pos
    )
  ),
  by = Gene
]


if (
  any(
    truth_consistency$n_truth != 1L
  )
) {

  stop(
    "truth_pos is inconsistent within at least one gene."
  )
}


# ---------------------------------------------------------------------------
# Truth distribution
# ---------------------------------------------------------------------------

gene_truth <- unique(
  ode_states[
    ,
    .(
      Gene,
      truth_pos
    )
  ]
)


cat(
  "\nTruth distribution:\n"
)


print(
  gene_truth[
    ,
    .N,
    by = truth_pos
  ]
)


# ---------------------------------------------------------------------------
# Null sigma_c
# ---------------------------------------------------------------------------

null_check <- unique(
  ode_states[
    truth_pos == 0,
    .(
      Gene,
      Base_sigma_c
    )
  ]
)


if (
  any(
    null_check$Base_sigma_c != 0
  )
) {

  stop(
    "A true-null gene has Base_sigma_c != 0."
  )
}


# ---------------------------------------------------------------------------
# Positive sigma_c
# ---------------------------------------------------------------------------

positive_check <- unique(
  ode_states[
    truth_pos == 1,
    .(
      Gene,
      Base_sigma_c
    )
  ]
)


if (
  any(
    positive_check$Base_sigma_c <= 0
  )
) {

  stop(
    "A positive gene has Base_sigma_c <= 0."
  )
}


# ---------------------------------------------------------------------------
# Perturbation / rho summary
# ---------------------------------------------------------------------------

cat(
  "\nPerturbation summary:\n"
)


perturbation_summary <- ode_states[
  ,
  .N,
  by = .(
    Perturbation,
    Post_R_fraction,
    Residual_transcription_pct
  )
][
  order(
    Post_R_fraction
  )
]


print(
  perturbation_summary
)


# =============================================================================
# 12. Save using ORIGINAL filename
# =============================================================================

output_file <- "ode_states_2k_20p_corrected_onset.rdata"


cat(
  "\nSaving complete synthetic dataset to:\n",
  output_file,
  "\n"
)


save(
  ode_states,
  file = output_file
)


# =============================================================================
# 13. Save simulation design separately
# =============================================================================

simulation_design <- list(

  ngenes =
    ngenes,

  frac_pos =
    frac_pos,

  range_n_replicates =
    range_n_replicates,

  range_time_samples =
    range_time_samples,

  range_tstep =
    range_tstep,

  range_post_R_fraction =
    range_post_R_fraction,

  T_star =
    T_star,

  max_onset_shift =
    MAX_ONSET_SHIFT,

  onset_semantics =
    "one gene-specific transcription onset; common pharmacological T_star",

  tmax =
    tmax,

  generated_at =
    Sys.time()
)


save(
  simulation_design,
  file =
    "ode_states_2k_20p_corrected_onset_design.rdata"
)


# =============================================================================
# 14. Session information
# =============================================================================

sink(
  "ode_states_2k_20p_corrected_onset_sessionInfo.txt"
)

print(
  sessionInfo()
)

sink()


# =============================================================================
# 15. Final report
# =============================================================================

cat(
  "\n============================================================\n"
)

cat(
  "SYNTHETIC DATA GENERATION COMPLETE\n"
)

cat(
  "============================================================\n"
)


cat(
  "Main dataset:",
  output_file,
  "\n"
)


cat(
  "Object name: ode_states\n"
)


cat(
  "Contains:\n"
)


cat(
  "  NONE\n"
)

cat(
  "  SHUTOFF\n"
)

cat(
  "  PSEUDO_SHUTOFF\n"
)

cat(
  "  SteadyState\n"
)


cat(
  "\nRows:",
  format(
    nrow(
      ode_states
    ),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Genes:",
  uniqueN(
    ode_states$Gene
  ),
  "\n"
)


cat(
  "Elapsed:",
  round(
    as.numeric(
      generation_elapsed
    ),
    2
  ),
  "hours\n"
)


cat(
  "\nAdditional files:\n"
)

cat(
  "  ode_states_2k_20p_design.rdata\n"
)

cat(
  "  ode_states_2k_20p_sessionInfo.txt\n"
)