# =============================================================================
# Cytoplasmic-only kinetic baseline
#
# Purpose
# -------
# Benchmark a reduced kinetic model that uses only the cytoplasmic states
# C(t) and C_s(t), without access to N(t) or N_s(t).
#
# This provides a stronger comparator than a Delta-PSI test because it:
#   - retains explicit time-course modeling;
#   - retains a non-negative post-export conversion rate sigma_c;
#   - performs a nested H0: sigma_c = 0 vs H1: sigma_c >= 0 comparison;
#   - uses bootstrap calibration.
#
# Reduced model
# -------------
# Because nuclear export is unobserved in the cytoplasmic-only baseline, the
# two nuclear inputs are represented by non-negative nuisance influx terms:
#
#   dC/dt   = beta_C  - (sigma_c + alpha) C
#   dC_s/dt = beta_Cs + sigma_c C - alpha_s C_s
#
# beta_C and beta_Cs are assumed constant over the analyzed post-shutoff
# sampling window. They are nuisance parameters, not biological estimates.
#
# Full cytoplasmic model:
#   theta = (beta_C, beta_Cs, sigma_c, alpha, alpha_s) >= 0
#
# Null model:
#   sigma_c = 0
#
# The interval-balance estimator uses the same trapezoidal idea as the main
# method:
#
#   (x_{i+1} - x_i) / Delta t
#       ~= f( (x_i + x_{i+1}) / 2 )
#
# Bootstrap
# ---------
# The null model is propagated exactly as a two-state affine linear system.
# Replicate-level residual pairs (C, C_s) are resampled within each time point,
# preserving the empirical covariance between the two cytoplasmic states.
#
# Design
# ------
# Uses the corrected synthetic SHUTOFF dataset and EXACTLY reproduces the
# benchmark measurement-noise realization through the saved measurement_seed.
#
# Representative designs:
#   5 time points x 3 replicates x T_step = 10
#   5 time points x 10 replicates x T_step = 10
#
# Platforms:
#   RT-qPCR, GAUSS, RNA-seq
#
# Noise:
#   Very low, Low, Medium, High
#
# Genes:
#   500 null + 250 alternative
#
# Bootstrap:
#   499 per test by default
#
# Outputs
# -------
# cytoplasmic_only_baseline/
#   cytoplasmic_only_raw.tsv
#   cytoplasmic_only_summary.tsv
#   cytoplasmic_only_comparison.tsv
#   cytoplasmic_only_progress.tsv
#   cytoplasmic_only_baseline.rdata
#   sessionInfo.txt
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

if (!requireNamespace("nnls", quietly = TRUE)) {
  stop(
    paste0(
      "Package 'nnls' is required. Install it with:\n",
      "  install.packages(\"nnls\")"
    )
  )
}

if (!requireNamespace("expm", quietly = TRUE)) {
  stop(
    paste0(
      "Package 'expm' is required. Install it with:\n",
      "  install.packages(\"expm\")"
    )
  )
}

data.table::setDTthreads(1)

source("../commons/platforms.r")


# =============================================================================
# 1. Files and settings
# =============================================================================

ODE_FILE <- "ode_states_2k_20p_corrected_onset.rdata"

FULL_RAW_FILE <- "benchmark_main_corrected_onset_raw.tsv"

OUTPUT_DIR <- "cytoplasmic_only_baseline"

CHECKPOINT_DIR <- file.path(
  OUTPUT_DIR,
  "checkpoints"
)

RAW_FILE <- file.path(
  OUTPUT_DIR,
  "cytoplasmic_only_raw.tsv"
)

SUMMARY_FILE <- file.path(
  OUTPUT_DIR,
  "cytoplasmic_only_summary.tsv"
)

COMPARISON_FILE <- file.path(
  OUTPUT_DIR,
  "cytoplasmic_only_comparison.tsv"
)

PROGRESS_FILE <- file.path(
  OUTPUT_DIR,
  "cytoplasmic_only_progress.tsv"
)

RDATA_FILE <- file.path(
  OUTPUT_DIR,
  "cytoplasmic_only_baseline.rdata"
)

SESSION_FILE <- file.path(
  OUTPUT_DIR,
  "sessionInfo.txt"
)

BENCHMARK_VERSION <- "cytoplasmic_only_v1"


# -----------------------------------------------------------------------------
# Inferential settings
# -----------------------------------------------------------------------------

N_BOOT <- 499L

REL_FLOOR <- 1e-8

TRUNCATE_NONNEGATIVE_BOOT <- TRUE

N_WORKERS <- 60L
TASK_BATCH_SIZE <- 50L
PSOCK_TIMEOUT_SECONDS <- 12L * 60L * 60L


# -----------------------------------------------------------------------------
# Gene counts
# -----------------------------------------------------------------------------

N_NULL <- 500L
N_ALT <- 250L

SELECTION_SEED <- 20260923L


# -----------------------------------------------------------------------------
# Representative designs
# -----------------------------------------------------------------------------

DESIGNS <- data.table(
  Design = c(
    "5tp_3rep_dt10",
    "5tp_10rep_dt10"
  ),
  N_tsamples = c(
    5L,
    5L
  ),
  N_replicates = c(
    3L,
    10L
  ),
  Tsteps = c(
    10L,
    10L
  )
)


# -----------------------------------------------------------------------------
# Platforms and noise levels
# -----------------------------------------------------------------------------

PLATFORMS <- c(
  "RT-qPCR",
  "GAUSS",
  "RNA-seq"
)

NOISE_LEVELS <- c(
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
# 2. Directories
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
# 3. Load corrected synthetic states and full benchmark
# =============================================================================

load(
  ODE_FILE
)

if (
  !exists("ode_states") ||
    !is.data.table(ode_states)
) {
  stop(
    "Object 'ode_states' was not found as a data.table."
  )
}

full_raw <- fread(
  FULL_RAW_FILE
)


# =============================================================================
# 4. Validate required columns
# =============================================================================

required_ode <- c(
  "Gene",
  "Perturbation",
  "truth_pos",
  "N_time_samples",
  "N_replicates",
  "T_step",
  "Post_R_fraction",
  "time",
  "replicate",
  "C",
  "C_s"
)

required_full <- c(
  "Gene",
  "Positive",
  "Perturbation",
  "Platform",
  "Exprs_noise",
  "N_tsamples",
  "N_replicates",
  "Tsteps",
  "measurement_seed",
  "p.value"
)

missing_ode <- setdiff(
  required_ode,
  names(
    ode_states
  )
)

missing_full <- setdiff(
  required_full,
  names(
    full_raw
  )
)

if (
  length(
    missing_ode
  ) > 0L
) {
  stop(
    paste0(
      "Missing columns in ode_states:\n  ",
      paste(
        missing_ode,
        collapse = "\n  "
      )
    )
  )
}

if (
  length(
    missing_full
  ) > 0L
) {
  stop(
    paste0(
      "Missing columns in full benchmark raw table:\n  ",
      paste(
        missing_full,
        collapse = "\n  "
      )
    )
  )
}


# =============================================================================
# 5. Restrict to complete SHUTOFF and selected designs
# =============================================================================

ode_sub <- ode_states[
  Perturbation == "SHUTOFF" &
    abs(
      Post_R_fraction
    ) <
      1e-12
]

ode_sub <- merge(
  ode_sub,
  DESIGNS[
    ,
    .(
      Design,
      N_time_samples =
        N_tsamples,
      N_replicates,
      T_step =
        Tsteps
    )
  ],
  by = c(
    "N_time_samples",
    "N_replicates",
    "T_step"
  ),
  all = FALSE
)

if (
  nrow(
    ode_sub
  ) ==
    0L
) {
  stop(
    "No selected complete-SHUTOFF trajectories were found."
  )
}


# =============================================================================
# 6. Select genes reproducibly
# =============================================================================

gene_truth <- unique(
  ode_sub[
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
      "null genes, but only",
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
      "alternative genes, but only",
      length(
        alt_genes_all
      ),
      "are available."
    )
  )
}

set.seed(
  SELECTION_SEED
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

ode_sub <- ode_sub[
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
# 7. Platform-noise helper
#
# Reproduces the observation model of the full benchmark.
# =============================================================================

add_platform_noise_main <- function(
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

  # Some reduced datasets passed to this helper still contain all four states.
  # The baseline later discards N and N_s.

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
# 8. Stable bootstrap seed
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

  for (
    ii in ints
  ) {
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


make_boot_seed <- function(
  gene,
  design,
  platform,
  noise
) {

  stable_seed_from_string(
    paste(
      BENCHMARK_VERSION,
      gene,
      design,
      platform,
      noise,
      "bootstrap",
      sep = "|"
    )
  )
}


# =============================================================================
# 9. Weighted trapezoidal design matrix
# =============================================================================

build_cyto_design <- function(
  dt,
  include_sigma = TRUE
) {

  means <- dt[
    ,
    .(
      C =
        mean(
          C,
          na.rm = TRUE
        ),
      C_s =
        mean(
          C_s,
          na.rm = TRUE
        )
    ),
    by =
      time
  ]

  setorder(
    means,
    time
  )

  if (
    nrow(
      means
    ) <
      3L
  ) {
    stop(
      "At least three sampled time points are required."
    )
  }


  # Mean-level measurement variance for each state/time.
  vars <- dt[
    ,
    .(
      var_C =
        if (
          .N >
            1L
        ) {
          var(
            C,
            na.rm = TRUE
          ) /
            .N
        } else {
          0
        },

      var_Cs =
        if (
          .N >
            1L
        ) {
          var(
            C_s,
            na.rm = TRUE
          ) /
            .N
        } else {
          0
        }
    ),
    by =
      time
  ]

  setorder(
    vars,
    time
  )


  X_list <- list()
  y_list <- list()
  var_list <- list()

  rr <- 1L


  for (
    ii in seq_len(
      nrow(
        means
      ) -
        1L
    )
  ) {

    dt_i <- means$time[ii + 1L] -
      means$time[ii]

    if (
      !is.finite(
        dt_i
      ) ||
        dt_i <= 0
    ) {
      stop(
        "Non-positive time interval."
      )
    }


    C0 <- means$C[ii]
    C1 <- means$C[ii + 1L]

    Cs0 <- means$C_s[ii]
    Cs1 <- means$C_s[ii + 1L]

    Cbar <- 0.5 *
      (
        C0 +
          C1
      )

    Csbar <- 0.5 *
      (
        Cs0 +
          Cs1
      )

    dC <- (
      C1 -
        C0
    ) /
      dt_i

    dCs <- (
      Cs1 -
        Cs0
    ) /
      dt_i


    # Full parameter order:
    # beta_C, beta_Cs, sigma_c, alpha, alpha_s

    X_C <- c(
      beta_C = 1,
      beta_Cs = 0,
      sigma_c = -Cbar,
      alpha = -Cbar,
      alpha_s = 0
    )

    X_Cs <- c(
      beta_C = 0,
      beta_Cs = 1,
      sigma_c = Cbar,
      alpha = 0,
      alpha_s = -Csbar
    )


    var_dC <- (
      vars$var_C[ii] +
        vars$var_C[ii + 1L]
    ) /
      dt_i^2

    var_dCs <- (
      vars$var_Cs[ii] +
        vars$var_Cs[ii + 1L]
    ) /
      dt_i^2


    X_list[[rr]] <- X_C
    y_list[[rr]] <- dC
    var_list[[rr]] <- var_dC
    rr <- rr + 1L

    X_list[[rr]] <- X_Cs
    y_list[[rr]] <- dCs
    var_list[[rr]] <- var_dCs
    rr <- rr + 1L
  }


  X <- do.call(
    rbind,
    X_list
  )

  y <- unlist(
    y_list,
    use.names = FALSE
  )

  vv <- unlist(
    var_list,
    use.names = FALSE
  )


  positive_var <- vv[
    is.finite(
      vv
    ) &
      vv >
        0
  ]

  floor_var <- if (
    length(
      positive_var
    ) >
      0L
  ) {
    max(
      quantile(
        positive_var,
        0.10,
        names = FALSE
      ) *
        0.1,
      REL_FLOOR
    )
  } else {
    1
  }

  vv[
    !is.finite(
      vv
    ) |
      vv <
        floor_var
  ] <- floor_var


  w <- 1 /
    sqrt(
      vv
    )

  Xw <- X *
    w

  yw <- y *
    w


  if (
    !include_sigma
  ) {

    Xw <- Xw[
      ,
      colnames(
        Xw
      ) !=
        "sigma_c",
      drop = FALSE
    ]
  }


  list(
    Xw =
      Xw,
    yw =
      yw,
    means =
      means
  )
}


# =============================================================================
# 10. Fit one model
# =============================================================================

fit_cyto_model <- function(
  dt,
  include_sigma = TRUE
) {

  design <- build_cyto_design(
    dt,
    include_sigma =
      include_sigma
  )

  fit <- nnls::nnls(
    A =
      design$Xw,
    b =
      design$yw
  )

  coef <- as.numeric(
    fit$x
  )

  names(
    coef
  ) <- colnames(
    design$Xw
  )

  residual <- design$yw -
    as.numeric(
      design$Xw %*%
        coef
    )

  rss <- sum(
    residual^2
  )

  list(
    coef =
      coef,
    rss =
      rss,
    means =
      design$means
  )
}


# =============================================================================
# 11. Exact null propagation
#
# Null model:
#   dC/dt   = beta_C  - alpha C
#   dC_s/dt = beta_Cs - alpha_s C_s
#
# Because sigma_c = 0, the two states are independent affine first-order ODEs.
# =============================================================================

propagate_null_mean <- function(
  times,
  initial_state,
  coef_null
) {

  beta_C <- unname(
    coef_null[
      "beta_C"
    ]
  )

  beta_Cs <- unname(
    coef_null[
      "beta_Cs"
    ]
  )

  alpha <- unname(
    coef_null[
      "alpha"
    ]
  )

  alpha_s <- unname(
    coef_null[
      "alpha_s"
    ]
  )


  t0 <- min(
    times
  )

  dt <- times -
    t0


  propagate_affine <- function(
    y0,
    beta,
    rate,
    dt
  ) {

    if (
      rate >
        1e-12
    ) {

      steady <- beta /
        rate

      return(
        steady +
          (
            y0 -
              steady
          ) *
          exp(
            -rate *
              dt
          )
      )
    }

    y0 +
      beta *
        dt
  }


  data.table(
    time =
      times,

    C =
      propagate_affine(
        y0 =
          initial_state[
            "C"
          ],
        beta =
          beta_C,
        rate =
          alpha,
        dt =
          dt
      ),

    C_s =
      propagate_affine(
        y0 =
          initial_state[
            "C_s"
          ],
        beta =
          beta_Cs,
        rate =
          alpha_s,
        dt =
          dt
      )
  )
}


# =============================================================================
# 12. Residual bootstrap under H0
# =============================================================================

bootstrap_cyto_test <- function(
  dt,
  B,
  seed
) {

  fit_full <- fit_cyto_model(
    dt,
    include_sigma =
      TRUE
  )

  fit_null <- fit_cyto_model(
    dt,
    include_sigma =
      FALSE
  )


  T_obs <- max(
    0,
    fit_null$rss -
      fit_full$rss
  )


  sigma_hat <- if (
    "sigma_c" %in%
      names(
        fit_full$coef
      )
  ) {
    unname(
      fit_full$coef[
        "sigma_c"
      ]
    )
  } else {
    NA_real_
  }


  times <- sort(
    unique(
      dt$time
    )
  )


  first_time <- min(
    times
  )


  initial_state <- c(
    C =
      mean(
        dt[
          time ==
            first_time,
          C
        ],
        na.rm =
          TRUE
      ),

    C_s =
      mean(
        dt[
          time ==
            first_time,
          C_s
        ],
        na.rm =
          TRUE
      )
  )


  null_mean <- propagate_null_mean(
    times =
      times,
    initial_state =
      initial_state,
    coef_null =
      fit_null$coef
  )


  # Residual pairs around observed time-specific means.
  obs_means <- dt[
    ,
    .(
      C_mean =
        mean(
          C,
          na.rm =
            TRUE
        ),
      C_s_mean =
        mean(
          C_s,
          na.rm =
            TRUE
        )
    ),
    by =
      time
  ]


  resid_dt <- merge(
    dt[
      ,
      .(
        time,
        replicate,
        C,
        C_s
      )
    ],
    obs_means,
    by =
      "time",
    all.x =
      TRUE
  )


  resid_dt[
    ,
    `:=`(
      resid_C =
        C -
          C_mean,
      resid_Cs =
        C_s -
          C_s_mean
    )
  ]


  set.seed(
    seed
  )


  T_boot <- numeric(
    B
  )

  failures <- 0L


  for (
    bb in seq_len(
      B
    )
  ) {

    boot_rows <- list()
    rr <- 1L


    for (
      tt in times
    ) {

      res_t <- resid_dt[
        time ==
          tt
      ]

      n_rep_t <- nrow(
        res_t
      )

      if (
        n_rep_t <
          1L
      ) {
        stop(
          "No residuals available at a sampled time point."
        )
      }


      idx <- sample(
        seq_len(
          n_rep_t
        ),
        size =
          n_rep_t,
        replace =
          TRUE
      )


      mu_t <- null_mean[
        time ==
          tt
      ]


      C_boot <- mu_t$C +
        res_t$resid_C[
          idx
        ]

      Cs_boot <- mu_t$C_s +
        res_t$resid_Cs[
          idx
        ]


      if (
        TRUNCATE_NONNEGATIVE_BOOT
      ) {

        C_boot <- pmax(
          0,
          C_boot
        )

        Cs_boot <- pmax(
          0,
          Cs_boot
        )
      }


      boot_rows[[rr]] <- data.table(
        time =
          tt,
        replicate =
          seq_len(
            n_rep_t
          ),
        C =
          C_boot,
        C_s =
          Cs_boot
      )

      rr <- rr + 1L
    }


    boot_dt <- rbindlist(
      boot_rows
    )


    one <- tryCatch(
      {

        ff <- fit_cyto_model(
          boot_dt,
          include_sigma =
            TRUE
        )

        fn <- fit_cyto_model(
          boot_dt,
          include_sigma =
            FALSE
        )

        max(
          0,
          fn$rss -
            ff$rss
        )

      },
      error =
        function(e) NA_real_
    )


    if (
      !is.finite(
        one
      )
    ) {
      failures <- failures + 1L
      T_boot[bb] <- NA_real_
    } else {
      T_boot[bb] <- one
    }
  }


  valid_boot <- T_boot[
    is.finite(
      T_boot
    )
  ]


  if (
    length(
      valid_boot
    ) <
      max(
        20L,
        floor(
          0.8 *
            B
        )
      )
  ) {

    return(
      list(
        p.value =
          NA_real_,
        T.obs =
          T_obs,
        sigma_c_hat =
          sigma_hat,
        rss_full =
          fit_full$rss,
        rss_null =
          fit_null$rss,
        bootstrap_failure_rate =
          failures /
            B,
        status =
          "bootstrap_failure"
      )
    )
  }


  p_value <- (
    1 +
      sum(
        valid_boot >=
          T_obs
      )
  ) /
    (
      1 +
        length(
          valid_boot
        )
    )


  list(
    p.value =
      p_value,
    T.obs =
      T_obs,
    sigma_c_hat =
      sigma_hat,
    rss_full =
      fit_full$rss,
    rss_null =
      fit_null$rss,
    bootstrap_failure_rate =
      failures /
        B,
    status =
      "ok"
  )
}


# =============================================================================
# 13. Run one scenario
# =============================================================================

run_one_scenario <- function(
  gene_i,
  design_i,
  platform_i,
  noise_i,
  B
) {

  latent <- ode_sub[
    Gene ==
      gene_i &
      Design ==
        design_i
  ]


  if (
    nrow(
      latent
    ) ==
      0L
  ) {
    stop(
      paste(
        "Missing latent data for gene",
        gene_i,
        "design",
        design_i
      )
    )
  }


  config_row <- DESIGNS[
    Design ==
      design_i
  ]


  full_match <- full_raw[
    Gene ==
      gene_i &
      Perturbation ==
        "SHUTOFF" &
      Platform ==
        platform_i &
      Exprs_noise ==
        noise_i &
      N_tsamples ==
        config_row$N_tsamples &
      N_replicates ==
        config_row$N_replicates &
      Tsteps ==
        config_row$Tsteps
  ]


  if (
    nrow(
      full_match
    ) !=
      1L
  ) {
    stop(
      paste(
        "Expected one full-benchmark row for gene",
        gene_i,
        design_i,
        platform_i,
        noise_i,
        "but found",
        nrow(
          full_match
        )
      )
    )
  }


  set.seed(
    as.integer(
      full_match$measurement_seed
    )
  )


  noisy <- add_platform_noise_main(
    dt =
      latent,
    platform =
      platform_i,
    noise =
      noise_i
  )


  boot_seed <- make_boot_seed(
    gene =
      gene_i,
    design =
      design_i,
    platform =
      platform_i,
    noise =
      noise_i
  )


  t0 <- Sys.time()


  test <- tryCatch(
    bootstrap_cyto_test(
      dt =
        noisy[
          ,
          .(
            time,
            replicate,
            C,
            C_s
          )
        ],
      B =
        B,
      seed =
        boot_seed
    ),
    error =
      function(e) {

        list(
          p.value =
            NA_real_,
          T.obs =
            NA_real_,
          sigma_c_hat =
            NA_real_,
          rss_full =
            NA_real_,
          rss_null =
            NA_real_,
          bootstrap_failure_rate =
            NA_real_,
          status =
            "test_error",
          error_message =
            conditionMessage(
              e
            )
        )
      }
  )


  data.table(
    Benchmark_version =
      BENCHMARK_VERSION,

    Gene =
      gene_i,

    Positive =
      full_match$Positive,

    Design =
      design_i,

    N_tsamples =
      config_row$N_tsamples,

    N_replicates =
      config_row$N_replicates,

    Tsteps =
      config_row$Tsteps,

    Platform =
      platform_i,

    Exprs_noise =
      noise_i,

    measurement_seed =
      full_match$measurement_seed,

    full_model_p_value =
      full_match$p.value,

    cyto_p_value =
      test$p.value,

    cyto_T_obs =
      test$T.obs,

    cyto_sigma_c_hat =
      test$sigma_c_hat,

    cyto_rss_full =
      test$rss_full,

    cyto_rss_null =
      test$rss_null,

    bootstrap_failure_rate =
      test$bootstrap_failure_rate,

    status =
      test$status,

    error_message =
      if (
        !is.null(
          test$error_message
        )
      ) {
        test$error_message
      } else {
        NA_character_
      },

    runtime_seconds =
      as.numeric(
        difftime(
          Sys.time(),
          t0,
          units =
            "secs"
        )
      )
  )
}


# =============================================================================
# 14. Run one gene
# =============================================================================

run_one_gene <- function(
  gene_i,
  B
) {

  out <- list()
  kk <- 1L


  for (
    design_i in DESIGNS$Design
  ) {

    for (
      platform_i in PLATFORMS
    ) {

      for (
        noise_i in NOISE_LEVELS
      ) {

        out[[kk]] <- run_one_scenario(
          gene_i =
            gene_i,
          design_i =
            design_i,
          platform_i =
            platform_i,
          noise_i =
            noise_i,
          B =
            B
        )

        kk <- kk + 1L
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
# 15. Checkpoints
# =============================================================================

EXPECTED_TESTS_PER_GENE <- (
  nrow(
    DESIGNS
  ) *
    length(
      PLATFORMS
    ) *
    length(
      NOISE_LEVELS
    )
)


checkpoint_file <- function(
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

  TRUE
}


run_gene_and_checkpoint <- function(
  gene_i,
  B
) {

  ff <- checkpoint_file(
    gene_i
  )


  if (
    file.exists(
      ff
    )
  ) {

    old <- tryCatch(
      readRDS(
        ff
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


  ans <- tryCatch(
    run_one_gene(
      gene_i =
        gene_i,
      B =
        B
    ),
    error =
      function(e) {

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
      ans,
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
          ans$error_message
      )
    )
  }


  if (
    !valid_gene_result(
      ans,
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
            ans
          ),
        error_message =
          "Gene result failed integrity validation."
      )
    )
  }


  tmp <- paste0(
    ff,
    ".tmp.",
    Sys.getpid()
  )


  saveRDS(
    ans,
    tmp,
    compress =
      "gzip"
  )


  if (
    !file.rename(
      tmp,
      ff
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
            ans
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
        ans
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
# 16. Smoke test
# =============================================================================

SMOKE_BOOT <- 19L
smoke_gene <- GENES_KEEP[1]


cat(
  "\n============================================================\n",
  "CYTOPLASMIC-ONLY BASELINE SMOKE TEST\n",
  "============================================================\n",
  "Gene: ",
  smoke_gene,
  "\n",
  "Expected tests/gene: ",
  EXPECTED_TESTS_PER_GENE,
  "\n",
  sep = ""
)


smoke <- run_one_gene(
  gene_i =
    smoke_gene,
  B =
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


if (
  any(
    smoke$status ==
      "test_error"
  )
) {

  cat(
    "\nSmoke errors:\n"
  )

  print(
    unique(
      smoke[
        status ==
          "test_error",
        error_message
      ]
    )
  )

  stop(
    "Smoke test failed."
  )
}


cat(
  "\nSMOKE TEST PASSED.\n"
)


rm(
  smoke
)

gc()


# =============================================================================
# 17. Determine remaining genes
# =============================================================================

done <- vapply(
  GENES_KEEP,
  function(
    gene_i
  ) {

    ff <- checkpoint_file(
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
# 18. Parallel run
# =============================================================================

if (
  length(
    genes_remaining
  ) >
    0L
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

          source("../commons/platforms.r")

          NULL
        }
      )


      export_names <- c(
        "ode_sub",
        "full_raw",
        "DESIGNS",
        "PLATFORMS",
        "NOISE_LEVELS",
        "RANGE_GAUSS_NOISE",
        "BENCHMARK_VERSION",
        "REL_FLOOR",
        "TRUNCATE_NONNEGATIVE_BOOT",
        "CHECKPOINT_DIR",
        "EXPECTED_TESTS_PER_GENE",
        "add_platform_noise_main",
        "stable_seed_from_string",
        "make_boot_seed",
        "build_cyto_design",
        "fit_cyto_model",
        "propagate_null_mean",
        "bootstrap_cyto_test",
        "run_one_scenario",
        "run_one_gene",
        "checkpoint_file",
        "valid_gene_result",
        "run_gene_and_checkpoint",
        "add_gaussian_noise",
        "simulate_rt_qpcr",
        "simulate_rnaseq",
        "sample_dispersion_gamma"
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
          B =
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
          ) >
            0L
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
# 19. Collect checkpoints
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

  ff <- checkpoint_file(
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
  ) ==
    0L
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
  Exprs_noise
)


# =============================================================================
# 20. Wilson interval helper
# =============================================================================

wilson_interval <- function(
  successes,
  n,
  conf =
    0.95
) {

  if (
    n <=
      0L
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
# 21. Summary by design/platform/noise
# =============================================================================

summary_dt <- results[
  status ==
    "ok",
  {

    null <- .SD[
      Positive ==
        0
    ]

    alt <- .SD[
      Positive ==
        1
    ]

    n0 <- nrow(
      null
    )

    n1 <- nrow(
      alt
    )

    k0_cyto <- sum(
      null$cyto_p_value <=
        0.05,
      na.rm =
        TRUE
    )

    k0_full <- sum(
      null$full_model_p_value <=
        0.05,
      na.rm =
        TRUE
    )

    ci_cyto <- wilson_interval(
      k0_cyto,
      n0
    )

    ci_full <- wilson_interval(
      k0_full,
      n0
    )

    .(
      N_null =
        n0,

      N_alt =
        n1,

      Cyto_TypeI_005 =
        k0_cyto /
          n0,

      Cyto_Wilson_low =
        ci_cyto[
          "low"
        ],

      Cyto_Wilson_high =
        ci_cyto[
          "high"
        ],

      Cyto_CI_contains_005 =
        ci_cyto[
          "low"
        ] <=
          0.05 &
        ci_cyto[
          "high"
        ] >=
          0.05,

      Cyto_Power_005 =
        mean(
          alt$cyto_p_value <=
            0.05,
          na.rm =
            TRUE
        ),

      Full_TypeI_005 =
        k0_full /
          n0,

      Full_Wilson_low =
        ci_full[
          "low"
        ],

      Full_Wilson_high =
        ci_full[
          "high"
        ],

      Full_CI_contains_005 =
        ci_full[
          "low"
        ] <=
          0.05 &
        ci_full[
          "high"
        ] >=
          0.05,

      Full_Power_005 =
        mean(
          alt$full_model_p_value <=
            0.05,
          na.rm =
            TRUE
        ),

      Cyto_sigma_boundary_null =
        mean(
          null$cyto_sigma_c_hat <=
            1e-12,
          na.rm =
            TRUE
        ),

      Cyto_sigma_boundary_alt =
        mean(
          alt$cyto_sigma_c_hat <=
            1e-12,
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
    Exprs_noise
  )
]


setorder(
  summary_dt,
  Design,
  Platform,
  Exprs_noise
)


# =============================================================================
# 22. Comparison restricted to scenarios calibrated for each method
# =============================================================================

comparison_dt <- summary_dt[
  ,
  .(
    Design,
    Platform,
    Exprs_noise,

    Cyto_TypeI_005,
    Cyto_CI_contains_005,
    Cyto_Power_005,

    Full_TypeI_005,
    Full_CI_contains_005,
    Full_Power_005,

    Power_difference_full_minus_cyto =
      Full_Power_005 -
        Cyto_Power_005
  )
]


# =============================================================================
# 23. Save
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
  comparison_dt,
  COMPARISON_FILE,
  sep =
    "\t"
)


cytoplasmic_only_baseline <- list(
  settings =
    list(
      benchmark_version =
        BENCHMARK_VERSION,
      n_boot =
        N_BOOT,
      n_null =
        N_NULL,
      n_alt =
        N_ALT,
      designs =
        DESIGNS,
      platforms =
        PLATFORMS,
      noise_levels =
        NOISE_LEVELS
    ),

  selected_null_genes =
    NULL_GENES,

  selected_alt_genes =
    ALT_GENES,

  raw =
    results,

  summary =
    summary_dt,

  comparison =
    comparison_dt
)


save(
  cytoplasmic_only_baseline,
  file =
    RDATA_FILE
)


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
  "CYTOPLASMIC-ONLY BASELINE COMPLETE\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nSummary:\n"
)

print(
  summary_dt
)


cat(
  "\nFull model minus cytoplasmic-only power:\n"
)

print(
  comparison_dt
)


cat(
  "\nGenerated in:\n  ",
  OUTPUT_DIR,
  "\n",
  sep = ""
)
