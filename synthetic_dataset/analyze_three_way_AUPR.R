# =============================================================================
# Three-way PR-AUC benchmark:
#   Full compartment model vs Cytoplasmic-only kinetic model vs Delta-PSI
#
# Purpose
# -------
# Compare ranking/discrimination performance using PR-AUC in the imbalanced
# synthetic benchmark while keeping null calibration as a separate diagnostic.
#
# The three comparators are:
#   1. Full compartment-resolved kinetic model
#   2. Cytoplasmic-only kinetic model
#   3. Delta-PSI non-kinetic baseline
#
# For a fair three-way ranking comparison, all three methods are ranked by the
# SAME type of evidence score:
#
#     score = -log10(p)
#
# Thus this analysis compares the discriminatory content of the statistical
# evidence itself and does not mix method-specific composite ranking scores.
#
# Delta-PSI baseline
# ------------------
# For each gene and replicate:
#   PSI(t) = C(t) / (C(t) + C_s(t) + eps)
#
# The early PSI is the mean over the first two sampled time points and the late
# PSI is the mean over the last two sampled time points. A paired one-sample
# t-test is then performed on replicate-level differences in logit(PSI):
#
#   logit(PSI_late) - logit(PSI_early)
#
# This implements explicitly the "aggregated early versus late windows"
# definition used for the PSI baseline in the manuscript.
#
# IMPORTANT:
# If the historical PSI code used a different early/late-window rule, change
# N_EARLY_TIMES and N_LATE_TIMES below before interpreting the comparison.
#
# Inputs
# ------
# ode_states_2k_20p_corrected_onset.rdata
# cytoplasmic_only_baseline/cytoplasmic_only_raw.tsv
# cytoplasmic_only_baseline/cytoplasmic_only_summary.tsv
#
# Outputs
# -------
# cytoplasmic_only_baseline/
#   three_way_AUPR_raw.tsv
#   three_way_AUPR_summary.tsv
#   three_way_AUPR_both_kinetic_calibrated.tsv
#   three_way_AUPR_overview.tsv
#   FigS_three_way_AUPR.pdf
#   FigS_three_way_AUPR.png
#   three_way_AUPR.rdata
#
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/synthetic_dataset")

library(data.table)
library(ggplot2)

source("../commons/platforms.r")


# =============================================================================
# 1. Files
# =============================================================================

ODE_FILE <- "ode_states_2k_20p_corrected_onset.rdata"

CYTO_RAW_FILE <- file.path(
  "cytoplasmic_only_baseline",
  "cytoplasmic_only_raw.tsv"
)

CYTO_SUMMARY_FILE <- file.path(
  "cytoplasmic_only_baseline",
  "cytoplasmic_only_summary.tsv"
)

OUT_DIR <- "cytoplasmic_only_baseline"

OUT_RAW <- file.path(
  OUT_DIR,
  "three_way_AUPR_raw.tsv"
)

OUT_SUMMARY <- file.path(
  OUT_DIR,
  "three_way_AUPR_summary.tsv"
)

OUT_CAL <- file.path(
  OUT_DIR,
  "three_way_AUPR_both_kinetic_calibrated.tsv"
)

OUT_OVERVIEW <- file.path(
  OUT_DIR,
  "three_way_AUPR_overview.tsv"
)

OUT_PDF <- file.path(
  OUT_DIR,
  "FigS_three_way_AUPR.pdf"
)

OUT_PNG <- file.path(
  OUT_DIR,
  "FigS_three_way_AUPR.png"
)

OUT_RDATA <- file.path(
  OUT_DIR,
  "three_way_AUPR.rdata"
)


# =============================================================================
# 2. Settings
# =============================================================================

N_BOOT_AUPR <- 2000L

BASE_SEED <- 20260923L

PSI_EPS <- 1e-6

N_EARLY_TIMES <- 2L
N_LATE_TIMES <- 2L

P_FLOOR <- 1e-300

# Parallel/cache settings. On Linux, mclapply uses forked workers and avoids
# repeatedly copying the large ode_states object.
N_WORKERS <- min(60L, max(1L, parallel::detectCores(logical = TRUE) - 2L))
USE_CACHE <- TRUE



NOISE_LEVELS <- c(
  "Very low",
  "Low",
  "Medium",
  "High"
)

PLATFORM_LEVELS <- c(
  "RT-qPCR",
  "GAUSS",
  "RNA-seq"
)

DESIGN_LEVELS <- c(
  "5tp_3rep_dt10",
  "5tp_10rep_dt10"
)

DESIGN_LABELS <- c(
  "5tp_3rep_dt10" =
    "5 time points, 3 replicates",
  "5tp_10rep_dt10" =
    "5 time points, 10 replicates"
)

RANGE_GAUSS_NOISE <- c(
  "Very low" = 0.02,
  "Low"      = 0.05,
  "Medium"   = 0.10,
  "High"     = 0.20
)


# =============================================================================
# 3. Load
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

cyto_raw <- fread(
  CYTO_RAW_FILE
)

cyto_summary <- fread(
  CYTO_SUMMARY_FILE
)


# =============================================================================
# 4. Validate
# =============================================================================

required_raw <- c(
  "Gene",
  "Positive",
  "Design",
  "N_tsamples",
  "N_replicates",
  "Tsteps",
  "Platform",
  "Exprs_noise",
  "measurement_seed",
  "status",
  "full_model_p_value",
  "cyto_p_value"
)

required_summary <- c(
  "Design",
  "Platform",
  "Exprs_noise",
  "Cyto_CI_contains_005",
  "Full_CI_contains_005",
  "Cyto_TypeI_005",
  "Full_TypeI_005"
)

required_ode <- c(
  "Gene",
  "Perturbation",
  "N_time_samples",
  "N_replicates",
  "T_step",
  "Post_R_fraction",
  "time",
  "replicate",
  "N",
  "N_s",
  "C",
  "C_s"
)

missing_raw <- setdiff(
  required_raw,
  names(
    cyto_raw
  )
)

missing_summary <- setdiff(
  required_summary,
  names(
    cyto_summary
  )
)

missing_ode <- setdiff(
  required_ode,
  names(
    ode_states
  )
)

if (
  length(
    missing_raw
  ) > 0L
) {
  stop(
    paste0(
      "Missing columns in cytoplasmic_only_raw.tsv:\n  ",
      paste(
        missing_raw,
        collapse = "\n  "
      )
    )
  )
}

if (
  length(
    missing_summary
  ) > 0L
) {
  stop(
    paste0(
      "Missing columns in cytoplasmic_only_summary.tsv:\n  ",
      paste(
        missing_summary,
        collapse = "\n  "
      )
    )
  )
}

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


# =============================================================================
# 5. Reproduce benchmark observation noise
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
# 6. Delta-PSI p-value
# =============================================================================

logit_safe <- function(
  x,
  eps =
    PSI_EPS
) {

  x <- pmin(
    1 -
      eps,
    pmax(
      eps,
      x
    )
  )

  log(
    x /
      (
        1 -
          x
      )
  )
}


compute_delta_psi_test <- function(
  noisy_dt,
  n_early =
    N_EARLY_TIMES,
  n_late =
    N_LATE_TIMES
) {

  tt <- sort(
    unique(
      noisy_dt$time
    )
  )

  if (
    length(
      tt
    ) <
      (
        n_early +
          n_late
      )
  ) {
    stop(
      "Not enough sampled time points for the requested early/late PSI windows."
    )
  }

  early_times <- head(
    tt,
    n_early
  )

  late_times <- tail(
    tt,
    n_late
  )


  d <- copy(
    noisy_dt
  )


  d[
    ,
    PSI :=
      C /
        (
          C +
            C_s +
            PSI_EPS
        )
  ]


  early <- d[
    time %in%
      early_times,
    .(
      PSI_early =
        mean(
          PSI,
          na.rm =
            TRUE
        )
    ),
    by =
      replicate
  ]


  late <- d[
    time %in%
      late_times,
    .(
      PSI_late =
        mean(
          PSI,
          na.rm =
            TRUE
        )
    ),
    by =
      replicate
  ]


  paired <- merge(
    early,
    late,
    by =
      "replicate",
    all =
      FALSE
  )


  paired[
    ,
    logit_difference :=
      logit_safe(
        PSI_late
      ) -
      logit_safe(
        PSI_early
      )
  ]


  delta_psi <- mean(
    paired$PSI_late -
      paired$PSI_early,
    na.rm =
      TRUE
  )


  if (
    nrow(
      paired
    ) <
      2L ||
      !is.finite(
        sd(
          paired$logit_difference
        )
      )
  ) {

    return(
      list(
        p.value =
          1,
        delta_psi =
          delta_psi
      )
    )
  }


  if (
    sd(
      paired$logit_difference
    ) <
      1e-14
  ) {

    p_value <- if (
      abs(
        mean(
          paired$logit_difference
        )
      ) >
        1e-14
    ) {
      0
    } else {
      1
    }

  } else {

    p_value <- t.test(
      paired$logit_difference,
      mu =
        0,
      alternative =
        "two.sided"
    )$p.value
  }


  list(
    p.value =
      p_value,
    delta_psi =
      delta_psi
  )
}


# =============================================================================
# 7. Reconstruct Delta-PSI results for every paired benchmark row
#
# The expensive reconstruction is cached in OUT_RAW. If that file already
# exists and contains the expected columns/number of rows, it is reused.
# Otherwise the reconstruction is parallelized across forked Linux workers.
# =============================================================================

cache_raw_ok <- FALSE

if (
  USE_CACHE &&
    file.exists(OUT_RAW)
) {

  cached_raw <- tryCatch(
    fread(OUT_RAW),
    error = function(e) NULL
  )

  cache_raw_ok <- !is.null(cached_raw) &&
    nrow(cached_raw) == nrow(cyto_raw) &&
    all(
      c(
        "psi_p_value",
        "delta_psi",
        "Full_score",
        "Cyto_score",
        "PSI_score"
      ) %in% names(cached_raw)
    )
}

if (cache_raw_ok) {

  cat(
    "\n============================================================\n",
    "USING CACHED DELTA-PSI RECONSTRUCTION\n",
    "============================================================\n",
    "Rows loaded: ",
    nrow(cached_raw),
    "\n",
    sep = ""
  )

  analysis_raw <- cached_raw
  rm(cached_raw)

} else {

  cat(
    "\n============================================================\n",
    "RECONSTRUCTING DELTA-PSI BASELINE IN PARALLEL\n",
    "============================================================\n",
    "Workers: ",
    N_WORKERS,
    "\n",
    sep = ""
  )

  # Keying materially accelerates repeated latent-trajectory lookup.
  setkeyv(
    ode_states,
    c(
      "Gene",
      "Perturbation",
      "N_time_samples",
      "N_replicates",
      "T_step"
    )
  )

  reconstruct_one_psi <- function(ii) {

    row_i <- cyto_raw[ii]

    latent <- ode_states[
      .(
        row_i$Gene,
        "SHUTOFF",
        row_i$N_tsamples,
        row_i$N_replicates,
        row_i$Tsteps
      )
    ][
      abs(Post_R_fraction) < 1e-12
    ]

    if (nrow(latent) == 0L) {
      stop(
        paste(
          "Missing latent data for gene",
          row_i$Gene,
          row_i$Design
        )
      )
    }

    set.seed(
      as.integer(row_i$measurement_seed)
    )

    noisy <- add_platform_noise_main(
      dt = latent,
      platform = row_i$Platform,
      noise = row_i$Exprs_noise
    )

    psi_test <- compute_delta_psi_test(
      noisy_dt = noisy
    )

    data.table(
      Gene = row_i$Gene,
      Design = row_i$Design,
      Platform = row_i$Platform,
      Exprs_noise = row_i$Exprs_noise,
      psi_p_value = psi_test$p.value,
      delta_psi = psi_test$delta_psi
    )
  }

  # Use a small number of large chunks rather than 18,000 individual fork jobs.
  n_chunks <- min(
    N_WORKERS,
    nrow(cyto_raw)
  )

  chunk_id <- cut(
    seq_len(nrow(cyto_raw)),
    breaks = n_chunks,
    labels = FALSE
  )

  index_chunks <- split(
    seq_len(nrow(cyto_raw)),
    chunk_id
  )

  worker_chunk <- function(idx) {
    rbindlist(
      lapply(
        idx,
        reconstruct_one_psi
      ),
      use.names = TRUE,
      fill = TRUE
    )
  }

  if (.Platform$OS.type == "unix" && N_WORKERS > 1L) {

    psi_parts <- parallel::mclapply(
      index_chunks,
      worker_chunk,
      mc.cores = N_WORKERS,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )

  } else {

    # Portable fallback. The user's cluster is Linux, so this branch is not
    # normally used there.
    psi_parts <- lapply(
      index_chunks,
      worker_chunk
    )
  }

  psi_dt <- rbindlist(
    psi_parts,
    use.names = TRUE,
    fill = TRUE
  )

  if (nrow(psi_dt) != nrow(cyto_raw)) {
    stop(
      paste(
        "Delta-PSI reconstruction returned",
        nrow(psi_dt),
        "rows; expected",
        nrow(cyto_raw)
      )
    )
  }

  analysis_raw <- merge(
    cyto_raw,
    psi_dt,
    by = c(
      "Gene",
      "Design",
      "Platform",
      "Exprs_noise"
    ),
    all.x = TRUE
  )

  analysis_raw[
    ,
    `:=`(
      Full_score = -log10(
        pmax(full_model_p_value, P_FLOOR)
      ),
      Cyto_score = -log10(
        pmax(cyto_p_value, P_FLOOR)
      ),
      PSI_score = -log10(
        pmax(psi_p_value, P_FLOOR)
      )
    )
  ]

  # Save immediately so an error in a later plotting/reporting stage can never
  # force the expensive reconstruction to be repeated.
  fwrite(
    analysis_raw,
    OUT_RAW,
    sep = "\t"
  )

  cat(
    "Delta-PSI reconstruction completed and cached in:\n  ",
    OUT_RAW,
    "\n",
    sep = ""
  )
}


# =============================================================================
# 8. PR-AUC helper
#
# Thresholds are evaluated at unique score values to avoid arbitrary ordering
# among tied scores. Trapezoidal integration is then used over recall.
# =============================================================================

pr_auc <- function(
  truth,
  score
) {

  good <- !is.na(
    truth
  ) &
    is.finite(
      score
    )

  truth <- as.integer(
    truth[
      good
    ]
  )

  score <- score[
    good
  ]


  n_pos <- sum(
    truth ==
      1L
  )

  if (
    n_pos <
      1L
  ) {
    return(
      NA_real_
    )
  }


  thresholds <- sort(
    unique(
      score
    ),
    decreasing =
      TRUE
  )


  recall <- numeric(
    length(
      thresholds
    )
  )

  precision <- numeric(
    length(
      thresholds
    )
  )


  for (
    ii in seq_along(
      thresholds
    )
  ) {

    pred <- score >=
      thresholds[ii]

    tp <- sum(
      pred &
        truth ==
          1L
    )

    fp <- sum(
      pred &
        truth ==
          0L
    )

    recall[ii] <- tp /
      n_pos

    precision[ii] <- if (
      tp +
        fp >
        0L
    ) {
      tp /
        (
          tp +
            fp
        )
    } else {
      1
    }
  }


  recall <- c(
    0,
    recall
  )

  precision <- c(
    1,
    precision
  )


  ord <- order(
    recall,
    precision
  )

  recall <- recall[
    ord
  ]

  precision <- precision[
    ord
  ]


  sum(
    diff(
      recall
    ) *
      (
        head(
          precision,
          -1L
        ) +
        tail(
          precision,
          -1L
        )
      ) /
      2
  )
}


# =============================================================================
# 9. Stable bootstrap seed
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
      out <=
        0L
  ) {
    out <- 1L
  }

  out
}


# =============================================================================
# 10. One paired three-way PR-AUC comparison
# =============================================================================

analyze_one_configuration <- function(
  dt,
  B,
  seed
) {

  dt <- copy(
    dt
  )


  dt <- dt[
    status ==
      "ok" &
      is.finite(
        Full_score
      ) &
      is.finite(
        Cyto_score
      ) &
      is.finite(
        PSI_score
      )
  ]


  null_dt <- dt[
    Positive ==
      0
  ]

  alt_dt <- dt[
    Positive ==
      1
  ]


  if (
    nrow(
      null_dt
    ) <
      2L ||
      nrow(
        alt_dt
      ) <
        2L
  ) {
    stop(
      "Not enough paired null/alternative observations for PR-AUC."
    )
  }


  auc_full <- pr_auc(
    dt$Positive,
    dt$Full_score
  )

  auc_cyto <- pr_auc(
    dt$Positive,
    dt$Cyto_score
  )

  auc_psi <- pr_auc(
    dt$Positive,
    dt$PSI_score
  )


  set.seed(
    seed
  )


  boot_full <- numeric(
    B
  )

  boot_cyto <- numeric(
    B
  )

  boot_psi <- numeric(
    B
  )


  for (
    bb in seq_len(
      B
    )
  ) {

    idx0 <- sample(
      seq_len(
        nrow(
          null_dt
        )
      ),
      nrow(
        null_dt
      ),
      replace =
        TRUE
    )

    idx1 <- sample(
      seq_len(
        nrow(
          alt_dt
        )
      ),
      nrow(
        alt_dt
      ),
      replace =
        TRUE
    )


    boot_dt <- rbind(
      null_dt[
        idx0
      ],
      alt_dt[
        idx1
      ]
    )


    boot_full[bb] <- pr_auc(
      boot_dt$Positive,
      boot_dt$Full_score
    )

    boot_cyto[bb] <- pr_auc(
      boot_dt$Positive,
      boot_dt$Cyto_score
    )

    boot_psi[bb] <- pr_auc(
      boot_dt$Positive,
      boot_dt$PSI_score
    )
  }


  diff_full_cyto <- boot_full -
    boot_cyto

  diff_full_psi <- boot_full -
    boot_psi

  diff_cyto_psi <- boot_cyto -
    boot_psi


  data.table(
    N_null =
      nrow(
        null_dt
      ),

    N_alt =
      nrow(
        alt_dt
      ),

    Positive_prevalence =
      nrow(
        alt_dt
      ) /
      (
        nrow(
          null_dt
        ) +
          nrow(
            alt_dt
          )
      ),

    Full_AUPR =
      auc_full,

    Full_AUPR_low =
      quantile(
        boot_full,
        0.025,
        names =
          FALSE
      ),

    Full_AUPR_high =
      quantile(
        boot_full,
        0.975,
        names =
          FALSE
      ),

    Cyto_AUPR =
      auc_cyto,

    Cyto_AUPR_low =
      quantile(
        boot_cyto,
        0.025,
        names =
          FALSE
      ),

    Cyto_AUPR_high =
      quantile(
        boot_cyto,
        0.975,
        names =
          FALSE
      ),

    PSI_AUPR =
      auc_psi,

    PSI_AUPR_low =
      quantile(
        boot_psi,
        0.025,
        names =
          FALSE
      ),

    PSI_AUPR_high =
      quantile(
        boot_psi,
        0.975,
        names =
          FALSE
      ),

    Full_minus_Cyto =
      auc_full -
        auc_cyto,

    Full_minus_Cyto_low =
      quantile(
        diff_full_cyto,
        0.025,
        names =
          FALSE
      ),

    Full_minus_Cyto_high =
      quantile(
        diff_full_cyto,
        0.975,
        names =
          FALSE
      ),

    Full_minus_PSI =
      auc_full -
        auc_psi,

    Full_minus_PSI_low =
      quantile(
        diff_full_psi,
        0.025,
        names =
          FALSE
      ),

    Full_minus_PSI_high =
      quantile(
        diff_full_psi,
        0.975,
        names =
          FALSE
      ),

    Cyto_minus_PSI =
      auc_cyto -
        auc_psi,

    Cyto_minus_PSI_low =
      quantile(
        diff_cyto_psi,
        0.025,
        names =
          FALSE
      ),

    Cyto_minus_PSI_high =
      quantile(
        diff_cyto_psi,
        0.975,
        names =
          FALSE
      )
  )
}


# =============================================================================
# 11. Run all configurations (cached or parallel)
# =============================================================================

summary_cache_ok <- FALSE

if (
  USE_CACHE &&
    file.exists(OUT_SUMMARY) &&
    file.exists(OUT_CAL) &&
    file.exists(OUT_OVERVIEW)
) {

  cached_summary <- tryCatch(
    fread(OUT_SUMMARY),
    error = function(e) NULL
  )

  cached_cal <- tryCatch(
    fread(OUT_CAL),
    error = function(e) NULL
  )

  cached_overview <- tryCatch(
    fread(OUT_OVERVIEW),
    error = function(e) NULL
  )

  summary_cache_ok <- !is.null(cached_summary) &&
    nrow(cached_summary) == 24L &&
    all(
      c(
        "Full_AUPR",
        "Cyto_AUPR",
        "PSI_AUPR",
        "Positive_prevalence",
        "Both_kinetic_calibrated"
      ) %in% names(cached_summary)
    ) &&
    !is.null(cached_cal) &&
    !is.null(cached_overview)
}

if (summary_cache_ok) {

  cat(
    "\n============================================================\n",
    "USING CACHED PR-AUC RESULTS\n",
    "============================================================\n",
    sep = ""
  )

  aupr_dt <- cached_summary
  both_calibrated <- cached_cal
  overview <- cached_overview

  rm(
    cached_summary,
    cached_cal,
    cached_overview
  )

} else {

  config_cols <- c(
    "Design",
    "N_tsamples",
    "N_replicates",
    "Tsteps",
    "Platform",
    "Exprs_noise"
  )

  configs <- unique(
    analysis_raw[
      ,
      ..config_cols
    ]
  )

  setorder(
    configs,
    Design,
    Platform,
    Exprs_noise
  )

  analyze_config_index <- function(ii) {

    cfg <- configs[ii]

    dt_i <- analysis_raw[
      Design == cfg$Design &
        N_tsamples == cfg$N_tsamples &
        N_replicates == cfg$N_replicates &
        Tsteps == cfg$Tsteps &
        Platform == cfg$Platform &
        Exprs_noise == cfg$Exprs_noise
    ]

    seed_i <- stable_seed_from_string(
      paste(
        BASE_SEED,
        cfg$Design,
        cfg$Platform,
        cfg$Exprs_noise,
        "AUPR",
        sep = "|"
      )
    )

    one <- analyze_one_configuration(
      dt = dt_i,
      B = N_BOOT_AUPR,
      seed = seed_i
    )

    cbind(
      cfg,
      one
    )
  }

  workers_auc <- min(
    N_WORKERS,
    nrow(configs)
  )

  cat(
    "\nComputing PR-AUC for ",
    nrow(configs),
    " configurations using ",
    workers_auc,
    " workers...\n",
    sep = ""
  )

  if (.Platform$OS.type == "unix" && workers_auc > 1L) {

    result_list <- parallel::mclapply(
      seq_len(nrow(configs)),
      analyze_config_index,
      mc.cores = workers_auc,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )

  } else {

    result_list <- lapply(
      seq_len(nrow(configs)),
      analyze_config_index
    )
  }

  aupr_dt <- rbindlist(
    result_list,
    use.names = TRUE,
    fill = TRUE
  )


  # ===========================================================================
  # 12. Attach kinetic-model calibration
  # ===========================================================================

  calibration_info <- cyto_summary[
    ,
    .(
      Design,
      Platform,
      Exprs_noise,
      Cyto_TypeI_005,
      Cyto_CI_contains_005,
      Full_TypeI_005,
      Full_CI_contains_005
    )
  ]

  aupr_dt <- merge(
    aupr_dt,
    calibration_info,
    by = c(
      "Design",
      "Platform",
      "Exprs_noise"
    ),
    all.x = TRUE
  )

  aupr_dt[
    ,
    Both_kinetic_calibrated :=
      Cyto_CI_contains_005 &
        Full_CI_contains_005
  ]

  setorder(
    aupr_dt,
    Design,
    Platform,
    Exprs_noise
  )

  both_calibrated <- aupr_dt[
    Both_kinetic_calibrated == TRUE
  ]


  # ===========================================================================
  # 13. Overview
  # ===========================================================================

  overview <- rbindlist(
    list(

      aupr_dt[
        ,
        .(
          Scope = "All configurations",
          N_configurations = .N,
          Full_AUPR_median = median(Full_AUPR),
          Cyto_AUPR_median = median(Cyto_AUPR),
          PSI_AUPR_median = median(PSI_AUPR),
          Full_minus_Cyto_median = median(Full_minus_Cyto),
          Full_minus_PSI_median = median(Full_minus_PSI)
        )
      ],

      both_calibrated[
        ,
        .(
          Scope = "Both kinetic models calibrated",
          N_configurations = .N,
          Full_AUPR_median = if (.N > 0L) median(Full_AUPR) else NA_real_,
          Cyto_AUPR_median = if (.N > 0L) median(Cyto_AUPR) else NA_real_,
          PSI_AUPR_median = if (.N > 0L) median(PSI_AUPR) else NA_real_,
          Full_minus_Cyto_median = if (.N > 0L) median(Full_minus_Cyto) else NA_real_,
          Full_minus_PSI_median = if (.N > 0L) median(Full_minus_PSI) else NA_real_
        )
      ]
    ),
    fill = TRUE
  )


  # ===========================================================================
  # 14. Save tables BEFORE plotting
  # ===========================================================================

  fwrite(
    analysis_raw,
    OUT_RAW,
    sep = "\t"
  )

  fwrite(
    aupr_dt,
    OUT_SUMMARY,
    sep = "\t"
  )

  fwrite(
    both_calibrated,
    OUT_CAL,
    sep = "\t"
  )

  fwrite(
    overview,
    OUT_OVERVIEW,
    sep = "\t"
  )
}


# =============================================================================
# 15. Figure
# =============================================================================

plot_long <- rbindlist(
  list(

    aupr_dt[
      ,
      .(
        Design,
        Platform,
        Exprs_noise,
        Positive_prevalence,
        Method =
          "Full compartment model",
        AUPR =
          Full_AUPR,
        Low =
          Full_AUPR_low,
        High =
          Full_AUPR_high
      )
    ],

    aupr_dt[
      ,
      .(
        Design,
        Platform,
        Exprs_noise,
        Positive_prevalence,
        Method =
          "Cytoplasmic-only model",
        AUPR =
          Cyto_AUPR,
        Low =
          Cyto_AUPR_low,
        High =
          Cyto_AUPR_high
      )
    ],

    aupr_dt[
      ,
      .(
        Design,
        Platform,
        Exprs_noise,
        Positive_prevalence,
        Method =
          "Delta-PSI",
        AUPR =
          PSI_AUPR,
        Low =
          PSI_AUPR_low,
        High =
          PSI_AUPR_high
      )
    ]
  ),
  fill =
    TRUE
)


plot_long[
  ,
  Method :=
    factor(
      as.character(
        Method
      ),
      levels = c(
        "Full compartment model",
        "Cytoplasmic-only model",
        "Delta-PSI"
      )
    )
]

plot_long[
  ,
  Exprs_noise :=
    factor(
      Exprs_noise,
      levels =
        NOISE_LEVELS
    )
]

plot_long[
  ,
  Platform :=
    factor(
      Platform,
      levels =
        PLATFORM_LEVELS
    )
]

plot_long[
  ,
  Design :=
    factor(
      Design,
      levels =
        DESIGN_LEVELS
    )
]


random_baseline <- unique(
  plot_long$Positive_prevalence
)

if (
  length(
    random_baseline
  ) !=
    1L
) {
  random_baseline <- median(
    plot_long$Positive_prevalence
  )
}


p <- ggplot(
  plot_long,
  aes(
    x =
      Exprs_noise,
    y =
      AUPR,
    group =
      Method,
    linetype =
      Method,
    shape =
      Method
  )
) +

  geom_hline(
    yintercept =
      random_baseline,
    linetype =
      "dotted",
    linewidth =
      0.5
  ) +

  geom_errorbar(
    aes(
      ymin =
        Low,
      ymax =
        High
    ),
    width =
      0.08,
    position =
      position_dodge(
        width =
          0.22
      ),
    linewidth =
      0.4
  ) +

  geom_line(
    linewidth =
      0.75,
    position =
      position_dodge(
        width =
          0.22
      )
  ) +

  geom_point(
    size =
      2,
    position =
      position_dodge(
        width =
          0.22
      )
  ) +

  facet_grid(
    Platform ~ Design,
    labeller =
      labeller(
        Design =
          DESIGN_LABELS
      )
  ) +

  scale_linetype_manual(
    values = c(
      "Full compartment model" =
        "solid",
      "Cytoplasmic-only model" =
        "longdash",
      "Delta-PSI" =
        "dotdash"
    )
  ) +

  scale_shape_manual(
    values = c(
      "Full compartment model" =
        16,
      "Cytoplasmic-only model" =
        17,
      "Delta-PSI" =
        15
    )
  ) +

  coord_cartesian(
    ylim =
      c(
        0,
        1
      )
  ) +

  labs(
    x =
      "Measurement-noise regime",
    y =
      "PR-AUC",
    linetype =
      NULL,
    shape =
      NULL
  ) +

  theme_classic(
    base_size =
      10.5
  ) +

  theme(
    panel.grid.major.y =
      element_line(
        colour =
          "grey92",
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
    axis.text.x =
      element_text(
        angle =
          30,
        hjust =
          1
      ),
    legend.position =
      "bottom",
    legend.key.width =
      grid::unit(
        1.5,
        "cm"
      )
  )


ggsave(
  filename =
    OUT_PDF,
  plot =
    p,
  width =
    10.5,
  height =
    7.8,
  units =
    "in",
  device =
    cairo_pdf
)


ggsave(
  filename =
    OUT_PNG,
  plot =
    p,
  width =
    10.5,
  height =
    7.8,
  units =
    "in",
  dpi =
    400,
  bg =
    "white"
)


# =============================================================================
# 16. Save R object and report
# =============================================================================

three_way_AUPR <- list(
  settings =
    list(
      n_boot =
        N_BOOT_AUPR,
      psi_epsilon =
        PSI_EPS,
      early_time_points =
        N_EARLY_TIMES,
      late_time_points =
        N_LATE_TIMES,
      score =
        "-log10(p)",
      aupr =
        "threshold-aggregated trapezoidal PR-AUC"
    ),
  raw =
    analysis_raw,
  summary =
    aupr_dt,
  both_kinetic_calibrated =
    both_calibrated,
  overview =
    overview
)


save(
  three_way_AUPR,
  file =
    OUT_RDATA
)


cat(
  "\n============================================================\n",
  "THREE-WAY PR-AUC COMPLETE\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nOverview:\n"
)

print(
  overview
)


cat(
  "\nBoth kinetic models calibrated:\n"
)

if (
  nrow(
    both_calibrated
  ) >
    0L
) {

  print(
    both_calibrated[
      ,
      .(
        Design,
        Platform,
        Exprs_noise,
        Full_TypeI_005,
        Cyto_TypeI_005,
        Full_AUPR,
        Cyto_AUPR,
        PSI_AUPR,
        Full_minus_Cyto,
        Full_minus_Cyto_low,
        Full_minus_Cyto_high,
        Full_minus_PSI,
        Full_minus_PSI_low,
        Full_minus_PSI_high
      )
    ]
  )

} else {

  cat(
    "None.\n"
  )
}


cat(
  "\nGenerated:\n",
  "  ",
  OUT_RAW,
  "\n",
  "  ",
  OUT_SUMMARY,
  "\n",
  "  ",
  OUT_CAL,
  "\n",
  "  ",
  OUT_OVERVIEW,
  "\n",
  "  ",
  OUT_PDF,
  "\n",
  "  ",
  OUT_PNG,
  "\n",
  "  ",
  OUT_RDATA,
  "\n",
  sep = ""
)
