# =============================================================================
# Title:
# Final Analysis of the Corrected-Onset Factorial Synthetic Benchmark
#
# Purpose:
#   Analyze the completed factorial benchmark generated after correcting the
#   synthetic ODE simulator so that:
#
#     - transcriptional onset may vary between genes;
#     - pharmacological shutoff time t_star is common and correctly aligned;
#     - no post-hoc observation-time shift is applied.
#
# This script treats:
#
#   NONE
#     as a reference / low-identifiability regime;
#
#   SHUTOFF
#     as the principal operational inferential regime.
#
# Primary inferential priority:
#   Type-I calibration is evaluated BEFORE power or parameter recovery.
#
# Calibration classification:
#   - "CI includes 0.05":
#       Wilson 95% CI for empirical Type-I contains the nominal 0.05 level.
#
#   - "Anti-conservative":
#       Wilson lower bound > 0.05.
#
#   - "Conservative":
#       Wilson upper bound < 0.05.
#
# A descriptive practical band 0.025-0.075 is also reported, but it is NOT
# treated as an inferential calibration criterion.
#
# Input:
#   benchmark_main_corrected_onset_summary.tsv
#
# Optional input:
#   benchmark_main_corrected_onset_raw.tsv
#
# Outputs:
#   benchmark_corrected_analysis/
#
# Main tables:
#   benchmark_global_summary.tsv
#   benchmark_shutoff_noise_summary.tsv
#   benchmark_design_summary.tsv
#   benchmark_shutoff_operational_domain.tsv
#   benchmark_calibration_classes.tsv
#   benchmark_flagged_scenarios.tsv
#   benchmark_runtime_summary.tsv
#
# Main figures:
#   Fig1_calibration_overview.pdf
#   Fig2_shutoff_design_dependence.pdf
#   Fig3_power_calibration_tradeoff.pdf
#   Fig4_parameter_recovery.pdf
#
# Supplementary:
#   FigS1_typeI_heatmaps_SHUTOFF.pdf
#   FigS2_power_heatmaps_SHUTOFF.pdf
#   FigS3_sigma_c_RMSE_heatmaps_SHUTOFF.pdf
#   FigS4_NONE_reference.pdf
#
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/synthetic_dataset")

library(data.table)
library(ggplot2)
library(patchwork)
library(grid)


# =============================================================================
# 1. Files and analysis settings
# =============================================================================

SUMMARY_FILE <- "benchmark_main_corrected_onset_summary.tsv"
RAW_FILE <- "benchmark_main_corrected_onset_raw.tsv"

OUTPUT_DIR <- "benchmark_corrected_analysis"

NOMINAL_ALPHA <- 0.05

PRACTICAL_LOW <- 0.025
PRACTICAL_HIGH <- 0.075

EXPECTED_PLATFORMS <- c(
  "RT-qPCR",
  "GAUSS",
  "RNA-seq"
)

EXPECTED_PERTURBATIONS <- c(
  "NONE",
  "SHUTOFF"
)

EXPECTED_NOISE <- c(
  "Very low",
  "Low",
  "Medium",
  "High"
)

EXPECTED_N_TSAMPLES <- c(
  3L,
  5L,
  10L,
  20L
)

EXPECTED_N_REPLICATES <- c(
  3L,
  5L,
  10L
)

EXPECTED_TSTEPS <- c(
  5L,
  10L,
  20L,
  50L
)


# =============================================================================
# 2. Output directory
# =============================================================================

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE
  )
}


# =============================================================================
# 3. Load benchmark summary
# =============================================================================

if (!file.exists(SUMMARY_FILE)) {
  stop(
    paste(
      "Summary file not found:",
      SUMMARY_FILE
    )
  )
}

dt <- fread(
  SUMMARY_FILE
)

cat(
  "\n============================================================\n",
  "CORRECTED-ONSET BENCHMARK LOADED\n",
  "============================================================\n",
  "Rows:    ",
  nrow(dt),
  "\n",
  "Columns: ",
  ncol(dt),
  "\n",
  sep = ""
)


# =============================================================================
# 4. Required columns
# =============================================================================

required_columns <- c(
  "Platform",
  "Perturbation",
  "Exprs_noise",
  "N_tsamples",
  "N_replicates",
  "Tsteps",
  "N_total",
  "N_valid",
  "Valid_fraction",
  "N_null",
  "N_alt",
  "TypeI_001",
  "TypeI_005",
  "TypeI_010",
  "Inflation_005",
  "TypeI_005_Wilson_low",
  "TypeI_005_Wilson_high",
  "TypeI_005_CI_contains_005",
  "Power_001",
  "Power_005",
  "Power_010",
  "Null_q005",
  "Alt_q005",
  "sigma_c_bias_alt",
  "sigma_c_MAE_alt",
  "sigma_c_RMSE_alt",
  "tau_MAE",
  "alpha_MAE",
  "Finite_condition_fraction",
  "Median_condition",
  "Condition_q95",
  "Median_atom_zero_null",
  "Mean_boot_failure",
  "Mean_boot_rank_deficient",
  "Median_test_seconds"
)

missing_columns <- setdiff(
  required_columns,
  names(dt)
)

if (length(missing_columns) > 0L) {
  stop(
    paste0(
      "Missing required columns:\n  ",
      paste(
        missing_columns,
        collapse = "\n  "
      )
    )
  )
}


# =============================================================================
# 5. Exact factorial integrity checks
# =============================================================================

expected_n_rows <- (
  length(EXPECTED_PLATFORMS) *
    length(EXPECTED_PERTURBATIONS) *
    length(EXPECTED_NOISE) *
    length(EXPECTED_N_TSAMPLES) *
    length(EXPECTED_N_REPLICATES) *
    length(EXPECTED_TSTEPS)
)

if (nrow(dt) != expected_n_rows) {
  stop(
    paste(
      "Unexpected number of benchmark rows.",
      "Expected:",
      expected_n_rows,
      "Observed:",
      nrow(dt)
    )
  )
}

check_set <- function(
  observed,
  expected,
  label
) {

  observed <- sort(
    unique(
      observed
    )
  )

  expected <- sort(
    unique(
      expected
    )
  )

  if (!identical(observed, expected)) {
    stop(
      paste0(
        "Unexpected values for ",
        label,
        ". Observed: ",
        paste(
          observed,
          collapse = ", "
        )
      )
    )
  }
}

check_set(
  dt$Platform,
  EXPECTED_PLATFORMS,
  "Platform"
)

check_set(
  dt$Perturbation,
  EXPECTED_PERTURBATIONS,
  "Perturbation"
)

check_set(
  dt$Exprs_noise,
  EXPECTED_NOISE,
  "Exprs_noise"
)

check_set(
  as.integer(dt$N_tsamples),
  EXPECTED_N_TSAMPLES,
  "N_tsamples"
)

check_set(
  as.integer(dt$N_replicates),
  EXPECTED_N_REPLICATES,
  "N_replicates"
)

check_set(
  as.integer(dt$Tsteps),
  EXPECTED_TSTEPS,
  "Tsteps"
)

duplicate_key <- dt[
  ,
  .N,
  by = .(
    Platform,
    Perturbation,
    Exprs_noise,
    N_tsamples,
    N_replicates,
    Tsteps
  )
][
  N != 1L
]

if (nrow(duplicate_key) > 0L) {
  stop(
    "Duplicate or missing factorial cells detected."
  )
}

if (any(!is.finite(dt$TypeI_005))) {
  stop(
    "Non-finite Type-I values found."
  )
}

if (any(!is.finite(dt$Power_005))) {
  stop(
    "Non-finite power values found."
  )
}

cat(
  "\nFactorial integrity checks PASSED.\n"
)


# =============================================================================
# 6. Ordered factors
# =============================================================================

dt[
  ,
  Platform :=
    factor(
      Platform,
      levels = EXPECTED_PLATFORMS
    )
]

dt[
  ,
  Perturbation :=
    factor(
      Perturbation,
      levels = EXPECTED_PERTURBATIONS
    )
]

dt[
  ,
  Exprs_noise :=
    factor(
      Exprs_noise,
      levels = EXPECTED_NOISE
    )
]


# =============================================================================
# 7. Calibration classes
# =============================================================================

dt[
  ,
  calibration_class :=
    fcase(

      TypeI_005_Wilson_low >
        NOMINAL_ALPHA,

      "Anti-conservative",

      TypeI_005_Wilson_high <
        NOMINAL_ALPHA,

      "Conservative",

      default =
        "CI includes 0.05"
    )
]

dt[
  ,
  practical_near_nominal :=
    TypeI_005 >=
      PRACTICAL_LOW &
    TypeI_005 <=
      PRACTICAL_HIGH
]

dt[
  ,
  calibration_distance :=
    abs(
      TypeI_005 -
        NOMINAL_ALPHA
    )
]


# =============================================================================
# 8. Helper functions
# =============================================================================

safe_quantile <- function(
  x,
  p
) {

  x <- x[
    is.finite(
      x
    )
  ]

  if (length(x) == 0L) {
    return(
      NA_real_
    )
  }

  unname(
    quantile(
      x,
      p,
      na.rm = TRUE
    )
  )
}


safe_median <- function(
  x
) {

  x <- x[
    is.finite(
      x
    )
  ]

  if (length(x) == 0L) {
    return(
      NA_real_
    )
  }

  median(
    x,
    na.rm = TRUE
  )
}


safe_mean <- function(
  x
) {

  x <- x[
    is.finite(
      x
    )
  ]

  if (length(x) == 0L) {
    return(
      NA_real_
    )
  }

  mean(
    x,
    na.rm = TRUE
  )
}


# =============================================================================
# 9. Global benchmark summary
# =============================================================================

global_summary <- dt[
  ,
  .(
    N_configurations =
      .N,

    TypeI_median =
      median(
        TypeI_005
      ),

    TypeI_mean =
      mean(
        TypeI_005
      ),

    TypeI_q05 =
      safe_quantile(
        TypeI_005,
        0.05
      ),

    TypeI_q25 =
      safe_quantile(
        TypeI_005,
        0.25
      ),

    TypeI_q75 =
      safe_quantile(
        TypeI_005,
        0.75
      ),

    TypeI_q95 =
      safe_quantile(
        TypeI_005,
        0.95
      ),

    Inflation_median =
      median(
        Inflation_005
      ),

    Fraction_CI_contains_005 =
      mean(
        TypeI_005_CI_contains_005
      ),

    Fraction_practical_near_nominal =
      mean(
        practical_near_nominal
      ),

    Fraction_anti_conservative =
      mean(
        calibration_class ==
          "Anti-conservative"
      ),

    Fraction_conservative =
      mean(
        calibration_class ==
          "Conservative"
      ),

    Power_median =
      median(
        Power_005
      ),

    Power_q25 =
      safe_quantile(
        Power_005,
        0.25
      ),

    Power_q75 =
      safe_quantile(
        Power_005,
        0.75
      ),

    sigma_c_bias_median =
      safe_median(
        sigma_c_bias_alt
      ),

    sigma_c_MAE_median =
      safe_median(
        sigma_c_MAE_alt
      ),

    sigma_c_RMSE_median =
      safe_median(
        sigma_c_RMSE_alt
      ),

    Median_atom_zero_null =
      safe_median(
        Median_atom_zero_null
      ),

    Finite_condition_median =
      median(
        Finite_condition_fraction
      ),

    Mean_boot_failure =
      mean(
        Mean_boot_failure
      ),

    Median_test_seconds =
      median(
        Median_test_seconds
      )
  ),
  by = .(
    Perturbation,
    Platform
  )
]

fwrite(
  global_summary,
  file.path(
    OUTPUT_DIR,
    "benchmark_global_summary.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 10. SHUTOFF calibration by platform and noise
# =============================================================================

shutoff_noise_summary <- dt[
  Perturbation ==
    "SHUTOFF",
  .(
    N_configurations =
      .N,

    TypeI_median =
      median(
        TypeI_005
      ),

    TypeI_mean =
      mean(
        TypeI_005
      ),

    TypeI_q05 =
      safe_quantile(
        TypeI_005,
        0.05
      ),

    TypeI_q95 =
      safe_quantile(
        TypeI_005,
        0.95
      ),

    Fraction_CI_contains_005 =
      mean(
        TypeI_005_CI_contains_005
      ),

    Fraction_practical_near_nominal =
      mean(
        practical_near_nominal
      ),

    Power_median =
      median(
        Power_005
      ),

    Power_q25 =
      safe_quantile(
        Power_005,
        0.25
      ),

    Power_q75 =
      safe_quantile(
        Power_005,
        0.75
      ),

    sigma_c_RMSE_median =
      safe_median(
        sigma_c_RMSE_alt
      ),

    Median_atom_zero_null =
      safe_median(
        Median_atom_zero_null
      ),

    Finite_condition_median =
      median(
        Finite_condition_fraction
      ),

    Mean_boot_failure =
      mean(
        Mean_boot_failure
      )
  ),
  by = .(
    Platform,
    Exprs_noise
  )
]

fwrite(
  shutoff_noise_summary,
  file.path(
    OUTPUT_DIR,
    "benchmark_shutoff_noise_summary.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 11. Design-level summaries
#
# Aggregate across platform x noise for each sampling design.
# This table is descriptive. It does not declare a "best" experimental design.
# =============================================================================

design_summary <- dt[
  ,
  .(
    N_scenarios =
      .N,

    Median_TypeI =
      median(
        TypeI_005
      ),

    Mean_TypeI =
      mean(
        TypeI_005
      ),

    Worst_TypeI =
      max(
        TypeI_005
      ),

    Minimum_TypeI =
      min(
        TypeI_005
      ),

    Median_calibration_distance =
      median(
        calibration_distance
      ),

    Fraction_CI_contains_005 =
      mean(
        TypeI_005_CI_contains_005
      ),

    Fraction_practical_near_nominal =
      mean(
        practical_near_nominal
      ),

    Median_power =
      median(
        Power_005
      ),

    Minimum_power =
      min(
        Power_005
      ),

    Maximum_power =
      max(
        Power_005
      ),

    Median_sigma_c_RMSE =
      safe_median(
        sigma_c_RMSE_alt
      ),

    Median_sigma_c_MAE =
      safe_median(
        sigma_c_MAE_alt
      ),

    Median_atom_zero_null =
      safe_median(
        Median_atom_zero_null
      ),

    Median_finite_condition_fraction =
      median(
        Finite_condition_fraction
      ),

    Median_test_seconds =
      median(
        Median_test_seconds
      )
  ),
  by = .(
    Perturbation,
    N_tsamples,
    N_replicates,
    Tsteps
  )
]

setorder(
  design_summary,
  Perturbation,
  Median_calibration_distance,
  -Fraction_CI_contains_005,
  -Median_power
)

fwrite(
  design_summary,
  file.path(
    OUTPUT_DIR,
    "benchmark_design_summary.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 12. SHUTOFF operational-domain table
#
# Primary operational-domain flag:
#   Wilson CI includes nominal alpha.
#
# Secondary descriptive flag:
#   empirical Type-I in 0.025-0.075 practical band.
#
# We preserve ALL scenarios so the table can be filtered transparently.
# =============================================================================

shutoff_operational <- copy(
  dt[
    Perturbation ==
      "SHUTOFF"
  ]
)

shutoff_operational[
  ,
  operational_CI :=
    TypeI_005_CI_contains_005
]

shutoff_operational[
  ,
  operational_practical :=
    practical_near_nominal
]

shutoff_operational[
  ,
  operational_both :=
    operational_CI &
    operational_practical
]

setorder(
  shutoff_operational,
  -operational_both,
  calibration_distance,
  -Power_005
)

fwrite(
  shutoff_operational,
  file.path(
    OUTPUT_DIR,
    "benchmark_shutoff_operational_domain.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 13. Calibration class counts
# =============================================================================

calibration_classes <- dt[
  ,
  .N,
  by = .(
    Perturbation,
    calibration_class
  )
]

fwrite(
  calibration_classes,
  file.path(
    OUTPUT_DIR,
    "benchmark_calibration_classes.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 14. Flagged scenarios
#
# Flag substantial departures from nominal calibration and numerical failures.
# =============================================================================

flagged_scenarios <- dt[
  calibration_class !=
    "CI includes 0.05" |
    Mean_boot_failure >
      0 |
    Valid_fraction <
      0.99 |
    Finite_condition_fraction <
      0.95
]

setorder(
  flagged_scenarios,
  Perturbation,
  -TypeI_005
)

fwrite(
  flagged_scenarios,
  file.path(
    OUTPUT_DIR,
    "benchmark_flagged_scenarios.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 15. Runtime summary
# =============================================================================

runtime_summary <- dt[
  ,
  .(
    N_configurations =
      .N,

    Median_test_seconds =
      median(
        Median_test_seconds
      ),

    Q25_test_seconds =
      safe_quantile(
        Median_test_seconds,
        0.25
      ),

    Q75_test_seconds =
      safe_quantile(
        Median_test_seconds,
        0.75
      ),

    Q95_test_seconds =
      safe_quantile(
        Median_test_seconds,
        0.95
      )
  ),
  by = .(
    Platform,
    Perturbation
  )
]

fwrite(
  runtime_summary,
  file.path(
    OUTPUT_DIR,
    "benchmark_runtime_summary.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 16. Optional raw-result diagnostics
# =============================================================================

raw_available <- file.exists(
  RAW_FILE
)

raw_summary <- NULL

if (raw_available) {

  raw <- fread(
    RAW_FILE
  )

  raw_summary <- data.table(
    N_rows =
      nrow(
        raw
      ),

    N_genes =
      if (
        "Gene" %in%
          names(raw)
      ) {
        uniqueN(
          raw$Gene
        )
      } else {
        NA_integer_
      },

    N_onset_values =
      if (
        "Onset_time" %in%
          names(raw)
      ) {
        uniqueN(
          raw$Onset_time
        )
      } else {
        NA_integer_
      },

    Fraction_finite_p =
      if (
        "p.value" %in%
          names(raw)
      ) {
        mean(
          is.finite(
            raw$p.value
          )
        )
      } else {
        NA_real_
      }
  )

  fwrite(
    raw_summary,
    file.path(
      OUTPUT_DIR,
      "benchmark_raw_integrity_summary.tsv"
    ),
    sep = "\t"
  )
}


# =============================================================================
# 17. Plot theme
# =============================================================================

theme_benchmark <- theme_classic(
  base_size = 10.5
) +
  theme(
    panel.grid.major =
      element_line(
        colour = "grey93",
        linewidth = 0.25
      ),

    panel.grid.minor =
      element_blank(),

    strip.background =
      element_blank(),

    strip.text =
      element_text(
        face = "bold",
        size = 10
      ),

    legend.position =
      "bottom",

    legend.title =
      element_text(
        size = 9
      ),

    legend.text =
      element_text(
        size = 9
      )
  )


# =============================================================================
# 18. Figure 1: calibration overview
# =============================================================================

p1a <- ggplot(
  dt,
  aes(
    x = Platform,
    y = TypeI_005,
    shape = Perturbation
  )
) +

  geom_hline(
    yintercept =
      NOMINAL_ALPHA,
    linetype =
      "dashed",
    linewidth =
      0.55,
    colour =
      "grey35"
  ) +

  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = PRACTICAL_LOW,
    ymax = PRACTICAL_HIGH,
    alpha = 0.05
  ) +

  geom_boxplot(
    aes(
      group =
        interaction(
          Platform,
          Perturbation
        )
    ),
    position =
      position_dodge(
        width = 0.7
      ),
    width =
      0.55,
    outlier.shape =
      NA,
    linewidth =
      0.45
  ) +

  geom_point(
    position =
      position_jitterdodge(
        jitter.width = 0.12,
        dodge.width = 0.7
      ),
    alpha =
      0.17,
    size =
      0.8
  ) +

  labs(
    x =
      NULL,

    y =
      "Empirical Type-I error",

    shape =
      NULL
  ) +

  coord_cartesian(
    ylim =
      c(
        0,
        max(
          0.60,
          max(
            dt$TypeI_005,
            na.rm = TRUE
          )
        )
      )
  ) +

  theme_benchmark


p1b_data <- copy(
  shutoff_noise_summary
)

p1b_data[
  ,
  Exprs_noise :=
    factor(
      Exprs_noise,
      levels =
        EXPECTED_NOISE
    )
]

p1b <- ggplot(
  p1b_data,
  aes(
    x =
      Exprs_noise,
    y =
      TypeI_median,
    group =
      Platform,
    linetype =
      Platform,
    shape =
      Platform
  )
) +

  annotate(
    "rect",
    xmin =
      -Inf,
    xmax =
      Inf,
    ymin =
      PRACTICAL_LOW,
    ymax =
      PRACTICAL_HIGH,
    alpha =
      0.05
  ) +

  geom_hline(
    yintercept =
      NOMINAL_ALPHA,
    linetype =
      "dashed",
    linewidth =
      0.55,
    colour =
      "grey35"
  ) +

  geom_line(
    linewidth =
      0.8
  ) +

  geom_point(
    size =
      2
  ) +

  labs(
    x =
      "Expression-noise regime",

    y =
      "Median Type-I error\n(SHUTOFF)",

    linetype =
      NULL,

    shape =
      NULL
  ) +

  theme_benchmark


p1c <- ggplot(
  dt[
    Perturbation ==
      "SHUTOFF"
  ],
  aes(
    x =
      TypeI_005,
    y =
      Power_005,
    shape =
      Platform
  )
) +

  geom_vline(
    xintercept =
      NOMINAL_ALPHA,
    linetype =
      "dashed",
    linewidth =
      0.55,
    colour =
      "grey35"
  ) +

  geom_point(
    alpha =
      0.35,
    size =
      1.3
  ) +

  labs(
    x =
      "Empirical Type-I error",

    y =
      "Empirical power",

    shape =
      NULL
  ) +

  theme_benchmark


fig1 <- (
  p1a |
    p1b |
    p1c
) +
  plot_layout(
    widths =
      c(
        1.05,
        1,
        1
      )
  )


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Fig1_calibration_overview.pdf"
    ),
  plot =
    fig1,
  width =
    14,
  height =
    4.7,
  units =
    "in",
  device =
    cairo_pdf
)


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Fig1_calibration_overview.png"
    ),
  plot =
    fig1,
  width =
    14,
  height =
    4.7,
  units =
    "in",
  dpi =
    400,
  bg =
    "white"
)


# =============================================================================
# 19. Figure 2: SHUTOFF design dependence
#
# Aggregate over platform and noise to show large-scale design dependence
# without claiming that platform/noise effects are absent.
# =============================================================================

design_plot_data <- dt[
  Perturbation ==
    "SHUTOFF",
  .(
    TypeI_median =
      median(
        TypeI_005
      ),

    TypeI_q25 =
      safe_quantile(
        TypeI_005,
        0.25
      ),

    TypeI_q75 =
      safe_quantile(
        TypeI_005,
        0.75
      ),

    Power_median =
      median(
        Power_005
      )
  ),
  by = .(
    N_tsamples,
    N_replicates,
    Tsteps
  )
]

design_plot_data[
  ,
  N_replicates :=
    factor(
      N_replicates,
      levels =
        EXPECTED_N_REPLICATES
    )
]


p2a <- ggplot(
  design_plot_data,
  aes(
    x =
      N_tsamples,
    y =
      TypeI_median,
    group =
      N_replicates,
    linetype =
      N_replicates,
    shape =
      N_replicates
  )
) +

  annotate(
    "rect",
    xmin =
      -Inf,
    xmax =
      Inf,
    ymin =
      PRACTICAL_LOW,
    ymax =
      PRACTICAL_HIGH,
    alpha =
      0.05
  ) +

  geom_hline(
    yintercept =
      NOMINAL_ALPHA,
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
      2
  ) +

  facet_wrap(
    ~Tsteps,
    nrow =
      1,
    labeller =
      label_both
  ) +

  scale_x_continuous(
    breaks =
      EXPECTED_N_TSAMPLES
  ) +

  labs(
    x =
      "Number of sampled time points",

    y =
      "Median Type-I error",

    linetype =
      "Replicates",

    shape =
      "Replicates"
  ) +

  theme_benchmark


p2b <- ggplot(
  design_plot_data,
  aes(
    x =
      N_tsamples,
    y =
      Power_median,
    group =
      N_replicates,
    linetype =
      N_replicates,
    shape =
      N_replicates
  )
) +

  geom_line(
    linewidth =
      0.8
  ) +

  geom_point(
    size =
      2
  ) +

  facet_wrap(
    ~Tsteps,
    nrow =
      1,
    labeller =
      label_both
  ) +

  scale_x_continuous(
    breaks =
      EXPECTED_N_TSAMPLES
  ) +

  labs(
    x =
      "Number of sampled time points",

    y =
      "Median empirical power",

    linetype =
      "Replicates",

    shape =
      "Replicates"
  ) +

  theme_benchmark


fig2 <- p2a / p2b


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Fig2_shutoff_design_dependence.pdf"
    ),
  plot =
    fig2,
  width =
    12,
  height =
    7.3,
  units =
    "in",
  device =
    cairo_pdf
)


# =============================================================================
# 20. Figure 3: calibration-power tradeoff
# =============================================================================

tradeoff_data <- copy(
  dt[
    Perturbation ==
      "SHUTOFF"
  ]
)

tradeoff_data[
  ,
  calibration_status :=
    factor(
      calibration_class,
      levels = c(
        "CI includes 0.05",
        "Conservative",
        "Anti-conservative"
      )
    )
]

p3 <- ggplot(
  tradeoff_data,
  aes(
    x =
      TypeI_005,
    y =
      Power_005,
    shape =
      calibration_status
  )
) +

  annotate(
    "rect",
    xmin =
      PRACTICAL_LOW,
    xmax =
      PRACTICAL_HIGH,
    ymin =
      -Inf,
    ymax =
      Inf,
    alpha =
      0.05
  ) +

  geom_vline(
    xintercept =
      NOMINAL_ALPHA,
    linetype =
      "dashed",
    linewidth =
      0.5,
    colour =
      "grey35"
  ) +

  geom_point(
    alpha =
      0.45,
    size =
      1.6
  ) +

  facet_grid(
    Platform ~ Exprs_noise
  ) +

  labs(
    x =
      "Empirical Type-I error",

    y =
      "Empirical power",

    shape =
      "Calibration class"
  ) +

  theme_benchmark


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Fig3_power_calibration_tradeoff.pdf"
    ),
  plot =
    p3,
  width =
    12,
  height =
    8.5,
  units =
    "in",
  device =
    cairo_pdf
)


# =============================================================================
# 21. Figure 4: parameter recovery
# =============================================================================

recovery_long <- melt(
  dt[
    Perturbation ==
      "SHUTOFF"
  ],
  id.vars = c(
    "Platform",
    "Exprs_noise",
    "N_tsamples",
    "N_replicates",
    "Tsteps"
  ),
  measure.vars = c(
    "sigma_c_RMSE_alt",
    "sigma_c_MAE_alt",
    "tau_MAE",
    "alpha_MAE"
  ),
  variable.name =
    "Metric",
  value.name =
    "Value"
)

metric_labels <- c(
  sigma_c_RMSE_alt =
    "sigma[c] RMSE",
  sigma_c_MAE_alt =
    "sigma[c] MAE",
  tau_MAE =
    "tau MAE",
  alpha_MAE =
    "alpha MAE"
)

recovery_long[
  ,
  Metric :=
    factor(
      Metric,
      levels =
        names(
          metric_labels
        ),
      labels =
        unname(
          metric_labels
        )
    )
]


p4 <- ggplot(
  recovery_long[
    is.finite(
      Value
    )
  ],
  aes(
    x =
      Exprs_noise,
    y =
      Value,
    shape =
      Platform
  )
) +

  geom_boxplot(
    aes(
      group =
        interaction(
          Exprs_noise,
          Platform
        )
    ),
    outlier.shape =
      NA,
    linewidth =
      0.4
  ) +

  geom_point(
    position =
      position_jitter(
        width =
          0.15
      ),
    alpha =
      0.16,
    size =
      0.7
  ) +

  facet_wrap(
    ~Metric,
    scales =
      "free_y",
    ncol =
      2
  ) +

  labs(
    x =
      "Expression-noise regime",

    y =
      NULL,

    shape =
      NULL
  ) +

  theme_benchmark


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Fig4_parameter_recovery.pdf"
    ),
  plot =
    p4,
  width =
    9,
  height =
    7,
  units =
    "in",
  device =
    cairo_pdf
)


# =============================================================================
# 22. Supplementary heatmap helper
# =============================================================================

make_heatmap_pdf <- function(
  value_column,
  output_filename,
  fill_label
) {

  plot_data <- copy(
    dt[
      Perturbation ==
        "SHUTOFF"
    ]
  )

  plot_data[
    ,
    value_to_plot :=
      get(
        value_column
      )
  ]

  p <- ggplot(
    plot_data,
    aes(
      x =
        factor(
          N_tsamples,
          levels =
            EXPECTED_N_TSAMPLES
        ),

      y =
        factor(
          N_replicates,
          levels =
            EXPECTED_N_REPLICATES
        ),

      fill =
        value_to_plot
    )
  ) +

    geom_tile(
      colour =
        "white",
      linewidth =
        0.4
    ) +

    geom_text(
      aes(
        label =
          sprintf(
            "%.3f",
            value_to_plot
          )
      ),
      size =
        2.5
    ) +

    facet_grid(
      Platform + Exprs_noise ~ Tsteps,
      labeller =
        label_both
    ) +

    labs(
      x =
        "Number of sampled time points",

      y =
        "Replicates",

      fill =
        fill_label
    ) +

    theme_minimal(
      base_size =
        8
    ) +

    theme(
      panel.grid =
        element_blank(),

      strip.text =
        element_text(
          size =
            7
        ),

      legend.position =
        "bottom"
    )

  ggsave(
    filename =
      file.path(
        OUTPUT_DIR,
        output_filename
      ),
    plot =
      p,
    width =
      12,
    height =
      18,
    units =
      "in",
    device =
      cairo_pdf,
    limitsize =
      FALSE
  )
}


make_heatmap_pdf(
  value_column =
    "TypeI_005",
  output_filename =
    "FigS1_typeI_heatmaps_SHUTOFF.pdf",
  fill_label =
    "Type-I"
)

make_heatmap_pdf(
  value_column =
    "Power_005",
  output_filename =
    "FigS2_power_heatmaps_SHUTOFF.pdf",
  fill_label =
    "Power"
)

make_heatmap_pdf(
  value_column =
    "sigma_c_RMSE_alt",
  output_filename =
    "FigS3_sigma_c_RMSE_heatmaps_SHUTOFF.pdf",
  fill_label =
    "sigma[c] RMSE"
)


# =============================================================================
# 23. Supplementary NONE reference figure
# =============================================================================

none_data <- dt[
  Perturbation ==
    "NONE"
]

p_none <- ggplot(
  none_data,
  aes(
    x =
      TypeI_005,
    y =
      Power_005,
    shape =
      Platform
  )
) +

  geom_vline(
    xintercept =
      NOMINAL_ALPHA,
    linetype =
      "dashed",
    colour =
      "grey35",
    linewidth =
      0.5
  ) +

  geom_point(
    alpha =
      0.4,
    size =
      1.4
  ) +

  facet_wrap(
    ~Exprs_noise,
    nrow =
      1
  ) +

  labs(
    x =
      "Empirical Type-I error",

    y =
      "Empirical positive rate under alternatives",

    shape =
      NULL
  ) +

  theme_benchmark


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "FigS4_NONE_reference.pdf"
    ),
  plot =
    p_none,
  width =
    11,
  height =
    3.8,
  units =
    "in",
  device =
    cairo_pdf
)


# =============================================================================
# 24. Save complete analysis object
# =============================================================================

benchmark_corrected_analysis <- list(

  settings =
    list(
      nominal_alpha =
        NOMINAL_ALPHA,

      practical_band =
        c(
          PRACTICAL_LOW,
          PRACTICAL_HIGH
        ),

      expected_rows =
        expected_n_rows
    ),

  benchmark_summary =
    dt,

  global_summary =
    global_summary,

  shutoff_noise_summary =
    shutoff_noise_summary,

  design_summary =
    design_summary,

  shutoff_operational =
    shutoff_operational,

  calibration_classes =
    calibration_classes,

  flagged_scenarios =
    flagged_scenarios,

  runtime_summary =
    runtime_summary,

  raw_integrity_summary =
    raw_summary
)


save(
  benchmark_corrected_analysis,
  file =
    file.path(
      OUTPUT_DIR,
      "benchmark_corrected_analysis.rdata"
    )
)


# =============================================================================
# 25. Session information
# =============================================================================

sink(
  file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()


# =============================================================================
# 26. Console report
# =============================================================================

cat(
  "\n============================================================\n",
  "FINAL CORRECTED-ONSET BENCHMARK ANALYSIS\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nGLOBAL NONE vs SHUTOFF SUMMARY\n"
)

print(
  global_summary
)


cat(
  "\nSHUTOFF BY PLATFORM AND NOISE\n"
)

print(
  shutoff_noise_summary[
    order(
      Platform,
      Exprs_noise
    )
  ]
)


cat(
  "\nCALIBRATION CLASS COUNTS\n"
)

print(
  calibration_classes
)


cat(
  "\nSHUTOFF OPERATIONAL-DOMAIN COUNTS\n"
)

print(
  shutoff_operational[
    ,
    .(
      N =
        .N
    ),
    by = .(
      operational_CI,
      operational_practical,
      operational_both
    )
  ][
    order(
      -operational_both,
      -operational_CI,
      -operational_practical
    )
  ]
)


cat(
  "\nTOP 20 SHUTOFF DESIGNS BY CALIBRATION DISTANCE\n"
)

print(
  design_summary[
    Perturbation ==
      "SHUTOFF"
  ][
    1:min(
      20L,
      .N
    )
  ]
)


cat(
  "\nGenerated in:\n  ",
  OUTPUT_DIR,
  "\n",
  sep = ""
)


cat(
  "\nInterpretation framework:\n",
  "  1. Evaluate Type-I calibration before power.\n",
  "  2. Treat SHUTOFF as the principal inferential regime.\n",
  "  3. Treat NONE as a reference / low-identifiability regime.\n",
  "  4. Wilson-CI inclusion of 0.05 is the primary calibration diagnostic.\n",
  "  5. The 0.025-0.075 interval is descriptive only.\n",
  "  6. Parameter recovery must be interpreted jointly with calibration and\n",
  "     numerical identifiability.\n",
  sep = ""
)
