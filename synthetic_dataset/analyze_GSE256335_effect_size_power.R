# =============================================================================
# Title:
# GSE256335 Matched Benchmark: Power and Effect-Size Analysis
#
# Purpose:
#   Analyze the corrected exact-matched mESC benchmark to quantify how
#   detection performance depends on the true post-export conversion rate
#   sigma_c and its corresponding half-time:
#
#       t1/2 = log(2) / sigma_c
#
#   The analysis focuses on the empirically relevant Medium and High
#   RNA-seq noise regimes identified by the mESC noise-matching analysis.
#
# Main outputs:
#   - power by sigma_c bin;
#   - power by sigma_c half-time bin;
#   - estimation bias / MAE / RMSE by effect-size bin;
#   - fraction of fits estimated on the sigma_c = 0 boundary;
#   - null Type-I summary for Medium and High noise;
#   - publication-ready power curves.
#
# IMPORTANT:
#   This script DOES NOT rerun the benchmark.
#   It uses the completed raw matched-benchmark output only.
#
# Expected input:
#   synthetic_dataset/GSE256335_matched_corrected_onset_raw.tsv
#
# Outputs:
#   synthetic_dataset/GSE256335_effect_size_analysis/
#       power_by_sigma_c_bin.tsv
#       power_by_half_time_bin.tsv
#       continuous_power_curve.tsv
#       estimation_by_sigma_c_bin.tsv
#       null_calibration_medium_high.tsv
#       Fig_power_vs_sigma_c.pdf
#       Fig_power_vs_half_time.pdf
#       Fig_sigma_c_recovery.pdf
#       effect_size_analysis.rdata
#       sessionInfo.txt
#
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics")

library(data.table)
library(ggplot2)


# =============================================================================
# 1. Settings
# =============================================================================

INPUT_FILE <-
  "synthetic_dataset/GSE256335_matched_corrected_onset_raw.tsv"

OUTPUT_DIR <-
  "synthetic_dataset/GSE256335_effect_size_analysis"

NOISE_KEEP <- c(
  "Medium",
  "High"
)

ALPHA <- 0.05

BOUNDARY_EPS <- 1e-12

# Number of quantile bins for true sigma_c.
N_SIGMA_BINS <- 6L

# Fixed half-time bins, in minutes.
HALF_TIME_BREAKS <- c(
  -Inf,
  5,
  10,
  20,
  40,
  80,
  Inf
)

# Sliding-window size for the continuous power curve.
POWER_WINDOW_N <- 150L

# Step between adjacent sliding windows.
POWER_WINDOW_STEP <- 25L


# =============================================================================
# 2. Output directory
# =============================================================================

if (
  !dir.exists(
    OUTPUT_DIR
  )
) {
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE
  )
}


# =============================================================================
# 3. Load data
# =============================================================================

dt <- fread(
  INPUT_FILE
)

cat(
  "\n============================================================\n",
  "MATCHED EFFECT-SIZE ANALYSIS\n",
  "============================================================\n",
  "Rows: ",
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
  "Gene",
  "truth",
  "Exprs_noise",
  "p.value",
  "Sigma",
  "sigma_true",
  "status",
  "condition_number",
  "bootstrap_failure_rate"
)

missing_columns <- setdiff(
  required_columns,
  names(dt)
)

if (
  length(
    missing_columns
  ) > 0L
) {
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
# 5. Harmonize / validate
# =============================================================================

dt <- dt[
  Exprs_noise %in%
    NOISE_KEEP
]

dt[
  ,
  valid :=
    is.finite(
      p.value
    ) &
    p.value >= 0 &
    p.value <= 1 &
    status !=
      "test_error"
]

dt[
  ,
  detected :=
    valid &
    p.value <=
      ALPHA
]

dt[
  ,
  sigma_hat :=
    Sigma
]

dt[
  ,
  sigma_hat_boundary :=
    is.finite(
      sigma_hat
    ) &
    sigma_hat <=
      BOUNDARY_EPS
]


# =============================================================================
# 6. Alternative set
# =============================================================================

alt <- dt[
  truth ==
    "alternative" &
    valid &
    is.finite(
      sigma_true
    ) &
    sigma_true >
      0
]

if (
  nrow(
    alt
  ) == 0L
) {
  stop(
    "No valid alternative rows found."
  )
}

alt[
  ,
  true_half_time :=
    log(2) /
      sigma_true
]

alt[
  ,
  estimated_half_time :=
    fifelse(
      is.finite(
        sigma_hat
      ) &
        sigma_hat >
          0,
      log(2) /
        sigma_hat,
      Inf
    )
]

alt[
  ,
  sigma_error :=
    sigma_hat -
      sigma_true
]

alt[
  ,
  sigma_abs_error :=
    abs(
      sigma_error
    )
]

alt[
  ,
  sigma_sq_error :=
    sigma_error^2
]


# =============================================================================
# 7. Quantile bins of true sigma_c
#
# Common bin boundaries are constructed from BOTH Medium and High combined,
# so the two noise regimes are compared on identical effect-size bins.
# =============================================================================

sigma_breaks <- unique(
  as.numeric(
    quantile(
      alt$sigma_true,
      probs =
        seq(
          0,
          1,
          length.out =
            N_SIGMA_BINS +
            1L
        ),
      na.rm =
        TRUE,
      names =
        FALSE,
      type =
        8
    )
  )
)

if (
  length(
    sigma_breaks
  ) <
    3L
) {
  stop(
    "Insufficient unique sigma_c values to construct effect-size bins."
  )
}

sigma_breaks[1] <- -Inf
sigma_breaks[
  length(
    sigma_breaks
  )
] <- Inf

alt[
  ,
  sigma_bin :=
    cut(
      sigma_true,
      breaks =
        sigma_breaks,
      include.lowest =
        TRUE,
      ordered_result =
        TRUE
    )
]


# =============================================================================
# 8. Power by sigma_c bin
# =============================================================================

power_by_sigma <- alt[
  ,
  .(
    N =
      .N,

    sigma_true_min =
      min(
        sigma_true
      ),

    sigma_true_median =
      median(
        sigma_true
      ),

    sigma_true_max =
      max(
        sigma_true
      ),

    half_time_median =
      median(
        true_half_time
      ),

    power =
      mean(
        detected
      ),

    boundary_fraction =
      mean(
        sigma_hat_boundary
      ),

    sigma_hat_median =
      median(
        sigma_hat,
        na.rm =
          TRUE
      ),

    sigma_bias =
      mean(
        sigma_error,
        na.rm =
          TRUE
      ),

    sigma_MAE =
      mean(
        sigma_abs_error,
        na.rm =
          TRUE
      ),

    sigma_RMSE =
      sqrt(
        mean(
          sigma_sq_error,
          na.rm =
            TRUE
        )
      )
  ),
  by =
    .(
      Exprs_noise,
      sigma_bin
    )
]

setorder(
  power_by_sigma,
  Exprs_noise,
  sigma_true_median
)

fwrite(
  power_by_sigma,
  file.path(
    OUTPUT_DIR,
    "power_by_sigma_c_bin.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 9. Power by true half-time bin
# =============================================================================

half_labels <- c(
  "<5",
  "5-10",
  "10-20",
  "20-40",
  "40-80",
  ">80"
)

alt[
  ,
  half_time_bin :=
    cut(
      true_half_time,
      breaks =
        HALF_TIME_BREAKS,
      labels =
        half_labels,
      include.lowest =
        TRUE,
      right =
        FALSE,
      ordered_result =
        TRUE
    )
]

power_by_half_time <- alt[
  !is.na(
    half_time_bin
  ),
  .(
    N =
      .N,

    true_half_time_median =
      median(
        true_half_time
      ),

    sigma_true_median =
      median(
        sigma_true
      ),

    power =
      mean(
        detected
      ),

    boundary_fraction =
      mean(
        sigma_hat_boundary
      ),

    sigma_bias =
      mean(
        sigma_error,
        na.rm =
          TRUE
      ),

    sigma_MAE =
      mean(
        sigma_abs_error,
        na.rm =
          TRUE
      ),

    sigma_RMSE =
      sqrt(
        mean(
          sigma_sq_error,
          na.rm =
            TRUE
        )
      )
  ),
  by =
    .(
      Exprs_noise,
      half_time_bin
    )
]

setorder(
  power_by_half_time,
  Exprs_noise,
  half_time_bin
)

fwrite(
  power_by_half_time,
  file.path(
    OUTPUT_DIR,
    "power_by_half_time_bin.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 10. Continuous sliding-window power curve
#
# This avoids making the visual result depend entirely on arbitrary bin edges.
# Each point is a local empirical detection fraction across neighboring
# alternatives ordered by true sigma_c.
# =============================================================================

make_sliding_power <- function(
  x,
  window_n =
    POWER_WINDOW_N,
  step_n =
    POWER_WINDOW_STEP
) {

  x <- copy(
    x
  )

  setorder(
    x,
    sigma_true
  )

  n <- nrow(
    x
  )

  if (
    n <
      window_n
  ) {
    return(
      data.table()
    )
  }

  starts <- seq(
    1L,
    n -
      window_n +
      1L,
    by =
      step_n
  )

  out <- lapply(
    starts,
    function(ss) {

      ii <- ss:
        (
          ss +
            window_n -
            1L
        )

      z <- x[ii]

      data.table(
        N =
          nrow(
            z
          ),

        sigma_true_median =
          median(
            z$sigma_true
          ),

        sigma_true_q25 =
          quantile(
            z$sigma_true,
            0.25,
            names =
              FALSE
          ),

        sigma_true_q75 =
          quantile(
            z$sigma_true,
            0.75,
            names =
              FALSE
          ),

        half_time_median =
          median(
            z$true_half_time
          ),

        power =
          mean(
            z$detected
          ),

        boundary_fraction =
          mean(
            z$sigma_hat_boundary
          )
      )
    }
  )

  rbindlist(
    out
  )
}


continuous_power <- rbindlist(
  lapply(
    NOISE_KEEP,
    function(noise_i) {

      xx <- make_sliding_power(
        alt[
          Exprs_noise ==
            noise_i
        ]
      )

      if (
        nrow(
          xx
        ) >
          0L
      ) {
        xx[
          ,
          Exprs_noise :=
            noise_i
        ]
      }

      xx
    }
  ),
  use.names =
    TRUE,
  fill =
    TRUE
)

fwrite(
  continuous_power,
  file.path(
    OUTPUT_DIR,
    "continuous_power_curve.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 11. Estimation metrics by sigma_c bin
# =============================================================================

estimation_by_sigma <- alt[
  ,
  .(
    N =
      .N,

    sigma_true_median =
      median(
        sigma_true
      ),

    sigma_hat_median =
      median(
        sigma_hat,
        na.rm =
          TRUE
      ),

    sigma_bias =
      mean(
        sigma_error,
        na.rm =
          TRUE
      ),

    sigma_MAE =
      mean(
        sigma_abs_error,
        na.rm =
          TRUE
      ),

    sigma_RMSE =
      sqrt(
        mean(
          sigma_sq_error,
          na.rm =
            TRUE
        )
      ),

    boundary_fraction =
      mean(
        sigma_hat_boundary
      ),

    detection_fraction =
      mean(
        detected
      )
  ),
  by =
    .(
      Exprs_noise,
      sigma_bin
    )
]

fwrite(
  estimation_by_sigma,
  file.path(
    OUTPUT_DIR,
    "estimation_by_sigma_c_bin.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 12. Null calibration: Medium and High
# =============================================================================

null <- dt[
  truth ==
    "null" &
    valid
]

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

  p <- successes /
    n

  den <- 1 +
    z^2 /
      n

  ctr <- (
    p +
      z^2 /
        (
          2 *
            n
        )
  ) /
    den

  half <- z *
    sqrt(
      p *
        (
          1 -
            p
        ) /
        n +
        z^2 /
          (
            4 *
              n^2
          )
    ) /
    den

  c(
    low =
      max(
        0,
        ctr -
          half
      ),

    high =
      min(
        1,
        ctr +
          half
      )
  )
}


null_calibration <- null[
  ,
  {

    n0 <- .N

    k <- sum(
      p.value <=
        ALPHA
    )

    ci <- wilson_interval(
      k,
      n0
    )

    .(
      N =
        n0,

      TypeI_005 =
        k /
          n0,

      Wilson_low =
        unname(
          ci[
            "low"
          ]
        ),

      Wilson_high =
        unname(
          ci[
            "high"
          ]
        ),

      CI_contains_005 =
        ci[
          "low"
        ] <=
          ALPHA &&
        ci[
          "high"
        ] >=
          ALPHA,

      fraction_p1 =
        mean(
          p.value ==
            1
        ),

      fraction_sigma0 =
        mean(
          sigma_hat_boundary
        ),

      mean_bootstrap_failure_rate =
        mean(
          bootstrap_failure_rate,
          na.rm =
            TRUE
        )
    )
  },
  by =
    Exprs_noise
]

fwrite(
  null_calibration,
  file.path(
    OUTPUT_DIR,
    "null_calibration_medium_high.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 13. Figure: power vs sigma_c
# =============================================================================

continuous_power[
  ,
  Exprs_noise :=
    factor(
      Exprs_noise,
      levels =
        NOISE_KEEP
    )
]

p_sigma <- ggplot(
  continuous_power,
  aes(
    x =
      sigma_true_median,
    y =
      power,
    linetype =
      Exprs_noise,
    shape =
      Exprs_noise
  )
) +

  geom_hline(
    yintercept =
      ALPHA,
    linewidth =
      0.4,
    linetype =
      "dotted",
    colour =
      "grey55"
  ) +

  geom_line(
    linewidth =
      0.9
  ) +

  geom_point(
    size =
      2.1
  ) +

  scale_y_continuous(
    limits =
      c(
        0,
        1
      ),
    breaks =
      seq(
        0,
        1,
        by =
          0.2
      )
  ) +

  labs(
    x =
      expression(
        "True " *
          sigma[c] *
          " (min"^{-1} * ")"
      ),

    y =
      "Empirical detection probability",

    linetype =
      "RNA-seq noise",

    shape =
      "RNA-seq noise"
  ) +

  theme_classic(
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

    legend.position =
      "bottom"
  )


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Fig_power_vs_sigma_c.pdf"
    ),
  plot =
    p_sigma,
  width =
    6.6,
  height =
    4.8,
  units =
    "in",
  device =
    cairo_pdf
)


# =============================================================================
# 14. Figure: power vs conversion half-time
#
# Smaller half-time = faster post-export conversion.
# Reverse the x-axis so faster conversion is visually on the right.
# =============================================================================

power_half_plot <- copy(
  power_by_half_time
)

power_half_plot[
  ,
  Exprs_noise :=
    factor(
      Exprs_noise,
      levels =
        NOISE_KEEP
    )
]

power_half_plot[
  ,
  half_time_bin :=
    factor(
      half_time_bin,
      levels =
        half_labels
    )
]


p_half <- ggplot(
  power_half_plot,
  aes(
    x =
      half_time_bin,
    y =
      power,
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
      0.9
  ) +

  geom_point(
    size =
      2.2
  ) +

  scale_y_continuous(
    limits =
      c(
        0,
        1
      ),
    breaks =
      seq(
        0,
        1,
        by =
          0.2
      )
  ) +

  labs(
    x =
      expression(
        "True conversion half-time " *
          t[1/2] *
          " (min)"
      ),

    y =
      "Empirical detection probability",

    linetype =
      "RNA-seq noise",

    shape =
      "RNA-seq noise"
  ) +

  theme_classic(
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

    legend.position =
      "bottom"
  )


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Fig_power_vs_half_time.pdf"
    ),
  plot =
    p_half,
  width =
    7,
  height =
    4.8,
  units =
    "in",
  device =
    cairo_pdf
)


# =============================================================================
# 15. Figure: sigma_c recovery
# =============================================================================

recovery_plot <- alt[
  is.finite(
    sigma_hat
  )
]

recovery_plot[
  ,
  Exprs_noise :=
    factor(
      Exprs_noise,
      levels =
        NOISE_KEEP
    )
]


p_recovery <- ggplot(
  recovery_plot,
  aes(
    x =
      sigma_true,
    y =
      sigma_hat
  )
) +

  geom_abline(
    slope =
      1,
    intercept =
      0,
    linewidth =
      0.55,
    linetype =
      "dashed",
    colour =
      "grey45"
  ) +

  geom_point(
    alpha =
      0.18,
    size =
      1
  ) +

  facet_wrap(
    ~Exprs_noise,
    nrow =
      1
  ) +

  labs(
    x =
      expression(
        "True " *
          sigma[c] *
          " (min"^{-1} * ")"
      ),

    y =
      expression(
        "Estimated " *
          hat(
            sigma
          )[c] *
          " (min"^{-1} * ")"
      )
  ) +

  theme_classic(
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
      )
  )


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "Fig_sigma_c_recovery.pdf"
    ),
  plot =
    p_recovery,
  width =
    8.2,
  height =
    4.2,
  units =
    "in",
  device =
    cairo_pdf
)


# =============================================================================
# 16. Useful thresholds
#
# Report the first sliding-window location at which empirical power reaches
# selected targets. If not reached, return NA.
# =============================================================================

POWER_TARGETS <- c(
  0.20,
  0.50,
  0.80
)


power_thresholds <- rbindlist(
  lapply(
    NOISE_KEEP,
    function(noise_i) {

      z <- continuous_power[
        Exprs_noise ==
          noise_i
      ]

      rbindlist(
        lapply(
          POWER_TARGETS,
          function(target_i) {

            zz <- z[
              power >=
                target_i
            ][
              order(
                sigma_true_median
              )
            ]

            if (
              nrow(
                zz
              ) == 0L
            ) {

              return(
                data.table(
                  Exprs_noise =
                    noise_i,
                  target_power =
                    target_i,
                  sigma_c_threshold =
                    NA_real_,
                  half_time_threshold =
                    NA_real_
                )
              )
            }

            data.table(
              Exprs_noise =
                noise_i,
              target_power =
                target_i,
              sigma_c_threshold =
                zz$sigma_true_median[1],
              half_time_threshold =
                log(2) /
                  zz$sigma_true_median[1]
            )
          }
        )
      )
    }
  )
)


fwrite(
  power_thresholds,
  file.path(
    OUTPUT_DIR,
    "power_detection_thresholds.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 17. Save complete analysis
# =============================================================================

effect_size_analysis <- list(

  settings =
    list(
      input_file =
        INPUT_FILE,
      noise_keep =
        NOISE_KEEP,
      alpha =
        ALPHA,
      N_sigma_bins =
        N_SIGMA_BINS,
      half_time_breaks =
        HALF_TIME_BREAKS,
      sliding_window_n =
        POWER_WINDOW_N,
      sliding_window_step =
        POWER_WINDOW_STEP
    ),

  alternative =
    alt,

  null =
    null,

  sigma_breaks =
    sigma_breaks,

  power_by_sigma =
    power_by_sigma,

  power_by_half_time =
    power_by_half_time,

  continuous_power =
    continuous_power,

  estimation_by_sigma =
    estimation_by_sigma,

  null_calibration =
    null_calibration,

  power_thresholds =
    power_thresholds
)


save(
  effect_size_analysis,
  file =
    file.path(
      OUTPUT_DIR,
      "effect_size_analysis.rdata"
    )
)


# =============================================================================
# 18. Session info
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
# 19. Console report
# =============================================================================

cat(
  "\n============================================================\n",
  "NULL CALIBRATION\n",
  "============================================================\n",
  sep = ""
)

print(
  null_calibration
)


cat(
  "\n============================================================\n",
  "POWER BY TRUE sigma_c BIN\n",
  "============================================================\n",
  sep = ""
)

print(
  power_by_sigma
)


cat(
  "\n============================================================\n",
  "POWER BY TRUE HALF-TIME BIN\n",
  "============================================================\n",
  sep = ""
)

print(
  power_by_half_time
)


cat(
  "\n============================================================\n",
  "POWER DETECTION THRESHOLDS\n",
  "============================================================\n",
  sep = ""
)

print(
  power_thresholds
)


cat(
  "\nSaved in:\n  ",
  OUTPUT_DIR,
  "\n",
  sep = ""
)
