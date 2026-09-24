# =============================================================================
# Practical-identifiability summary in the calibrated SHUTOFF domain
#
# Purpose:
#   Summarize true-versus-estimated kinetic parameters using the completed
#   corrected-onset factorial benchmark, restricting the analysis to SHUTOFF
#   configurations whose Wilson 95% CI for empirical Type-I error contains
#   the nominal 0.05 level.
#
# This restriction is deliberate: parameter recovery should not be interpreted
# independently of null calibration.
#
# Inputs:
#   benchmark_main_corrected_onset_raw.tsv
#   benchmark_main_corrected_onset_summary.tsv
#
# Outputs:
#   practical_identifiability/
#       parameter_recovery_operational_domain.tsv
#       parameter_recovery_metrics.tsv
#       parameter_error_correlations.tsv
#       FigS_parameter_recovery_operational_domain.pdf
#       FigS_parameter_recovery_operational_domain.png
#
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/synthetic_dataset")

library(data.table)
library(ggplot2)


# =============================================================================
# 1. Files
# =============================================================================

RAW_FILE <- "benchmark_main_corrected_onset_raw.tsv"
SUMMARY_FILE <- "benchmark_main_corrected_onset_summary.tsv"

OUTPUT_DIR <- "practical_identifiability"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE
  )
}


# =============================================================================
# 2. Load
# =============================================================================

raw <- fread(
  RAW_FILE
)

summary_dt <- fread(
  SUMMARY_FILE
)


# =============================================================================
# 3. Required columns
# =============================================================================

config_cols <- c(
  "Platform",
  "Perturbation",
  "Exprs_noise",
  "N_tsamples",
  "N_replicates",
  "Tsteps"
)

required_raw <- c(
  config_cols,
  "Positive",
  "status",
  "sigma_true",
  "sigma_c_hat",
  "tau_true",
  "tau_hat",
  "tau_s_true",
  "tau_s_hat",
  "sigma_n_true",
  "sigma_n_hat",
  "alpha_true",
  "alpha_hat",
  "alpha_s_true",
  "alpha_s_hat",
  "condition_number",
  "bootstrap_failure_rate"
)

required_summary <- c(
  config_cols,
  "TypeI_005",
  "TypeI_005_CI_contains_005"
)

missing_raw <- setdiff(
  required_raw,
  names(raw)
)

missing_summary <- setdiff(
  required_summary,
  names(summary_dt)
)

if (length(missing_raw) > 0L) {
  stop(
    paste0(
      "Missing columns in raw benchmark:\n  ",
      paste(
        missing_raw,
        collapse = "\n  "
      )
    )
  )
}

if (length(missing_summary) > 0L) {
  stop(
    paste0(
      "Missing columns in summary benchmark:\n  ",
      paste(
        missing_summary,
        collapse = "\n  "
      )
    )
  )
}


# =============================================================================
# 4. Restrict to calibrated SHUTOFF configurations
# =============================================================================

operational <- summary_dt[
  Perturbation == "SHUTOFF" &
    TypeI_005_CI_contains_005 == TRUE,
  c(
    config_cols,
    "TypeI_005"
  ),
  with = FALSE
]

if (nrow(operational) == 0L) {
  stop(
    "No CI-calibrated SHUTOFF configurations found."
  )
}


dt <- merge(
  raw,
  operational,
  by = config_cols,
  all = FALSE,
  suffixes = c(
    "",
    "_config"
  )
)

dt <- dt[
  status == "ok" &
    Positive == 1
]

cat(
  "\n============================================================\n",
  "PRACTICAL IDENTIFIABILITY\n",
  "============================================================\n",
  "Operational SHUTOFF configurations: ",
  nrow(operational),
  "\n",
  "Alternative fitted rows retained:   ",
  nrow(dt),
  "\n",
  sep = ""
)


# =============================================================================
# 5. Convert to long parameter table
# =============================================================================

param_map <- data.table(
  Parameter = c(
    "sigma_c",
    "sigma_n",
    "tau",
    "tau_s",
    "alpha",
    "alpha_s"
  ),
  True_col = c(
    "sigma_true",
    "sigma_n_true",
    "tau_true",
    "tau_s_true",
    "alpha_true",
    "alpha_s_true"
  ),
  Hat_col = c(
    "sigma_c_hat",
    "sigma_n_hat",
    "tau_hat",
    "tau_s_hat",
    "alpha_hat",
    "alpha_s_hat"
  )
)


long_list <- lapply(
  seq_len(
    nrow(
      param_map
    )
  ),
  function(ii) {

    pp <- param_map[ii]

    data.table(
      Parameter =
        pp$Parameter,

      Platform =
        dt$Platform,

      Exprs_noise =
        dt$Exprs_noise,

      N_tsamples =
        dt$N_tsamples,

      N_replicates =
        dt$N_replicates,

      Tsteps =
        dt$Tsteps,

      True =
        dt[[pp$True_col]],

      Estimated =
        dt[[pp$Hat_col]],

      condition_number =
        dt$condition_number,

      bootstrap_failure_rate =
        dt$bootstrap_failure_rate
    )
  }
)


recovery <- rbindlist(
  long_list
)

recovery <- recovery[
  is.finite(True) &
    is.finite(Estimated)
]

recovery[
  ,
  Error :=
    Estimated -
      True
]

recovery[
  ,
  Abs_error :=
    abs(
      Error
    )
]

recovery[
  ,
  Sq_error :=
    Error^2
]


# =============================================================================
# 6. Parameter-recovery metrics
# =============================================================================

metrics <- recovery[
  ,
  .(
    N =
      .N,

    Bias =
      mean(
        Error
      ),

    MAE =
      mean(
        Abs_error
      ),

    RMSE =
      sqrt(
        mean(
          Sq_error
        )
      ),

    Pearson =
      if (
        sd(True) > 0 &
          sd(Estimated) > 0
      ) {
        cor(
          True,
          Estimated,
          method = "pearson"
        )
      } else {
        NA_real_
      },

    Spearman =
      if (
        length(
          unique(
            True
          )
        ) > 1L &
          length(
            unique(
              Estimated
            )
          ) > 1L
      ) {
        cor(
          True,
          Estimated,
          method = "spearman"
        )
      } else {
        NA_real_
      },

    Boundary_fraction =
      mean(
        Estimated <=
          1e-12
      )
  ),
  by =
    .(
      Parameter
    )
]

setorder(
  metrics,
  Parameter
)


# =============================================================================
# 7. Error correlations relevant to confounding with sigma_c
#
# Reviewer concern:
#   sigma_c may be confounded particularly with export and cytoplasmic decay.
#
# We therefore correlate sigma_c estimation error with the estimation errors
# of tau, tau_s, alpha, alpha_s and sigma_n within the same benchmark rows.
# =============================================================================

corr_source <- dt[
  is.finite(
    sigma_c_hat
  ) &
    is.finite(
      sigma_true
    )
]


corr_source[
  ,
  sigma_c_error :=
    sigma_c_hat -
      sigma_true
]

corr_source[
  ,
  tau_error_local :=
    tau_hat -
      tau_true
]

corr_source[
  ,
  tau_s_error_local :=
    tau_s_hat -
      tau_s_true
]

corr_source[
  ,
  sigma_n_error_local :=
    sigma_n_hat -
      sigma_n_true
]

corr_source[
  ,
  alpha_error_local :=
    alpha_hat -
      alpha_true
]

corr_source[
  ,
  alpha_s_error_local :=
    alpha_s_hat -
      alpha_s_true
]


other_errors <- c(
  tau =
    "tau_error_local",
  tau_s =
    "tau_s_error_local",
  sigma_n =
    "sigma_n_error_local",
  alpha =
    "alpha_error_local",
  alpha_s =
    "alpha_s_error_local"
)


correlation_table <- rbindlist(
  lapply(
    names(
      other_errors
    ),
    function(parameter_i) {

      yy <- corr_source[[other_errors[[parameter_i]]]]

      xx <- corr_source$sigma_c_error

      good <- is.finite(
        xx
      ) &
        is.finite(
          yy
        )

      xx <- xx[
        good
      ]

      yy <- yy[
        good
      ]

      data.table(
        Parameter =
          parameter_i,

        N =
          length(
            xx
          ),

        Pearson_error_correlation =
          if (
            length(
              xx
            ) > 2L &
              sd(
                xx
              ) > 0 &
              sd(
                yy
              ) > 0
          ) {
            cor(
              xx,
              yy,
              method = "pearson"
            )
          } else {
            NA_real_
          },

        Spearman_error_correlation =
          if (
            length(
              unique(
                xx
              )
            ) > 1L &
              length(
                unique(
                  yy
                )
              ) > 1L
          ) {
            cor(
              xx,
              yy,
              method = "spearman"
            )
          } else {
            NA_real_
          }
      )
    }
  )
)


# =============================================================================
# 8. Save tables
# =============================================================================

fwrite(
  recovery,
  file.path(
    OUTPUT_DIR,
    "parameter_recovery_operational_domain.tsv"
  ),
  sep = "\t"
)

fwrite(
  metrics,
  file.path(
    OUTPUT_DIR,
    "parameter_recovery_metrics.tsv"
  ),
  sep = "\t"
)

fwrite(
  correlation_table,
  file.path(
    OUTPUT_DIR,
    "parameter_error_correlations.tsv"
  ),
  sep = "\t"
)


# =============================================================================
# 9. Plot preparation
#
# Plotting every raw point can create a very large PDF. We therefore draw a
# deterministic random subset from each parameter solely for visualization.
# Metrics above always use ALL retained benchmark rows.
# =============================================================================

MAX_POINTS_PER_PARAMETER <- 15000L

set.seed(
  20260923L
)

plot_dt <- recovery[
  ,
  if (
    .N >
      MAX_POINTS_PER_PARAMETER
  ) {
    .SD[
      sample(
        .N,
        MAX_POINTS_PER_PARAMETER
      )
    ]
  } else {
    .SD
  },
  by =
    Parameter
]


parameter_labels <- c(
  sigma_c =
    expression(
      sigma[c]
    ),
  sigma_n =
    expression(
      sigma[n]
    ),
  tau =
    expression(
      tau
    ),
  tau_s =
    expression(
      tau[s]
    ),
  alpha =
    expression(
      alpha
    ),
  alpha_s =
    expression(
      alpha[s]
    )
)


# =============================================================================
# 10. True-versus-estimated figure
# =============================================================================

p <- ggplot(
  plot_dt,
  aes(
    x =
      True,
    y =
      Estimated
  )
) +

  geom_abline(
    slope =
      1,
    intercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.55,
    colour =
      "grey35"
  ) +

  geom_point(
    alpha =
      0.10,
    size =
      0.7
  ) +

  facet_wrap(
    ~Parameter,
    scales =
      "free",
    ncol =
      3,
    labeller =
      as_labeller(
        parameter_labels,
        default =
          label_parsed
      )
  ) +

  labs(
    x =
      "True kinetic rate",
    y =
      "Estimated kinetic rate"
  ) +

  theme_classic(
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
      )
  )


ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "FigS_parameter_recovery_operational_domain.pdf"
    ),
  plot =
    p,
  width =
    9.2,
  height =
    6.4,
  units =
    "in",
  device =
    cairo_pdf
)

ggsave(
  filename =
    file.path(
      OUTPUT_DIR,
      "FigS_parameter_recovery_operational_domain.png"
    ),
  plot =
    p,
  width =
    9.2,
  height =
    6.4,
  units =
    "in",
  dpi =
    400,
  bg =
    "white"
)


# =============================================================================
# 11. Report
# =============================================================================

cat(
  "\nParameter-recovery metrics:\n"
)

print(
  metrics
)


cat(
  "\nCorrelation of sigma_c estimation error with other parameter errors:\n"
)

print(
  correlation_table
)


cat(
  "\nGenerated in:\n  ",
  OUTPUT_DIR,
  "\n",
  sep = ""
)
