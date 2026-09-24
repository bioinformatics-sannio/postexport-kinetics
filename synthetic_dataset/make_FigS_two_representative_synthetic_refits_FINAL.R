# =============================================================================
# Supplementary Figure: two representative noisy synthetic refits
#
# Final version for Supplementary Material.
#
# Shows, for two objectively selected alternatives from the SAME calibrated
# SHUTOFF configuration:
#
#   A-B. Representative detected alternative
#   C-D. Representative non-detected alternative
#
# Trajectory panels include:
#   - generating trajectory
#   - full-model fit
#   - null-model fit
#   - noisy replicate-level observations
#
# Parameter panels include:
#   - truth
#   - full fit
#   - null fit
#
# Inputs:
#   ode_states_2k_20p_corrected_onset.rdata
#   benchmark_main_corrected_onset_raw.tsv
#   benchmark_main_corrected_onset_summary.tsv
#
# Outputs:
#   practical_identifiability/FigS_two_representative_synthetic_refits_final.pdf
#   practical_identifiability/FigS_two_representative_synthetic_refits_final.png
#   practical_identifiability/two_representative_synthetic_refits_selection.tsv
#
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/synthetic_dataset")

library(data.table)
library(ggplot2)
library(patchwork)
library(deSolve)

source("../commons/platforms.r")


# =============================================================================
# 1. Files
# =============================================================================

ODE_FILE <- "ode_states_2k_20p_corrected_onset.rdata"
RAW_FILE <- "benchmark_main_corrected_onset_raw.tsv"
SUMMARY_FILE <- "benchmark_main_corrected_onset_summary.tsv"

OUTPUT_DIR <- "practical_identifiability"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE
  )
}

OUT_PDF <- file.path(
  OUTPUT_DIR,
  "FigS_two_representative_synthetic_refits_final.pdf"
)

OUT_PNG <- file.path(
  OUTPUT_DIR,
  "FigS_two_representative_synthetic_refits_final.png"
)

OUT_SELECTION <- file.path(
  OUTPUT_DIR,
  "two_representative_synthetic_refits_selection.tsv"
)


# =============================================================================
# 2. Load data
# =============================================================================

load(ODE_FILE)

if (
  !exists("ode_states") ||
    !data.table::is.data.table(ode_states)
) {
  stop("Object 'ode_states' was not found as a data.table.")
}

raw <- fread(RAW_FILE)
summary_dt <- fread(SUMMARY_FILE)


# =============================================================================
# 3. Validate required columns
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
  "Gene",
  "Positive",
  "status",
  "p.value",
  "T_star",
  "measurement_seed",
  "sigma_true",
  "sigma_c_hat",
  "sigma_n_true",
  "sigma_n_hat",
  "tau_true",
  "tau_hat",
  "tau_s_true",
  "tau_s_hat",
  "alpha_true",
  "alpha_hat",
  "alpha_s_true",
  "alpha_s_hat",
  "R_true",
  "R_hat",
  "R_hat_null",
  "sigma_n_hat_null",
  "tau_hat_null",
  "tau_s_hat_null",
  "alpha_hat_null",
  "alpha_s_hat_null"
)

required_summary <- c(
  config_cols,
  "TypeI_005",
  "TypeI_005_CI_contains_005"
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

missing_raw <- setdiff(required_raw, names(raw))
missing_summary <- setdiff(required_summary, names(summary_dt))
missing_ode <- setdiff(required_ode, names(ode_states))

if (length(missing_raw) > 0L) {
  stop(
    paste0(
      "Missing columns in raw benchmark:\n  ",
      paste(missing_raw, collapse = "\n  ")
    )
  )
}

if (length(missing_summary) > 0L) {
  stop(
    paste0(
      "Missing columns in benchmark summary:\n  ",
      paste(missing_summary, collapse = "\n  ")
    )
  )
}

if (length(missing_ode) > 0L) {
  stop(
    paste0(
      "Missing columns in ode_states:\n  ",
      paste(missing_ode, collapse = "\n  ")
    )
  )
}


# =============================================================================
# 4. Observation model used in the principal benchmark
# =============================================================================

RANGE_GAUSS_NOISE <- c(
  "Very low" = 0.02,
  "Low"      = 0.05,
  "Medium"   = 0.10,
  "High"     = 0.20
)

add_platform_noise_main <- function(
  dt,
  platform,
  noise
) {

  dt <- copy(as.data.table(dt))

  targets <- c(
    "N",
    "C",
    "C_s",
    "N_s"
  )

  if (platform == "GAUSS") {

    noise_sd <- unname(
      RANGE_GAUSS_NOISE[noise]
    )

    return(
      add_gaussian_noise(
        dt,
        cols = targets,
        noise_sd = noise_sd
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

  stop(
    paste(
      "Unknown platform:",
      platform
    )
  )
}


# =============================================================================
# 5. Choose one calibrated SHUTOFF configuration
# =============================================================================

operational <- summary_dt[
  Perturbation == "SHUTOFF" &
    TypeI_005_CI_contains_005 == TRUE
]

if (nrow(operational) == 0L) {
  stop("No CI-calibrated SHUTOFF configuration was found.")
}

operational[
  ,
  config_distance :=
    ifelse(
      Platform == "RNA-seq",
      0,
      10
    ) +
    ifelse(
      Exprs_noise == "Medium",
      0,
      5
    ) +
    abs(
      N_tsamples - 5
    ) +
    0.5 * abs(
      N_replicates - 10
    ) +
    0.1 * abs(
      Tsteps - 10
    )
]

operational[
  ,
  calibration_distance :=
    abs(
      TypeI_005 - 0.05
    )
]

setorder(
  operational,
  config_distance,
  calibration_distance
)

chosen_config <- operational[1]


cat(
  "\n============================================================\n",
  "TWO REPRESENTATIVE SYNTHETIC REFITS\n",
  "============================================================\n",
  "Chosen calibrated configuration:\n",
  sep = ""
)

print(
  chosen_config[
    ,
    .(
      Platform,
      Exprs_noise,
      N_tsamples,
      N_replicates,
      Tsteps,
      TypeI_005,
      TypeI_005_CI_contains_005
    )
  ]
)


# =============================================================================
# 6. Candidate positive genes
# =============================================================================

candidate <- raw[
  Perturbation == "SHUTOFF" &
    Positive == 1 &
    status == "ok" &
    Platform == chosen_config$Platform &
    Exprs_noise == chosen_config$Exprs_noise &
    N_tsamples == chosen_config$N_tsamples &
    N_replicates == chosen_config$N_replicates &
    Tsteps == chosen_config$Tsteps &
    is.finite(p.value) &
    is.finite(sigma_true) &
    is.finite(sigma_c_hat)
]

if (nrow(candidate) == 0L) {
  stop(
    "No positive fitted genes found in the selected configuration."
  )
}

candidate[
  ,
  abs_sigma_error :=
    abs(
      sigma_c_hat - sigma_true
    )
]


# =============================================================================
# 7. Objective selection of detected and non-detected examples
# =============================================================================

candidate_detected <- candidate[
  p.value <= 0.05 &
    sigma_c_hat > 1e-12
]

candidate_nondetected <- candidate[
  p.value > 0.05 |
    sigma_c_hat <= 1e-12
]

if (nrow(candidate_detected) == 0L) {
  stop(
    "No detected alternative is available in the chosen configuration."
  )
}

if (nrow(candidate_nondetected) == 0L) {
  stop(
    "No non-detected alternative is available in the chosen configuration."
  )
}

select_median_error_case <- function(
  dt
) {

  out <- copy(dt)

  med <- median(
    out$abs_sigma_error,
    na.rm = TRUE
  )

  out[
    ,
    distance_from_group_median :=
      abs(
        abs_sigma_error - med
      )
  ]

  setorder(
    out,
    distance_from_group_median,
    Gene
  )

  out[1]
}

selected_detected <- select_median_error_case(
  candidate_detected
)

selected_nondetected <- select_median_error_case(
  candidate_nondetected
)

selected_detected[
  ,
  Example :=
    "Detected alternative"
]

selected_nondetected[
  ,
  Example :=
    "Non-detected alternative"
]

selected_cases <- rbind(
  selected_detected,
  selected_nondetected,
  fill = TRUE
)


cat(
  "\nSelected examples:\n"
)

print(
  selected_cases[
    ,
    .(
      Example,
      Gene,
      p.value,
      sigma_true,
      sigma_c_hat,
      abs_sigma_error,
      measurement_seed
    )
  ]
)


# =============================================================================
# 8. ODE propagation helper
# =============================================================================

state_cols <- c(
  "N",
  "N_s",
  "C",
  "C_s"
)

ode_rhs <- function(
  t,
  y,
  pars
) {

  R_t <- if (
    is.finite(
      pars["t_star"]
    ) &&
      t >= pars["t_star"]
  ) {
    0
  } else {
    pars["R"]
  }

  N <- y["N"]
  N_s <- y["N_s"]
  C <- y["C"]
  C_s <- y["C_s"]

  dN <-
    R_t -
      pars["sigma_n"] * N -
      pars["tau"] * N

  dN_s <-
    pars["sigma_n"] * N -
      pars["tau_s"] * N_s

  dC <-
    pars["tau"] * N -
      pars["sigma_c"] * C -
      pars["alpha"] * C

  dC_s <-
    pars["tau_s"] * N_s +
      pars["sigma_c"] * C -
      pars["alpha_s"] * C_s

  list(
    c(
      dN,
      dN_s,
      dC,
      dC_s
    )
  )
}

propagate_fit <- function(
  pars,
  initial_state,
  time_grid
) {

  fit <- deSolve::ode(
    y = initial_state,
    times = time_grid,
    func = ode_rhs,
    parms = pars,
    method = "lsoda",
    rtol = 1e-9,
    atol = 1e-11
  )

  fit <- as.data.table(fit)

  melt(
    fit,
    id.vars = "time",
    measure.vars = state_cols,
    variable.name = "State",
    value.name = "Value"
  )
}


# =============================================================================
# 9. Build plotting data for one selected case
# =============================================================================

build_case_data <- function(
  selected_row
) {

  latent <- ode_states[
    Gene == selected_row$Gene &
      Perturbation == "SHUTOFF" &
      N_time_samples == selected_row$N_tsamples &
      N_replicates == selected_row$N_replicates &
      T_step == selected_row$Tsteps &
      abs(
        Post_R_fraction
      ) < 1e-12
  ]

  if (nrow(latent) == 0L) {
    stop(
      paste(
        "Could not recover latent trajectory for gene",
        selected_row$Gene
      )
    )
  }

  set.seed(
    as.integer(
      selected_row$measurement_seed
    )
  )

  noisy <- add_platform_noise_main(
    dt = latent,
    platform = selected_row$Platform,
    noise = selected_row$Exprs_noise
  )


  latent_long <- melt(
    copy(latent),
    id.vars = c(
      "time",
      "replicate"
    ),
    measure.vars = state_cols,
    variable.name = "State",
    value.name = "Value"
  )

  latent_mean <- latent_long[
    ,
    .(
      Value =
        mean(
          Value,
          na.rm = TRUE
        )
    ),
    by = .(
      time,
      State
    )
  ]

  latent_mean[
    ,
    Model :=
      "Generating trajectory"
  ]


  noisy_long <- melt(
    copy(noisy),
    id.vars = c(
      "time",
      "replicate"
    ),
    measure.vars = state_cols,
    variable.name = "State",
    value.name = "Value"
  )


  noisy_summary <- noisy_long[
    ,
    .(
      Mean =
        mean(
          Value,
          na.rm = TRUE
        ),
      SD =
        sd(
          Value,
          na.rm = TRUE
        )
    ),
    by = .(
      time,
      State
    )
  ]


  first_time <- min(
    noisy$time
  )

  initial_state_dt <- noisy[
    time == first_time,
    .(
      N =
        mean(
          N,
          na.rm = TRUE
        ),
      N_s =
        mean(
          N_s,
          na.rm = TRUE
        ),
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
    )
  ]

  initial_state <- c(
    N = initial_state_dt$N[1],
    N_s = initial_state_dt$N_s[1],
    C = initial_state_dt$C[1],
    C_s = initial_state_dt$C_s[1]
  )


  t_star_i <- as.numeric(
    selected_row$T_star
  )

  if (!is.finite(t_star_i)) {
    stop(
      paste(
        "Non-finite T_star for gene",
        selected_row$Gene
      )
    )
  }


  pars_full <- c(
    R = selected_row$R_hat,
    sigma_n = selected_row$sigma_n_hat,
    tau = selected_row$tau_hat,
    tau_s = selected_row$tau_s_hat,
    sigma_c = selected_row$sigma_c_hat,
    alpha = selected_row$alpha_hat,
    alpha_s = selected_row$alpha_s_hat,
    t_star = t_star_i
  )


  pars_null <- c(
    R = selected_row$R_hat_null,
    sigma_n = selected_row$sigma_n_hat_null,
    tau = selected_row$tau_hat_null,
    tau_s = selected_row$tau_s_hat_null,
    sigma_c = 0,
    alpha = selected_row$alpha_hat_null,
    alpha_s = selected_row$alpha_s_hat_null,
    t_star = t_star_i
  )


  if (
    any(
      !is.finite(
        pars_full
      )
    ) ||
      any(
        !is.finite(
          pars_null
        )
      )
  ) {
    stop(
      paste(
        "Non-finite fitted parameter for gene",
        selected_row$Gene
      )
    )
  }


  time_grid <- seq(
    min(
      noisy$time
    ),
    max(
      noisy$time
    ),
    length.out = 400
  )


  full_fit <- propagate_fit(
    pars = pars_full,
    initial_state = initial_state,
    time_grid = time_grid
  )

  full_fit[
    ,
    Model :=
      "Full fit"
  ]


  null_fit <- propagate_fit(
    pars = pars_null,
    initial_state = initial_state,
    time_grid = time_grid
  )

  null_fit[
    ,
    Model :=
      "Null fit"
  ]


  fit_long <- rbind(
    latent_mean,
    full_fit,
    null_fit,
    fill = TRUE
  )


  fit_long[
    ,
    Model :=
      factor(
        Model,
        levels = c(
          "Generating trajectory",
          "Full fit",
          "Null fit"
        )
      )
  ]


  parameter_dt <- rbindlist(
    list(

      data.table(
        Parameter = c(
          "sigma_n",
          "sigma_c",
          "tau",
          "tau_s",
          "alpha",
          "alpha_s"
        ),
        Estimate = c(
          selected_row$sigma_n_true,
          selected_row$sigma_true,
          selected_row$tau_true,
          selected_row$tau_s_true,
          selected_row$alpha_true,
          selected_row$alpha_s_true
        ),
        Model = "Truth"
      ),

      data.table(
        Parameter = c(
          "sigma_n",
          "sigma_c",
          "tau",
          "tau_s",
          "alpha",
          "alpha_s"
        ),
        Estimate = c(
          selected_row$sigma_n_hat,
          selected_row$sigma_c_hat,
          selected_row$tau_hat,
          selected_row$tau_s_hat,
          selected_row$alpha_hat,
          selected_row$alpha_s_hat
        ),
        Model = "Full fit"
      ),

      data.table(
        Parameter = c(
          "sigma_n",
          "sigma_c",
          "tau",
          "tau_s",
          "alpha",
          "alpha_s"
        ),
        Estimate = c(
          selected_row$sigma_n_hat_null,
          0,
          selected_row$tau_hat_null,
          selected_row$tau_s_hat_null,
          selected_row$alpha_hat_null,
          selected_row$alpha_s_hat_null
        ),
        Model = "Null fit"
      )
    )
  )


  parameter_dt[
    ,
    Parameter :=
      factor(
        Parameter,
        levels = c(
          "sigma_n",
          "sigma_c",
          "tau",
          "tau_s",
          "alpha",
          "alpha_s"
        )
      )
  ]

  parameter_dt[
    ,
    Model :=
      factor(
        Model,
        levels = c(
          "Truth",
          "Full fit",
          "Null fit"
        )
      )
  ]


  list(
    noisy_long = noisy_long,
    noisy_summary = noisy_summary,
    fit_long = fit_long,
    parameter_dt = parameter_dt,
    t_star = t_star_i
  )
}


# =============================================================================
# 10. Plot one selected case
# =============================================================================

plot_case <- function(
  selected_row,
  case_data
) {

  line_values <- c(
    "Generating trajectory" = "dashed",
    "Full fit" = "solid",
    "Null fit" = "dotdash"
  )


  p_traj <- ggplot() +

    geom_line(
      data = case_data$fit_long,
      aes(
        x = time,
        y = Value,
        linetype = Model,
        group = Model
      ),
      linewidth = 0.78
    ) +

    geom_errorbar(
      data = case_data$noisy_summary,
      aes(
        x = time,
        ymin = pmax(
          0,
          Mean - SD
        ),
        ymax = Mean + SD
      ),
      width = 0.7,
      linewidth = 0.3
    ) +

    geom_point(
      data = case_data$noisy_long,
      aes(
        x = time,
        y = Value
      ),
      alpha = 0.28,
      size = 1.0,
      position = position_jitter(
        width = 0.14,
        height = 0
      )
    ) +

    geom_point(
      data = case_data$noisy_summary,
      aes(
        x = time,
        y = Mean
      ),
      size = 1.8,
      shape = 21,
      fill = "white",
      stroke = 0.5
    ) +

    geom_vline(
      xintercept = case_data$t_star,
      linetype = "dotted",
      linewidth = 0.45
    ) +

    facet_wrap(
      ~State,
      scales = "free_y",
      ncol = 2
    ) +

    scale_linetype_manual(
      values = line_values,
      drop = FALSE
    ) +

    labs(
      x = "Time",
      y = "Abundance",
      linetype = NULL
    ) +

    theme_classic(
      base_size = 9.5
    ) +

    theme(
      panel.grid.major =
        element_line(
          colour = "grey93",
          linewidth = 0.22
        ),

      panel.grid.minor =
        element_blank(),

      strip.background =
        element_blank(),

      strip.text =
        element_text(
          face = "bold"
        ),

      legend.position =
        "bottom",

      legend.key.width =
        grid::unit(
          1.5,
          "cm"
        )
    )


  p_param <- ggplot(
    case_data$parameter_dt,
    aes(
      x = Model,
      y = Estimate
    )
  ) +

    geom_col(
      width = 0.68
    ) +

    facet_wrap(
      ~Parameter,
      scales = "free_y",
      ncol = 3,
      labeller =
        labeller(
          Parameter = c(
            sigma_n = "sigma[n]",
            sigma_c = "sigma[c]",
            tau = "tau",
            tau_s = "tau[s]",
            alpha = "alpha",
            alpha_s = "alpha[s]"
          ),
          .default =
            label_parsed
        )
    ) +

    labs(
      x = NULL,
      y =
        expression(
          "Kinetic rate (min"^{-1}*")"
        )
    ) +

    theme_classic(
      base_size = 9.5
    ) +

    theme(
      panel.grid.major.y =
        element_line(
          colour = "grey93",
          linewidth = 0.22
        ),

      panel.grid.minor =
        element_blank(),

      strip.background =
        element_blank(),

      strip.text =
        element_text(
          face = "bold"
        ),

      axis.text.x =
        element_text(
          angle = 35,
          hjust = 1
        )
    )


  (
    p_traj /
      p_param
  ) +
    plot_layout(
      heights = c(
        1.5,
        1
      )
    )
}


# =============================================================================
# 11. Build both examples
# =============================================================================

detected_data <- build_case_data(
  selected_detected
)

nondetected_data <- build_case_data(
  selected_nondetected
)

p_detected <- plot_case(
  selected_detected,
  detected_data
)

p_nondetected <- plot_case(
  selected_nondetected,
  nondetected_data
)


# =============================================================================
# 12. Final figure
# =============================================================================

left_title <- paste0(
  "Detected alternative: gene ",
  selected_detected$Gene,
  "; p=",
  format(
    selected_detected$p.value,
    digits = 3
  ),
  "; true sigma_c=",
  format(
    selected_detected$sigma_true,
    digits = 3
  ),
  "; fitted sigma_c=",
  format(
    selected_detected$sigma_c_hat,
    digits = 3
  )
)

right_title <- paste0(
  "Non-detected alternative: gene ",
  selected_nondetected$Gene,
  "; p=",
  format(
    selected_nondetected$p.value,
    digits = 3
  ),
  "; true sigma_c=",
  format(
    selected_nondetected$sigma_true,
    digits = 3
  ),
  "; fitted sigma_c=",
  format(
    selected_nondetected$sigma_c_hat,
    digits = 3
  )
)

p_detected <- p_detected +
  plot_annotation(
    title = left_title,
    tag_levels = "A"
  )

p_nondetected <- p_nondetected +
  plot_annotation(
    title = right_title,
    tag_levels = "C"
  )


config_subtitle <- paste0(
  chosen_config$Platform,
  ", ",
  chosen_config$Exprs_noise,
  " noise; ",
  chosen_config$N_tsamples,
  " time points, ",
  chosen_config$N_replicates,
  " replicates, T_step=",
  chosen_config$Tsteps,
  "; configuration Type-I=",
  format(
    chosen_config$TypeI_005,
    digits = 3
  )
)


final_plot <- (
  p_detected |
    p_nondetected
) +
  plot_annotation(
    subtitle = config_subtitle
  )


# =============================================================================
# 13. Save
# =============================================================================

ggsave(
  filename = OUT_PDF,
  plot = final_plot,
  width = 16.5,
  height = 10.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = OUT_PNG,
  plot = final_plot,
  width = 16.5,
  height = 10.5,
  units = "in",
  dpi = 400,
  bg = "white"
)


# =============================================================================
# 14. Save selection metadata
# =============================================================================

selection_out <- selected_cases[
  ,
  .(
    Example,
    Gene,
    Platform,
    Exprs_noise,
    N_tsamples,
    N_replicates,
    Tsteps,
    TypeI_configuration =
      chosen_config$TypeI_005,
    p.value,
    sigma_true,
    sigma_c_hat,
    abs_sigma_error,
    measurement_seed,
    T_star
  )
]

fwrite(
  selection_out,
  OUT_SELECTION,
  sep = "\t"
)


# =============================================================================
# 15. Report
# =============================================================================

cat(
  "\nGenerated:\n",
  "  ",
  OUT_PDF,
  "\n",
  "  ",
  OUT_PNG,
  "\n",
  "  ",
  OUT_SELECTION,
  "\n",
  sep = ""
)

cat(
  "\nFinal selected examples:\n"
)

print(
  selection_out
)
