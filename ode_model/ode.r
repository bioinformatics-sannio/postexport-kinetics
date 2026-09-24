# =============================================================================
# Title: ODE-Based Simulation Utilities for RNA Kinetics
#
# Description:
#   Four-state RNA kinetics model with support for:
#     - continuous transcription,
#     - complete transcriptional shutoff,
#     - partial / pseudo-shutoff with residual transcription,
#     - replicate-level biological variability,
#     - destructive sampling at selected time points,
#     - gene-specific transcription-onset heterogeneity.
#
# IMPORTANT SEMANTICS
# -------------------
#   The pharmacological intervention time t_star is COMMON to all genes and
#   replicates.
#
#   Gene-specific temporal heterogeneity is implemented by changing the
#   transcriptional HISTORY before t_star, NOT by shifting observations after
#   the intervention has been simulated.
#
#   A positive onset_shift means later transcriptional onset.
#   A negative onset_shift means earlier transcriptional onset; this is
#   implemented by a pre-run before the nominal simulation origin.
#
#   Therefore, for every gene:
#
#       observed time t_star == true intervention time t_star
#
#   and after t_star:
#
#       R_post = post_R_fraction * R_pre
#
#   post_R_fraction = 0   -> complete shutoff
#   post_R_fraction = 0.1 -> 10% residual transcription
#   post_R_fraction = 1   -> no effective shutoff
#
# Backward compatibility:
#   - generate_ODE_states() retains the old argument use_time_shift.
#   - use_time_shift is now interpreted as transcription-onset heterogeneity.
#   - the returned field time_shift is retained as an alias of onset_shift.
#   - apply_observation_shift() is retained only as a legacy helper and is NOT
#     used by generate_ODE_states().
# =============================================================================


library(deSolve)


# =============================================================================
# 1. Four-state ODE
# =============================================================================

rna_kinetics <- function(t, y, params) {

  with(
    as.list(c(y, params)),
    {

      dN <-
        R -
        sigma_n * N -
        tau * N

      dN_s <-
        sigma_n * N -
        tau_s * N_s

      dC <-
        tau * N -
        sigma_c * C -
        alpha * C

      dC_s <-
        tau_s * N_s +
        sigma_c * C -
        alpha_s * C_s

      list(
        c(
          dN,
          dN_s,
          dC,
          dC_s
        )
      )
    }
  )
}


# =============================================================================
# 2. Parameter ranges
# =============================================================================

r_tau_min <- 0.006
r_tau_max <- 0.06

r_tau_s_min <- 0.003
r_tau_s_max <- 0.03

r_alpha_min <- 0.03
r_alpha_max <- 0.69

r_alpha_s_min <- 0.01
r_alpha_s_max <- 0.23

r_sigma_n_min <- 0.05
r_sigma_n_max <- 0.2

r_sigma_c_min <- 0.05
r_sigma_c_max <- 0.2


# =============================================================================
# 3. Closed-form steady state
# =============================================================================

steady_states <- function(params) {

  N <-
    params$R /
    (
      params$tau +
      params$sigma_n
    )

  N_s <-
    params$R *
    params$sigma_n /
    (
      (
        params$tau +
        params$sigma_n
      ) *
      params$tau_s
    )

  C <-
    params$R *
    params$tau /
    (
      (
        params$tau +
        params$sigma_n
      ) *
      (
        params$alpha +
        params$sigma_c
      )
    )

  C_s <-
    (
      params$sigma_n +
      params$sigma_c *
      params$tau /
      (
        params$sigma_c +
        params$alpha
      )
    ) *
    params$R /
    (
      (
        params$tau +
        params$sigma_n
      ) *
      params$alpha_s
    )

  list(
    N = N,
    N_s = N_s,
    C = C,
    C_s = C_s
  )
}


# =============================================================================
# 4. Random parameter generator
# =============================================================================

random_params <- function(
  rtau_min = r_tau_min,
  rtau_max = r_tau_max,

  rtau_s_min = r_tau_s_min,
  rtau_s_max = r_tau_s_max,

  ralpha_min = r_alpha_min,
  ralpha_max = r_alpha_max,

  ralpha_s_min = r_alpha_s_min,
  ralpha_s_max = r_alpha_s_max,

  rsigma_n_min = r_sigma_n_min,
  rsigma_n_max = r_sigma_n_max,

  rsigma_c_min = r_sigma_c_min,
  rsigma_c_max = r_sigma_c_max
) {

  list(
    tau =
      runif(
        1,
        rtau_min,
        rtau_max
      ),

    tau_s =
      runif(
        1,
        rtau_s_min,
        rtau_s_max
      ),

    alpha =
      runif(
        1,
        ralpha_min,
        ralpha_max
      ),

    alpha_s =
      runif(
        1,
        ralpha_s_min,
        ralpha_s_max
      ),

    sigma_n =
      runif(
        1,
        rsigma_n_min,
        rsigma_n_max
      ),

    sigma_c =
      runif(
        1,
        rsigma_c_min,
        rsigma_c_max
      )
  )
}


# =============================================================================
# 5. Internal helper: simulate one constant-parameter latent trajectory
# =============================================================================

simulate_latent_trajectory <- function(
  y0,
  times,
  params,
  model_kinetics = rna_kinetics
) {

  times <- sort(
    unique(
      as.numeric(times)
    )
  )

  if (length(times) < 1L) {
    stop("times must contain at least one value.")
  }

  # deSolve needs at least two output times for actual integration.
  if (length(times) == 1L) {

    out <- data.frame(
      time = times[1],
      N = unname(y0["N"]),
      N_s = unname(y0["N_s"]),
      C = unname(y0["C"]),
      C_s = unname(y0["C_s"])
    )

    return(out)
  }

  out <- deSolve::ode(
    y = y0,
    times = times,
    func = model_kinetics,
    parms = params
  )

  as.data.frame(out)
}


# =============================================================================
# 6. Helper: integrate over one interval with a fixed transcription rate
# =============================================================================

integrate_interval_fixed_R <- function(
  y_start,
  t_start,
  t_end,
  params,
  R_interval,
  model_kinetics = rna_kinetics
) {

  if (t_end < t_start) {
    stop("t_end must be >= t_start.")
  }

  if (t_end == t_start) {
    return(y_start)
  }

  p <- params
  p$R <- R_interval

  out <- deSolve::ode(
    y = y_start,
    times = c(t_start, t_end),
    func = model_kinetics,
    parms = p
  )

  y_end <- as.numeric(
    out[
      nrow(out),
      c(
        "N",
        "N_s",
        "C",
        "C_s"
      )
    ]
  )

  names(y_end) <- c(
    "N",
    "N_s",
    "C",
    "C_s"
  )

  y_end
}


# =============================================================================
# 7. Helper: scheduled transcription trajectory
#
# This is the core correction.
#
# R(t) is defined on the OBSERVED / EXPERIMENTAL time axis:
#
#   before onset_time:
#       R = 0
#
#   onset_time <= t < t_star:
#       R = R_pre
#
#   t >= t_star:
#       R = post_R_fraction * R_pre
#
# If t_star is NULL, transcription remains at R_pre after onset.
#
# Negative onset_time is handled through a pre-run before the simulation
# origin so that the state at time min(times) already contains RNA accumulated
# during the earlier transcriptional history.
# =============================================================================

simulate_scheduled_trajectory <- function(
  y0,
  times,
  params,
  onset_time = 0,
  t_star = NULL,
  post_R_fraction = 0,
  model_kinetics = rna_kinetics
) {

  times <- sort(
    unique(
      as.numeric(times)
    )
  )

  if (length(times) < 1L) {
    stop("times must contain at least one value.")
  }

  if (
    length(onset_time) != 1L ||
    !is.finite(onset_time)
  ) {
    stop("onset_time must be one finite number.")
  }

  if (
    !is.null(t_star) &&
    (
      length(t_star) != 1L ||
      !is.finite(t_star)
    )
  ) {
    stop("t_star must be NULL or one finite number.")
  }

  if (
    length(post_R_fraction) != 1L ||
    !is.finite(post_R_fraction) ||
    post_R_fraction < 0 ||
    post_R_fraction > 1
  ) {
    stop("post_R_fraction must be in [0,1].")
  }

  if (
    is.null(params$R) ||
    length(params$R) != 1L ||
    !is.finite(params$R)
  ) {
    stop("params$R must be one finite transcription rate.")
  }

  t_min <- min(times)
  t_max <- max(times)

  R_pre <- params$R
  R_post <- post_R_fraction * R_pre

  # ---------------------------------------------------------------------------
  # Earlier-than-origin onset:
  #
  # If onset_time < t_min, accumulate a pre-history from onset_time to t_min
  # with normal transcription. This changes the state at the beginning of the
  # observable simulation without changing t_star.
  # ---------------------------------------------------------------------------

  y_start <- y0

  if (onset_time < t_min) {

    # If the intervention itself happened before t_min, the pre-run must be
    # split at t_star. This is uncommon in the intended benchmark but handled
    # correctly for completeness.
    if (
      !is.null(t_star) &&
      t_star > onset_time &&
      t_star < t_min
    ) {

      y_at_star <- integrate_interval_fixed_R(
        y_start = y_start,
        t_start = onset_time,
        t_end = t_star,
        params = params,
        R_interval = R_pre,
        model_kinetics = model_kinetics
      )

      y_start <- integrate_interval_fixed_R(
        y_start = y_at_star,
        t_start = t_star,
        t_end = t_min,
        params = params,
        R_interval = R_post,
        model_kinetics = model_kinetics
      )

    } else {

      R_pre_history <- if (
        !is.null(t_star) &&
        t_min >= t_star
      ) {
        R_post
      } else {
        R_pre
      }

      y_start <- integrate_interval_fixed_R(
        y_start = y_start,
        t_start = onset_time,
        t_end = t_min,
        params = params,
        R_interval = R_pre_history,
        model_kinetics = model_kinetics
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Breakpoints include every requested time plus all rate-change times that
  # fall inside the observable simulation interval.
  # ---------------------------------------------------------------------------

  breakpoints <- times

  if (
    onset_time > t_min &&
    onset_time < t_max
  ) {
    breakpoints <- c(
      breakpoints,
      onset_time
    )
  }

  if (
    !is.null(t_star) &&
    t_star > t_min &&
    t_star < t_max
  ) {
    breakpoints <- c(
      breakpoints,
      t_star
    )
  }

  breakpoints <- sort(
    unique(
      breakpoints
    )
  )

  # If a requested grid does not start at t_min after sorting, this should
  # never happen, but retain an explicit guard.
  if (breakpoints[1] != t_min) {
    breakpoints <- sort(
      unique(
        c(
          t_min,
          breakpoints
        )
      )
    )
  }

  states <- matrix(
    NA_real_,
    nrow = length(breakpoints),
    ncol = 4L
  )

  colnames(states) <- c(
    "N",
    "N_s",
    "C",
    "C_s"
  )

  states[1, ] <- y_start

  current_y <- y_start

  # ---------------------------------------------------------------------------
  # Integrate each interval with the correct constant R.
  #
  # Rate is evaluated at the interval midpoint; because onset_time and t_star
  # are inserted as explicit breakpoints, no interval can cross a rate change.
  # ---------------------------------------------------------------------------

  if (length(breakpoints) >= 2L) {

    for (
      ii in 2:length(breakpoints)
    ) {

      t_left <- breakpoints[ii - 1L]
      t_right <- breakpoints[ii]
      t_mid <- (t_left + t_right) / 2

      R_interval <- if (
        t_mid < onset_time
      ) {

        0

      } else if (
        !is.null(t_star) &&
        t_mid >= t_star
      ) {

        R_post

      } else {

        R_pre
      }

      current_y <- integrate_interval_fixed_R(
        y_start = current_y,
        t_start = t_left,
        t_end = t_right,
        params = params,
        R_interval = R_interval,
        model_kinetics = model_kinetics
      )

      states[ii, ] <- current_y
    }
  }

  out <- data.frame(
    time = breakpoints,
    N = states[, "N"],
    N_s = states[, "N_s"],
    C = states[, "C"],
    C_s = states[, "C_s"]
  )

  # Right-continuous transcription metadata:
  # at t == onset_time transcription is ON;
  # at t == t_star the intervention is already active.
  out$R <- ifelse(
    out$time < onset_time,
    0,
    R_pre
  )

  if (!is.null(t_star)) {
    out$R[
      out$time >= t_star
    ] <- R_post
  }

  # Return only originally requested times.
  out <- out[
    out$time %in% times,
    ,
    drop = FALSE
  ]

  rownames(out) <- NULL

  out
}


# =============================================================================
# 8. Backward-compatible intervention helper
#
# No onset heterogeneity is introduced here. This function now correctly
# handles t_star == min(times), which means the intervention is active from
# the first observed time onward.
# =============================================================================

simulate_piecewise_trajectory <- function(
  y0,
  times,
  params,
  t_star = NULL,
  post_R_fraction = 0,
  model_kinetics = rna_kinetics
) {

  simulate_scheduled_trajectory(
    y0 = y0,
    times = times,
    params = params,
    onset_time = min(
      as.numeric(times)
    ),
    t_star = t_star,
    post_R_fraction = post_R_fraction,
    model_kinetics = model_kinetics
  )
}


# =============================================================================
# 9. Legacy observational time-shift helper
#
# RETAINED ONLY FOR REPRODUCING OLD SIMULATIONS.
# generate_ODE_states() below DOES NOT use this helper.
# =============================================================================

apply_observation_shift <- function(
  latent_df,
  actual_times,
  shift
) {

  latent_df <- as.data.frame(
    latent_df
  )

  actual_times <- sort(
    unique(
      as.numeric(
        actual_times
      )
    )
  )

  lookup_times <-
    actual_times -
    shift

  out <- data.frame(
    time = actual_times,
    N = 0,
    N_s = 0,
    C = 0,
    C_s = 0
  )

  for (
    ii in seq_along(
      actual_times
    )
  ) {

    tt <- lookup_times[ii]

    if (
      tt <
      min(
        latent_df$time
      )
    ) {
      next
    }

    jj <- match(
      tt,
      latent_df$time
    )

    if (is.na(jj)) {
      stop(
        paste(
          "Requested latent time",
          tt,
          "is missing from simulation grid."
        )
      )
    }

    out[
      ii,
      c(
        "N",
        "N_s",
        "C",
        "C_s"
      )
    ] <-
      latent_df[
        jj,
        c(
          "N",
          "N_s",
          "C",
          "C_s"
        )
      ]
  }

  out
}


# =============================================================================
# 10. Main replicate-level simulation
# =============================================================================

generate_ODE_states <- function(
  base_params,
  y0,
  times,

  n_replicates = 3,

  model_kinetics = rna_kinetics,

  param_cv = 0.05,

  stimes = c(
    10,
    40,
    50,
    100
  ),

  shutofftimes = c(
    10,
    40,
    50,
    100
  ),

  max_shift = NULL,

  t_star = NULL,

  post_R_fraction = 0,

  # ---------------------------------------------------------------------------
  # Backward-compatible name.
  #
  # IMPORTANT:
  #   This now controls GENE-SPECIFIC TRANSCRIPTION-ONSET heterogeneity.
  #   It no longer shifts the observed intervention trajectory.
  # ---------------------------------------------------------------------------
  use_time_shift = TRUE,

  # Preferred explicit name. If non-NULL, overrides use_time_shift.
  use_onset_shift = NULL,

  # Nominal transcription onset on the simulation time axis.
  nominal_onset_time = 0
) {


  # ---------------------------------------------------------------------------
  # Validation
  # ---------------------------------------------------------------------------

  if (
    length(post_R_fraction) != 1L ||
    !is.finite(post_R_fraction) ||
    post_R_fraction < 0 ||
    post_R_fraction > 1
  ) {
    stop(
      "post_R_fraction must be in [0,1]."
    )
  }

  if (
    length(n_replicates) != 1L ||
    !is.finite(n_replicates) ||
    n_replicates < 1L
  ) {
    stop(
      "n_replicates must be >= 1."
    )
  }

  if (
    length(nominal_onset_time) != 1L ||
    !is.finite(nominal_onset_time)
  ) {
    stop(
      "nominal_onset_time must be one finite number."
    )
  }

  actual_times <- sort(
    unique(
      as.numeric(
        times
      )
    )
  )

  if (
    length(actual_times) < 2L
  ) {
    stop(
      "At least two simulation time points are required."
    )
  }

  if (
    !is.null(t_star) &&
    (
      length(t_star) != 1L ||
      !is.finite(t_star)
    )
  ) {
    stop(
      "t_star must be NULL or one finite number."
    )
  }


  # ---------------------------------------------------------------------------
  # Biological-variability scale.
  # ---------------------------------------------------------------------------

  sdlog <- if (
    param_cv > 0
  ) {

    sqrt(
      log(
        1 +
        param_cv^2
      )
    )

  } else {

    0
  }


  # ---------------------------------------------------------------------------
  # Onset-shift semantics.
  # ---------------------------------------------------------------------------

  shift_enabled <- if (
    is.null(
      use_onset_shift
    )
  ) {
    isTRUE(
      use_time_shift
    )
  } else {
    isTRUE(
      use_onset_shift
    )
  }

  if (
    is.null(max_shift)
  ) {

    max_shift <- floor(
      0.1 *
      max(
        actual_times
      )
    )
  }

  if (
    length(max_shift) != 1L ||
    !is.finite(max_shift) ||
    max_shift < 0
  ) {
    stop(
      "max_shift must be one finite value >= 0."
    )
  }

  max_shift <- as.integer(
    floor(
      max_shift
    )
  )

  onset_shift <- if (
    shift_enabled &&
    max_shift > 0L
  ) {

    sample(
      -max_shift:max_shift,
      1
    )

  } else {

    0L
  }

  onset_time <-
    nominal_onset_time +
    onset_shift


  # ---------------------------------------------------------------------------
  # The actual output grid is NOT shifted.
  #
  # Add requested sampling times and t_star to the integration grid so that
  # all requested output times and intervention boundaries are represented
  # exactly.
  # ---------------------------------------------------------------------------

  simulation_times <- sort(
    unique(
      c(
        actual_times,
        stimes,
        shutofftimes,
        if (
          !is.null(
            t_star
          )
        ) {
          t_star
        } else {
          numeric()
        }
      )
    )
  )

  # Only integrate across the requested simulation interval.
  simulation_times <- simulation_times[
    simulation_times >=
      min(
        actual_times
      ) &
    simulation_times <=
      max(
        actual_times
      )
  ]

  simulation_times <- sort(
    unique(
      c(
        min(
          actual_times
        ),
        simulation_times,
        max(
          actual_times
        )
      )
    )
  )


  # ---------------------------------------------------------------------------
  # Containers.
  # ---------------------------------------------------------------------------

  data_list <- vector(
    "list",
    n_replicates
  )

  intervention_list <- vector(
    "list",
    n_replicates
  )

  parameters_list <- vector(
    "list",
    n_replicates
  )

  steady_state_list <- vector(
    "list",
    n_replicates
  )


  # ===========================================================================
  # Replicate loop
  # ===========================================================================

  for (
    i in seq_len(
      n_replicates
    )
  ) {


    # -------------------------------------------------------------------------
    # Biological parameter variation.
    # -------------------------------------------------------------------------

    perturbed_params <- lapply(
      base_params,
      function(param) {

        if (
          sdlog > 0
        ) {

          param *
          exp(
            stats::rnorm(
              1,
              0,
              sdlog
            )
          )

        } else {

          param
        }
      }
    )


    # Keep transcription rate fixed across replicates, as in the old code.
    perturbed_params$R <-
      base_params$R


    # -------------------------------------------------------------------------
    # Parameter bookkeeping.
    # -------------------------------------------------------------------------

    parameters_list[[i]] <-
      as.data.frame(
        c(
          perturbed_params,
          replicate = i
        )
      )


    # -------------------------------------------------------------------------
    # A. Baseline / continuous-transcription trajectory.
    #
    # Uses the SAME gene-specific onset_time as the intervention trajectory.
    # No post-hoc observation shift is applied.
    # -------------------------------------------------------------------------

    observed_none <- simulate_scheduled_trajectory(
      y0 = y0,
      times = simulation_times,
      params = perturbed_params,
      onset_time = onset_time,
      t_star = NULL,
      post_R_fraction = 1,
      model_kinetics = model_kinetics
    )


    # -------------------------------------------------------------------------
    # B. Intervention trajectory.
    #
    # SAME pre-shutoff onset history as baseline.
    # The intervention ALWAYS occurs at the common t_star.
    # -------------------------------------------------------------------------

    observed_intervention <- simulate_scheduled_trajectory(
      y0 = y0,
      times = simulation_times,
      params = perturbed_params,
      onset_time = onset_time,
      t_star = t_star,
      post_R_fraction = post_R_fraction,
      model_kinetics = model_kinetics
    )


    # -------------------------------------------------------------------------
    # Replicate identifiers.
    # -------------------------------------------------------------------------

    observed_none$replicate <- i

    observed_intervention$replicate <- i


    # -------------------------------------------------------------------------
    # Steady state before intervention.
    # -------------------------------------------------------------------------

    ssi <- steady_states(
      perturbed_params
    )


    steady_state_list[[i]] <-
      data.frame(
        replicate = i,
        R = perturbed_params$R,
        N = ssi$N,
        N_s = ssi$N_s,
        C = ssi$C,
        C_s = ssi$C_s
      )


    # -------------------------------------------------------------------------
    # Store.
    # -------------------------------------------------------------------------

    data_list[[i]] <-
      observed_none

    intervention_list[[i]] <-
      observed_intervention
  }


  # =============================================================================
  # 11. Combine replicates
  # =============================================================================

  data <- do.call(
    rbind,
    data_list
  )

  intervention_data <- do.call(
    rbind,
    intervention_list
  )

  parameters <- do.call(
    rbind,
    parameters_list
  )

  steady_state <- do.call(
    rbind,
    steady_state_list
  )

  rownames(data) <- NULL
  rownames(intervention_data) <- NULL
  rownames(parameters) <- NULL
  rownames(steady_state) <- NULL


  # =============================================================================
  # 12. Sample requested observation times
  # =============================================================================

  df_tsampled <- data[
    data$time %in%
      stimes,
    ,
    drop = FALSE
  ]

  df_intervention_tsampled <-
    intervention_data[
      intervention_data$time %in%
        shutofftimes,
      ,
      drop = FALSE
    ]


  # =============================================================================
  # 13. Return
  # =============================================================================

  list(

    # Original names retained for backward compatibility.
    data =
      data,

    shutoff_data =
      intervention_data,

    tsampled_data =
      df_tsampled,

    shutoff_tsampled_data =
      df_intervention_tsampled,


    # Clearer aliases.
    intervention_data =
      intervention_data,

    intervention_tsampled_data =
      df_intervention_tsampled,


    ss_data =
      steady_state,

    # Backward-compatible alias.
    time_shift =
      onset_shift,

    # Preferred explicit metadata.
    onset_shift =
      onset_shift,

    onset_time =
      onset_time,

    nominal_onset_time =
      nominal_onset_time,

    use_onset_shift =
      shift_enabled,

    parameters =
      parameters,

    post_R_fraction =
      post_R_fraction,

    t_star =
      t_star
  )
}
