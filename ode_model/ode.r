# =============================================================================
# Title: ODE-Based Simulation Utilities for RNA Kinetics
# Description:
#   This script defines the compartmental RNA kinetics model and helper
#   functions to:
#     - compute its ODE dynamics,
#     - derive steady-state values in closed form,
#     - sample random parameter sets,
#     - generate replicate-level simulated trajectories,
#     - optionally simulate transcriptional shutoff,
#     - extract sampled observations at selected time points.
#
# Included functionality:
#   - four-state ODE system,
#   - parameter-range definitions,
#   - steady-state formulas,
#   - random parameter generation,
#   - replicate simulation with biological variability,
#   - optional intervention at time t_star,
#   - time subsampling for downstream inference.
#
# Intended use:
#   Research code accompanying the paper and shared for transparency/review.
#
# Copyright (c) 2026 Luigi Cerulo
#
# Permission is hereby granted to use, copy, modify, and distribute this
# software for academic and research purposes, provided that this notice is
# retained in all copies.
#
# This software is provided "as is", without warranty of any kind, express
# or implied, including but not limited to the warranties of merchantability,
# fitness for a particular purpose, and noninfringement. In no event shall
# the authors be liable for any claim, damages, or other liability arising
# from, out of, or in connection with the software or its use.
# =============================================================================


library(deSolve)


# -----------------------------------------------------------------------------
# Four-state ODE system for RNA kinetics.
#
# State variables:
#   - N   : nuclear unprocessed RNA
#   - N_s : nuclear processed RNA
#   - C   : cytoplasmic unprocessed RNA
#   - C_s : cytoplasmic processed RNA
#
# Parameters:
#   - R        : transcription/input rate
#   - sigma_n  : nuclear processing rate
#   - tau      : export/transition rate from N to C
#   - tau_s    : export/transition rate from N_s to C_s
#   - sigma_c  : cytoplasmic processing rate
#   - alpha    : degradation rate of C
#   - alpha_s  : degradation rate of C_s
#
# Output:
#   A list containing the derivatives in the order expected by deSolve::ode().
# -----------------------------------------------------------------------------
rna_kinetics <- function(t, y, params) {
  with(as.list(c(y, params)), {
    dN <- R - sigma_n * N - tau * N
    dN_s <- sigma_n * N - tau_s * N_s
    dC <- tau * N - sigma_c * C - alpha * C
    dC_s <- tau_s * N_s + sigma_c * C - alpha_s * C_s

    # Optional debug line:
    # print(paste("Time:", t, "dN:", dN, "dN_s:", dN_s, "dC:", dC, "dC_s:", dC_s))

    list(c(dN, dN_s, dC, dC_s))
  })
}


# -----------------------------------------------------------------------------
# Default parameter ranges used for simulation.
#
# These ranges define the support of the uniform distributions sampled by
# random_params(). They can be overridden when needed.
# -----------------------------------------------------------------------------
r_tau_min = 0.006
r_tau_max = 0.06
r_tau_s_min = 0.003
r_tau_s_max = 0.03

r_alpha_min = 0.03
r_alpha_max = 0.69
r_alpha_s_min = 0.01
r_alpha_s_max = 0.23

r_sigma_n_min = 0.05
r_sigma_n_max = 0.2
r_sigma_c_min = 0.05
r_sigma_c_max = 0.2


# -----------------------------------------------------------------------------
# Closed-form steady states of the ODE system.
#
# Given a parameter list containing:
#   R, tau, tau_s, alpha, alpha_s, sigma_n, sigma_c
# this function returns the steady-state values of:
#   N, N_s, C, C_s
#
# These formulas are useful for diagnostics and for checking that simulated
# trajectories are consistent with the chosen parameter regime.
# -----------------------------------------------------------------------------
steady_states <- function(params) {
  N = params$R / (params$tau + params$sigma_n)

  N_s = params$R * params$sigma_n /
    ((params$tau + params$sigma_n) * params$tau_s)

  C = params$R * params$tau /
    ((params$tau + params$sigma_n) * (params$alpha + params$sigma_c))

  C_s = (params$sigma_n + params$sigma_c * params$tau / (params$sigma_c + params$alpha)) *
    params$R / ((params$tau + params$sigma_n) * params$alpha_s)

  return(list(N = N, N_s = N_s, C = C, C_s = C_s))
}


# -----------------------------------------------------------------------------
# Sample one random parameter set from uniform ranges.
#
# The function samples:
#   - tau, tau_s,
#   - alpha, alpha_s,
#   - sigma_n, sigma_c
#
# It does not sample R, which is expected to be specified separately when
# constructing the full parameter list.
# -----------------------------------------------------------------------------
random_params = function(
  rtau_min = r_tau_min, rtau_max = r_tau_max,
  rtau_s_min = r_tau_s_min, rtau_s_max = r_tau_s_max,
  ralpha_min = r_alpha_min, ralpha_max = r_alpha_max,
  ralpha_s_min = r_alpha_s_min, ralpha_s_max = r_alpha_s_max,
  rsigma_n_min = r_sigma_n_min, rsigma_n_max = r_sigma_n_max,
  rsigma_c_min = r_sigma_c_min, rsigma_c_max = r_sigma_c_max
) {

  r_tau = runif(1, rtau_min, rtau_max)
  r_tau_s = runif(1, rtau_s_min, rtau_s_max)
  r_alpha = runif(1, ralpha_min, ralpha_max)
  r_alpha_s = runif(1, ralpha_s_min, ralpha_s_max)

  # Optional alternative coupling between parameters:
  # r_rho = runif(1, 0.1, 1)
  # r_alpha_s = r_rho * r_alpha

  r_sigma_n = runif(1, rsigma_n_min, rsigma_n_max)
  r_sigma_c = runif(1, rsigma_c_min, rsigma_c_max)

  retval = list(
    tau = r_tau,
    tau_s = r_tau_s,
    alpha = r_alpha,
    alpha_s = r_alpha_s,
    sigma_n = r_sigma_n,
    sigma_c = r_sigma_c
  )

  return(retval)
}


# -----------------------------------------------------------------------------
# Generate replicate-level ODE simulations with optional transcriptional shutoff.
#
# Main features:
#   - simulates multiple biological replicates,
#   - perturbs parameters across replicates via log-normal noise,
#   - keeps R fixed across replicates by design,
#   - optionally introduces a random temporal shift,
#   - optionally applies a transcriptional shutoff at time t_star,
#   - returns full trajectories, subsampled trajectories, steady states,
#     and replicate-specific parameters.
#
# Arguments:
#   - base_params: named parameter list including R and kinetic rates
#   - y0: initial state vector
#   - times: intended observation times
#   - n_replicates: number of biological replicates
#   - model_kinetics: ODE right-hand side function
#   - param_cv: coefficient of variation for replicate-to-replicate variability
#   - stimes: sampling times extracted from the full simulation
#   - shutofftimes: sampling times extracted from shutoff trajectories
#   - max_shift: maximum temporal shift applied to the simulation
#   - t_star: intervention time; after this point, R is set to zero
#
# Output:
#   A list containing:
#     - data: full trajectories without shutoff
#     - shutoff_data: trajectories with shutoff, if enabled
#     - tsampled_data: data sampled at stimes
#     - shutoff_tsampled_data: shutoff data sampled at shutofftimes
#     - ss_data: closed-form steady-state values per replicate
#     - time_shift: global random shift applied to simulated times
#     - parameters: replicate-specific perturbed parameter values
# -----------------------------------------------------------------------------
generate_ODE_states <- function(
  base_params, y0, times,
  n_replicates = 3,
  model_kinetics = rna_kinetics,
  param_cv = 0.05,
  stimes = c(10, 40, 50, 100),
  shutofftimes = c(10, 40, 50, 100),
  max_shift = NULL,
  t_star = NULL
){

  # Containers for outputs across replicates.
  data <- c()
  shutoff_data <- c()
  parameters <- c()
  steady_state <- c()

  # Convert coefficient of variation into log-normal sd.
  sdlog <- if (param_cv > 0) sqrt(log(1 + param_cv^2)) else 0

  # ---------------------------------------------------------------------------
  # Prepare simulation time grid.
  #
  # Times are sorted and deduplicated. Extra time points are appended to allow
  # a random temporal shift while still retaining values at the intended
  # observed times after shifting.
  # ---------------------------------------------------------------------------
  all_times <- sort(unique(times))

  # Default shift window: 10% of the maximum simulated time.
  if (is.null(max_shift)) {
    max_shift = floor(0.1 * max(all_times))
  }

  actual_times = all_times

  # Extend the simulation beyond the requested times to accommodate shifting.
  all_times = c(all_times, (max(all_times) + 1):(max(all_times) + max_shift))

  # Draw one global integer shift applied to all replicates in this simulation.
  rnd_shift = sample(-max_shift:max_shift, 1)

  # Determine whether the shutoff intervention is active and falls inside
  # the simulated time window.
  do_perturb <- (!is.null(t_star) &&
                   t_star > min(all_times) && t_star < max(all_times))

  # Ensure t_star is present in the integration grid when needed.
  if (do_perturb) {
    all_times <- sort(unique(c(all_times, t_star)))
  }

  # ---------------------------------------------------------------------------
  # Simulate each biological replicate.
  # ---------------------------------------------------------------------------
  for (i in seq_len(n_replicates)) {

    # Apply replicate-specific biological variability to parameters
    # via multiplicative log-normal perturbation.
    perturbed_params <- lapply(base_params, function(param) {
      if (sdlog > 0) param * exp(stats::rnorm(1, 0, sdlog)) else param
    })

    # Keep transcription rate R fixed across replicates by design.
    perturbed_params$R <- base_params$R

    # Store perturbed parameters.
    parameters <- rbind(parameters, as.data.frame(c(perturbed_params, replicate = i)))

    # Simulate the baseline trajectory without shutoff.
    out1 <- deSolve::ode(y = y0, times = all_times, func = model_kinetics, parms = perturbed_params)
    out_df <- as.data.frame(out1)

    # Apply the random time shift.
    out_df$time = out_df$time + rnd_shift

    # Keep only rows corresponding to the original intended times.
    out_df = out_df[out_df$time %in% actual_times, ]

    # If the shift moves the trajectory to the right, prepend zeros so that
    # the early time points remain represented.
    if (rnd_shift > 0) {
      df_zero = data.frame(
        time = actual_times[1:rnd_shift],
        N = rep(0, rnd_shift),
        N_s = rep(0, rnd_shift),
        C = rep(0, rnd_shift),
        C_s = rep(0, rnd_shift)
      )
      out_df = rbind(df_zero, out_df)
    }

    # Initialize shutoff trajectory as the unperturbed one.
    out_df_shutoff = out_df

    # -------------------------------------------------------------------------
    # If t_star is active, build a piecewise trajectory:
    #   - before t_star: original parameters,
    #   - from t_star onward: same parameters but R = 0.
    # -------------------------------------------------------------------------
    if (do_perturb) {
      times_post <- all_times[all_times >= t_star]

      # State at the intervention time.
      y_star <- as.numeric(out_df[out_df$time == t_star, c("N","N_s","C","C_s")])
      names(y_star) <- c("N","N_s","C","C_s")

      # Copy parameters and impose transcriptional shutoff.
      perturbed_params_post <- perturbed_params
      perturbed_params_post$R <- 0

      # Simulate post-intervention segment starting from y_star.
      out2 <- deSolve::ode(y = y_star, times = times_post, func = model_kinetics, parms = perturbed_params_post)

      # Keep pre-intervention rows from the original simulation.
      out_df_shutoff = out_df[out_df$time < t_star, ]

      # Append post-intervention trajectory.
      # This may include the t_star row again from the second segment.
      out_df_shutoff <- rbind(out_df_shutoff, out2)

      # Annotate the active value of R over time.
      out_df_shutoff$R <- ifelse(out_df_shutoff$time >= t_star, 0, perturbed_params$R)
    }

    # Add replicate identifiers.
    out_df$replicate <- i
    out_df$R = perturbed_params$R
    out_df_shutoff$replicate <- i

    # Compute and store closed-form steady states for this replicate.
    ssi <- steady_states(perturbed_params)
    steady_state <- rbind(
      steady_state,
      data.frame(replicate = i, R = perturbed_params$R, ssi)
    )

    # Append replicate data to global outputs.
    data <- rbind(data, out_df)
    shutoff_data <- rbind(shutoff_data, out_df_shutoff)
  }

  # ---------------------------------------------------------------------------
  # Extract sampled observations at selected time points.
  # ---------------------------------------------------------------------------
  df_tsampled = subset(data, time %in% stimes)
  df_shutoff_tsampled = subset(shutoff_data, time %in% shutofftimes)

  # Return all outputs in a structured list.
  list(
    data = data,
    shutoff_data = shutoff_data,
    tsampled_data = df_tsampled,
    shutoff_tsampled_data = df_shutoff_tsampled,
    ss_data = steady_state,
    time_shift = rnd_shift,
    parameters = parameters
  )
}