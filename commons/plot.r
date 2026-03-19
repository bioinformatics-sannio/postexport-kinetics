# =============================================================================
# Title: ODE Simulation and Visualization Utilities for Model Comparison
# Description:
#   This script provides utilities to:
#     - simulate the compartmental ODE model,
#     - estimate initial conditions for plotting,
#     - compare full and null model fits against replicate data,
#     - visualize simulation parameters and trajectories.
#
# Included functionality:
#   - numerical simulation of the four-state ODE system,
#   - plotting of observed data together with full/null model predictions,
#   - heuristic adjustment of initial conditions,
#   - diagnostic plotting for simulation studies.
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


library(data.table)
library(deSolve)
library(ggplot2)


# -----------------------------------------------------------------------------
# Plot observed time-course summaries together with predictions from
# the full and null models.
#
# Main steps:
#   1. reshape replicate-level data into long format,
#   2. compute mean and standard deviation at each time point,
#   3. simulate dense trajectories for the full and null fitted models,
#   4. overlay model predictions, replicate points, and mean ± SD summaries.
#
# Arguments:
#   - ts_gene: observed replicate-level time-course data
#   - res_nested: output of the nested fit/testing procedure
#   - t_star: optional shutoff time used in simulation
#   - step: time spacing for dense simulation
#   - show_shutoff: whether to draw a vertical line at t_star
#   - state_labels: optional mapping of state names
#   - model_labels: labels for legend entries
#
# Output:
#   A ggplot object with one facet per state variable.
# -----------------------------------------------------------------------------
plot_fit_two_models_paper <- function(ts_gene, res_nested,
                                      t_star = NULL, step = 1,
                                      show_shutoff = FALSE,
                                      state_labels = c(N="N", N_s="N_s", C="C", C_s="C_s"),
                                      model_labels = c(full = "Full ($\\sigma_c>0$)",
                                                       null = "Null ($\\sigma_c=0$)")) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  stopifnot(requireNamespace("latex2exp", quietly = TRUE))
  library(data.table)
  library(ggplot2)

  dt <- as.data.table(ts_gene)

  # Reshape replicate data into long format for plotting raw observations.
  d_long <- melt(
    dt,
    id.vars = c("time","replicate"),
    measure.vars = c("N","N_s","C","C_s"),
    variable.name = "State",
    value.name = "Value"
  )

  # Compute time-wise mean and standard deviation across replicates.
  d_sum <- d_long[, .(
    mean = mean(Value, na.rm = TRUE),
    sd   = sd(Value, na.rm = TRUE)
  ), by = .(time, State)]

  # Dense time grid used to draw smooth model trajectories.
  times_dense <- seq(min(dt$time), max(dt$time), by = step)

  # Simulate the full model using optimized initial conditions.
  y0_full <- optimal_y0(ts_gene, res_nested$coef_full)
  pred_full <- simulate_fit(res_nested$coef_full, y0_full, times_dense, t_star = t_star)
  pred_full[, Model := "full"]

  # Simulate the null model using its own optimized initial conditions.
  y0_null <- optimal_y0(ts_gene, res_nested$coef_null)
  pred_null <- simulate_fit(res_nested$coef_null, y0_null, times_dense, t_star = t_star)
  pred_null[, Model := "null"]

  # Combine full and null predictions.
  pred_all <- rbindlist(list(pred_full, pred_null), use.names = TRUE)

  # Human-readable facet labels.
  state_labels_txt <- c(
    N="N(t)", N_s="Ns(t)", C="C(t)", C_s="Cs(t)"
  )

  # Long-format predicted trajectories.
  p_long <- melt(
    pred_all,
    id.vars = c("time","Model"),
    measure.vars = c("N","N_s","C","C_s"),
    variable.name = "State",
    value.name = "Value"
  )

  # Harmonize facet ordering.
  p_long[, State := factor(State, levels = c("N","N_s","C","C_s"))]
  d_sum[,  State := factor(State,  levels = c("N","N_s","C","C_s"))]

  # Use the minimum spacing between observed time points to size error bars.
  dx <- min(diff(sort(unique(d_sum$time))))

  # Build the comparison plot:
  #   - model fits as lines,
  #   - replicate observations as small grey points,
  #   - means as larger square points,
  #   - SD bars around the means.
  g <- ggplot() +
    geom_line(
      data = p_long,
      aes(x = time, y = Value, linetype = Model, colour = Model),
      linewidth = 1
    ) +
    geom_errorbar(
      data = d_sum,
      aes(x = time, ymin = pmax(mean - sd, 0), ymax = mean + sd),
      width = 0.2 * dx,
      linewidth = 0.35,
      colour = "grey55"
    ) +
    geom_point(
      data = d_long,
      aes(x = time, y = Value),
      size = 0.5,
      colour = "grey40"
    ) +
    geom_point(
      data = d_sum,
      aes(x = time, y = mean),
      shape = 22, size = 5,
      stroke = 0.35, colour = "grey25", fill = "white"
    ) +
    facet_wrap(
      ~State, scales="free_y", ncol=2,
      labeller = labeller(State = state_labels_txt)
    ) +
    scale_colour_manual(
      values = c(full = "#2a24e4", null = "#d20910"),
      breaks = c("full","null"),
      labels = c(model_labels["full"], model_labels["null"])
    ) +
    scale_linetype_manual(
      values = c(full = "solid", null = "22"),
      breaks = c("full","null"),
      labels = c(model_labels["full"], model_labels["null"])
    ) +
    labs(
      x = "Time (min)",
      colour = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      strip.text = element_text(size = 14),
      legend.key.width = unit(1.4, "cm"),
      strip.background = element_blank(),
      panel.spacing = unit(0.6, "lines")
    )

  # Optionally indicate the transcriptional shutoff time.
  if (isTRUE(show_shutoff) && !is.null(t_star)) {
    g <- g + geom_vline(
      xintercept = t_star,
      linewidth = 0.35,
      linetype = "dotted",
      colour = "grey40"
    )
  }

  g
}


# -----------------------------------------------------------------------------
# Heuristic adjustment of the initial condition y0 for plotting.
#
# Starting from the mean observed state at the first time point, this function
# adjusts each component separately by a small grid search. For each state:
#   - a range around the initial observed value is explored,
#   - trajectories are simulated,
#   - the value minimizing the squared error on that state is retained.
#
# This is a plotting-oriented heuristic rather than a formal optimization.
# -----------------------------------------------------------------------------
optimal_y0 = function(ts_gene, coeff) {
  # Initial time point.
  t0 <- min(ts_gene$time)

  # Initial state estimated as the mean across replicates at the first time.
  y0 <- ts_gene[time == t0, .(
    N   = mean(N,   na.rm=TRUE),
    N_s = mean(N_s, na.rm=TRUE),
    C   = mean(C,   na.rm=TRUE),
    C_s = mean(C_s, na.rm=TRUE)
  )]

  y0 <- as.numeric(y0[1,])
  names(y0) <- c("N","N_s","C","C_s")

  new_y0 = y0

  # Dense time grid for simulation during the search.
  times_dense <- seq(min(ts_gene$time), max(ts_gene$time), by=1)

  # Adjust each state independently.
  for (X in names(y0)) {
    # Search radius based on the observed range of the state.
    drg = diff(range(ts_gene[[X]], na.rm=TRUE))
    rr = seq(-drg, drg, by=drg/50)

    s_mse = c()

    for (j in 1:length(rr)) {
      y1 = y0
      y1[X] = y1[X] + rr[j]

      # Simulate candidate trajectory.
      sim_d = simulate_fit(coeff, y1, times_dense, 0)

      # Retain only simulated times matching observed times.
      sim_d = sim_d[sim_d$time %in% ts_gene$time,]

      # Compare simulated values to observed time-wise means for the same state.
      s_mse = c(
        s_mse,
        sum((sim_d[[X]] - tapply(ts_gene[[X]], ts_gene$time, mean))^2)
      )
    }

    # Keep the shift that minimizes the state-specific squared error.
    new_y0[X] = y0[X] + rr[which.min(s_mse)]
  }

  return(new_y0)
}


# -----------------------------------------------------------------------------
# Right-hand side of the four-state ODE system.
#
# State variables:
#   - N   : nuclear unprocessed
#   - N_s : nuclear processed
#   - C   : cytoplasmic unprocessed
#   - C_s : cytoplasmic processed
#
# Parameters:
#   - R, sigma_n, tau, tau_s, sigma_c, alpha, alpha_s
#
# Optional transcriptional shutoff:
#   If t_star is provided, the production term R is set to zero for t > t_star.
# -----------------------------------------------------------------------------
ode_rhs <- function(t, y, p) {
  N  <- y["N"]
  Ns <- y["N_s"]
  C  <- y["C"]
  Cs <- y["C_s"]

  # Shutoff mechanism: production is active before t_star and zero afterward.
  R0 <- p[["R"]]
  t_star <- p[["t_star"]]
  R <- if (!is.null(t_star) && !is.na(t_star) && t > t_star) 0 else R0

  # ODE system.
  dN  <- R - p[["sigma_n"]] * N - p[["tau"]] * N
  dNs <- p[["sigma_n"]] * N - p[["tau_s"]] * Ns
  dC  <- p[["tau"]] * N - p[["sigma_c"]] * C - p[["alpha"]] * C
  dCs <- p[["tau_s"]] * Ns + p[["sigma_c"]] * C - p[["alpha_s"]] * Cs

  list(c(dN, dNs, dC, dCs))
}


# -----------------------------------------------------------------------------
# Simulate the ODE model for a given parameter vector and initial state.
#
# Arguments:
#   - coef_vec: named vector of model parameters
#   - y0: named vector of initial conditions
#   - times: time grid for numerical integration
#   - t_star: optional shutoff time
#
# Output:
#   data.table with columns: time, N, N_s, C, C_s
# -----------------------------------------------------------------------------
simulate_fit <- function(coef_vec, y0, times, t_star=0) {
  # Convert coefficients to parameter list and append t_star.
  p <- as.list(coef_vec)
  p$t_star <- t_star

  # Solve the ODE system numerically.
  out <- ode(y = y0, times = times, func = ode_rhs, parms = p, method = "lsoda")

  # Convert output to data.table and set standard column names.
  out <- as.data.table(out)
  setnames(out, c("time","N","N_s","C","C_s"))
  out
}


# -----------------------------------------------------------------------------
# Boxplot of simulated parameter values across replicates.
#
# This is a compact diagnostic plot used to inspect the distribution of
# parameter values generated in simulation studies.
# -----------------------------------------------------------------------------
plot_simulation_parameters = function(params) {
  require(tidyr)

  # Reshape parameter table into long format.
  df_long <- params |>
    pivot_longer(
      cols = -c(replicate, R),
      names_to = "Par",
      values_to = "Value"
    )

  # Boxplot of parameters across replicates.
  ggplot(df_long, aes(x = Par, y = Value)) +
    geom_boxplot(fill = "lightblue") +
    labs(title = "", x = NULL, y = NULL) +
    scale_x_discrete(
      labels = c(
        alpha   = TeX("$\\alpha$"),
        alpha_s = TeX("$\\alpha_{s}$"),
        sigma_c = TeX("$\\sigma_{c}$"),
        sigma_n = TeX("$\\sigma_{n}$"),
        tau_s   = TeX("$\\tau_{s}$"),
        tau     = TeX("$\\tau$")
      )
    ) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18)
    ) +
    theme(legend.position = "none")
}


# -----------------------------------------------------------------------------
# Plot simulated trajectories for one replicate or for one selected variable.
#
# Two modes:
#   - variable = "all":
#       shows all four trajectories for a chosen replicate, with optional
#       reference steady-state lines and parameter annotations.
#   - variable != "all":
#       shows the selected variable across replicates.
#
# Arguments:
#   - data: simulated time-course data
#   - ssdata: steady-state reference values
#   - r_param: parameter values used in simulation
#   - replicate: replicate to highlight in "all" mode
#   - variable: "all" or the name of one variable to plot
#   - stimes: optional vertical reference times
# -----------------------------------------------------------------------------
plot_simulation = function(data, ssdata, r_param, replicate=1,
                           variable="all", stimes=c(10,40,50,100)) {
  library(ggplot2)
  require(latex2exp)
  
  # Format parameter annotations for the plot.
  param_text1 = sprintf(
    "$\\tau = %.3f, \\tau_s = %.3f, \\alpha = %.3f",
    r_param$tau, r_param$tau_s, r_param$alpha
  )
  param_text2 = sprintf(
    "$\\alpha_s = %.3f, \\sigma_n = %.3f, \\sigma_c = %.3f$",
    r_param$alpha_s, r_param$sigma_n, r_param$sigma_c
  )

  if (variable == "all") {
    # Select steady-state values for the chosen replicate.
    ss_row = ssdata[ssdata$replicate == replicate, ]

    # Plot all four state trajectories for the selected replicate.
    sim_plot = ggplot(data[data$replicate == replicate, ], aes(x = time)) +
      geom_line(aes(y = N,   color = "N"),   linewidth = 1.5) +
      geom_line(aes(y = N_s, color = "N_s"), linewidth = 1.5) +
      geom_line(aes(y = C,   color = "C"),   linewidth = 1.5) +
      geom_line(aes(y = C_s, color = "C_s"), linewidth = 1.5) +

      # Add steady-state reference levels.
      geom_hline(yintercept = ss_row$N,   linetype = "dotted", color = "black") +
      geom_hline(yintercept = ss_row$N_s, linetype = "dotted", color = "black") +
      geom_hline(yintercept = ss_row$C,   linetype = "dotted", color = "black") +
      geom_hline(yintercept = ss_row$C_s, linetype = "dotted", color = "black") +

      # Label the steady-state reference lines.
      annotate("text", x = max(data$time), y = ss_row$N,   label = "N",
               hjust = -0.1, vjust = 0, size = 4) +
      annotate("text", x = max(data$time), y = ss_row$N_s, label = "Ns",
               hjust = -0.1, vjust = 0, size = 4) +
      annotate("text", x = max(data$time), y = ss_row$C,   label = "C",
               hjust = -0.1, vjust = 0, size = 4) +
      annotate("text", x = max(data$time), y = ss_row$C_s, label = "Cs",
               hjust = -0.1, vjust = 0, size = 4) +

      # Add vertical markers at selected times.
      geom_vline(xintercept = stimes, linetype = "dashed", color = "grey40") +
      labs(color = "Molecule") +

      # Annotate the parameter values used in the simulation.
      annotate(
        "text",
        x = max(data$time) * 0.60,
        y = max(data$N, data$N_s, data$C, data$C_s) * 0.7,
        label = TeX(param_text1, output = "character"),
        hjust = 0,
        vjust = 1,
        parse = TRUE,
        size = 3
      ) +
      annotate(
        "text",
        x = max(data$time) * 0.60,
        y = max(data$N, data$N_s, data$C, data$C_s) * 0.65,
        label = TeX(param_text2, output = "character"),
        hjust = 0,
        vjust = 1,
        parse = TRUE,
        size = 3
      )
    
  } else {
    # Plot a single variable across replicates.
    sim_plot = ggplot(data, aes(x = time, color = as.factor(replicate))) +
      geom_line(aes_string(y = variable), linewidth = 1) +
      geom_vline(xintercept = stimes, linetype = "dashed", color = "grey40") +
      labs(color = "Molecule")
  }

  # Common minimalist styling.
  sim_plot = sim_plot +
    labs(title = "", x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18)
    ) +
    theme(legend.position = "none")

  return(sim_plot)
}