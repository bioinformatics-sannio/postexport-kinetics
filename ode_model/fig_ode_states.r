# =============================================================================
# Title: Example Script for Simulating and Plotting ODE Trajectories
# Description:
#   This script demonstrates how to:
#     - initialize the simulation environment,
#     - define a time horizon and sampling scheme,
#     - generate ODE trajectories under different parameter settings,
#     - visualize baseline and shutoff dynamics,
#     - inspect replicate-level kinetic parameter variability,
#     - assemble publication-style figures.
#
# Workflow:
#   1. load model and plotting utilities,
#   2. define simulation times and sampling times,
#   3. simulate one or more scenarios,
#   4. summarize sampled shutoff trajectories,
#   5. build multi-panel comparison plots,
#   6. export the final figure to PDF.
#
# Intended use:
#   Reproducible example / figure-generation script accompanying the paper.
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


# -----------------------------------------------------------------------------
# Set working directory and load required project scripts.
#
# ode.r:
#   defines the ODE system and simulation utilities.
#
# ../commons/plot.r:
#   provides helper functions for trajectory and parameter visualization.
# -----------------------------------------------------------------------------
setwd("~/postexport-kinetics/ode_model")
source("ode.r")
source("../commons/plot.r")


# -----------------------------------------------------------------------------
# Define the simulation time horizon.
#
# The upper bound tmax is chosen from the slowest kinetic timescale so that
# the simulation covers a sufficiently long interval relative to the model
# dynamics.
# -----------------------------------------------------------------------------
tmax = 2 / min(
  r_sigma_n_min + r_tau_min,
  r_tau_s_min,
  r_sigma_c_min + r_alpha_min,
  r_alpha_s_min
)

# Integer-valued simulation grid.
times = seq(0, tmax, by = 1)

# Define the shutoff/intervention time as approximately one third
# of the full simulated time grid.
T_star = times[floor(length(times) / 3)]


# -----------------------------------------------------------------------------
# Global simulation settings.
#
# ALPHA and BETA parameterize the Gamma distribution used to sample R.
# -----------------------------------------------------------------------------
ALPHA = 5
BETA = 20

# Sampling design controls.
tstep_i = 20
n_replicates_i = 5
time_samples_i = 5

# -----------------------------------------------------------------------------
# Build a non-uniform baseline sampling design.
#
# The vector u spans [0, 1], then is transformed by a power p > 1
# to obtain denser early sampling. The resulting indices are mapped
# into the time grid.
# -----------------------------------------------------------------------------
u <- seq(0, 1, length.out = time_samples_i + 2)
p <- 3
idx <- floor(max(times) / 3 * (u ^ p))

# Retain unique internal sample times, excluding boundaries.
sampled_times = sort(unique(idx))[-1]
sampled_times = sampled_times[-length(sampled_times)]

# -----------------------------------------------------------------------------
# Build shutoff sampling times centered around T_star.
#
# These include:
#   - one point before shutoff,
#   - the shutoff time itself,
#   - multiple points after shutoff spaced by tstep_i.
# -----------------------------------------------------------------------------
step_adders = cumsum(rep(tstep_i, time_samples_i - 2))
shutoff_sampled_times = sort(c(T_star - tstep_i, T_star, T_star + step_adders))


# -----------------------------------------------------------------------------
# Scenario 1: example parameter configuration with sigma_c > 0.
#
# Start from a random parameter set, then overwrite selected entries to obtain
# a controlled illustrative scenario.
# -----------------------------------------------------------------------------
r_param = random_params()
r_param$sigma_n = 0.04
r_param$sigma_c = 0.03
r_param$alpha_s = 0.02
r_param$alpha = 0.04
r_param$tau = 0.02
r_param$tau_s = 0.07
r_param$R = rgamma(1, shape = ALPHA, scale = BETA)

# Start from zero abundance in all four compartments.
y0 = c(N = 0, N_s = 0, C = 0, C_s = 0)

# Generate replicate trajectories, including a shutoff condition at T_star.
dd = generate_ODE_states(
  base_params = r_param,
  y0 = y0,
  times = times,
  t_star = T_star,
  n_replicates = n_replicates_i,
  stimes = sampled_times,
  shutofftimes = shutoff_sampled_times
)


# -----------------------------------------------------------------------------
# Summarize the sampled shutoff data across replicates.
#
# This section reshapes the sampled shutoff data, computes mean and SD
# at each time point for each state, and plots a quick diagnostic summary.
# -----------------------------------------------------------------------------
require(data.table)

d_long <- melt(
  as.data.table(dd$shutoff_tsampled_data),
  id.vars = c("time", "replicate"),
  measure.vars = c("N", "N_s", "C", "C_s"),
  variable.name = "State",
  value.name = "Value"
)

d_sum <- d_long[, .(
  mean = mean(Value),
  sd   = sd(Value)
), by = .(time, State)]

# Quick faceted summary plot of sampled shutoff trajectories.
ggplot(d_sum, aes(x = time, y = mean)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 5
  ) +
  facet_wrap(~State, scales = "free_y") +
  labs(
    x = "Time (min)",
    y = "Mean ± SD"
  ) +
  theme_bw()


# -----------------------------------------------------------------------------
# Build plots for Scenario 1.
#
# s1:
#   baseline trajectory for one replicate.
#
# sh1:
#   shutoff trajectory for one replicate.
#
# p1:
#   replicate-level variability of perturbed kinetic parameters.
# -----------------------------------------------------------------------------
s1 = plot_simulation(
  data = dd$data,
  ssdata = dd$ss_data,
  r_param = r_param,
  replicate = 1,
  variable = "all",
  stimes = sampled_times
)

sh1 = plot_simulation(
  data = dd$shutoff_data,
  ssdata = dd$ss_data,
  r_param = r_param,
  replicate = 3,
  variable = "all",
  stimes = shutoff_sampled_times
)

p1 = plot_simulation_parameters(dd$parameters)


# -----------------------------------------------------------------------------
# Scenario 2: null-like configuration with sigma_c = 0.
#
# This scenario is generated but not included in the final multi-panel figure
# below (the corresponding patchwork row is commented out).
# -----------------------------------------------------------------------------
r_param = random_params()
r_param$sigma_c = 0
r_param$alpha_s = runif(1, r_alpha_s_min, r_param$alpha)
r_param$R = rgamma(1, shape = ALPHA, scale = BETA)

dd = generate_ODE_states(
  base_params = r_param,
  y0 = y0,
  times = times,
  t_star = T_star,
  n_replicates = n_replicates_i,
  stimes = sampled_times,
  shutofftimes = shutoff_sampled_times
)

s2 = plot_simulation(
  data = dd$data,
  ssdata = dd$ss_data,
  r_param = r_param,
  replicate = 1,
  variable = "all",
  stimes = sampled_times
)

sh2 = plot_simulation(
  data = dd$shutoff_data,
  ssdata = dd$ss_data,
  r_param = r_param,
  replicate = 3,
  variable = "all",
  stimes = shutoff_sampled_times
)

p2 = plot_simulation_parameters(dd$parameters)


# -----------------------------------------------------------------------------
# Scenario 3: alternative parameter regime with sigma_c = 0.
#
# This scenario is included in the final figure together with Scenario 1,
# allowing visual comparison of two different kinetic regimes.
# -----------------------------------------------------------------------------
r_param$sigma_n = 0.03
r_param$sigma_c = 0
r_param$alpha_s = 0.04
r_param$alpha = 0.02
r_param$tau = 0.03
r_param$tau_s = 0.1
r_param$R = rgamma(1, shape = ALPHA, scale = BETA)

dd = generate_ODE_states(
  base_params = r_param,
  y0 = y0,
  times = times,
  t_star = T_star,
  n_replicates = n_replicates_i,
  stimes = sampled_times,
  shutofftimes = shutoff_sampled_times
)

s3 = plot_simulation(
  data = dd$data,
  ssdata = dd$ss_data,
  r_param = r_param,
  replicate = 1,
  variable = "all",
  stimes = sampled_times
)

sh3 = plot_simulation(
  data = dd$shutoff_data,
  ssdata = dd$ss_data,
  r_param = r_param,
  replicate = 3,
  variable = "all",
  stimes = shutoff_sampled_times
)

p3 = plot_simulation_parameters(dd$parameters)


# -----------------------------------------------------------------------------
# Assemble final multi-panel figure using patchwork.
#
# Top row:
#   Scenario 1 (baseline, shutoff, parameter boxplots)
#
# Bottom row:
#   Scenario 3 (baseline, shutoff, parameter boxplots)
#
# Scenario 2 is currently omitted from the final layout but can be restored
# by uncommenting the corresponding row below.
# -----------------------------------------------------------------------------
library(patchwork)

s1 = s1 + ggtitle("No perturbation")
sh1 = sh1 + ggtitle("Shutoff")
p1 = p1 + ggtitle("Kinetic parameters")

final_plot <- (s1 | sh1 | p1) /
              #(s2 | sh2 | p2) /
              (s3 | sh3 | p3)

# Display the final figure.
final_plot

# Export to PDF.
ggsave("ode-states.pdf", final_plot, width = 16, height = 12)