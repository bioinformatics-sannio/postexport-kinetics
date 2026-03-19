# =============================================================================
# Title: Synthetic Dataset Generation for ODE-Based RNA Kinetics
# Description:
#   This script generates a large synthetic benchmark dataset by simulating
#   replicate-level trajectories from the ODE model under multiple sampling
#   designs.
#
#   For each simulated gene:
#     - a random kinetic parameter set is sampled,
#     - the gene is assigned as positive or negative according to frac_pos,
#     - sigma_c is set accordingly,
#     - trajectories are generated under both baseline and shutoff settings,
#     - steady-state summaries are appended,
#     - metadata describing the design and ground truth are stored.
#
#   Simulations are repeated over a grid of:
#     - number of replicates,
#     - number of sampled time points,
#     - spacing of shutoff samples.
#
#   Parallel execution is used to scale generation across many genes.
#
# Intended use:
#   Synthetic benchmark generation for method evaluation, power analysis,
#   and ROC-style comparisons.
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
library(parallel)

# -----------------------------------------------------------------------------
# Set working directory and load the ODE simulation utilities.
# -----------------------------------------------------------------------------
setwd("~/postexport-kinetics/synthetic_dataset")
source("../ode_model/ode.r")


# -----------------------------------------------------------------------------
# Define the sampling-design grid explored in the benchmark.
#
# N_replicates:
#   number of biological replicates per gene.
#
# N_time_samples:
#   number of sampled time points used in the baseline design.
#
# T_step:
#   spacing between shutoff follow-up samples.
# -----------------------------------------------------------------------------
range_n_replicates <- c(3, 5, 10)
range_time_samples <- c(3, 5, 10, 20)
range_tstep        <- c(5, 10, 20, 50)


# -----------------------------------------------------------------------------
# Global benchmark settings.
#
# ngenes:
#   total number of synthetic genes to simulate.
#
# frac_pos:
#   fraction of genes assigned to the positive class, i.e. genes with
#   sigma_c > 0.
# -----------------------------------------------------------------------------
ngenes   <- 2000
frac_pos <- 0.2


# -----------------------------------------------------------------------------
# Define the global simulation time horizon.
#
# The upper bound is chosen from the slowest characteristic timescale so that
# the simulation covers a broad dynamic range.
# -----------------------------------------------------------------------------
tmax  <- 3 / min(
  r_sigma_n_min + r_tau_min,
  r_tau_s_min,
  r_sigma_c_min + r_alpha_min,
  r_alpha_s_min
)

times <- seq(0, tmax, by = 1)

# Intervention/shutoff time used throughout the benchmark.
T_star <- times[floor(length(times) / 3)]


# -----------------------------------------------------------------------------
# Build the full Cartesian product of design settings.
#
# Each row of grid corresponds to one sampling configuration that will be
# simulated for every synthetic gene.
# -----------------------------------------------------------------------------
grid <- CJ(
  N_replicates   = range_n_replicates,
  N_time_samples = range_time_samples,
  T_step         = range_tstep
)


# -----------------------------------------------------------------------------
# Simulate one synthetic gene across all sampling designs in the grid.
#
# Workflow for each gene:
#   1. set a reproducible seed,
#   2. sample a random parameter set,
#   3. assign the gene to the positive or negative class,
#   4. generate trajectories for every sampling design,
#   5. append baseline, shutoff, and steady-state outputs,
#   6. add truth labels and parameter metadata.
#
# Arguments:
#   - gene_id: integer identifier of the synthetic gene
#   - grid: data.table of sampling-design combinations
#   - frac_pos: probability that a gene belongs to the positive class
#   - times: dense simulation time grid
#   - T_star: shutoff/intervention time
#
# Output:
#   A data.table containing all simulated rows for that gene.
# -----------------------------------------------------------------------------
simulate_one_gene <- function(gene_id, grid, frac_pos, times, T_star) {

  # Reproducible gene-specific seed.
  # This is useful for benchmarking and ROC reproducibility.
  set.seed(1000000 + gene_id)

  # Sample a random kinetic parameter set.
  r_param <- random_params()

  # Assign gene class:
  #   truth_pos = 1 -> positive gene, sigma_c sampled from its range
  #   truth_pos = 0 -> negative gene, sigma_c fixed at zero
  truth_pos <- as.integer(runif(1) <= frac_pos)
  if (truth_pos == 1L) {
    r_param$sigma_c <- runif(1, r_sigma_c_min, r_sigma_c_max)
  } else {
    r_param$sigma_c <- 0
  }

  # Sample transcription rate R from a Gamma distribution.
  ALPHA <- 5
  BETA  <- 20
  r_param$R <- rgamma(1, shape = ALPHA, scale = BETA)

  # Constrain alpha_s relative to alpha.
  r_param$alpha_s <- runif(1, r_alpha_s_min, r_param$alpha)

  # Store the baseline parameter set as metadata.
  base_param_vec <- unlist(r_param)
  names(base_param_vec) <- paste0("Base_", names(base_param_vec))

  # One output block per row of the sampling-design grid.
  out_list <- vector("list", nrow(grid))

  # ---------------------------------------------------------------------------
  # Loop over sampling designs.
  # ---------------------------------------------------------------------------
  for (gi in seq_len(nrow(grid))) {

    n_replicates_i <- grid[gi, N_replicates]
    time_samples_i <- grid[gi, N_time_samples]
    tstep_i        <- grid[gi, T_step]

    # -------------------------------------------------------------------------
    # Build baseline sampled time points.
    #
    # A power transform is used to favor earlier sampling times.
    # -------------------------------------------------------------------------
    u <- seq(0, 1, length.out = time_samples_i + 2)
    p <- 3
    idx <- floor(max(times) / 3 * (u ^ p))

    sampled_times <- sort(unique(idx))[-1]
    sampled_times <- sampled_times[-length(sampled_times)]

    # -------------------------------------------------------------------------
    # Build shutoff sampling times.
    #
    # These include:
    #   - one time point before T_star,
    #   - T_star itself,
    #   - a sequence of post-shutoff time points spaced by tstep_i.
    # -------------------------------------------------------------------------
    step_adders <- cumsum(rep(tstep_i, time_samples_i - 2))
    shutoff_sampled_times <- sort(c(T_star - tstep_i, T_star, T_star + step_adders))

    # Zero initial condition for all compartments.
    y0 <- c(N = 0, N_s = 0, C = 0, C_s = 0)

    # Generate ODE trajectories under both baseline and shutoff conditions.
    dd <- generate_ODE_states(
      base_params   = r_param,
      y0            = y0,
      times         = times,
      t_star        = T_star,
      n_replicates  = n_replicates_i,
      stimes        = sampled_times,
      shutofftimes  = shutoff_sampled_times
    )

    # -------------------------------------------------------------------------
    # Convert outputs into tagged tables.
    #
    # Perturbation labels:
    #   NONE       -> baseline trajectory
    #   SHUTOFF    -> post-intervention trajectory
    #   SteadyState -> closed-form steady-state summary
    # -------------------------------------------------------------------------
    dt_none <- as.data.table(dd$tsampled_data)[, Perturbation := "NONE"]
    dt_shut <- as.data.table(dd$shutoff_tsampled_data)[, Perturbation := "SHUTOFF"]

    dt_ss <- as.data.table(dd$ss_data)
    dt_ss[, `:=`(Perturbation = "SteadyState", time = 0)]

    # Merge all outputs for the current design.
    dt_one <- rbindlist(list(dt_none, dt_shut, dt_ss), use.names = TRUE, fill = TRUE)

    # Add gene-level truth labels and design metadata.
    dt_one[, `:=`(
      Gene           = gene_id,
      truth_pos      = truth_pos,
      N_replicates   = n_replicates_i,
      N_time_samples = time_samples_i,
      T_step         = tstep_i,
      T_star         = T_star
    )]

    # Add baseline parameter values as columns.
    dt_one[, (names(base_param_vec)) := as.list(base_param_vec)]

    out_list[[gi]] <- dt_one
  }

  # Return the full per-gene simulation table.
  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}


# -----------------------------------------------------------------------------
# Determine the number of CPU cores to use.
#
# The expression below leaves a large number of cores free. This may be useful
# on shared systems, though on smaller machines it may reduce to a single core.
# -----------------------------------------------------------------------------
ncores <- max(1L, detectCores() - 80L)


# -----------------------------------------------------------------------------
# Run the benchmark in parallel across genes.
#
# Each gene is simulated independently, so the outer gene loop is naturally
# parallelizable.
# -----------------------------------------------------------------------------
res_list <- mclapply(
  X = seq_len(ngenes),
  FUN = simulate_one_gene,
  grid = grid,
  frac_pos = frac_pos,
  times = times,
  T_star = T_star,
  mc.cores = ncores,
  mc.set.seed = TRUE
)


# -----------------------------------------------------------------------------
# Combine all per-gene outputs into one large table.
# -----------------------------------------------------------------------------
ode_states <- rbindlist(res_list, use.names = TRUE, fill = TRUE)


# -----------------------------------------------------------------------------
# Save the benchmark dataset to disk.
#
# File name indicates:
#   - 2k genes
#   - 20% positives
# -----------------------------------------------------------------------------
save(ode_states, file = "ode_states_2k_20p.rdata")