library(data.table)
library(parallel)

setwd("~/postexport-kinetics/synthetic_dataset")
source("../ode_model/ode.r")

range_n_replicates <- c(3,5,10)
range_time_samples <- c(3,5,10,20)
range_tstep        <- c(5,10,20,50)

ngenes   <- 2000
frac_pos <- 0.2

tmax  <- 3/min(r_sigma_n_min+r_tau_min, r_tau_s_min, r_sigma_c_min+r_alpha_min, r_alpha_s_min)
times <- seq(0, tmax, by = 1)
T_star <- times[floor(length(times)/3)]

grid <- CJ(
  N_replicates   = range_n_replicates,
  N_time_samples = range_time_samples,
  T_step         = range_tstep
)

simulate_one_gene <- function(gene_id, grid, frac_pos, times, T_star) {

  # riproducibile per gene (utile per benchmark/ROC)
  set.seed(1000000 + gene_id)

  r_param <- random_params()

  truth_pos <- as.integer(runif(1) <= frac_pos)
  if (truth_pos == 1L) {
    r_param$sigma_c <- runif(1, r_sigma_c_min, r_sigma_c_max)
  } else {
    r_param$sigma_c <- 0
  }

  ALPHA <- 5
  BETA  <- 20
  r_param$R <- rgamma(1, shape = ALPHA, scale = BETA)
  r_param$alpha_s <- runif(1, r_alpha_s_min, r_param$alpha)

  base_param_vec <- unlist(r_param)
  names(base_param_vec) <- paste0("Base_", names(base_param_vec))

  out_list <- vector("list", nrow(grid))

  for (gi in seq_len(nrow(grid))) {

    n_replicates_i <- grid[gi, N_replicates]
    time_samples_i <- grid[gi, N_time_samples]
    tstep_i        <- grid[gi, T_step]

    # sampled_times
    u <- seq(0, 1, length.out = time_samples_i + 2)
    p <- 3
    idx <- floor(max(times)/3 * (u ^ p))
    sampled_times <- sort(unique(idx))[-1]
    sampled_times <- sampled_times[-length(sampled_times)]

    # shutoff times
    step_adders <- cumsum(rep(tstep_i, time_samples_i - 2))
    shutoff_sampled_times <- sort(c(T_star - tstep_i, T_star, T_star + step_adders))

    y0 <- c(N=0, N_s=0, C=0, C_s=0)

    dd <- generate_ODE_states(
      base_params   = r_param,
      y0            = y0,
      times         = times,
      t_star        = T_star,
      n_replicates  = n_replicates_i,
      stimes        = sampled_times,
      shutofftimes  = shutoff_sampled_times
    )

    dt_none <- as.data.table(dd$tsampled_data)[, Perturbation := "NONE"]
    dt_shut <- as.data.table(dd$shutoff_tsampled_data)[, Perturbation := "SHUTOFF"]

    dt_ss <- as.data.table(dd$ss_data)
    dt_ss[, `:=`(Perturbation="SteadyState", time=0)]

    dt_one <- rbindlist(list(dt_none, dt_shut, dt_ss), use.names=TRUE, fill=TRUE)

    dt_one[, `:=`(
      Gene           = gene_id,
      truth_pos      = truth_pos,
      N_replicates   = n_replicates_i,
      N_time_samples = time_samples_i,
      T_step         = tstep_i,
      T_star         = T_star
    )]

    dt_one[, (names(base_param_vec)) := as.list(base_param_vec)]

    out_list[[gi]] <- dt_one
  }

  rbindlist(out_list, use.names=TRUE, fill=TRUE)
}

ncores <- max(1L, detectCores() - 80L)

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

ode_states <- rbindlist(res_list, use.names=TRUE, fill=TRUE)
save(ode_states, file = "ode_states_2k_20p.rdata")
