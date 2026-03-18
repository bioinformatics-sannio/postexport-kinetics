library(deSolve)

# Define the ODE system
rna_kinetics <- function(t, y, params) {
  with(as.list(c(y, params)), {
    dN <- R - sigma_n * N - tau * N
    dN_s <- sigma_n * N - tau_s * N_s
    dC <- tau * N - sigma_c * C - alpha * C
    dC_s <- tau_s * N_s + sigma_c * C - alpha_s * C_s
#    print(paste("Time:", t, "dN:", dN, "dN_s:", dN_s, "dC:", dC, "dC_s:", dC_s))
    list(c(dN, dN_s, dC, dC_s))
  })
}

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

# closed form of steady states
steady_states <- function(params) {
    N = params$R / (params$tau+params$sigma_n)
    N_s = params$R * params$sigma_n / ((params$tau+params$sigma_n)*params$tau_s)
    C = params$R * params$tau / ((params$tau+params$sigma_n)*(params$alpha+params$sigma_c))
    C_s =  (params$sigma_n + params$sigma_c*params$tau/(params$sigma_c+params$alpha)) * params$R / ((params$tau+params$sigma_n)*params$alpha_s)
    return(list(N = N, N_s = N_s, C = C, C_s = C_s))
}


random_params = function(rtau_min = r_tau_min, rtau_max = r_tau_max, rtau_s_min = r_tau_s_min, rtau_s_max = r_tau_s_max,
  ralpha_min = r_alpha_min, ralpha_max = r_alpha_max, ralpha_s_min = r_alpha_s_min, ralpha_s_max = r_alpha_s_max,
  rsigma_n_min = r_sigma_n_min, rsigma_n_max = r_sigma_n_max,
  rsigma_c_min = r_sigma_c_min, rsigma_c_max = r_sigma_c_max) {

    r_tau = runif(1, rtau_min, rtau_max)
    r_tau_s = runif(1, rtau_s_min, rtau_s_max)
    r_alpha = runif(1, ralpha_min,ralpha_max)
    r_alpha_s = runif(1, ralpha_s_min,ralpha_s_max)

    #r_rho = runif(1, 0.1,1)
    #r_alpha_s = r_rho * r_alpha
    r_sigma_n = runif(1, rsigma_n_min,rsigma_n_max)
    r_sigma_c = runif(1, rsigma_c_min,rsigma_c_max)

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



generate_ODE_states <- function(
  base_params, y0, times,
  n_replicates = 3,
  model_kinetics = rna_kinetics,
  param_cv = 0.05,             # biological CV across replicates (log-normal)
  stimes = c(10, 40, 50, 100),
  shutofftimes = c(10, 40, 50, 100),
  max_shift = NULL, # delta to use for random shift the simulation later
  t_star = NULL            # intervention time (in same units as 'times')
){

  data <- c(); shutoff_data=c(); parameters <- c(); steady_state <- c()
  sdlog <- if (param_cv > 0) sqrt(log(1 + param_cv^2)) else 0

  # ---------------- piecewise integration (pre/post t_star) ----------------
  # ensure times are sorted/unique and include t_star if applicable
  all_times <- sort(unique(times))

  # add a delta to the simulation to use for random shift the simulation later
  if (is.null(max_shift)) {
    max_shift = floor(0.1*max(all_times))
  }
  actual_times = all_times
  all_times = c(all_times, (max(all_times)+1):(max(all_times)+max_shift))
  rnd_shift = sample(-max_shift:max_shift,1)

  do_perturb <- (!is.null(t_star) &&
                   t_star > min(all_times) && t_star < max(all_times))
  if (do_perturb) {
    all_times <- sort(unique(c(all_times, t_star)))
  }

  for (i in seq_len(n_replicates)) {
    # biological replicate variability
    perturbed_params <- lapply(base_params, function(param) {
      if (sdlog > 0) param * exp(stats::rnorm(1, 0, sdlog)) else param
    })
    # keep transcription rate R constant across replicates by design
    perturbed_params$R <- base_params$R
    parameters <- rbind(parameters, as.data.frame(c(perturbed_params, replicate = i)))

    out1 <- deSolve::ode(y = y0, times = all_times, func = model_kinetics, parms = perturbed_params)
    out_df <- as.data.frame(out1)
    out_df$time = out_df$time + rnd_shift
    out_df = out_df[out_df$time %in% actual_times,]
    if (rnd_shift > 0) {
      df_zero = data.frame(time=actual_times[1:rnd_shift],
                           N=rep(0,rnd_shift),
                           N_s=rep(0,rnd_shift),
                           C=rep(0,rnd_shift), C_s=rep(0,rnd_shift))
      out_df = rbind(df_zero, out_df)
    }
    

    out_df_shutoff = out_df
    if (do_perturb) {
      times_post <- all_times[all_times >= t_star]
      # pre-perturbation
      y_star <- as.numeric(out_df[out_df$time==t_star, c("N","N_s","C","C_s")])           # state at t_star
      names(y_star) <- c("N","N_s","C","C_s")
      # post-perturbation parameters
      perturbed_params_post <- perturbed_params
      perturbed_params_post$R <- 0

      # post-perturbation (start at t_star)
      out2 <- deSolve::ode(y = y_star, times = times_post, func = model_kinetics, parms = perturbed_params_post)
      out_df_shutoff = out_df[out_df$time<t_star,]
      # merge, dropping duplicate t_star row from second segment
      out_df_shutoff <- rbind(out_df_shutoff, out2)
      out_df_shutoff$R <- ifelse(out_df_shutoff$time >= t_star, 0, perturbed_params$R)
    }

    out_df$replicate <- i
    out_df$R = perturbed_params$R
    out_df_shutoff$replicate <- i

    ssi <- steady_states(perturbed_params)
    steady_state <- rbind(steady_state, data.frame(replicate = i, R = perturbed_params$R, ssi))
    data <- rbind(data, out_df)
    shutoff_data <- rbind(shutoff_data, out_df_shutoff)
  }

  df_tsampled = subset(data, time %in% stimes) 
  df_shutoff_tsampled = subset(shutoff_data, time %in% shutofftimes)
  

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
