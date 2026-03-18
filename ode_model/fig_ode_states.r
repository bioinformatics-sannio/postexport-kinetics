setwd("~/postexport-kinetics/ode_model")
source("ode.r")
source("../commons/plot.r")

tmax = 2/min(r_sigma_n_min+r_tau_min, r_tau_s_min, r_sigma_c_min+r_alpha_min,r_alpha_s_min)
times = seq(0, tmax, by = 1) # simulation time points
T_star = times[floor(length(times)/3)]
ALPHA = 5
BETA = 20
tstep_i = 20
n_replicates_i = 5
time_samples_i = 5
u <- seq(0, 1, length.out = time_samples_i+2)
p <- 3
idx <- floor(max(times)/3 * (u ^ p))
sampled_times = sort(unique(idx))[-1]
sampled_times = sampled_times[-length(sampled_times)]
step_adders = cumsum(rep(tstep_i,time_samples_i-2))
shutoff_sampled_times = sort(c(T_star-tstep_i,T_star, T_star+step_adders))



r_param = random_params()
r_param$sigma_n = 0.04 #runif(1, r_sigma_c_min,r_sigma_c_max)
r_param$sigma_c = 0.03 #0.01 #runif(1, r_sigma_c_min,r_sigma_c_max)
r_param$alpha_s = 0.02 #runif(1,r_alpha_s_min,r_param$alpha)
r_param$alpha = 0.04 #runif(1,r_alpha_s_min,r_param$alpha)
r_param$tau = 0.02
r_param$tau_s = 0.07
r_param$R = rgamma(1, shape = ALPHA, scale = BETA)

y0 = c(N = 0, N_s = 0, C = 0, C_s = 0) 
dd = generate_ODE_states(base_params = r_param, y0=y0, times=times,
                          t_star = T_star,
                          n_replicates = n_replicates_i,
                          stimes = sampled_times,
                          shutofftimes = shutoff_sampled_times
                          )


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


ggplot(d_sum, aes(x = time, y = mean)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 5) +
  facet_wrap(~State, scales = "free_y") +
  labs(
    x = "Time (min)",
    y = "Mean ± SD"
  ) +
  theme_bw()





s1 = plot_simulation(data=dd$data,ssdata=dd$ss_data, r_param=r_param, replicate=1,variable="all",stimes = sampled_times)
sh1 = plot_simulation(data=dd$shutoff_data,ssdata=dd$ss_data, r_param=r_param, replicate=3,variable="all",stimes = shutoff_sampled_times)
p1 = plot_simulation_parameters(dd$parameters)

r_param = random_params()
r_param$sigma_c = 0
r_param$alpha_s = runif(1,r_alpha_s_min,r_param$alpha)
r_param$R = rgamma(1, shape = ALPHA, scale = BETA)

dd = generate_ODE_states(base_params = r_param, y0=y0, times=times,
                          t_star = T_star,
                          n_replicates = n_replicates_i,
                          stimes = sampled_times,
                          shutofftimes = shutoff_sampled_times
                          )

s2 = plot_simulation(data=dd$data,ssdata=dd$ss_data, r_param=r_param, replicate=1,variable="all",stimes = sampled_times)
sh2 = plot_simulation(data=dd$shutoff_data,ssdata=dd$ss_data, r_param=r_param, replicate=3,variable="all",stimes = shutoff_sampled_times)
p2 = plot_simulation_parameters(dd$parameters)

r_param$sigma_n = 0.03 #runif(1, r_sigma_c_min,r_sigma_c_max)
r_param$sigma_c = 0 #0.01 #runif(1, r_sigma_c_min,r_sigma_c_max)
r_param$alpha_s = 0.04 #runif(1,r_alpha_s_min,r_param$alpha)
r_param$alpha = 0.02 #runif(1,r_alpha_s_min,r_param$alpha)
r_param$tau = 0.03
r_param$tau_s = 0.1
r_param$R = rgamma(1, shape = ALPHA, scale = BETA)

dd = generate_ODE_states(base_params = r_param, y0=y0, times=times,
                          t_star = T_star,
                          n_replicates = n_replicates_i,
                          stimes = sampled_times,
                          shutofftimes = shutoff_sampled_times
                          )

s3 = plot_simulation(data=dd$data,ssdata=dd$ss_data, r_param=r_param, replicate=1,variable="all",stimes = sampled_times)
sh3 = plot_simulation(data=dd$shutoff_data,ssdata=dd$ss_data, r_param=r_param, replicate=3,variable="all",stimes = shutoff_sampled_times)
p3 = plot_simulation_parameters(dd$parameters)

library(patchwork)

s1 = s1 + ggtitle("No perturbation")
sh1 = sh1 + ggtitle("Shutoff")
p1 = p1 + ggtitle("Kinetic parameters")
final_plot <- (s1 | sh1 | p1) /
              #(s2 | sh2 | p2) /
              (s3 | sh3 | p3) 

final_plot
ggsave("ode-states.pdf", final_plot, width = 16, height = 12)
