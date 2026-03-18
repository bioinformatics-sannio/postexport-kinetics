library(data.table)
library(deSolve)
library(ggplot2)

plot_fit_two_models_paper <- function(ts_gene, res_nested,
                                      t_star = NULL, step = 1,
                                      show_shutoff = FALSE,
                                      state_labels = c(N="N", N_s="N_s", C="C", C_s="C_s"),
                                      model_labels = c(full = "Full ($\\sigma_c>0$)",
                                                       null = "Null ($\\sigma_c=0$)")) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  stopifnot(requireNamespace("latex2exp", quietly = TRUE))
  library(data.table); library(ggplot2)

  dt <- as.data.table(ts_gene)

  # ---- DATA: mean ± sd across replicates ----
  d_long <- melt(
    dt,
    id.vars = c("time","replicate"),
    measure.vars = c("N","N_s","C","C_s"),
    variable.name = "State",
    value.name = "Value"
  )

  d_sum <- d_long[, .(
    mean = mean(Value, na.rm = TRUE),
    sd   = sd(Value, na.rm = TRUE)
  ), by = .(time, State)]

  # ---- simulate fits ----
  times_dense <- seq(min(dt$time), max(dt$time), by = step)

  y0_full <- optimal_y0(ts_gene, res_nested$coef_full)
  pred_full <- simulate_fit(res_nested$coef_full, y0_full, times_dense, t_star = t_star)
  pred_full[, Model := "full"]

  y0_null <- optimal_y0(ts_gene, res_nested$coef_null)
  pred_null <- simulate_fit(res_nested$coef_null, y0_null, times_dense, t_star = t_star)
  pred_null[, Model := "null"]

  pred_all <- rbindlist(list(pred_full, pred_null), use.names = TRUE)

state_labels_txt <- c(
  N="N(t)", N_s="Ns(t)", C="C(t)", C_s="Cs(t)"
)

p_long <- melt(
    pred_all,
    id.vars = c("time","Model"),
    measure.vars = c("N","N_s","C","C_s"),
    variable.name = "State",
    value.name = "Value"
  )

  # nicer facet labels
  p_long[, State := factor(State, levels = c("N","N_s","C","C_s"))]
  d_sum[,  State := factor(State,  levels = c("N","N_s","C","C_s"))]
  dx <- min(diff(sort(unique(d_sum$time))))
  # ---- plot ----
  g <- ggplot() +
    # uncertainty (thin)
        # model fits
    geom_line(
      data = p_long,
      aes(x = time, y = Value, linetype = Model, colour = Model),
      linewidth = 1
    ) +
    geom_errorbar(
      data = d_sum,
      aes(x = time, ymin = pmax(mean - sd, 0), ymax = mean + sd),
      width = 0.2*dx,
      linewidth = 0.35,
      colour = "grey55"
    ) +
    geom_point(
      data = d_long,
      aes(x=time, y=Value),
      size = 0.5,
      colour = "grey40"
      #position = position_jitter(width = 1, height = 0)
    ) +
    # data points (squares with light border)
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
      #y = "Abundance (mean \u00B1 SD)",
      colour = NULL,
      linetype = NULL
    ) +
  theme_bw(base_size = 11) +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    axis.title.y = element_blank(),
    strip.text = element_text(size = 14),
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    panel.spacing = unit(0.6, "lines")
  )

  if (isTRUE(show_shutoff) && !is.null(t_star)) {
    g <- g + geom_vline(xintercept = t_star, linewidth = 0.35,
                        linetype = "dotted", colour = "grey40")
  }

  g
}

optimal_y0 = function(ts_gene, coeff) {
  t0 <- min(ts_gene$time)
  y0 <- ts_gene[time==t0, .(
    N   = mean(N,   na.rm=TRUE),
    N_s = mean(N_s, na.rm=TRUE),
    C   = mean(C,   na.rm=TRUE),
    C_s = mean(C_s, na.rm=TRUE)
  )]
  y0 <- as.numeric(y0[1,])
  names(y0) <- c("N","N_s","C","C_s")
  new_y0 = y0
  times_dense <- seq(min(ts_gene$time), max(ts_gene$time), by=1)

  
  for (X in names(y0)) {
    drg = diff(range(ts_gene[[X]], na.rm=TRUE))
    rr = seq(-drg,drg,by=drg/50)
    s_mse=c()
    for (j in 1:length(rr)) {
      y1 = y0
      y1[X] = y1[X] + rr[j]
      sim_d = simulate_fit(coeff,y1,times_dense,0)
      sim_d = sim_d[sim_d$time %in% ts_gene$time,]      
      s_mse = c(s_mse, sum((sim_d[[X]]-tapply(ts_gene[[X]], ts_gene$time, mean))^2))
    }
    new_y0[X] = y0[X] + rr[which.min(s_mse)]
  }
  return(new_y0)
}





ode_rhs <- function(t, y, p) {
  N  <- y["N"];  Ns <- y["N_s"]; C <- y["C"]; Cs <- y["C_s"]

  # shutoff: R=0 after t_star
  R0 <- p[["R"]]
  t_star <- p[["t_star"]]
  R <- if (!is.null(t_star) && !is.na(t_star) && t > t_star) 0 else R0

  dN  <- R - p[["sigma_n"]]*N - p[["tau"]]*N
  dNs <- p[["sigma_n"]]*N - p[["tau_s"]]*Ns
  dC  <- p[["tau"]]*N - p[["sigma_c"]]*C - p[["alpha"]]*C
  dCs <- p[["tau_s"]]*Ns + p[["sigma_c"]]*C - p[["alpha_s"]]*Cs
  list(c(dN, dNs, dC, dCs))
}

simulate_fit <- function(coef_vec, y0, times, t_star=0) {
  p <- as.list(coef_vec)
  p$t_star <- t_star
  out <- ode(y=y0, times=times, func=ode_rhs, parms=p, method="lsoda")
  out <- as.data.table(out)
  setnames(out, c("time","N","N_s","C","C_s"))
  out
}



plot_simulation_parameters = function(params) {
require(tidyr)
  df_long <- params |>
  pivot_longer(
    cols = -c(replicate,R),             # keep Replicate as ID
    names_to = "Par",        # new column for parameter names
    values_to = "Value"            # new column for parameter values
  )
  # Boxplot comparing perturbed and estimated parameters
  ggplot(df_long, aes(x = Par, y = Value)) +
    geom_boxplot(fill = "lightblue") +
    labs(title = "", x = NULL, y = NULL) +
    scale_x_discrete(
      labels = c(
        alpha            = TeX("$\\alpha$"),
        alpha_s     = TeX("$\\alpha_{s}$"),
        sigma_c     = TeX("$\\sigma_{c}$"),
        sigma_n     = TeX("$\\sigma_{n}$"),
        tau_s     = TeX("$\\tau_{s}$"),
        tau     = TeX("$\\tau$")
      )
    ) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),   # remove x and y labels
      plot.title = element_text(hjust = 0.5, size = 18)
    ) +
    theme(legend.position = "none")
}


plot_simulation = function(data,ssdata,r_param,replicate=1,variable="all",stimes=c(10,40,50,100)) {
  library(ggplot2)
  require(latex2exp)
  
  param_text1 = sprintf("$\\tau = %.3f, \\tau_s = %.3f, \\alpha = %.3f", 
  r_param$tau, r_param$tau_s, r_param$alpha)
  param_text2 = sprintf("$\\alpha_s = %.3f, \\sigma_n = %.3f, \\sigma_c = %.3f$", 
  r_param$alpha_s, r_param$sigma_n, r_param$sigma_c)
  if (variable=="all") {
    ss_row = ssdata[ssdata$replicate == replicate, ]
    sim_plot = ggplot(data[data$replicate == replicate, ], aes(x = time)) +
    geom_line(aes(y = N, color = "N"),linewidth=1.5) +
    geom_line(aes(y = N_s, color = "N_s"),linewidth=1.5) +
    geom_line(aes(y = C, color = "C"),linewidth=1.5) +
    geom_line(aes(y = C_s, color = "C_s"),linewidth=1.5) +
    geom_hline(yintercept = ss_row$N,   linetype = "dotted", color = "black") +
    geom_hline(yintercept = ss_row$N_s, linetype = "dotted", color = "black") +
    geom_hline(yintercept = ss_row$C,   linetype = "dotted", color = "black") +
    geom_hline(yintercept = ss_row$C_s, linetype = "dotted", color = "black") +
    annotate("text", x = max(data$time), y = ss_row$N,   label = "N",   hjust = -0.1, vjust = 0, size=4) +
    annotate("text", x = max(data$time), y = ss_row$N_s, label = "Ns", hjust = -0.1, vjust = 0, size=4) +
    annotate("text", x = max(data$time), y = ss_row$C,   label = "C",   hjust = -0.1, vjust = 0, size=4) +
    annotate("text", x = max(data$time), y = ss_row$C_s, label = "Cs", hjust = -0.1, vjust = 0, size=4) +
    geom_vline(xintercept = stimes, linetype = "dashed", color = "grey40")+
    labs(color = "Molecule") +
    annotate("text",
        x = max(data$time) * 0.60,    # adjust horizontal position
        y = max(data$N,data$N_s,data$C,data$C_s) * 0.7,       # adjust vertical position
        label = TeX(param_text1, output = "character"),
        hjust = 0,
        vjust = 1,
        parse=TRUE, size=3
    ) +
    annotate("text",
        x = max(data$time) * 0.60,    # adjust horizontal position
        y = max(data$N,data$N_s,data$C,data$C_s) * 0.65,       # adjust vertical position
        label = TeX(param_text2, output = "character"),
        hjust = 0,
        vjust = 1,
        parse=TRUE, size=3
    )
    
  } else {
    sim_plot = ggplot(data, aes(x = time, color = as.factor(replicate))) +    
    geom_line(aes_string(y = variable),linewidth=1) +
    geom_vline(xintercept = stimes, linetype = "dashed", color = "grey40")+
    labs(color = "Molecule")
    
  }
  sim_plot = sim_plot + 
    labs(title = "", x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),   # remove x and y labels
      axis.text  = element_blank(),   # remove tick labels
      axis.ticks = element_blank(),    # remove tick marks
      plot.title = element_text(hjust = 0.5, size = 18)
    ) +
    theme(legend.position = "none")
  return(sim_plot)
}

