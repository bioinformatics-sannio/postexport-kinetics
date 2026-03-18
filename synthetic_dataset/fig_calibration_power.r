setwd("~/postexport-kinetics/synthetic_dataset")

require(data.table)
library(ggplot2)
library(scales)

load("test_sigma_2k_new.rdata")




# *************
# FIG1-A **** CALIBRATION PLOT del test sotto H0 (Type I error) ****
# *************

dt0 <- dt_test_results[
  Test_method %in% c("Nested-wild_whitened") & 
  Positive == 0 & 
  Lambda == 0.3 &
  Tsteps ==10 
  #N_tsamples ==10
]

dt0[, Test_method := factor(Test_method,
                              levels = c("Nested-wild-rademacher", "Nested-wild-mammen"),
                              labels = c("Nested test (r)", "Nested test (m)"))]

# opzionale: rimuovi NA/inf
dt0 <- dt0[is.finite(p.value)]

# controllo quick
dt0[, .N, by = .(Test_method)]
dt0[, .(pmin = min(p.value), pmax = max(p.value)), by = Test_method]

alpha_grid <- c(0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.001)

size_by_group <- function(DT, alphas = alpha_grid) {
  DT[, {
    pv <- p.value
    out <- data.table(alpha = alphas)
    out[, size := sapply(alpha, function(a) mean(pv <= a, na.rm = TRUE))]
    out
  }, by = .(
    Test_method,
    N_replicates, N_tsamples, Tsteps,
    Perturbation, Exprs_noise, Platform
  )]
}
dt_size <- size_by_group(dt0, alpha_grid)

# deviazione dalla linea ideale size=alpha
dt_size[, diff := size - alpha]

# peggiori scenari in assoluto a alpha=0.05
dt_size[alpha == 0.05][order(-abs(diff))][1:30]

calibration_curve <- function(DT, alphas = seq(0, 1, by = 0.01)) {
  DT[, {
    pv <- p.value
    data.table(alpha = alphas,
               size = sapply(alphas, function(a) mean(pv <= a)))
  }, by = .(
    Test_method,
    N_replicates, N_tsamples, Tsteps,
    Perturbation, Exprs_noise, Platform
  )]
}

dt_cal <- calibration_curve(dt0)


# filtro scenario principale (fisso Tsteps=20, RNA-seq)
dt_cal_sub <- dt_cal[
  #Exprs_noise %in% c("Medium") &
  N_replicates == 5 &
  Tsteps == 10  & Exprs_noise=="Medium" #& Platform %in% c("RNA-seq") 
]

# Se vuoi anche stratificare per N_replicates, decommenta e scegli livelli:
# dt_cal_sub <- dt_cal_sub[N_replicates %in% c(3,5,10)]

# ---- riassumi: mediana + bande across scenari ----
# (qui “scenari” = tutto ciò che resta variabile dentro al gruppo; tipicamente geni/replicate-setups)
dt_band <- dt_cal_sub[, .(
  N = .N,
  med = median(size, na.rm = TRUE),
  q25 = quantile(size, 0.25, na.rm = TRUE),
  q75 = quantile(size, 0.75, na.rm = TRUE),
  q05 = quantile(size, 0.05, na.rm = TRUE),
  q95 = quantile(size, 0.95, na.rm = TRUE)
), by = .(
  Test_method,
     N_tsamples, Tsteps,
    Perturbation, Exprs_noise, Platform, alpha)]

dt_band[, N_tsamples := factor(N_tsamples, levels = c(3,5,10,20))]
#dt_band[, N_replicates := factor(N_replicates, levels = c(3,5,10))]

require(latex2exp)
require(scales)
p <- ggplot(dt_band, aes(x = alpha, group = N_tsamples)) +
  # diagonale discreta
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",alpha = 0.3,linewidth = 0.8) +
  # banda “larga” 90%
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = N_tsamples),
              alpha = 0.08, colour = NA) +
  # banda “stretta” 50%
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = N_tsamples),
              alpha = 0.12, colour = NA) +
  # linea mediana
  geom_line(aes(y = med, color = N_tsamples), linewidth = 0.9) +
  scale_x_continuous(
  trans = scales::pseudo_log_trans(base = 10),
  limits = c(0, 1),
  breaks = c(0.001, 0.05, 0.1, 0.2, 0.4),
  labels = c("0.001", "0.05", "0.1", "0.2", "0.4")
) +
  scale_y_continuous(
  trans = scales::pseudo_log_trans(base = 10),
  limits = c(0, 1),
  breaks = c(0.001, 0.05, 0.1, 0.2, 0.4),
  labels = c("0.001", "0.05", "0.1", "0.2", "0.4")
) +
  facet_grid(Platform ~ Perturbation) +

  # assi: mostra bene la regione rilevante (0–0.1)
  #scale_x_continuous(breaks = c(0.01, 0.02, 0.05, 0.1,0.2)) +
  #scale_y_continuous(breaks = c(0.01, 0.02, 0.05, 0.1,0.2)) +
  #coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 0.2)) +
  labs(
    x = TeX("Nominal $\\; \\alpha$"), y = TeX("Empirical $\\; P(p \\leq \\alpha \\;|\\; H_0)$"),
    color = "# time points",
    fill  = "# time points"
  ) +
    theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.width = unit(1.6, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )
p
ggsave("calibration-curve.pdf", p, width = 8, height = 12)

# ********************
# TabS1 Calibration data for all noise and replicates regimes
# alpha = 0.05 only
# ********************

library(data.table)

alpha0 <- 0.05

dt_suppl <- dt_size[alpha == alpha0 & Tsteps==10,
  .(
    n = .N,
    Mean = round(mean(size, na.rm = TRUE), 3),
    SD = round(sd(size, na.rm = TRUE), 3),
    q05 = round(quantile(size, 0.05, na.rm = TRUE), 3),
    q95 = round(quantile(size, 0.95, na.rm = TRUE), 3),
    Worst = round(max(size, na.rm = TRUE), 3),
    Deviation = round(mean(size, na.rm = TRUE) - alpha0, 3),
    Inflation = round(mean(size, na.rm = TRUE) / alpha0, 2)
  ),
  by = .(Test_method, Exprs_noise, Platform, Perturbation)
][order(Test_method, Exprs_noise, Platform, Perturbation)]

dt_suppl

dt_suppl_min <- dt_suppl[, .(
  Platform,
  Exprs_noise,
  Perturbation,
  Mean,
  SD,
  Inflation
)]
kbl <- knitr::kable(
  dt_suppl_min,
  format = "latex",
  booktabs = TRUE,
  longtable = FALSE,
  escape = TRUE,
  digits = 3,
  caption = "Empirical type~I error rates at nominal level $\\alpha = 0.05$ across simulation settings. 
Mean size, standard deviation across design configurations, and inflation factor (Mean/$\\alpha$) are reported.",
  label = "tab:type1_supp",
  align = "llllrrr"
)

kbl



# ***********************
# FIG1-B **** POWER CURVE (P(p ≤ α | H1)) al variare di alpha e della vera sigma_c (effetto) ****
# ***********************

# --- Wilson CI helper (binomial proportion) ---
wilson_ci <- function(x, n, z = 1.96) {
  if (n == 0) return(c(lo = NA_real_, hi = NA_real_))
  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2*n)) / denom
  half <- (z * sqrt((p*(1-p) + z^2/(4*n)) / n)) / denom
  c(lo = max(0, center - half), hi = min(1, center + half))
}

# --- filter scenario ---
w_dt <- dt_test_results[
  #N_tsamples == 20 &
    #N_replicates == 10 &
    Tsteps == 10 &
    Lambda==0.7 &
    Test_method == "Nested-wild_whitened" &
    Perturbation == "SHUTOFF" &
    Positive == 1L &
    is.finite(p.value) &
    is.finite(sigma_true_median)
]

# --- quantile binning of sigma_true_median ---
K <- 10  # 8–12 typically good
brks <- unique(quantile(w_dt$sigma_true_median, probs = seq(0, 1, length.out = K + 1), na.rm = TRUE))
# safeguard in case of ties in quantiles
if (length(brks) < 4) stop("Not enough unique sigma_true_median quantiles to bin. Reduce K or check sigma_true_median.")

w_dt[, sigma_bin := cut(sigma_true_median, breaks = brks, include.lowest = TRUE)]

# representative x for each bin: median sigma in bin (per facet)
w_dt[, sigma_x := median(sigma_true_median), 
    by = .(N_tsamples, N_replicates, Platform, Exprs_noise, sigma_bin)]

# --- aggregate power per bin and alpha + Wilson CI ---
make_pow <- function(alpha0) {
  dt <- w_dt[, .(
    rej = sum(p.value <= alpha0),
    n = .N,
    sigma = unique(sigma_x)
  ), by = .(N_tsamples, N_replicates, Platform, Exprs_noise, sigma_bin)]
  dt[, power := rej / n]
  dt[, c("lo","hi") := {
    ci <- wilson_ci(rej, n)
    .(ci["lo"], ci["hi"])
  }, by = .(N_tsamples, N_replicates, Platform, Exprs_noise, sigma_bin)]
  dt[, alpha := factor(sprintf("%.2f", alpha0), levels = c("0.01", "0.05"))]
  dt[]
}

dt_pow <- rbindlist(list(make_pow(0.01), make_pow(0.05)), use.names = TRUE)

# order within each facet by sigma
setorder(dt_pow, Platform, Exprs_noise, alpha, sigma)

p_power <- ggplot(dt_pow[N_tsamples==10 & N_replicates==5], 
    aes(x = sigma, y = power, color=alpha, group = alpha)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.04, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.6) +
  facet_grid(Platform ~ Exprs_noise) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = TeX("True effect size ($\\sigma_c$)"), 
    y = TeX("Power $\\; P(p \\leq \\alpha \\;|\\; H_1)$"),
    color = expression(alpha)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_power
ggsave("power-curve.pdf", p_power, width = 16, height = 10)

library(patchwork)

p_combined <- p + p_power +
  plot_layout(ncol = 2, widths = c(1, 2)) &
  theme(
    legend.position = "top"
  )

p_combined
ggsave("calibration_power_combined.pdf", p_combined,
       width = 20, height = 10)
       
# ***********************
# FIG Supp power analysis fatta solo per RNA-seq Medium
# ***********************

dt_plot = dt_pow[alpha==0.05 & #Platform=="RNA-seq" & 
          Exprs_noise=="Medium" ]
require(latex2exp)
p_power <- ggplot(dt_plot, 
 aes(x = sigma, y = power, colour = Platform)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = 0.005,
    linewidth = 0.1
  ) +
  facet_grid(N_tsamples ~ N_replicates, labeller = labeller(
      N_tsamples = function(x) paste0("T=", x),
      N_replicates = function(x) paste0("R=", x)
    )) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    #x = "Platform",
    y = TeX("Power $\\; P(p \\leq \\alpha \\;|\\; H_1)$"),
    x = TeX("True effect size ($\\sigma_c$)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_power
ggsave("power_curve_suppl.pdf", p_power, width = 8, height = 10)





# =========================
# decision plot at different p-values
# =========================

# filtro scenario
w_dt <- dt_test_results[
  N_tsamples == 10 &
  N_replicates == 5 &
  Tsteps == 10 & 
  Lambda==0.5 &
  Test_method == "Nested-wild_whitened" 

  #Platform %in% c("RNA-seq") &
  #Exprs_noise %in% c("Medium") &
  #Test_method %in% c("Nested-wild-rademacher")   
  #Perturbation == "SHUTOFF" 
]

w_dt[, Test_method := factor(Test_method,
                              levels = c("Nested-wild-rademacher", "PSIBaseline"),
                              labels = c("Nested test", "PSI Baseline"))]

w_dt[, qvalue := p.value]
w_dt[, qvalue := p.adjust(p.value, "BH"),
                by = .(Platform, Exprs_noise, Perturbation,
                       Test_method, N_tsamples, Tsteps, N_replicates)]

#qs <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1)
qs <- c(0.01, 0.02, 0.05, 0.1, 0.2)

w_dt = w_dt[!is.na(qvalue)]
dt_qs_pred <- rbindlist(lapply(seq_along(qs), function(i) {

  qN <- qs[i]  
  w_dt[
    , .(
      alpha = qN,
      pred = qvalue <= qN),
    by = .(Platform, Exprs_noise, Perturbation,
                       Test_method, N_tsamples, Tsteps, N_replicates, Gene, Positive)
  ]
}))

pr_metrics <- function(pred, pos) {
  TP <- sum(pred==1 & pos==1L)
  FP <- sum(pred==1 & pos==0L)
  FN <- sum(pred==0 & pos==1L)

  P <- if ((TP + FP) == 0) NA_real_ else TP / (TP + FP)
  R <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)

  list(precision = P, recall = R)
}


B <- 1000L

pr_ci_dt <- dt_qs_pred[, {

  pr0 = pr_metrics(pred, Positive)

  n <- .N

  # bootstrap con replicate
  boot_mat <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    pr_b = pr_metrics(pred[idx], Positive[idx])
    c(
      recall = pr_b$recall,
      precision = pr_b$precision
    )
  })

  # boot_mat è 2 x B
  recall_ci <- quantile(boot_mat["recall", ], c(0.025, 0.975), na.rm = TRUE)
  precision_ci <- quantile(boot_mat["precision", ], c(0.025, 0.975), na.rm = TRUE)

  .(
    Recall = pr0$recall,
    Recall_low = recall_ci[1],
    Recall_high = recall_ci[2],
    Precision = pr0$precision,
    Precision_low = precision_ci[1],
    Precision_high = precision_ci[2]
  )

}, by = .(alpha, Platform, Exprs_noise, Perturbation,
                       Test_method, N_tsamples, Tsteps, N_replicates)]


library(grid)  
library(scales)
library(ggrepel)

pr_ci_dt$alpha_level <- factor(
  pr_ci_dt$alpha,
  levels = qs,
  ordered = TRUE
)
alpha_vals <- seq(1,0.4, -0.6/(length(qs)-1))
names(alpha_vals) <- as.character(qs)


p_decision <- ggplot(
  pr_ci_dt,
  aes(x = Recall, y = Precision,
      fill = Perturbation,
      alpha = alpha_level)
) +
  geom_errorbar(
    aes(xmin = Recall_low, xmax = Recall_high),
    width = 0, linewidth = 0.30, alpha = 0.25, colour = "grey25"
  ) +
  geom_errorbar(
    aes(ymin = Precision_low, ymax = Precision_high),
    width = 0, linewidth = 0.30, alpha = 0.25, colour = "grey25"
  ) +
  geom_point(shape = 21, size = 2.8, stroke = 0.4, colour = "grey25") +
  scale_alpha_manual(
    name = expression(alpha),
    values = alpha_vals
  ) +
  guides(
  fill = guide_legend(
    title = NULL,   # <-- niente titolo per methods
    order = 1
  ),
  alpha = guide_legend(
    order = 2,
    override.aes = list(
      fill = gray(seq(0.2, 0.85, length.out = length(alpha_vals))),
      colour = "grey25",
      shape = 21,
      size = 3
    )
  )
) +
  facet_grid(Platform ~ Exprs_noise) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

p_decision
ggsave("decision-plot.pdf", p_decision, width = 16, height = 10)

