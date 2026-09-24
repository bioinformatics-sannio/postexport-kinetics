# =============================================================================
# Numerical validation of trapezoidal / Crank-Nicolson interval balance
# against the exact linear ODE transition
# =============================================================================

setwd("~/postexport-kinetics/synthetic_dataset")

library(data.table)
library(ggplot2)

if (!requireNamespace("expm", quietly = TRUE)) {
  stop("Package 'expm' is required. Install it with install.packages('expm').")
}

if (!requireNamespace("nnls", quietly = TRUE)) {
  stop("Package 'nnls' is required. Install it with install.packages('nnls').")
}

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required. Install it with install.packages('patchwork').")
}

INPUT_FILE <- "benchmark_main_corrected_onset_raw.tsv"
OUTPUT_DIR <- "cn_exact_validation"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

RAW_FILE <- file.path(OUTPUT_DIR, "cn_exact_raw.tsv")
SUMMARY_FILE <- file.path(OUTPUT_DIR, "cn_exact_summary.tsv")
PARAM_SUMMARY_FILE <- file.path(OUTPUT_DIR, "cn_exact_parameter_summary.tsv")
FIG_PDF <- file.path(OUTPUT_DIR, "FigS_CN_vs_exact.pdf")
FIG_PNG <- file.path(OUTPUT_DIR, "FigS_CN_vs_exact.png")
SESSION_FILE <- file.path(OUTPUT_DIR, "sessionInfo.txt")
N_PARAM_SETS_REQUESTED <- 500L
N_TIME_POINTS <- 5L
DELTA_T_LEVELS <- c(1, 2, 5, 10, 20, 50)
SEED <- 20260923L
EPS <- 1e-12

raw <- fread(INPUT_FILE)

required <- c(
  "Gene", "Positive", "sigma_true", "sigma_n_true", "tau_true",
  "tau_s_true", "alpha_true", "alpha_s_true", "R_true"
)

missing <- setdiff(required, names(raw))

if (length(missing) > 0L) {
  stop(
    paste0(
      "Missing required columns in benchmark raw table:\n  ",
      paste(missing, collapse = "\n  ")
    )
  )
}

param_sets <- unique(
  raw[
    Positive == 1 &
      is.finite(sigma_true) &
      is.finite(sigma_n_true) &
      is.finite(tau_true) &
      is.finite(tau_s_true) &
      is.finite(alpha_true) &
      is.finite(alpha_s_true) &
      is.finite(R_true),
    .(
      Gene,
      sigma_c = sigma_true,
      sigma_n = sigma_n_true,
      tau = tau_true,
      tau_s = tau_s_true,
      alpha = alpha_true,
      alpha_s = alpha_s_true,
      R = R_true
    )
  ]
)

setorder(param_sets, Gene)
param_sets <- param_sets[!duplicated(Gene)]


N_PARAM_SETS <- min(
  N_PARAM_SETS_REQUESTED,
  nrow(param_sets)
)

if (
  N_PARAM_SETS <
    N_PARAM_SETS_REQUESTED
) {
  cat(
    "
Requested ",
    N_PARAM_SETS_REQUESTED,
    " parameter sets, but only ",
    nrow(param_sets),
    " unique alternative parameter sets are available.
",
    "Using all ",
    N_PARAM_SETS,
    " available sets.
",
    sep = ""
  )
}

set.seed(
  SEED
)

if (
  nrow(param_sets) >
    N_PARAM_SETS
) {
  param_sets <- param_sets[
    sample(
      .N,
      N_PARAM_SETS,
      replace = FALSE
    )
  ]
}
setorder(param_sets, Gene)

cat(
  "\n============================================================\n",
  "CRANK-NICOLSON VS EXACT LINEAR TRANSITION\n",
  "============================================================\n",
  "Parameter sets: ", nrow(param_sets), "\n",
  "Delta t values: ", paste(DELTA_T_LEVELS, collapse = ", "), "\n",
  "Post-shutoff time points per fit: ", N_TIME_POINTS, "\n",
  sep = ""
)

kinetic_matrix_local <- function(
  sigma_n,
  sigma_c,
  tau,
  tau_s,
  alpha,
  alpha_s
) {
  matrix(
    c(
      -(sigma_n + tau), 0, 0, 0,
      sigma_n, -tau_s, 0, 0,
      tau, 0, -(sigma_c + alpha), 0,
      0, tau_s, sigma_c, -alpha_s
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("N", "N_s", "C", "C_s"),
      c("N", "N_s", "C", "C_s")
    )
  )
}

exact_transition <- function(x, A, dt) {
  as.numeric(expm::expm(A * dt) %*% x)
}

cn_transition <- function(x, A, dt) {
  I4 <- diag(4)
  lhs <- I4 - 0.5 * dt * A
  rhs <- (I4 + 0.5 * dt * A) %*% x
  as.numeric(solve(lhs, rhs))
}

generate_exact_trajectory <- function(x0, A, dt, n_time_points) {
  times <- (0:(n_time_points - 1L)) * dt

  states <- t(
    vapply(
      times,
      function(tt) exact_transition(x0, A, tt),
      numeric(4)
    )
  )

  colnames(states) <- c("N", "N_s", "C", "C_s")

  data.table(
    time = times,
    N = states[, "N"],
    N_s = states[, "N_s"],
    C = states[, "C"],
    C_s = states[, "C_s"]
  )
}

fit_rates_from_trapezoidal_balance <- function(trajectory) {
  n_int <- nrow(trajectory) - 1L

  if (n_int < 1L) {
    stop("Trajectory must contain at least two time points.")
  }

  X_list <- vector("list", n_int)
  y_list <- vector("list", n_int)

  for (ii in seq_len(n_int)) {
    t0 <- trajectory$time[ii]
    t1 <- trajectory$time[ii + 1L]
    dt_i <- t1 - t0

    x0 <- c(
      N = trajectory$N[ii],
      N_s = trajectory$N_s[ii],
      C = trajectory$C[ii],
      C_s = trajectory$C_s[ii]
    )

    x1 <- c(
      N = trajectory$N[ii + 1L],
      N_s = trajectory$N_s[ii + 1L],
      C = trajectory$C[ii + 1L],
      C_s = trajectory$C_s[ii + 1L]
    )

    xm <- 0.5 * (x0 + x1)
    deriv <- (x1 - x0) / dt_i

    X_i <- matrix(0, nrow = 4, ncol = 6)
    colnames(X_i) <- c(
      "sigma_n", "sigma_c", "tau", "tau_s", "alpha", "alpha_s"
    )

    X_i[1, "sigma_n"] <- -xm["N"]
    X_i[1, "tau"] <- -xm["N"]

    X_i[2, "sigma_n"] <- xm["N"]
    X_i[2, "tau_s"] <- -xm["N_s"]

    X_i[3, "tau"] <- xm["N"]
    X_i[3, "sigma_c"] <- -xm["C"]
    X_i[3, "alpha"] <- -xm["C"]

    X_i[4, "tau_s"] <- xm["N_s"]
    X_i[4, "sigma_c"] <- xm["C"]
    X_i[4, "alpha_s"] <- -xm["C_s"]

    X_list[[ii]] <- X_i
    y_list[[ii]] <- as.numeric(deriv)
  }

  X <- do.call(rbind, X_list)
  y <- unlist(y_list, use.names = FALSE)

  state_id <- rep(1:4, times = n_int)

  scale_by_state <- vapply(
    1:4,
    function(ss) {
      max(abs(y[state_id == ss]), na.rm = TRUE)
    },
    numeric(1)
  )

  scale_by_state[
    !is.finite(scale_by_state) |
      scale_by_state < EPS
  ] <- 1

  w <- 1 / scale_by_state[state_id]
  Xw <- X * w
  yw <- y * w

  fit <- nnls::nnls(A = Xw, b = yw)

  coef <- as.numeric(fit$x)
  names(coef) <- colnames(X)

  fitted_y <- as.numeric(X %*% coef)

  residual_rel <- sqrt(sum((y - fitted_y)^2)) /
    max(sqrt(sum(y^2)), EPS)

  list(
    coef = coef,
    relative_balance_residual = residual_rel
  )
}

run_one_parameter_set <- function(row_i) {
  A <- kinetic_matrix_local(
    sigma_n = row_i$sigma_n,
    sigma_c = row_i$sigma_c,
    tau = row_i$tau,
    tau_s = row_i$tau_s,
    alpha = row_i$alpha,
    alpha_s = row_i$alpha_s
  )

  b_pre <- c(row_i$R, 0, 0, 0)

  x0 <- tryCatch(
    as.numeric(solve(A, -b_pre)),
    error = function(e) rep(NA_real_, 4)
  )

  if (
    any(!is.finite(x0)) ||
      any(x0 < 0)
  ) {
    return(NULL)
  }

  names(x0) <- c("N", "N_s", "C", "C_s")

  true_rates <- c(
    sigma_n = row_i$sigma_n,
    sigma_c = row_i$sigma_c,
    tau = row_i$tau,
    tau_s = row_i$tau_s,
    alpha = row_i$alpha,
    alpha_s = row_i$alpha_s
  )

  half_lives <- log(2) / true_rates

  t_half_min <- min(
    half_lives[
      is.finite(half_lives) &
        half_lives > 0
    ]
  )

  out <- list()

  for (jj in seq_along(DELTA_T_LEVELS)) {
    dt_i <- DELTA_T_LEVELS[jj]

    traj <- generate_exact_trajectory(
      x0 = x0,
      A = A,
      dt = dt_i,
      n_time_points = N_TIME_POINTS
    )

    transition_errors <- numeric(N_TIME_POINTS - 1L)

    for (kk in seq_len(N_TIME_POINTS - 1L)) {
      x_exact_start <- c(
        traj$N[kk],
        traj$N_s[kk],
        traj$C[kk],
        traj$C_s[kk]
      )

      x_exact_end <- c(
        traj$N[kk + 1L],
        traj$N_s[kk + 1L],
        traj$C[kk + 1L],
        traj$C_s[kk + 1L]
      )

      x_cn_end <- cn_transition(
        x = x_exact_start,
        A = A,
        dt = dt_i
      )

      transition_errors[kk] <- sqrt(
        sum((x_cn_end - x_exact_end)^2)
      ) / max(
        sqrt(sum(x_exact_end^2)),
        EPS
      )
    }

    fit <- fit_rates_from_trapezoidal_balance(traj)
    est <- fit$coef

    param_rel_error <- abs(est - true_rates) /
      pmax(abs(true_rates), EPS)

    out[[jj]] <- data.table(
      Gene = row_i$Gene,
      delta_t = dt_i,
      t_half_min = t_half_min,
      delta_t_over_t_half_min = dt_i / t_half_min,
      transition_error_median = median(transition_errors),
      transition_error_max = max(transition_errors),
      balance_residual_relative = fit$relative_balance_residual,

      sigma_n_true = true_rates["sigma_n"],
      sigma_n_hat = est["sigma_n"],
      sigma_n_rel_error = param_rel_error["sigma_n"],

      sigma_c_true = true_rates["sigma_c"],
      sigma_c_hat = est["sigma_c"],
      sigma_c_rel_error = param_rel_error["sigma_c"],

      tau_true = true_rates["tau"],
      tau_hat = est["tau"],
      tau_rel_error = param_rel_error["tau"],

      tau_s_true = true_rates["tau_s"],
      tau_s_hat = est["tau_s"],
      tau_s_rel_error = param_rel_error["tau_s"],

      alpha_true = true_rates["alpha"],
      alpha_hat = est["alpha"],
      alpha_rel_error = param_rel_error["alpha"],

      alpha_s_true = true_rates["alpha_s"],
      alpha_s_hat = est["alpha_s"],
      alpha_s_rel_error = param_rel_error["alpha_s"]
    )
  }

  rbindlist(out)
}

result_list <- vector("list", nrow(param_sets))

for (ii in seq_len(nrow(param_sets))) {
  if (ii %% 50L == 0L) {
    cat(
      "Completed ",
      ii,
      " / ",
      nrow(param_sets),
      " parameter sets\n",
      sep = ""
    )
  }

  result_list[[ii]] <- run_one_parameter_set(
    param_sets[ii]
  )
}

results <- rbindlist(
  result_list,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(results) == 0L) {
  stop("No valid numerical comparisons were generated.")
}

summary_dt <- results[
  ,
  .(
    N = .N,
    normalized_interval_median =
      median(delta_t_over_t_half_min),
    transition_error_median =
      median(transition_error_median),
    transition_error_q95 =
      quantile(
        transition_error_median,
        0.95,
        names = FALSE
      ),
    transition_error_max_median =
      median(transition_error_max),
    balance_residual_median =
      median(balance_residual_relative),
    sigma_c_rel_error_median =
      median(sigma_c_rel_error),
    sigma_c_rel_error_q95 =
      quantile(
        sigma_c_rel_error,
        0.95,
        names = FALSE
      )
  ),
  by = delta_t
]

setorder(summary_dt, delta_t)

param_long <- melt(
  results,
  id.vars = c(
    "Gene",
    "delta_t",
    "t_half_min",
    "delta_t_over_t_half_min"
  ),
  measure.vars = c(
    "sigma_n_rel_error",
    "sigma_c_rel_error",
    "tau_rel_error",
    "tau_s_rel_error",
    "alpha_rel_error",
    "alpha_s_rel_error"
  ),
  variable.name = "Parameter",
  value.name = "Relative_error"
)

param_long[
  ,
  Parameter :=
    sub(
      "_rel_error$",
      "",
      Parameter
    )
]

param_summary <- param_long[
  ,
  .(
    N = .N,
    Relative_error_median =
      median(Relative_error),
    Relative_error_q75 =
      quantile(
        Relative_error,
        0.75,
        names = FALSE
      ),
    Relative_error_q95 =
      quantile(
        Relative_error,
        0.95,
        names = FALSE
      )
  ),
  by = .(
    delta_t,
    Parameter
  )
]

setorder(
  param_summary,
  Parameter,
  delta_t
)

fwrite(
  results,
  RAW_FILE,
  sep = "\t"
)

fwrite(
  summary_dt,
  SUMMARY_FILE,
  sep = "\t"
)

fwrite(
  param_summary,
  PARAM_SUMMARY_FILE,
  sep = "\t"
)

plot_transition <- results[
  ,
  .(
    x =
      median(delta_t_over_t_half_min),
    Median =
      median(transition_error_median),
    Q25 =
      quantile(
        transition_error_median,
        0.25,
        names = FALSE
      ),
    Q75 =
      quantile(
        transition_error_median,
        0.75,
        names = FALSE
      )
  ),
  by = delta_t
]

p1 <- ggplot(
  plot_transition,
  aes(
    x = x,
    y = Median
  )
) +
  geom_ribbon(
    aes(
      ymin = Q25,
      ymax = Q75
    ),
    alpha = 0.15
  ) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x =
      expression(
        Delta * t / t[1/2 * "," * min]
      ),
    y =
      "Relative state-transition error"
  ) +
  theme_classic(
    base_size = 10.5
  ) +
  theme(
    panel.grid.major =
      element_line(
        colour = "grey93",
        linewidth = 0.25
      ),
    panel.grid.minor =
      element_blank()
  )

plot_param <- param_long[
  ,
  .(
    x =
      median(delta_t_over_t_half_min),
    Median =
      median(Relative_error)
  ),
  by = .(
    delta_t,
    Parameter
  )
]

plot_param[
  ,
  Parameter :=
    factor(
      Parameter,
      levels = c(
        "sigma_n",
        "sigma_c",
        "tau",
        "tau_s",
        "alpha",
        "alpha_s"
      )
    )
]

p2 <- ggplot(
  plot_param,
  aes(
    x = x,
    y = Median,
    linetype = Parameter,
    shape = Parameter,
    group = Parameter
  )
) +
  geom_line(
    linewidth = 0.75
  ) +
  geom_point(
    size = 1.8
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x =
      expression(
        Delta * t / t[1/2 * "," * min]
      ),
    y =
      "Median relative rate-estimation error",
    linetype =
      "Parameter",
    shape =
      "Parameter"
  ) +
  theme_classic(
    base_size = 10.5
  ) +
  theme(
    panel.grid.major =
      element_line(
        colour = "grey93",
        linewidth = 0.25
      ),
    panel.grid.minor =
      element_blank(),
    legend.position =
      "bottom"
  )

final_plot <- p1 /
  p2 +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

ggsave(
  filename = FIG_PDF,
  plot = final_plot,
  width = 8.5,
  height = 8.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = FIG_PNG,
  plot = final_plot,
  width = 8.5,
  height = 8.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

sink(SESSION_FILE)
print(sessionInfo())
sink()

cat(
  "\n============================================================\n",
  "CN VS EXACT VALIDATION COMPLETE\n",
  "============================================================\n",
  sep = ""
)

cat(
  "\nSummary by absolute sampling interval:\n"
)

print(summary_dt)

cat(
  "\nParameter-recovery error:\n"
)

print(param_summary)

cat(
  "\nGenerated in:\n  ",
  OUTPUT_DIR,
  "\n",
  sep = ""
)
