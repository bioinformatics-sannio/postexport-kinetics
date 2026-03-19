# =============================================================================
# Title: Weighted NNLS with Variance Shrinkage and Bootstrap Inference
# Description:
#   This script implements a weighted non-negative least squares (NNLS)
#   framework for time-aggregated data, combining:
#     - variance shrinkage for numerical stability,
#     - construction of a linear system from time summaries,
#     - whitening through inverse-variance weights,
#     - nested-model testing via bootstrap under non-negativity constraints.
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


# -----------------------------------------------------------------------------
# Compute time-specific means and stabilized variances.
#
# For each time point, this function:
#   1. computes the sample mean of each selected variable,
#   2. computes the sample variance within the time point,
#   3. shrinks that variance toward a pooled variance estimated from the
#      full dataset,
#   4. converts the shrunk variance into the variance of the sample mean.
#
# Shrinkage is useful when some time points have few observations, making
# empirical variances unstable or undefined.
# -----------------------------------------------------------------------------
time_summary_diag_shrink <- function(df,
                                    vars=c("N","N_s","C","C_s"),
                                    eps_var=1e-8,
                                    lambda_var=0.7) {
  # Enforce a valid shrinkage weight in [0, 1].
  stopifnot(lambda_var >= 0, lambda_var <= 1)

  # Compute pooled variances across the entire dataset for each variable.
  # These are used as shrinkage targets.
  pooled_var <- sapply(vars, function(v) stats::var(df[[v]], na.rm=TRUE))

  # Replace missing pooled variances with zero, then bound away from zero.
  pooled_var[is.na(pooled_var)] <- 0
  pooled_var <- pmax(pooled_var, eps_var)

  # Split data by time point.
  # This works both for data.frame and data.table inputs.
  spl <- split(df, df$time)

  # Extract and sort time points numerically.
  times <- sort(as.numeric(names(spl)))

  # Preallocate output list.
  out <- vector("list", length(times))
  names(out) <- times

  # Process each time point separately.
  for (ti in times) {
    DT <- spl[[as.character(ti)]]

    # Robust column selection, compatible with data.table syntax.
    sub <- as.data.frame(DT[, ..vars])

    # Number of observations at the current time point.
    n <- nrow(sub)

    # Sample means at the current time point.
    mu <- sapply(vars, function(v) mean(sub[[v]], na.rm=TRUE))

    # Sample variances at the current time point.
    v <- sapply(vars, function(vn) stats::var(sub[[vn]], na.rm=TRUE))

    # Replace undefined variances with zero, then bound away from zero.
    v[is.na(v)] <- 0
    v <- pmax(v, eps_var)

    # Shrink time-specific variances toward the pooled variances.
    # lambda_var = 0   -> no shrinkage
    # lambda_var = 1   -> full pooling
    v_shr <- (1 - lambda_var) * v + lambda_var * pooled_var

    # Convert variance into variance of the sample mean.
    # max(n, 1) prevents division by zero in degenerate cases.
    v_mean <- pmax(v_shr / max(n, 1), eps_var)

    # Store summaries for this time point.
    out[[as.character(ti)]] <- list(time=ti, n=n, mean=mu, var_mean=v_mean)
  }

  # Attach pooled variances as an attribute for later reuse.
  attr(out, "pooled_var") <- pooled_var
  out
}


# -----------------------------------------------------------------------------
# Build the weighted linear system A * theta ≈ b.
#
# The system is obtained from consecutive time-point summaries.
# For each interval:
#   - dX contains observed changes in the mean state variables,
#   - I contains midpoint-approximated integrals over the interval,
#   - VdX contains diagonal variance approximations for the increments.
#
# The function returns:
#   - A: design matrix,
#   - b: stacked increments,
#   - w_sqrt: whitening weights = 1 / sqrt(var_b),
#   - col_norms: optional scaling factors applied to A,
#   - times, var_b, pooled_var: diagnostic outputs.
# -----------------------------------------------------------------------------
build_Ab_w_diag_shrink <- function(tsampled_data,
                                   scaling_A=FALSE,
                                   t_star=NULL,
                                   eps_var=1e-8,
                                   lambda_var=0.7) {
  vars <- c("N","N_s","C","C_s")

  # Compute time-specific means and stabilized variances.
  S <- time_summary_diag_shrink(tsampled_data, vars,
                                eps_var=eps_var,
                                lambda_var=lambda_var)

  times <- sort(as.numeric(names(S)))
  K <- length(times)

  # At least two time points are required to define increments.
  if (K < 2) stop("Serve almeno 2 timepoint")

  # Matrix of time-specific means.
  Xbar <- do.call(rbind, lapply(S, `[[`, "mean"))
  colnames(Xbar) <- vars

  # Matrix of variances of the means.
  Vbar <- do.call(rbind, lapply(S, `[[`, "var_mean"))
  colnames(Vbar) <- vars

  # Time-step lengths.
  deltaT <- diff(times)

  # Increments between consecutive time points.
  dX <- Xbar[-1, , drop=FALSE] - Xbar[-K, , drop=FALSE]

  # Midpoint values used for interval integration.
  mX <- (Xbar[-1, , drop=FALSE] + Xbar[-K, , drop=FALSE]) / 2

  # Midpoint rule approximation of integral terms.
  I  <- mX * deltaT

  # Approximate diagonal variances of the increments:
  # Var(X_{t+1} - X_t) ≈ Var(X_{t+1}) + Var(X_t)
  VdX <- Vbar[-1, , drop=FALSE] + Vbar[-K, , drop=FALSE]

  # Contribution of R * Delta t in the first equation.
  # By default, the full interval length is used.
  R_dt <- deltaT

  # If t_star is provided, only the portion of each interval before t_star
  # contributes to R. The value is clipped to [0, Delta t].
  if (!is.null(t_star)) {
    t0 <- times[-K]
    t1 <- times[-1]

    # Cases:
    #   t_star <= t0  -> contribution = 0
    #   t_star >= t1  -> contribution = Delta t
    #   otherwise     -> contribution = t_star - t0
    R_dt <- pmax(0, pmin(deltaT, t_star - t0))
  }

  # Build the block design matrix A, one 4-row block per time interval.
  # Only the first row uses R_dt instead of Delta t.
  A <- do.call(rbind, lapply(1:(K-1), function(k){
    IN  <- I[k, "N"]
    INs <- I[k, "N_s"]
    IC  <- I[k, "C"]
    ICs <- I[k, "C_s"]

    rbind(
      c(R_dt[k], -IN,  0,      0,     -IN, 0,      0),
      c(0,       0,   -INs,    0,      IN, 0,      0),
      c(0,       IN,   0,     -IC,     0, -IC,     0),
      c(0,       0,    INs,    IC,     0,  0,     -ICs)
    )
  }))
  colnames(A) <- c("R","tau","tau_s","sigma_c","sigma_n","alpha","alpha_s")

  # Stack increments row-wise into the response vector b.
  b <- as.vector(t(dX))

  # Stack corresponding variances into var_b.
  var_b <- as.vector(t(VdX))
  var_b <- pmax(var_b, eps_var)

  # Whitening weights for weighted least squares.
  w_sqrt <- 1 / sqrt(var_b)

  # Optional column scaling to improve conditioning of A.
  col_norms <- rep(1, ncol(A))
  if (scaling_A) {
    col_norms <- sqrt(colSums(A^2))
    col_norms[col_norms == 0] <- 1
    A <- sweep(A, 2, col_norms, "/")
  }

  list(A=A, b=b, w_sqrt=w_sqrt, col_norms=col_norms, times=times,
       var_b=var_b, pooled_var=attr(S, "pooled_var"))
}


library(nnls)

# -----------------------------------------------------------------------------
# Nested-model bootstrap test under weighted NNLS.
#
# The test compares:
#   - a null model excluding one parameter (column col_test),
#   - a full model including all parameters.
#
# The fitting is performed on whitened data:
#   A_w = A * w_sqrt
#   b_w = b * w_sqrt
#
# The test statistic is the reduction in residual sum of squares:
#   T = RSS_null - RSS_full
#
# Its null distribution is approximated by bootstrap:
#   - parametric_whitened: Gaussian perturbations around the null fit,
#   - wild_whitened: sign/weight perturbations of null residuals.
# -----------------------------------------------------------------------------
nnls_nested_boot_whitened <- function(
  A, b, w_sqrt, col_test=4, B=1999, seed=NULL,
  bootstrap=c("parametric_whitened","wild_whitened"),
  wild_dist=c("rademacher","mammen","normal")
){
  bootstrap <- match.arg(bootstrap)
  wild_dist <- match.arg(wild_dist)

  # Optional reproducibility.
  if (!is.null(seed)) set.seed(seed)

  # Ensure standard matrix/vector types.
  A <- as.matrix(A)
  b <- as.numeric(b)
  w_sqrt <- as.numeric(w_sqrt)

  stopifnot(length(b) == nrow(A), length(w_sqrt) == length(b))

  # Sanitize non-finite inputs.
  w_sqrt[!is.finite(w_sqrt)] <- 0
  A[!is.finite(A)] <- 0
  b[!is.finite(b)] <- 0

  # Whitening transformation.
  A_w <- A * w_sqrt
  b_w <- b * w_sqrt

  # Null model excludes the tested column; full model keeps all columns.
  A_null <- A_w[, -col_test, drop=FALSE]
  A_full <- A_w

  # Fit null and full weighted NNLS models.
  fit_null <- tryCatch(nnls(A_null, b_w), error=function(e) NULL)
  fit_full <- tryCatch(nnls(A_full, b_w), error=function(e) NULL)

  # Return NA-based output if fitting fails.
  if (is.null(fit_null) || is.null(fit_full)) {
    return(list(
      p_value=NA_real_, T_obs=NA_real_, T_boot=numeric(),
      rss_null=NA_real_, rss_full=NA_real_,
      coef_full=rep(NA_real_, ncol(A_full)),
      coef_null=rep(NA_real_, ncol(A_null))
    ))
  }

  # Fitted mean under the null model.
  mu0 <- as.vector(A_null %*% coef(fit_null))

  # Residual sums of squares under null and full models.
  rss_null <- sum((b_w - mu0)^2)
  rss_full <- sum((b_w - as.vector(A_full %*% coef(fit_full)))^2)

  # Observed improvement in fit.
  T_obs <- max(0, rss_null - rss_full)

  # Null residuals used by the bootstrap.
  resid_null <- b_w - mu0
  resid_null[!is.finite(resid_null)] <- 0

  # Generator for wild-bootstrap multipliers.
  draw_w <- switch(
    wild_dist,
    "rademacher" = function(n) sample(c(-1,1), n, TRUE),
    "normal"     = function(n) rnorm(n),
    "mammen"     = function(n){
      p <- (sqrt(5)+1)/(2*sqrt(5))
      w1 <- -(sqrt(5)-1)/2
      w2 <- +(sqrt(5)+1)/2
      ifelse(runif(n) < p, w1, w2)
    }
  )

  n <- length(b_w)
  T_boot <- rep(NA_real_, B)

  # Bootstrap loop.
  for (i in seq_len(B)) {

    # Generate bootstrap pseudo-data under the null.
    b_star <- if (bootstrap == "parametric_whitened") {
      # Parametric bootstrap assumes unit-variance Gaussian noise after whitening.
      mu0 + rnorm(n)
    } else {
      # Wild bootstrap reweights the null residuals.
      mu0 + resid_null * draw_w(n)
    }

    # Conservative fallback if numerical issues appear.
    if (any(!is.finite(b_star))) b_star <- mu0

    # Refit null and full models to bootstrap data.
    fit0 <- tryCatch(nnls(A_null, b_star), error=function(e) NULL)
    fit1 <- tryCatch(nnls(A_full, b_star), error=function(e) NULL)

    # If fitting fails, record zero improvement as a conservative choice.
    if (is.null(fit0) || is.null(fit1)) {
      T_boot[i] <- 0
      next
    }

    # Bootstrap test statistic.
    r0 <- sum((b_star - as.vector(A_null %*% coef(fit0)))^2)
    r1 <- sum((b_star - as.vector(A_full %*% coef(fit1)))^2)
    T_boot[i] <- max(0, r0 - r1)
  }

  # Ignore any remaining non-finite bootstrap values.
  T_boot_ok <- T_boot[is.finite(T_boot)]

  if (length(T_boot_ok) == 0) {
    p_val <- NA_real_
  } else {
    # Randomized p-value to break ties in a numerically stable way.
    tol <- 1e-12
    gt <- sum(T_boot_ok >  T_obs + tol)
    eq <- sum(abs(T_boot_ok - T_obs) <= tol)
    p_val <- (1 + gt + runif(1) * eq) / (length(T_boot_ok) + 1)

    # Deterministic alternative:
    # p_val <- (1 + sum(T_boot_ok >= T_obs)) / (length(T_boot_ok) + 1)
  }

  list(
    p_value=p_val,
    T_obs=T_obs,
    T_boot=T_boot_ok,
    rss_null=rss_null,
    rss_full=rss_full,
    coef_full=fit_full$x,
    coef_null=fit_null$x
  )
}


# -----------------------------------------------------------------------------
# High-level wrapper to test sigma_c.
#
# Workflow:
#   1. build the weighted system (A, b, w_sqrt),
#   2. test the significance of column 4 = sigma_c via nested NNLS bootstrap,
#   3. rescale coefficients back if A was column-normalized,
#   4. return estimates and diagnostic summaries.
# -----------------------------------------------------------------------------
test_sigma_nested <- function(
  tsampled_data,
  scaling_A=FALSE,
  t_star=NULL,
  B_n=1999,
  seed=NULL,
  eps_var=1e-8,
  boot_mode="parametric_whitened",
  lambda_var=0.7
){
  # Build regression system and associated weights.
  built <- build_Ab_w_diag_shrink(
    tsampled_data=tsampled_data,
    scaling_A=scaling_A,
    t_star=t_star,
    eps_var=eps_var,
    lambda_var=lambda_var
  )

  A <- built$A
  b <- built$b
  w_sqrt <- built$w_sqrt
  col_norms <- built$col_norms

  # Run the nested bootstrap test for sigma_c (column 4).
  res <- nnls_nested_boot_whitened(
    A,
    b,
    w_sqrt,
    col_test = 4,
    B = B_n,
    bootstrap = boot_mode,
    wild_dist = "rademacher"
  )

  # Rescale full-model coefficients back to the original parameterization.
  coef_full <- res$coef_full / col_norms
  names(coef_full) <- colnames(A)

  # Rescale null-model coefficients back as well.
  cn_null <- colnames(A)[-4]
  coef_null <- res$coef_null / col_norms[-4]
  names(coef_null) <- cn_null

  # Reconstruct a full-length null coefficient vector with sigma_c fixed at zero.
  coef_null_full <- coef_full
  coef_null_full[] <- NA_real_
  coef_null_full[names(coef_null)] <- coef_null
  coef_null_full["sigma_c"] <- 0

  # Number of observations in the stacked system.
  n_obs <- length(b)

  # Information-style and variance-explained diagnostics.
  deltaAIC <- n_obs * log(res$rss_null / res$rss_full) - 2
  DeltaRSS = pmax(res$rss_null - res$rss_full, 0)

  list(
    p.value = res$p_value,
    Sigma = coef_full["sigma_c"],
    Alpha = coef_full["alpha"],
    T.obs = res$T_obs,
    RSS0 = res$rss_null,
    RSS1 = res$rss_full,
    IR = fifelse(res$rss_null > 0, DeltaRSS / res$rss_null, 0),
    # N_obs = n_obs,
    # deltaAIC = deltaAIC,
    coef_full = coef_full,
    coef_null = coef_null_full,
    pooled_var = built$pooled_var
  )
}