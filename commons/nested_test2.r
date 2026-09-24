# =============================================================================
# Bootstrap-calibrated weighted NNLS test
# for post-export conversion
#
# Revised statistical formulation
# -----------------------------------------------------------------------------
# - destructive biological sampling
# - no pairing of replicates across time points
# - within-time covariance estimation
# - scale-equivariant covariance shrinkage
# - full covariance propagation to interval differences
# - GLS whitening
# - NNLS null/full comparison
# - replicate-level generative bootstrap under H0
# - A*, b*, Sigma_b* reconstructed at every bootstrap iteration
# - deterministic add-one bootstrap p-value
# - explicit NNLS boundary diagnostics
#
# =============================================================================

library(nnls)
library(MASS)


# =============================================================================
# Constants
# =============================================================================

KINETIC_VARS <- c(
  "N",
  "N_s",
  "C",
  "C_s"
)

PARAM_NAMES <- c(
  "R",
  "tau",
  "tau_s",
  "sigma_c",
  "sigma_n",
  "alpha",
  "alpha_s"
)


# =============================================================================
# Matrix utilities
# =============================================================================

make_spd <- function(
  S,
  rel_floor = 1e-8
) {

  S <- as.matrix(S)

  if (
    nrow(S) != ncol(S) ||
    nrow(S) == 0L
  ) {
    return(NULL)
  }

  # numerical symmetry
  S <- (S + t(S)) / 2

  if (any(!is.finite(S))) {
    return(NULL)
  }

  # Scale reference from the diagonal.
  d <- diag(S)

  positive_d <- d[
    is.finite(d) &
      d > 0
  ]

  if (length(positive_d) > 0L) {

    scale_ref <- median(positive_d)

  } else {

    positive_entries <- abs(
      S[
        is.finite(S) &
          S != 0
      ]
    )

    if (length(positive_entries) == 0L) {
      return(NULL)
    }

    scale_ref <- median(positive_entries)
  }

  if (
    !is.finite(scale_ref) ||
      scale_ref <= 0
  ) {
    return(NULL)
  }

  ee <- tryCatch(
    eigen(
      S,
      symmetric = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(ee)) {
    return(NULL)
  }

  # Relative, scale-equivariant floor.
  eig_floor <- max(
    scale_ref * rel_floor,
    .Machine$double.xmin
  )

  eig_values <- pmax(
    ee$values,
    eig_floor
  )

  S_spd <- ee$vectors %*%
    diag(
      eig_values,
      nrow = length(eig_values)
    ) %*%
    t(ee$vectors)

  S_spd <- (
    S_spd +
      t(S_spd)
  ) / 2

  dimnames(S_spd) <- dimnames(S)

  S_spd
}


inverse_sqrt_matrix <- function(
  Sigma,
  rel_floor = 1e-8
) {

  Sigma <- make_spd(
    Sigma,
    rel_floor = rel_floor
  )

  if (is.null(Sigma)) {
    return(NULL)
  }

  ee <- tryCatch(
    eigen(
      Sigma,
      symmetric = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(ee)) {
    return(NULL)
  }

  if (
    any(!is.finite(ee$values)) ||
      any(ee$values <= 0)
  ) {
    return(NULL)
  }

  W12 <- ee$vectors %*%
    diag(
      1 / sqrt(ee$values),
      nrow = length(ee$values)
    ) %*%
    t(ee$vectors)

  (
    W12 +
      t(W12)
  ) / 2
}


# =============================================================================
# Destructive-sampling summaries
# =============================================================================

time_summary_cov_shrink <- function(
  df,
  vars = KINETIC_VARS,
  lambda_time = 0.5,
  lambda_diag = 0.1,
  rel_floor = 1e-8
) {

  stopifnot(
    lambda_time >= 0,
    lambda_time <= 1,
    lambda_diag >= 0,
    lambda_diag <= 1
  )

  df <- as.data.frame(df)

  required_cols <- c(
    "time",
    vars
  )

  if (!all(required_cols %in% names(df))) {

    stop(
      paste(
        "Missing columns:",
        paste(
          setdiff(
            required_cols,
            names(df)
          ),
          collapse = ", "
        )
      )
    )
  }

  times <- sort(
    unique(df$time)
  )

  K <- length(times)
  P <- length(vars)

  if (K < 2L) {
    stop("At least two time points are required.")
  }

  means <- matrix(
    NA_real_,
    nrow = K,
    ncol = P,
    dimnames = list(
      as.character(times),
      vars
    )
  )

  n_rep <- integer(K)

  raw_cov <- vector(
    "list",
    K
  )

  # ---------------------------------------------------------------------------
  # Pooled WITHIN-TIME covariance.
  #
  # Mean changes between time points are explicitly excluded.
  # ---------------------------------------------------------------------------

  pooled_scatter <- matrix(
    0,
    nrow = P,
    ncol = P,
    dimnames = list(
      vars,
      vars
    )
  )

  pooled_df <- 0L

  for (k in seq_along(times)) {

    ti <- times[k]

    sub <- df[
      df$time == ti,
      vars,
      drop = FALSE
    ]

    sub <- sub[
      complete.cases(sub),
      ,
      drop = FALSE
    ]

    nk <- nrow(sub)

    if (nk == 0L) {

      stop(
        paste(
          "No complete observations at time",
          ti
        )
      )
    }

    n_rep[k] <- nk

    means[k, ] <- colMeans(sub)

    if (nk >= 2L) {

      Sk <- stats::cov(sub)

      if (
        all(dim(Sk) == c(P, P)) &&
          all(is.finite(Sk))
      ) {

        raw_cov[[k]] <- Sk

        pooled_scatter <-
          pooled_scatter +
          (nk - 1L) * Sk

        pooled_df <-
          pooled_df +
          (nk - 1L)
      }
    }
  }

  if (pooled_df <= 0L) {

    stop(
      "Insufficient replication to estimate within-time covariance."
    )
  }

  S_pool <-
    pooled_scatter /
    pooled_df

  # ---------------------------------------------------------------------------
  # Scale-equivariant shrinkage toward diagonal covariance.
  # ---------------------------------------------------------------------------

  S_diag <- diag(
    diag(S_pool)
  )

  dimnames(S_diag) <-
    dimnames(S_pool)

  S_pool_shr <-
    (1 - lambda_diag) *
      S_pool +
    lambda_diag *
      S_diag

  S_pool_shr <- make_spd(
    S_pool_shr,
    rel_floor = rel_floor
  )

  if (is.null(S_pool_shr)) {

    stop(
      "Unable to stabilize pooled covariance."
    )
  }

  # ---------------------------------------------------------------------------
  # Shrink individual time-point covariance matrices toward pooled covariance.
  # ---------------------------------------------------------------------------

  cov_obs <- vector(
    "list",
    K
  )

  cov_mean <- vector(
    "list",
    K
  )

  for (k in seq_along(times)) {

    nk <- n_rep[k]

    if (
      nk >= 2L &&
        !is.null(raw_cov[[k]])
    ) {

      Sk <-
        (1 - lambda_time) *
          raw_cov[[k]] +
        lambda_time *
          S_pool_shr

    } else {

      Sk <- S_pool_shr
    }

    Sk <- make_spd(
      Sk,
      rel_floor = rel_floor
    )

    if (is.null(Sk)) {

      stop(
        paste(
          "Unable to stabilize covariance at time",
          times[k]
        )
      )
    }

    # Covariance of individual destructive samples.
    cov_obs[[k]] <- Sk

    # Covariance of the sample mean.
    cov_mean[[k]] <-
      Sk / nk
  }

  names(cov_obs) <-
    as.character(times)

  names(cov_mean) <-
    as.character(times)

  list(
    times = times,
    means = means,
    n_rep = n_rep,
    cov_obs = cov_obs,
    cov_mean = cov_mean,
    pooled_cov = S_pool_shr
  )
}


# =============================================================================
# Covariance propagation
# =============================================================================

build_sigma_means <- function(
  summary_obj
) {

  K <- length(
    summary_obj$times
  )

  P <- length(
    KINETIC_VARS
  )

  # Species-major ordering:
  #
  # N(t1)...N(tK),
  # Ns(t1)...Ns(tK),
  # C(t1)...C(tK),
  # Cs(t1)...Cs(tK)

  Sigma_m <- matrix(
    0,
    nrow = P * K,
    ncol = P * K
  )

  idx <- function(
    species,
    time_index
  ) {

    (species - 1L) *
      K +
      time_index
  }

  for (k in seq_len(K)) {

    Sk <- summary_obj$cov_mean[[k]]

    if (
      is.null(Sk) ||
        any(!is.finite(Sk))
    ) {

      stop(
        paste(
          "Invalid covariance at time",
          summary_obj$times[k]
        )
      )
    }

    for (a in seq_len(P)) {

      for (bb in seq_len(P)) {

        ia <- idx(
          a,
          k
        )

        ib <- idx(
          bb,
          k
        )

        Sigma_m[
          ia,
          ib
        ] <- Sk[
          a,
          bb
        ]
      }
    }
  }

  Sigma_m
}


build_difference_matrix <- function(
  K
) {

  if (K < 2L) {
    stop("K must be at least 2.")
  }

  D1 <- matrix(
    0,
    nrow = K - 1L,
    ncol = K
  )

  for (
    k in seq_len(
      K - 1L
    )
  ) {

    D1[
      k,
      k
    ] <- -1

    D1[
      k,
      k + 1L
    ] <- 1
  }

  kronecker(
    diag(
      length(KINETIC_VARS)
    ),
    D1
  )
}


# =============================================================================
# Construct A, b and Sigma_b
# =============================================================================

build_Ab_fullcov <- function(
  tsampled_data,
  scaling_A = TRUE,
  t_star = NULL,
  lambda_time = 0.5,
  lambda_diag = 0.1,
  rel_floor = 1e-8
) {

  S <- time_summary_cov_shrink(
    df = tsampled_data,
    lambda_time = lambda_time,
    lambda_diag = lambda_diag,
    rel_floor = rel_floor
  )

  times <- S$times

  K <- length(times)

  X <- S$means

  deltaT <- diff(times)

  if (
    any(!is.finite(deltaT)) ||
      any(deltaT <= 0)
  ) {

    stop(
      "Times must be strictly increasing."
    )
  }

  # ---------------------------------------------------------------------------
  # Observed interval changes.
  # ---------------------------------------------------------------------------

  deltaN <-
    diff(
      X[, "N"]
    )

  deltaNs <-
    diff(
      X[, "N_s"]
    )

  deltaC <-
    diff(
      X[, "C"]
    )

  deltaCs <-
    diff(
      X[, "C_s"]
    )

  # ---------------------------------------------------------------------------
  # Trapezoidal means.
  # ---------------------------------------------------------------------------

  meanN <-
    (
      X[-1, "N"] +
        X[-K, "N"]
    ) / 2

  meanNs <-
    (
      X[-1, "N_s"] +
        X[-K, "N_s"]
    ) / 2

  meanC <-
    (
      X[-1, "C"] +
        X[-K, "C"]
    ) / 2

  meanCs <-
    (
      X[-1, "C_s"] +
        X[-K, "C_s"]
    ) / 2

  IN <-
    deltaT *
    meanN

  INs <-
    deltaT *
    meanNs

  IC <-
    deltaT *
    meanC

  ICs <-
    deltaT *
    meanCs

  # ---------------------------------------------------------------------------
  # Active transcription duration.
  #
  # Only the R term is truncated by shutoff.
  # All other kinetic processes continue for the full interval.
  # ---------------------------------------------------------------------------

  R_dt <- deltaT

  if (!is.null(t_star)) {

    if (
      length(t_star) != 1L ||
        !is.finite(t_star)
    ) {

      stop(
        "t_star must be NULL or a single finite number."
      )
    }

    t0 <- times[-K]

    R_dt <- pmax(
      0,
      pmin(
        deltaT,
        t_star - t0
      )
    )
  }

  zeros <- rep(
    0,
    K - 1L
  )

  # ---------------------------------------------------------------------------
  # A matrix.
  #
  # Row ordering:
  #
  # N intervals
  # Ns intervals
  # C intervals
  # Cs intervals
  # ---------------------------------------------------------------------------

  A_N <- cbind(
    R_dt,
    -IN,
    zeros,
    zeros,
    -IN,
    zeros,
    zeros
  )

  A_Ns <- cbind(
    zeros,
    zeros,
    -INs,
    zeros,
    IN,
    zeros,
    zeros
  )

  A_C <- cbind(
    zeros,
    IN,
    zeros,
    -IC,
    zeros,
    -IC,
    zeros
  )

  A_Cs <- cbind(
    zeros,
    zeros,
    INs,
    IC,
    zeros,
    zeros,
    -ICs
  )

  A <- rbind(
    A_N,
    A_Ns,
    A_C,
    A_Cs
  )

  colnames(A) <-
    PARAM_NAMES

  b <- c(
    deltaN,
    deltaNs,
    deltaC,
    deltaCs
  )

  # ---------------------------------------------------------------------------
  # Covariance propagation.
  # ---------------------------------------------------------------------------

  Sigma_m <-
    build_sigma_means(S)

  D <-
    build_difference_matrix(K)

  if (
    ncol(D) !=
      nrow(Sigma_m)
  ) {

    stop(
      "Dimension mismatch between D and Sigma_m."
    )
  }

  Sigma_b <-
    D %*%
    Sigma_m %*%
    t(D)

  if (
    nrow(Sigma_b) !=
      length(b)
  ) {

    stop(
      "Dimension mismatch between Sigma_b and b."
    )
  }

  Sigma_b <- make_spd(
    Sigma_b,
    rel_floor = rel_floor
  )

  if (is.null(Sigma_b)) {

    stop(
      "Unable to stabilize Sigma_b."
    )
  }

  # ---------------------------------------------------------------------------
  # Optional column scaling.
  # ---------------------------------------------------------------------------

  col_norms <- rep(
    1,
    ncol(A)
  )

  if (scaling_A) {

    col_norms <-
      sqrt(
        colSums(
          A^2
        )
      )

    invalid <- (
      !is.finite(col_norms) |
        col_norms <= 0
    )

    col_norms[
      invalid
    ] <- 1

    A <- sweep(
      A,
      2,
      col_norms,
      "/"
    )
  }

  if (
    nrow(A) != length(b) ||
      nrow(Sigma_b) != length(b) ||
      ncol(Sigma_b) != length(b)
  ) {

    stop(
      "Internal A/b/Sigma_b dimension mismatch."
    )
  }

  list(
    A = A,
    b = b,
    Sigma_b = Sigma_b,
    Sigma_m = Sigma_m,
    D = D,
    col_norms = col_norms,
    summary = S
  )
}


# =============================================================================
# NNLS null/full fit
# =============================================================================

fit_nnls_nested_once <- function(
  A,
  b,
  Sigma_b,
  col_test = 4L,
  rel_floor = 1e-8
) {

  A <-
    as.matrix(A)

  b <-
    as.numeric(b)

  Sigma_b <-
    as.matrix(Sigma_b)

  if (
    nrow(A) != length(b) ||
      nrow(Sigma_b) != length(b) ||
      ncol(Sigma_b) != length(b)
  ) {

    return(NULL)
  }

  if (
    any(!is.finite(A)) ||
      any(!is.finite(b)) ||
      any(!is.finite(Sigma_b))
  ) {

    return(NULL)
  }

  W12 <- inverse_sqrt_matrix(
    Sigma_b,
    rel_floor = rel_floor
  )

  if (is.null(W12)) {
    return(NULL)
  }

  A_w <-
    W12 %*%
    A

  b_w <-
    as.vector(
      W12 %*%
        b
    )

  A0 <- A_w[
    ,
    -col_test,
    drop = FALSE
  ]

  fit0 <- tryCatch(
    nnls::nnls(
      A0,
      b_w
    ),
    error = function(e) NULL
  )

  fit1 <- tryCatch(
    nnls::nnls(
      A_w,
      b_w
    ),
    error = function(e) NULL
  )

  if (
    is.null(fit0) ||
      is.null(fit1)
  ) {

    return(NULL)
  }

  pred0 <-
    as.vector(
      A0 %*%
        fit0$x
    )

  pred1 <-
    as.vector(
      A_w %*%
        fit1$x
    )

  rss0 <-
    sum(
      (b_w - pred0)^2
    )

  rss1 <-
    sum(
      (b_w - pred1)^2
    )

  if (
    !is.finite(rss0) ||
      !is.finite(rss1)
  ) {

    return(NULL)
  }

  Tstat <- max(
    0,
    rss0 - rss1
  )

  # ---------------------------------------------------------------------------
  # Identifiability diagnostics.
  # ---------------------------------------------------------------------------

  singular_values <- tryCatch(
    svd(A_w)$d,
    error = function(e) numeric()
  )

  if (
    length(singular_values) > 0L &&
      all(is.finite(singular_values))
  ) {

    max_sv <-
      max(singular_values)

    min_sv <-
      min(singular_values)

    condition_number <- if (
      min_sv > 0
    ) {

      max_sv / min_sv

    } else {

      Inf
    }

  } else {

    min_sv <-
      NA_real_

    condition_number <-
      Inf
  }

  list(
    T = Tstat,

    RSS0 = rss0,
    RSS1 = rss1,

    coef_null_scaled =
      fit0$x,

    coef_full_scaled =
      fit1$x,

    condition_number =
      condition_number,

    min_singular_value =
      min_sv,

    rank_full =
      qr(A_w)$rank,

    rank_null =
      qr(A0)$rank
  )
}


# =============================================================================
# Crank-Nicolson null trajectory
# =============================================================================

# -----------------------------------------------------------------------------
# Continuous-time kinetic matrix.
#
# State order:
# N, N_s, C, C_s
# -----------------------------------------------------------------------------

kinetic_matrix <- function(
  theta
) {

  matrix(
    c(
      -(theta["sigma_n"] + theta["tau"]),
      0,
      0,
      0,

      theta["sigma_n"],
      -theta["tau_s"],
      0,
      0,

      theta["tau"],
      0,
      -(theta["sigma_c"] + theta["alpha"]),
      0,

      0,
      theta["tau_s"],
      theta["sigma_c"],
      -theta["alpha_s"]
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      KINETIC_VARS,
      KINETIC_VARS
    )
  )
}


# -----------------------------------------------------------------------------
# One Crank-Nicolson interval.
#
# The equation is:
#
# x1 - x0 =
#     dt/2 * K * (x0 + x1)
#     + input
#
# therefore:
#
# (I - dt/2 K)x1 =
#     (I + dt/2 K)x0 + input
#
# R_active_dt is the duration over which transcription remains active
# inside the interval.
# -----------------------------------------------------------------------------

cn_interval <- function(
  x0,
  dt,
  theta,
  R_active_dt
) {

  if (
    !is.finite(dt) ||
      dt <= 0
  ) {

    stop(
      "Invalid interval length."
    )
  }

  Kmat <-
    kinetic_matrix(theta)

  I4 <-
    diag(4)

  lhs <-
    I4 -
    (dt / 2) *
      Kmat

  rhs <-
    (
      I4 +
      (dt / 2) *
        Kmat
    ) %*%
    x0

  rhs <-
    as.numeric(rhs)

  rhs[1] <-
    rhs[1] +
    theta["R"] *
      R_active_dt

  x1 <- tryCatch(
    solve(
      lhs,
      rhs
    ),
    error = function(e) NULL
  )

  if (
    is.null(x1) ||
      any(!is.finite(x1))
  ) {

    stop(
      "Crank-Nicolson propagation failed."
    )
  }

  names(x1) <-
    KINETIC_VARS

  x1
}


predict_null_cn <- function(
  times,
  x0,
  theta0,
  t_star = NULL
) {

  times <-
    sort(
      unique(times)
    )

  theta0 <-
    as.numeric(theta0)

  names(theta0) <-
    PARAM_NAMES

  theta0["sigma_c"] <- 0

  if (
    any(!is.finite(theta0)) ||
      any(theta0 < 0)
  ) {

    stop(
      "Invalid null parameter estimates."
    )
  }

  x0 <-
    as.numeric(x0)

  names(x0) <-
    KINETIC_VARS

  X <- matrix(
    NA_real_,
    nrow = length(times),
    ncol = 4,
    dimnames = list(
      as.character(times),
      KINETIC_VARS
    )
  )

  X[1, ] <-
    x0

  if (length(times) == 1L) {
    return(X)
  }

  for (
    k in seq_len(
      length(times) - 1L
    )
  ) {

    t0 <-
      times[k]

    t1 <-
      times[k + 1L]

    dt <-
      t1 - t0

    R_active_dt <-
      dt

    if (!is.null(t_star)) {

      R_active_dt <- max(
        0,
        min(
          dt,
          t_star - t0
        )
      )
    }

    X[
      k + 1L,
    ] <- cn_interval(
      x0 = X[k, ],
      dt = dt,
      theta = theta0,
      R_active_dt =
        R_active_dt
    )
  }

  X
}


# =============================================================================
# Destructive-sampling bootstrap generator
# =============================================================================

simulate_destructive_null <- function(
  null_means,
  summary_obj,
  truncate_nonnegative = FALSE
) {

  times <-
    summary_obj$times

  if (
    nrow(null_means) !=
      length(times) ||
      ncol(null_means) !=
      length(KINETIC_VARS)
  ) {

    stop(
      "null_means dimensions are inconsistent."
    )
  }

  out <- vector(
    "list",
    length(times)
  )

  for (
    k in seq_along(times)
  ) {

    nk <-
      summary_obj$n_rep[k]

    Sigma_k <-
      summary_obj$cov_obs[[k]]

    if (
      nk < 1L ||
        is.null(Sigma_k) ||
        any(!is.finite(Sigma_k))
    ) {

      stop(
        paste(
          "Invalid bootstrap parameters at time",
          times[k]
        )
      )
    }

    # Biological samples from DIFFERENT time points
    # are generated independently.

    Y <- MASS::mvrnorm(
      n = nk,
      mu = null_means[k, ],
      Sigma = Sigma_k
    )

    if (nk == 1L) {

      Y <- matrix(
        Y,
        nrow = 1L
      )
    }

    colnames(Y) <-
      KINETIC_VARS

    if (truncate_nonnegative) {

      Y[
        Y < 0
      ] <- 0
    }

    if (any(!is.finite(Y))) {

      stop(
        "Non-finite bootstrap observations."
      )
    }

    out[[k]] <- data.frame(
      time =
        rep(
          times[k],
          nk
        ),

      # Local label only.
      #
      # There is NO relation between t1_r1 and t2_r1.
      replicate =
        paste0(
          "t",
          k,
          "_r",
          seq_len(nk)
        ),

      N =
        Y[, "N"],

      N_s =
        Y[, "N_s"],

      C =
        Y[, "C"],

      C_s =
        Y[, "C_s"],

      stringsAsFactors =
        FALSE
    )
  }

  do.call(
    rbind,
    out
  )
}


# =============================================================================
# Main bootstrap test
# =============================================================================

test_sigma_nested <- function(
  tsampled_data,

  scaling_A = TRUE,

  t_star = NULL,

  B_n = 1999,

  seed = NULL,

  # New covariance shrinkage parameters
  lambda_time = 0.5,
  lambda_diag = 0.1,

  # Backward-compatible alias:
  # if provided, lambda_var overrides lambda_time.
  lambda_var = NULL,

  rel_floor = 1e-8,

  truncate_nonnegative_boot = FALSE,

  max_failure_rate = 0.05,

  return_boot = FALSE,

  verbose = FALSE
) {

  # ---------------------------------------------------------------------------
  # Compatibility with old experimental script.
  # ---------------------------------------------------------------------------

  if (!is.null(lambda_var)) {

    lambda_time <-
      lambda_var
  }

  if (
    length(B_n) != 1L ||
      !is.finite(B_n) ||
      B_n < 1
  ) {

    stop(
      "B_n must be a positive integer."
    )
  }

  B_n <-
    as.integer(B_n)

  if (!is.null(seed)) {

    set.seed(seed)
  }

  # ===========================================================================
  # Observed system
  # ===========================================================================

  obs_error <-
    NULL

  built_obs <- tryCatch(

    build_Ab_fullcov(
      tsampled_data =
        tsampled_data,

      scaling_A =
        scaling_A,

      t_star =
        t_star,

      lambda_time =
        lambda_time,

      lambda_diag =
        lambda_diag,

      rel_floor =
        rel_floor
    ),

    error = function(e) {

      obs_error <<-
        conditionMessage(e)

      NULL
    }
  )

  if (is.null(built_obs)) {

    return(
      list(
        p.value =
          NA_real_,

        status =
          "observed_system_failed",

        error.message =
          obs_error
      )
    )
  }

  # ===========================================================================
  # Observed NNLS fit
  # ===========================================================================

  fit_obs <- fit_nnls_nested_once(
    A =
      built_obs$A,

    b =
      built_obs$b,

    Sigma_b =
      built_obs$Sigma_b,

    col_test =
      4L,

    rel_floor =
      rel_floor
  )

  if (is.null(fit_obs)) {

    return(
      list(
        p.value =
          NA_real_,

        status =
          "observed_fit_failed",

        error.message =
          "Observed GLS-NNLS fit failed."
      )
    )
  }

  # ===========================================================================
  # Rescale coefficients
  # ===========================================================================

  coef_full <-
    fit_obs$coef_full_scaled /
    built_obs$col_norms

  names(coef_full) <-
    PARAM_NAMES

  null_names <-
    PARAM_NAMES[-4L]

  coef_null_short <-
    fit_obs$coef_null_scaled /
    built_obs$col_norms[-4L]

  names(coef_null_short) <-
    null_names

  coef_null <- setNames(
    rep(
      0,
      length(PARAM_NAMES)
    ),
    PARAM_NAMES
  )

  coef_null[
    null_names
  ] <- coef_null_short

  coef_null[
    "sigma_c"
  ] <- 0

  # ===========================================================================
  # Null mean trajectory
  # ===========================================================================

  times <-
    built_obs$summary$times

  x0 <-
    built_obs$summary$means[
      1,
      KINETIC_VARS
    ]

  null_error <-
    NULL

  null_means <- tryCatch(

    predict_null_cn(
      times =
        times,

      x0 =
        x0,

      theta0 =
        coef_null,

      t_star =
        t_star
    ),

    error = function(e) {

      null_error <<-
        conditionMessage(e)

      NULL
    }
  )

  if (
    is.null(null_means) ||
      any(!is.finite(null_means))
  ) {

    return(
      list(
        p.value =
          NA_real_,

        Sigma =
          coef_full["sigma_c"],

        Alpha =
          coef_full["alpha"],

        T.obs =
          fit_obs$T,

        RSS0 =
          fit_obs$RSS0,

        RSS1 =
          fit_obs$RSS1,

        status =
          "null_trajectory_failed",

        error.message =
          null_error
      )
    )
  }

  # ===========================================================================
  # Bootstrap
  # ===========================================================================

  T_boot <-
    rep(
      NA_real_,
      B_n
    )

  condition_boot <-
    rep(
      NA_real_,
      B_n
    )

  rank_boot <-
    rep(
      NA_integer_,
      B_n
    )

  failed <-
    logical(
      B_n
    )

  for (
    bb in seq_len(B_n)
  ) {

    # -------------------------------------------------------------------------
    # New independent destructive samples.
    # -------------------------------------------------------------------------

    data_star <- tryCatch(

      simulate_destructive_null(
        null_means =
          null_means,

        summary_obj =
          built_obs$summary,

        truncate_nonnegative =
          truncate_nonnegative_boot
      ),

      error = function(e)
        NULL
    )

    if (is.null(data_star)) {

      failed[bb] <-
        TRUE

      next
    }

    # -------------------------------------------------------------------------
    # Reconstruct A*, b*, Sigma_b*.
    # -------------------------------------------------------------------------

    built_star <- tryCatch(

      build_Ab_fullcov(
        tsampled_data =
          data_star,

        scaling_A =
          scaling_A,

        t_star =
          t_star,

        lambda_time =
          lambda_time,

        lambda_diag =
          lambda_diag,

        rel_floor =
          rel_floor
      ),

      error = function(e)
        NULL
    )

    if (is.null(built_star)) {

      failed[bb] <-
        TRUE

      next
    }

    # -------------------------------------------------------------------------
    # Null/full NNLS on bootstrap data.
    # -------------------------------------------------------------------------

    fit_star <-
      fit_nnls_nested_once(
        A =
          built_star$A,

        b =
          built_star$b,

        Sigma_b =
          built_star$Sigma_b,

        col_test =
          4L,

        rel_floor =
          rel_floor
      )

    if (
      is.null(fit_star) ||
        !is.finite(fit_star$T)
    ) {

      failed[bb] <-
        TRUE

      next
    }

    T_boot[bb] <-
      fit_star$T

    condition_boot[bb] <-
      fit_star$condition_number

    rank_boot[bb] <-
      fit_star$rank_full
  }

  valid <-
    !failed &
    is.finite(T_boot)

  T_ok <-
    T_boot[
      valid
    ]

  condition_ok <-
    condition_boot[
      valid
    ]

  rank_ok <-
    rank_boot[
      valid
    ]

  failure_rate <-
    mean(
      !valid
    )

  if (
    length(T_ok) == 0L ||
      failure_rate >
      max_failure_rate
  ) {

    return(
      list(
        p.value =
          NA_real_,

        Sigma =
          coef_full["sigma_c"],

        Alpha =
          coef_full["alpha"],

        T.obs =
          fit_obs$T,

        RSS0 =
          fit_obs$RSS0,

        RSS1 =
          fit_obs$RSS1,

        coef_full =
          coef_full,

        coef_null =
          coef_null,

        bootstrap.failure.rate =
          failure_rate,

        n.bootstrap.valid =
          length(T_ok),

        condition.number =
          fit_obs$condition_number,

        rank.full =
          fit_obs$rank_full,

        rank.null =
          fit_obs$rank_null,

        status =
          "bootstrap_unstable",

        error.message =
          paste(
            "Bootstrap failure rate =",
            signif(
              failure_rate,
              4
            )
          )
      )
    )
  }

  # ===========================================================================
  # Boundary tolerance
  # ===========================================================================

  # T = RSS0 - RSS1 is computed on the whitened scale.
  #
  # Values much smaller than numerical precision relative to the RSS scale
  # are interpreted as the NNLS boundary T = 0.

  rss_scale <- max(
    1,
    abs(
      fit_obs$RSS0
    ),
    abs(
      fit_obs$RSS1
    )
  )

  tol_zero <-
    1e-10 *
    rss_scale

  # ===========================================================================
  # Bootstrap p-value
  # ===========================================================================

  if (
    fit_obs$T <=
      tol_zero
  ) {

    p_val <-
      1

  } else {

    p_val <-
      (
        1 +
        sum(
          T_ok >=
            fit_obs$T
        )
      ) /
      (
        length(T_ok) +
        1
      )
  }

  # Empirical point mass at boundary.
  atom_zero <-
    mean(
      T_ok <=
        tol_zero
    )

  # ===========================================================================
  # Fit improvement
  # ===========================================================================

  DeltaRSS <- max(
    fit_obs$RSS0 -
      fit_obs$RSS1,
    0
  )

  IR <- if (
    is.finite(
      fit_obs$RSS0
    ) &&
      fit_obs$RSS0 >
      0
  ) {

    DeltaRSS /
      fit_obs$RSS0

  } else {

    NA_real_
  }

  # ===========================================================================
  # Bootstrap diagnostics
  # ===========================================================================

  condition_median <-
    median(
      condition_ok,
      na.rm = TRUE
    )

  condition_q95 <-
    as.numeric(
      quantile(
        condition_ok,
        probs = 0.95,
        na.rm = TRUE,
        names = FALSE
      )
    )

  condition_max <-
    max(
      condition_ok,
      na.rm = TRUE
    )

  rank_deficient_fraction <-
    mean(
      rank_ok <
        ncol(
          built_obs$A
        ),
      na.rm = TRUE
    )

  # ===========================================================================
  # Output
  # ===========================================================================

  retval <- list(
    p.value =
      p_val,

    Sigma =
      coef_full["sigma_c"],

    Alpha =
      coef_full["alpha"],

    T.obs =
      fit_obs$T,

    RSS0 =
      fit_obs$RSS0,

    RSS1 =
      fit_obs$RSS1,

    IR =
      IR,

    coef_full =
      coef_full,

    coef_null =
      coef_null,

    # Boundary diagnostic
    atom.zero =
      atom_zero,

    boundary.tolerance =
      tol_zero,

    # Bootstrap diagnostics
    bootstrap.failure.rate =
      failure_rate,

    n.bootstrap.valid =
      length(T_ok),

    bootstrap.condition.median =
      condition_median,

    bootstrap.condition.q95 =
      condition_q95,

    bootstrap.condition.max =
      condition_max,

    bootstrap.rank.deficient.fraction =
      rank_deficient_fraction,

    # Observed identifiability diagnostics
    condition.number =
      fit_obs$condition_number,

    min.singular.value =
      fit_obs$min_singular_value,

    rank.full =
      fit_obs$rank_full,

    rank.null =
      fit_obs$rank_null,

    # Shrinkage
    lambda.time =
      lambda_time,

    lambda.diag =
      lambda_diag,

    # Other useful information
    null.means =
      null_means,

    pooled.covariance =
      built_obs$summary$pooled_cov,

    n.replicates.by.time =
      built_obs$summary$n_rep,

    times =
      times,

    status =
      "ok",

    error.message =
      NULL
  )

  if (return_boot) {

    retval$T.boot <-
      T_ok

    retval$bootstrap.condition <-
      condition_ok

    retval$bootstrap.rank <-
      rank_ok
  }

  retval
}