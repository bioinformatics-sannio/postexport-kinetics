time_summary_diag_shrink <- function(df,
                                    vars=c("N","N_s","C","C_s"),
                                    eps_var=1e-8,
                                    lambda_var=0.7) {
  stopifnot(lambda_var >= 0, lambda_var <= 1)

  # pooled variance per specie (usa tutto il dataset)
  pooled_var <- sapply(vars, function(v) stats::var(df[[v]], na.rm=TRUE))
  pooled_var[is.na(pooled_var)] <- 0
  pooled_var <- pmax(pooled_var, eps_var)

  # split per time: se df è data.table ok lo stesso
  spl <- split(df, df$time)
  times <- sort(as.numeric(names(spl)))
  out <- vector("list", length(times))
  names(out) <- times

  for (ti in times) {
    DT <- spl[[as.character(ti)]]

    # selezione colonne robusta a data.table
    sub <- as.data.frame(DT[, ..vars])
    n <- nrow(sub)

    mu <- sapply(vars, function(v) mean(sub[[v]], na.rm=TRUE))

    v <- sapply(vars, function(vn) stats::var(sub[[vn]], na.rm=TRUE))
    v[is.na(v)] <- 0
    v <- pmax(v, eps_var)

    # shrinkage verso pooled_var
    v_shr <- (1 - lambda_var) * v + lambda_var * pooled_var
    v_mean <- pmax(v_shr / max(n, 1), eps_var)

    out[[as.character(ti)]] <- list(time=ti, n=n, mean=mu, var_mean=v_mean)
  }

  attr(out, "pooled_var") <- pooled_var
  out
}
build_Ab_w_diag_shrink <- function(tsampled_data,
                                   scaling_A=FALSE,
                                   t_star=NULL,
                                   eps_var=1e-8,
                                   lambda_var=0.7) {
  vars <- c("N","N_s","C","C_s")
  S <- time_summary_diag_shrink(tsampled_data, vars,
                                eps_var=eps_var,
                                lambda_var=lambda_var)
  times <- sort(as.numeric(names(S)))
  K <- length(times)
  if (K < 2) stop("Serve almeno 2 timepoint")

  Xbar <- do.call(rbind, lapply(S, `[[`, "mean"))
  colnames(Xbar) <- vars
  Vbar <- do.call(rbind, lapply(S, `[[`, "var_mean"))
  colnames(Vbar) <- vars

  deltaT <- diff(times)

  dX <- Xbar[-1, , drop=FALSE] - Xbar[-K, , drop=FALSE]
  mX <- (Xbar[-1, , drop=FALSE] + Xbar[-K, , drop=FALSE]) / 2
  I  <- mX * deltaT

  VdX <- Vbar[-1, , drop=FALSE] + Vbar[-K, , drop=FALSE]

  # Termini per R*Δt nella prima equazione (N) per ogni intervallo:
  # default: Δt
  R_dt <- deltaT

  if (!is.null(t_star)) {
    t0 <- times[-K]
    t1 <- times[-1]

    # porzione di intervallo prima di t_star (clippata tra 0 e Δt)
    # se t_star <= t0 -> 0
    # se t_star >= t1 -> Δt
    # altrimenti -> t_star - t0
    R_dt <- pmax(0, pmin(deltaT, t_star - t0))
  }

  # A per intervallo (solo la prima riga usa R_dt invece di deltaT)
  A <- do.call(rbind, lapply(1:(K-1), function(k){
    IN  <- I[k, "N"];   INs <- I[k, "N_s"]; IC <- I[k, "C"]; ICs <- I[k, "C_s"]
    rbind(
      c(R_dt[k], -IN,  0,      0,     -IN, 0,      0),  # <-- qui
      c(0,       0,   -INs,    0,      IN, 0,      0),
      c(0,       IN,   0,     -IC,     0, -IC,     0),
      c(0,       0,    INs,    IC,     0,  0,     -ICs)
    )
  }))
  colnames(A) <- c("R","tau","tau_s","sigma_c","sigma_n","alpha","alpha_s")

  b <- as.vector(t(dX))
  var_b <- as.vector(t(VdX))
  var_b <- pmax(var_b, eps_var)
  w_sqrt <- 1 / sqrt(var_b)

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
nnls_nested_boot_whitened <- function(
  A, b, w_sqrt, col_test=4, B=1999, seed=NULL,
  bootstrap=c("parametric_whitened","wild_whitened"),
  wild_dist=c("rademacher","mammen","normal")
){
  bootstrap <- match.arg(bootstrap)
  wild_dist <- match.arg(wild_dist)
  if (!is.null(seed)) set.seed(seed)

  A <- as.matrix(A); b <- as.numeric(b); w_sqrt <- as.numeric(w_sqrt)

  stopifnot(length(b) == nrow(A), length(w_sqrt) == length(b))

  # sanitize
  w_sqrt[!is.finite(w_sqrt)] <- 0
  A[!is.finite(A)] <- 0
  b[!is.finite(b)] <- 0

  A_w <- A * w_sqrt
  b_w <- b * w_sqrt

  A_null <- A_w[, -col_test, drop=FALSE]
  A_full <- A_w

  fit_null <- tryCatch(nnls(A_null, b_w), error=function(e) NULL)
  fit_full <- tryCatch(nnls(A_full, b_w), error=function(e) NULL)
  if (is.null(fit_null) || is.null(fit_full)) {
    return(list(
      p_value=NA_real_, T_obs=NA_real_, T_boot=numeric(),
      rss_null=NA_real_, rss_full=NA_real_,
      coef_full=rep(NA_real_, ncol(A_full)),
      coef_null=rep(NA_real_, ncol(A_null))
    ))
  }

  mu0 <- as.vector(A_null %*% coef(fit_null))
  rss_null <- sum((b_w - mu0)^2)
  rss_full <- sum((b_w - as.vector(A_full %*% coef(fit_full)))^2)
  T_obs <- max(0, rss_null - rss_full)

  resid_null <- b_w - mu0
  resid_null[!is.finite(resid_null)] <- 0

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

  for (i in seq_len(B)) {
    b_star <- if (bootstrap == "parametric_whitened") {
      mu0 + rnorm(n)
    } else {
      mu0 + resid_null * draw_w(n)
    }

    # se per qualsiasi motivo compaiono NA, forza a mu0 (conservativo)
    if (any(!is.finite(b_star))) b_star <- mu0

    fit0 <- tryCatch(nnls(A_null, b_star), error=function(e) NULL)
    fit1 <- tryCatch(nnls(A_full, b_star), error=function(e) NULL)

    if (is.null(fit0) || is.null(fit1)) {
      T_boot[i] <- 0  # conservativo: conta come "non supera"
      next
    }

    r0 <- sum((b_star - as.vector(A_null %*% coef(fit0)))^2)
    r1 <- sum((b_star - as.vector(A_full %*% coef(fit1)))^2)
    T_boot[i] <- max(0, r0 - r1)
  }

  # calcolo p-value ignorando eventuali NA (se ce ne fossero)
  T_boot_ok <- T_boot[is.finite(T_boot)]
  if (length(T_boot_ok) == 0) {
    p_val <- NA_real_
  } else {
    tol <- 1e-12
    gt <- sum(T_boot_ok >  T_obs + tol)
    eq <- sum(abs(T_boot_ok - T_obs) <= tol)
    p_val <- (1 + gt + runif(1) * eq) / (length(T_boot_ok) + 1) # randomized pvalue to remove ties
    #p_val <- (1 + sum(T_boot_ok >= T_obs)) / (length(T_boot_ok) + 1)
  }

  list(
    p_value=p_val, T_obs=T_obs, T_boot=T_boot_ok,
    rss_null=rss_null, rss_full=rss_full,
    coef_full=fit_full$x, coef_null=fit_null$x
  )
}

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
  built <- build_Ab_w_diag_shrink(
    tsampled_data=tsampled_data,
    scaling_A=scaling_A,
    t_star=t_star,
    eps_var=eps_var,
    lambda_var=lambda_var
  )

  A <- built$A; b <- built$b; w_sqrt <- built$w_sqrt
  col_norms <- built$col_norms
  
  

  res <- nnls_nested_boot_whitened(
    A,
    b,
    w_sqrt,
    col_test = 4,
    B = B_n,
    bootstrap = boot_mode,
    wild_dist = "rademacher"
  )

  coef_full <- res$coef_full / col_norms
  names(coef_full) <- colnames(A)

  cn_null <- colnames(A)[-4]
  coef_null <- res$coef_null / col_norms[-4]
  names(coef_null) <- cn_null

  coef_null_full <- coef_full; coef_null_full[] <- NA_real_
  coef_null_full[names(coef_null)] <- coef_null
  coef_null_full["sigma_c"] <- 0

  n_obs <- length(b)
  deltaAIC <- n_obs * log(res$rss_null / res$rss_full) - 2
  DeltaRSS = pmax(res$rss_null - res$rss_full, 0)
  
  list(
    p.value = res$p_value, 
    Sigma = coef_full["sigma_c"],
    Alpha      = coef_full["alpha"],
    T.obs      = res$T_obs,
    RSS0       = res$rss_null,
    RSS1       = res$rss_full,
    IR = fifelse(res$rss_null > 0, DeltaRSS / res$rss_null, 0),
    #N_obs      = n_obs,
    #deltaAIC   = deltaAIC,
    coef_full  = coef_full,
    coef_null  = coef_null_full,
    pooled_var = built$pooled_var
  )
}