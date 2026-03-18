
baseline_psi_logit_ttest_pre_vs_post <- function(
  tsampled_data, ss_data,
  scaling_A = FALSE,
  usa_W = TRUE,
  add_SS = FALSE,
  B_n = -1,
  t_star = NULL,          # separatore: pre <= t_star, post > t_star
  eps = 1e-6,
  min_total = 1L
){
  tmain <- tsampled_data[order(tsampled_data$replicate, tsampled_data$time), ]

  req <- c("time","replicate","C","C_s")
  miss <- setdiff(req, names(tmain))
  if (length(miss) > 0) stop("tsampled_data is missing columns: ", paste(miss, collapse=", "))

  times_all <- sort(unique(tmain$time))
  if (length(times_all) < 2) stop("Need at least 2 distinct time points.")

  # define pre/post
  if (!is.null(t_star)) {
    pre_times  <- times_all[times_all <= t_star]
    post_times <- times_all[times_all >  t_star]
    if (length(pre_times) < 1 || length(post_times) < 1) {
      stop("t_star yields empty pre or post window.")
    }
  } else {
    # if not provided: pre = first time only, post = rest
    pre_times  <- times_all[1]
    post_times <- times_all[-1]
  }

  # aggregate counts per replicate in pre and post
  dt_pre <- tmain[time %in% pre_times,
                  .(C_pre = sum(C, na.rm=TRUE), Cs_pre = sum(C_s, na.rm=TRUE)),
                  by = replicate]
  dt_post <- tmain[time %in% post_times,
                   .(C_post = sum(C, na.rm=TRUE), Cs_post = sum(C_s, na.rm=TRUE)),
                   by = replicate]

  dt <- merge(dt_pre, dt_post, by="replicate", all=FALSE)
  dt <- dt[((C_pre + Cs_pre) > min_total) & ((C_post + Cs_post) > min_total)]

  if (nrow(dt) < 2) {
    return(list(Tab_coeff=as.data.frame(dt), Sigma=NA_real_, p.value=NA_real_, T.obs=NA_real_))
  }

  psi_pre  <- dt$C_pre  / (dt$C_pre  + dt$Cs_pre  + eps)
  psi_post <- dt$C_post / (dt$C_post + dt$Cs_post + eps)

  logit <- function(x) log((x + eps) / (1 - x + eps))
  d <- logit(psi_post) - logit(psi_pre)

  # one-sided: PSI decreases => post < pre => d < 0
  tt <- t.test(d, mu = 0, alternative = "two.sided")  

  # effect size on PSI scale
  dt[, `:=`(
    PSI_pre = psi_pre,
    PSI_post = psi_post,
    dPSI = psi_post - psi_pre,
    score = -(psi_post - psi_pre)
  )]
  Sigma <- median(dt$score, na.rm=TRUE)

  list(
    Tab_coeff = as.data.frame(dt),
    Sigma = Sigma,
    Alpha = Sigma,
    p.value = tt$p.value,
    RSS0 = 1,
    RSS1 = 1,
    IR=1,
    T.obs = as.numeric(tt$statistic),
    deltaAIC = as.numeric(tt$statistic)
  )
}
