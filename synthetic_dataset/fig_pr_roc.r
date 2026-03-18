setwd("~/postexport-kinetics/synthetic_dataset")

require(data.table)
library(ggplot2)
library(scales)

load("test_sigma_2k_new.rdata")


# stack in long format
dt_test_results <- melt(
  dt_test_results,
  #measure.vars = c("score_unified_soft","score_aiclogp","score_siglogp","score_aic","score_neglogp","score_sigma"),
  #measure.vars = c("score_unified_soft","score_siglogp"),
  measure.vars = c("score_siglog_IR_logq","score_sig_IR_logq","score_siglogq"),
  variable.name = "Score_method",
  value.name    = "Score"
)


# ***********************
# FIG2 Analisi di Prioritizzazione
# PR-AUC ROC-AUC Top 50 Top 100
# only Nested
# per mostrare SHUTOFF migliore e dipendenza da noise e platform ****
# in supp ROC-AUC e tab completa
# ***********************


# ---- helpers ----
pr_auc <- function(score, y01) {
  s_pos <- score[y01 == 1L]
  s_neg <- score[y01 == 0L]
  if (length(s_pos) == 0L || length(s_neg) == 0L) return(NA_real_)
  PRROC::pr.curve(scores.class0 = s_pos,
                  scores.class1 = s_neg,
                  curve = FALSE)$auc.integral
}

roc_auc <- function(score, y01) {
  if (length(unique(y01)) < 2L) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(response = y01,
                                 predictor = score,
                                 quiet = TRUE)))
}

top_pr <- function(score, y01, K = 100L) {
  topk_idx = order(score, decreasing = TRUE)[1:K]
  TP = sum(y01[topk_idx] == 1L)
  FP = sum(y01[topk_idx] == 0L)
  TN = sum(y01[-topk_idx] == 0L)
  FN = sum(y01[-topk_idx] == 1L)
  P = TP / (TP + FP)
  R = TP / (TP + FN)
  list(precision = P, recall = R)
}


# filtro scenario
w_dt <- dt_test_results[
  #N_tsamples == 10 &
  N_replicates == 5 &
  Score_method=="score_siglogq" &
  Tsteps == 10 & 
  Lambda == 0.5 &
  #Platform %in% c("RNA-seq") &
  #Exprs_noise %in% c("Medium") &
  Test_method %in% c("Nested-wild_whitened")  #,"PSIBaseline")  &
  #Perturbation == "SHUTOFF" 
]

require(plotROC)
gg = ggplot(w_dt, aes(d = Positive, m = Score, color=Score_method)) +
    geom_roc(n.cuts = 0,size = 0.5) +
    facet_grid(Lambda+Platform ~ Perturbation+Exprs_noise, labeller = label_both)
a=as.data.table(calc_auc(gg))
a[Exprs_noise=="Medium"][order(-AUC)][1:10]
gg


w_dt[, Test_method := factor(Test_method,
                              levels = c("Nested-wild_whitened", "PSIBaseline"),
                              labels = c("Nested test", "PSI Baseline"))]


w_dt[is.na(Score)]
B <- 1000L
K = 100L

auc_ci_dt <- w_dt[, {

  # point estimates
  pr0  <- pr_auc(Score, Positive)
  roc0 <- roc_auc(Score, Positive)
  top_pr0 = top_pr(Score, Positive, K = K)

  n <- .N

  # bootstrap con replicate
  boot_mat <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    top_pr_b = top_pr(Score[idx], Positive[idx], K = K)
    c(
      pr  = pr_auc(Score[idx], Positive[idx]),
      roc = roc_auc(Score[idx], Positive[idx]),
      recall = top_pr_b$recall,
      precision = top_pr_b$precision
    )
  })

  # boot_mat è 2 x B
  pr_ci  <- quantile(boot_mat["pr", ],  c(0.025, 0.975), na.rm = TRUE)
  roc_ci <- quantile(boot_mat["roc", ], c(0.025, 0.975), na.rm = TRUE)
  recall_ci <- quantile(boot_mat["recall", ], c(0.025, 0.975), na.rm = TRUE)
  precision_ci <- quantile(boot_mat["precision", ], c(0.025, 0.975), na.rm = TRUE)

  .(
    PR_AUC = pr0,
    PR_AUC_low = pr_ci[1],
    PR_AUC_high = pr_ci[2],
    ROC_AUC = roc0,
    ROC_AUC_low = roc_ci[1],
    ROC_AUC_high = roc_ci[2],
    Recall = top_pr0$recall,
    Recall_low = recall_ci[1],
    Recall_high = recall_ci[2],
    Precision = top_pr0$precision,
    Precision_low = precision_ci[1],
    Precision_high = precision_ci[2]
  )

}, by = .(Score_method,N_replicates, Tsteps, Test_method, N_tsamples, Platform, Exprs_noise, Perturbation)]

save(auc_ci_dt,file="auc_ci_dt.rdata")

library(ggplot2)

auc_ci_dt[, N_tsamples_f := factor(N_tsamples, levels = c(3, 5, 10, 20))]

p_pr = ggplot(auc_ci_dt, 
    aes(x = N_tsamples, y = PR_AUC, colour = Perturbation)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  geom_errorbar(
    aes(ymin = PR_AUC_low, ymax = PR_AUC_high),
    width = 0.5,
    linewidth = 0.45
  ) +
  facet_grid(Platform ~ Exprs_noise) +
  scale_x_continuous(breaks = c(3,5,10,20)) +
  labs(
    x = "# time points",
    y = "PR-AUC"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
p_pr
ggsave("pr-auc_plot.pdf", p_pr, width = 16, height = 10, useDingbats = FALSE)



p_roc = ggplot(auc_ci_dt, 
    aes(x = N_tsamples, y = ROC_AUC, colour = Perturbation)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  geom_errorbar(
    aes(ymin = ROC_AUC_low, ymax = ROC_AUC_high),
    width = 0.5,
    linewidth = 0.45
  ) +
  facet_grid(Platform ~ Exprs_noise) +
  scale_x_continuous(breaks = c(3,5,10,20)) +
  labs(
    x = "# time points",
    y = "ROC-AUC"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
p_roc
ggsave("roc-auc_plot.pdf", p_roc, width = 16, height = 10, useDingbats = FALSE)









# ***********************
# FIG4 *** PR curves with ribbon CI and topK points (Nested vs Baseline) ****
# per mostrare confronto con baseline expression based
# fisso su SHUTOFF replicates time points 
# in supp ROC curve
# ***********************

# ---- helpers ----
pr_xy <- function(score, y01) {
  n_y01 = 1*(!y01)
  R = cumsum(y01[order(score, decreasing = TRUE)]) / sum(y01)
  P = cumsum(y01[order(score, decreasing = TRUE)]) / 1:length(score)
  FPR = cumsum(n_y01[order(score, decreasing = TRUE)]) / sum(n_y01)
  list(recall = R, precision = P, fpr = FPR)
}

top_pr <- function(score, y01, K = 100L) {
  n <- length(score)
  if (n == 0L) return(list(precision = NA_real_, recall = NA_real_))
  K <- min(K, n)
  ord <- order(score, decreasing = TRUE)
  topk_idx <- ord[seq_len(K)]
  TP <- sum(y01[topk_idx] == 1L)
  FP <- sum(y01[topk_idx] == 0L)
  FN <- sum(y01[-topk_idx] == 1L)
  P <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  R <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  list(precision = P, recall = R)
}

# filtro scenario
w_dt <- dt_test_results[
  N_tsamples == 10 &
  N_replicates == 5 &
  Score_method=="score_siglogq" &
  Tsteps == 10 & Perturbation == "SHUTOFF"  &
  ((Lambda == 0.5 & Test_method == "Nested-wild_whitened") |
  (Lambda == 0 & Test_method == "PSIBaseline"))
  #Platform %in% c("RNA-seq") &
  #Exprs_noise %in% c("Medium") &
  #Test_method %in% c("Nested-wild_whitened","PSIBaseline")  &
  
]

w_dt[Test_method == "PSIBaseline",Score:=-Score]
w_dt[, Test_method := factor(Test_method,
                              levels = c("Nested-parametric_whitened","Nested-wild_whitened", "PSIBaseline"),
                              labels = c("Nested test param","Nested test wild", "PSI Baseline"))]

B = 1000L
pr_curve_dt <- w_dt[, {
  pr0 <- pr_xy(Score, Positive)

  n = length(Score)
  boot_pr = matrix(NA_real_, nrow = length(pr0$recall), ncol = B)
  lapply(seq_len(B), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    pr_b <- pr_xy(Score[idx], Positive[idx])
    boot_pr[, b] <<- pr_b$precision
  })

  pr_hi = apply(boot_pr, 1, quantile, probs = 0.975, na.rm = TRUE)
  pr_lo = apply(boot_pr, 1, quantile, probs = 0.025, na.rm = TRUE)
  
  data.table(
    recall = pr0$recall,
    fpr = pr0$fpr,
    precision = pr0$precision,
    precision_low = pr_lo,
    precision_high = pr_hi
  )
}, by = .(N_replicates, Tsteps, Test_method, N_tsamples, Platform, Exprs_noise, Perturbation)]


topK_dt <- rbindlist(lapply(c(50L, 100L, 200L), function(K_top) {
  w_dt[, {
    tp <- top_pr(Score, Positive, K = K_top)
    data.table(recall = tp$recall, precision = tp$precision, K = K_top)
  }, by = .(N_replicates, Tsteps, Test_method, N_tsamples, Platform, Exprs_noise, Perturbation)]
}))

topK_dt <- topK_dt[is.finite(recall) & is.finite(precision)]
topK_dt[, K := factor(K, levels = c(50, 100, 200),
                      labels = c("Top-50", "Top-100", "Top-200"))]
          # palette coerente (usa i tuoi livelli/nomi esatti)
meth_cols <- c(
  "Nested test wild" = "#E76F61",
  "Nested test param" = "#2a19bd",
  "PSI Baseline" = "#2CA02C"
)

p_pr_curves_main <- ggplot(
  pr_curve_dt,
  aes(x = recall, y = precision, color = Test_method)
) +
  geom_ribbon(aes(ymin=precision_low, ymax=precision_high,
                  fill=Test_method),
              alpha=0.15, color=NA) +
  geom_line(linewidth=0.5) +
  # punti Top-K: fill=metodo, shape=K
  geom_point(
    data = topK_dt,
    aes(x = recall, y = precision, fill = Test_method, shape = K),
    colour = "grey55",
    size = 2.6,
    stroke = 0.8,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +

  scale_color_manual(values = meth_cols) +
  scale_fill_manual(values = meth_cols) +
  scale_shape_manual(values = c("Top-50" = 21, "Top-100" = 24, "Top-200" = 22)) +

  # qui forziamo la struttura della legenda
  guides(
    color = guide_legend(order = 1, override.aes = list(linetype = 1, linewidth = 1.2)),
    fill  = "none",  # evitiamo una seconda legenda duplicata per fill
    shape = guide_legend(order = 2, override.aes = list(fill = "white", colour = "grey40"))
  ) +

  facet_grid(Platform ~ Exprs_noise) +
  scale_x_continuous(limits = c(0, 1)) +
# scale_x_continuous(
#   trans = scales::pseudo_log_trans(base = 10),
#   limits = c(0, 1),
#   breaks = c(0.001, 0.05, 0.1, 0.2, 0.4),
#   labels = c("0.001", "0.05", "0.1", "0.2", "0.4")
# ) +
  labs(x = "Recall", y = "Precision", color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

p_pr_curves_main
ggsave("pr-curves_with_topK.pdf", p_pr_curves_main, width = 8.5, height = 6, useDingbats = FALSE)

# ROC curves
p_roc_curves_main <- ggplot(
  pr_curve_dt,
  aes(x = fpr, y = recall, color = Test_method)
) +
  geom_line(linewidth=0.5) +
  # punti Top-K: fill=metodo, shape=K

  facet_grid(Platform ~ Exprs_noise) +
  #scale_x_continuous(limits = c(0, 1)) +
# scale_x_continuous(
#   trans = scales::pseudo_log_trans(base = 10),
#   limits = c(0, 1),
#   breaks = c(0.001, 0.05, 0.1, 0.2, 0.4),
#   labels = c("0.001", "0.05", "0.1", "0.2", "0.4")
# ) +
  labs(x = "TPR", y = "FPR", color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

p_roc_curves_main
ggsave("roc-curves.pdf", p_roc_curves_main, width = 8.5, height = 6, useDingbats = FALSE)




