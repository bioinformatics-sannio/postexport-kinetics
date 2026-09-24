# =============================================================================
# REGENERATE DISCRIMINATION FIGURES FROM THE FINAL CORRECTED-ONSET BENCHMARK
#
# This script regenerates:
#   1. ROC-AUC vs number of sampled time points
#   2. PR-AUC vs number of sampled time points
#   3. ROC curves for SHUTOFF at 10 sampled time points
#   4. Precision@Top-K (K = 50,100,200) for SHUTOFF at 10 sampled time points
#
# All quantities are recomputed from the FINAL corrected-onset raw benchmark.
#
# Ranking score:
#   E  = min(-log10(BH q), 6)
#   IR = max(RSS0 - RSS1, 0) / RSS0
#   S  = sigma_c_hat * IR * E
#
# This is the exact exploratory composite score documented in the revised
# manuscript. Formal calibration/power remain based on bootstrap p-values.
#
# Default display design:
#   N_replicates = 5
#   Tsteps       = 10
#
# Change these two constants ONLY if the manuscript figure is intentionally
# defined for another design.
#
# Input:
#   benchmark_main_corrected_onset_raw.tsv
#
# Outputs:
#   benchmark_corrected_analysis/Fig4_AUROC_corrected.pdf/png
#   benchmark_corrected_analysis/FigS17_AUPR_corrected.pdf/png
#   benchmark_corrected_analysis/FigS18_ROC_curves_corrected.pdf/png
#   benchmark_corrected_analysis/FigS_topK_precision_corrected.pdf/png
#   benchmark_corrected_analysis/discrimination_metrics_corrected.tsv
#   benchmark_corrected_analysis/topK_precision_corrected.tsv
# =============================================================================

setwd("~/postexport-kinetics/synthetic_dataset")

library(data.table)
library(ggplot2)
library(parallel)

INFILE <- "benchmark_main_corrected_onset_raw.tsv"
OUTDIR <- "benchmark_corrected_analysis"

if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR, recursive = TRUE)
}

N_REPLICATES_KEEP <- 5L
TSTEPS_KEEP <- 10L
N_BOOT <- 1000L
N_WORKERS <- min(48L, max(1L, detectCores() - 2L))
SEED <- 20260924L
P_FLOOR <- 1e-300

NOISE_LEVELS <- c("Very low", "Low", "Medium", "High")
PLATFORM_LEVELS <- c("RT-qPCR", "GAUSS", "RNA-seq")
PERT_LEVELS <- c("NONE", "SHUTOFF")
TOP_K <- c(50L, 100L, 200L)

dt <- fread(INFILE)

pick_col <- function(candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) > 0L) return(hit[1])
  if (required) {
    stop("Cannot find ", label, ". Tried: ", paste(candidates, collapse = ", "))
  }
  NULL
}

p_col <- pick_col(
  c("p.value", "p_value", "pvalue"),
  "p-value column"
)

sigma_col <- pick_col(
  c("sigma_c_hat", "Sigma", "sigma_hat", "sigma_c_est"),
  "sigma_c estimate"
)

rss0_col <- pick_col(
  c("RSS0", "rss0", "RSS_null", "rss_null"),
  "null RSS"
)

rss1_col <- pick_col(
  c("RSS1", "rss1", "RSS_full", "rss_full"),
  "full RSS"
)

positive_col <- pick_col(
  c("Positive", "positive", "Truth", "truth", "Alternative"),
  "truth/Positive indicator"
)

required_design <- c(
  "Perturbation",
  "Platform",
  "Exprs_noise",
  "N_tsamples",
  "N_replicates",
  "Tsteps"
)

missing_design <- setdiff(required_design, names(dt))
if (length(missing_design) > 0L) {
  stop("Missing design columns: ", paste(missing_design, collapse = ", "))
}

if ("status" %in% names(dt)) {
  dt <- dt[status == "ok"]
}

dt[, truth := as.integer(get(positive_col))]
dt[, p_internal := as.numeric(get(p_col))]
dt[, sigma_internal := pmax(as.numeric(get(sigma_col)), 0)]
dt[, rss0_internal := as.numeric(get(rss0_col))]
dt[, rss1_internal := as.numeric(get(rss1_col))]

dt <- dt[
  truth %in% c(0L, 1L) &
    is.finite(p_internal) &
    is.finite(sigma_internal) &
    is.finite(rss0_internal) &
    is.finite(rss1_internal)
]

config_cols <- c(
  "Perturbation",
  "Platform",
  "Exprs_noise",
  "N_tsamples",
  "N_replicates",
  "Tsteps"
)

# BH correction is performed within each synthetic benchmark configuration.
dt[
  ,
  q_internal := p.adjust(p_internal, method = "BH"),
  by = config_cols
]

dt[, evidence := pmin(-log10(pmax(q_internal, P_FLOOR)), 6)]
dt[
  ,
  IR := fifelse(
    rss0_internal > 0,
    pmax(rss0_internal - rss1_internal, 0) / rss0_internal,
    0
  )
]
dt[, score := sigma_internal * IR * evidence]

# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------

auc_rank <- function(truth, score) {
  good <- !is.na(truth) & is.finite(score)
  truth <- truth[good]
  score <- score[good]

  n1 <- sum(truth == 1L)
  n0 <- sum(truth == 0L)
  if (n1 < 1L || n0 < 1L) return(NA_real_)

  rr <- rank(score, ties.method = "average")
  sr <- sum(rr[truth == 1L])

  (sr - n1 * (n1 + 1) / 2) / (n1 * n0)
}

pr_auc <- function(truth, score) {
  good <- !is.na(truth) & is.finite(score)
  truth <- truth[good]
  score <- score[good]

  npos <- sum(truth == 1L)
  nneg <- sum(truth == 0L)
  if (npos < 1L || nneg < 1L) return(NA_real_)

  th <- sort(unique(score), decreasing = TRUE)

  rec <- numeric(length(th))
  prec <- numeric(length(th))

  for (i in seq_along(th)) {
    pred <- score >= th[i]
    tp <- sum(pred & truth == 1L)
    fp <- sum(pred & truth == 0L)
    rec[i] <- tp / npos
    prec[i] <- if ((tp + fp) > 0L) tp / (tp + fp) else 1
  }

  rec <- c(0, rec)
  prec <- c(1, prec)

  sum(diff(rec) * (head(prec, -1L) + tail(prec, -1L)) / 2)
}

precision_at_k <- function(truth, score, k) {
  good <- !is.na(truth) & is.finite(score)
  truth <- truth[good]
  score <- score[good]

  k <- min(k, length(score))
  ord <- order(score, decreasing = TRUE)
  mean(truth[ord[seq_len(k)]] == 1L)
}

roc_points <- function(truth, score) {
  good <- !is.na(truth) & is.finite(score)
  truth <- truth[good]
  score <- score[good]

  npos <- sum(truth == 1L)
  nneg <- sum(truth == 0L)

  th <- c(Inf, sort(unique(score), decreasing = TRUE), -Inf)

  out <- rbindlist(
    lapply(
      th,
      function(z) {
        pred <- score >= z
        data.table(
          FPR = sum(pred & truth == 0L) / nneg,
          TPR = sum(pred & truth == 1L) / npos
        )
      }
    )
  )

  unique(out, by = c("FPR", "TPR"))[
    order(FPR, TPR)
  ]
}

# -----------------------------------------------------------------------------
# Analyze one configuration with paired stratified bootstrap
# -----------------------------------------------------------------------------

analyze_one <- function(d, seed) {

  null <- d[truth == 0L]
  alt <- d[truth == 1L]

  auc0 <- auc_rank(d$truth, d$score)
  aupr0 <- pr_auc(d$truth, d$score)

  set.seed(seed)

  b_auc <- numeric(N_BOOT)
  b_aupr <- numeric(N_BOOT)
  b_top <- matrix(
    NA_real_,
    nrow = N_BOOT,
    ncol = length(TOP_K),
    dimnames = list(NULL, paste0("K", TOP_K))
  )

  for (b in seq_len(N_BOOT)) {

    i0 <- sample(seq_len(nrow(null)), nrow(null), replace = TRUE)
    i1 <- sample(seq_len(nrow(alt)), nrow(alt), replace = TRUE)

    db <- rbind(
      null[i0],
      alt[i1]
    )

    b_auc[b] <- auc_rank(db$truth, db$score)
    b_aupr[b] <- pr_auc(db$truth, db$score)

    for (j in seq_along(TOP_K)) {
      b_top[b, j] <- precision_at_k(
        db$truth,
        db$score,
        TOP_K[j]
      )
    }
  }

  ans <- data.table(
    AUROC = auc0,
    AUROC_low = quantile(b_auc, 0.025, names = FALSE),
    AUROC_high = quantile(b_auc, 0.975, names = FALSE),
    AUPR = aupr0,
    AUPR_low = quantile(b_aupr, 0.025, names = FALSE),
    AUPR_high = quantile(b_aupr, 0.975, names = FALSE),
    Positive_prevalence = mean(d$truth == 1L)
  )

  top <- rbindlist(
    lapply(
      seq_along(TOP_K),
      function(j) {
        k <- TOP_K[j]
        data.table(
          K = k,
          Precision = precision_at_k(d$truth, d$score, k),
          Precision_low = quantile(b_top[, j], 0.025, names = FALSE),
          Precision_high = quantile(b_top[, j], 0.975, names = FALSE)
        )
      }
    )
  )

  list(metrics = ans, top = top)
}

cfg <- unique(dt[, ..config_cols])

run_cfg <- function(i) {

  c0 <- cfg[i]

  d <- dt[
    Perturbation == c0$Perturbation &
      Platform == c0$Platform &
      Exprs_noise == c0$Exprs_noise &
      N_tsamples == c0$N_tsamples &
      N_replicates == c0$N_replicates &
      Tsteps == c0$Tsteps
  ]

  z <- analyze_one(
    d,
    seed = as.integer((SEED + i * 104729L) %% 2000000000L)
  )

  list(
    metrics = cbind(c0, z$metrics),
    top = cbind(c0, z$top)
  )
}

cat("\nComputing corrected-benchmark discrimination metrics...\n")
cat("Configurations: ", nrow(cfg), "\n", sep = "")
cat("Workers: ", min(N_WORKERS, nrow(cfg)), "\n", sep = "")

if (.Platform$OS.type == "unix") {
  zz <- mclapply(
    seq_len(nrow(cfg)),
    run_cfg,
    mc.cores = min(N_WORKERS, nrow(cfg)),
    mc.preschedule = TRUE
  )
} else {
  zz <- lapply(seq_len(nrow(cfg)), run_cfg)
}

metrics <- rbindlist(lapply(zz, `[[`, "metrics"))
topk <- rbindlist(lapply(zz, `[[`, "top"))

fwrite(
  metrics,
  file.path(OUTDIR, "discrimination_metrics_corrected.tsv"),
  sep = "\t"
)

fwrite(
  topk,
  file.path(OUTDIR, "topK_precision_corrected.tsv"),
  sep = "\t"
)

# -----------------------------------------------------------------------------
# Display subset used for the manuscript trend figures
# -----------------------------------------------------------------------------

mplot <- metrics[
  N_replicates == N_REPLICATES_KEEP &
    Tsteps == TSTEPS_KEEP
]

tplot <- topk[
  N_replicates == N_REPLICATES_KEEP &
    Tsteps == TSTEPS_KEEP &
    Perturbation == "SHUTOFF" &
    N_tsamples == 10
]

mplot[, Platform := factor(Platform, levels = PLATFORM_LEVELS)]
mplot[, Exprs_noise := factor(Exprs_noise, levels = NOISE_LEVELS)]
mplot[, Perturbation := factor(Perturbation, levels = PERT_LEVELS)]

# -----------------------------------------------------------------------------
# 1. AUROC figure
# -----------------------------------------------------------------------------

p_auc <- ggplot(
  mplot,
  aes(
    x = N_tsamples,
    y = AUROC,
    group = Perturbation,
    linetype = Perturbation,
    shape = Perturbation
  )
) +
  geom_errorbar(
    aes(ymin = AUROC_low, ymax = AUROC_high),
    width = 0.25,
    linewidth = 0.35,
    position = position_dodge(width = 0.25)
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted",
    linewidth = 0.45
  ) +
  facet_grid(Platform ~ Exprs_noise) +
  scale_x_continuous(breaks = sort(unique(mplot$N_tsamples))) +
  coord_cartesian(ylim = c(0.45, 1)) +
  labs(
    x = "Number of sampled time points",
    y = "ROC-AUC",
    linetype = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25)
  )

ggsave(
  file.path(OUTDIR, "Fig4_AUROC_corrected.pdf"),
  p_auc,
  width = 11,
  height = 6.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(OUTDIR, "Fig4_AUROC_corrected.png"),
  p_auc,
  width = 11,
  height = 6.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 2. AUPR figure
# -----------------------------------------------------------------------------

p_aupr <- ggplot(
  mplot,
  aes(
    x = N_tsamples,
    y = AUPR,
    group = Perturbation,
    linetype = Perturbation,
    shape = Perturbation
  )
) +
  geom_errorbar(
    aes(ymin = AUPR_low, ymax = AUPR_high),
    width = 0.25,
    linewidth = 0.35,
    position = position_dodge(width = 0.25)
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  geom_hline(
    yintercept = unique(mplot$Positive_prevalence)[1],
    linetype = "dotted",
    linewidth = 0.45
  ) +
  facet_grid(Platform ~ Exprs_noise) +
  scale_x_continuous(breaks = sort(unique(mplot$N_tsamples))) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Number of sampled time points",
    y = "PR-AUC",
    linetype = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25)
  )

ggsave(
  file.path(OUTDIR, "FigS17_AUPR_corrected.pdf"),
  p_aupr,
  width = 11,
  height = 6.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(OUTDIR, "FigS17_AUPR_corrected.png"),
  p_aupr,
  width = 11,
  height = 6.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 3. ROC curves, SHUTOFF, 10 sampled time points
# -----------------------------------------------------------------------------

curve_dt <- dt[
  N_replicates == N_REPLICATES_KEEP &
    Tsteps == TSTEPS_KEEP &
    Perturbation == "SHUTOFF" &
    N_tsamples == 10
]

roc_list <- curve_dt[
  ,
  {
    rr <- roc_points(truth, score)
    rr
  },
  by = .(
    Platform,
    Exprs_noise
  )
]

roc_list[, Platform := factor(Platform, levels = PLATFORM_LEVELS)]
roc_list[, Exprs_noise := factor(Exprs_noise, levels = NOISE_LEVELS)]

p_roc <- ggplot(
  roc_list,
  aes(x = FPR, y = TPR)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dotted",
    linewidth = 0.45
  ) +
  geom_path(linewidth = 0.75) +
  facet_grid(Platform ~ Exprs_noise) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "False-positive rate",
    y = "True-positive rate"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.25)
  )

ggsave(
  file.path(OUTDIR, "FigS18_ROC_curves_corrected.pdf"),
  p_roc,
  width = 9.5,
  height = 7,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(OUTDIR, "FigS18_ROC_curves_corrected.png"),
  p_roc,
  width = 9.5,
  height = 7,
  units = "in",
  dpi = 400,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 4. Top-K precision
# -----------------------------------------------------------------------------

tplot[, Platform := factor(Platform, levels = PLATFORM_LEVELS)]
tplot[, Exprs_noise := factor(Exprs_noise, levels = NOISE_LEVELS)]
tplot[, K_label := factor(
  paste0("Top-", K),
  levels = paste0("Top-", TOP_K)
)]

p_top <- ggplot(
  tplot,
  aes(
    x = K_label,
    y = Precision,
    group = 1
  )
) +
  geom_errorbar(
    aes(
      ymin = Precision_low,
      ymax = Precision_high
    ),
    width = 0.12,
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.9) +
  facet_grid(Platform ~ Exprs_noise) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = NULL,
    y = "Precision among top-ranked events"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25)
  )

ggsave(
  file.path(OUTDIR, "FigS_topK_precision_corrected.pdf"),
  p_top,
  width = 9.5,
  height = 7,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(OUTDIR, "FigS_topK_precision_corrected.png"),
  p_top,
  width = 9.5,
  height = 7,
  units = "in",
  dpi = 400,
  bg = "white"
)

cat("\n============================================================\n")
cat("CORRECTED DISCRIMINATION FIGURES COMPLETE\n")
cat("============================================================\n")
cat("Design shown in figures:\n")
cat("  N_replicates = ", N_REPLICATES_KEEP, "\n", sep = "")
cat("  Tsteps       = ", TSTEPS_KEEP, "\n", sep = "")
cat("\nGenerated in:\n  ", OUTDIR, "\n", sep = "")
