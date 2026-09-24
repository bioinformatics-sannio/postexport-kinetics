# =============================================================================
# Composite ranking-score ablation on the corrected SHUTOFF benchmark
#
# Reproduces the score definition used in the original benchmark code:
#
#   E_i = min[-log10(q_i), 6]
#   IR_i = max(RSS0_i - RSS1_i, 0) / RSS0_i
#   S_full = sigma_hat_i * IR_i * E_i
#
# and compares it with:
#   significance only       E
#   effect only             sigma_hat
#   fit improvement only    IR
#   significance + effect   sigma_hat * E
#   significance + fit      IR * E
#   effect + fit            sigma_hat * IR
#   complete score          sigma_hat * IR * E
#
# PR-AUC is the primary metric.
#
# Input:
#   benchmark_main_corrected_onset_raw.tsv
#
# Output:
#   score_ablation/
#     score_ablation_raw.tsv
#     score_ablation_summary.tsv
#     score_ablation_overview.tsv
#     FigS_score_ablation.pdf
#     FigS_score_ablation.png
#     score_ablation.rdata
#
# =============================================================================

setwd("~/postexport-kinetics/synthetic_dataset")

library(data.table)
library(ggplot2)
library(parallel)

INPUT_FILE <- "benchmark_main_corrected_onset_raw.tsv"
OUTPUT_DIR <- "score_ablation"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

OUT_RAW <- file.path(OUTPUT_DIR, "score_ablation_raw.tsv")
OUT_SUMMARY <- file.path(OUTPUT_DIR, "score_ablation_summary.tsv")
OUT_OVERVIEW <- file.path(OUTPUT_DIR, "score_ablation_overview.tsv")
OUT_PDF <- file.path(OUTPUT_DIR, "FigS_score_ablation.pdf")
OUT_PNG <- file.path(OUTPUT_DIR, "FigS_score_ablation.png")
OUT_RDATA <- file.path(OUTPUT_DIR, "score_ablation.rdata")

N_BOOT <- 2000L
N_WORKERS <- min(60L, max(1L, parallel::detectCores() - 2L))
SEED <- 20260924L
EPS <- 1e-10
NEGLOGQ_CAP <- 6

# Restrict the ablation to SHUTOFF, the principal inferential regime.
PERTURBATION_KEEP <- "SHUTOFF"

dt <- fread(INPUT_FILE)

# -----------------------------------------------------------------------------
# Resolve column names robustly.
# -----------------------------------------------------------------------------

pick_col <- function(candidates, nm, label) {
  hit <- candidates[candidates %in% nm]
  if (length(hit) == 0L) {
    stop(
      paste0(
        "Cannot find ", label, ". Tried: ",
        paste(candidates, collapse = ", ")
      )
    )
  }
  hit[1]
}

p_col <- pick_col(
  c("p.value", "p_value", "pvalue"),
  names(dt),
  "p-value column"
)

sigma_col <- pick_col(
  c("sigma_c_hat", "Sigma", "sigma_hat", "sigma_c_est"),
  names(dt),
  "estimated sigma_c column"
)

rss0_col <- pick_col(
  c("RSS0", "rss0", "RSS_null", "rss_null"),
  names(dt),
  "null RSS column"
)

rss1_col <- pick_col(
  c("RSS1", "rss1", "RSS_full", "rss_full"),
  names(dt),
  "full RSS column"
)

required <- c(
  "Gene",
  "Positive",
  "Perturbation",
  "Platform",
  "Exprs_noise",
  "N_tsamples",
  "N_replicates",
  "Tsteps",
  p_col,
  sigma_col,
  rss0_col,
  rss1_col
)

missing <- setdiff(required, names(dt))

if (length(missing) > 0L) {
  stop(
    paste0(
      "Missing required columns:\n  ",
      paste(missing, collapse = "\n  ")
    )
  )
}

# Keep successful fits when a status column exists.
if ("status" %in% names(dt)) {
  dt <- dt[status == "ok"]
}

dt <- dt[
  Perturbation == PERTURBATION_KEEP &
    Positive %in% c(0, 1)
]

dt[, p_internal := get(p_col)]
dt[, sigma_internal := pmax(get(sigma_col), 0)]
dt[, rss0_internal := get(rss0_col)]
dt[, rss1_internal := get(rss1_col)]

dt <- dt[
  is.finite(p_internal) &
    is.finite(sigma_internal) &
    is.finite(rss0_internal) &
    is.finite(rss1_internal)
]

# -----------------------------------------------------------------------------
# q-values and exact original score components.
#
# The historical code adjusted p-values within each benchmark stratum.
# Lambda and Test_method are included if they exist, preserving the old code's
# grouping structure without requiring them in the corrected benchmark.
# -----------------------------------------------------------------------------

group_cols <- c(
  "N_replicates",
  "Tsteps",
  "N_tsamples",
  "Platform",
  "Exprs_noise",
  "Perturbation"
)

for (optional_col in c("Lambda", "Test_method")) {
  if (optional_col %in% names(dt)) {
    group_cols <- c(optional_col, group_cols)
  }
}

dt[
  ,
  q_internal := p.adjust(
    p_internal,
    method = "fdr"
  ),
  by = group_cols
]

dt[
  ,
  evidence := pmin(
    -log10(
      pmax(
        q_internal,
        EPS
      )
    ),
    NEGLOGQ_CAP
  )
]

dt[
  ,
  delta_rss := pmax(
    rss0_internal - rss1_internal,
    0
  )
]

dt[
  ,
  IR := fifelse(
    rss0_internal > 0,
    delta_rss / rss0_internal,
    0
  )
]

# Exact original recommended composite:
dt[, score_complete := sigma_internal * IR * evidence]

# Ablations:
dt[, score_p := evidence]
dt[, score_sigma := sigma_internal]
dt[, score_IR := IR]
dt[, score_p_sigma := evidence * sigma_internal]
dt[, score_p_IR := evidence * IR]
dt[, score_sigma_IR := sigma_internal * IR]

score_cols <- c(
  "score_p",
  "score_sigma",
  "score_IR",
  "score_p_sigma",
  "score_p_IR",
  "score_sigma_IR",
  "score_complete"
)

score_labels <- c(
  score_p = "Statistical evidence",
  score_sigma = "Effect size",
  score_IR = "Fit improvement",
  score_p_sigma = "Evidence x effect",
  score_p_IR = "Evidence x fit",
  score_sigma_IR = "Effect x fit",
  score_complete = "Complete score"
)

# -----------------------------------------------------------------------------
# PR-AUC using threshold aggregation; ties are handled at common thresholds.
# -----------------------------------------------------------------------------

pr_auc <- function(truth, score) {
  good <- !is.na(truth) & is.finite(score)
  truth <- as.integer(truth[good])
  score <- score[good]

  n_pos <- sum(truth == 1L)
  if (n_pos < 1L || sum(truth == 0L) < 1L) {
    return(NA_real_)
  }

  thresholds <- sort(unique(score), decreasing = TRUE)

  recall <- numeric(length(thresholds))
  precision <- numeric(length(thresholds))

  for (ii in seq_along(thresholds)) {
    pred <- score >= thresholds[ii]
    tp <- sum(pred & truth == 1L)
    fp <- sum(pred & truth == 0L)

    recall[ii] <- tp / n_pos
    precision[ii] <- if ((tp + fp) > 0L) tp / (tp + fp) else 1
  }

  recall <- c(0, recall)
  precision <- c(1, precision)

  sum(
    diff(recall) *
      (head(precision, -1L) + tail(precision, -1L)) / 2
  )
}

# -----------------------------------------------------------------------------
# Analyze one benchmark configuration with paired stratified bootstrap.
# -----------------------------------------------------------------------------

analyze_config <- function(d, B, seed) {
  null <- d[Positive == 0]
  alt <- d[Positive == 1]

  if (nrow(null) < 2L || nrow(alt) < 2L) {
    stop("Insufficient null or alternative rows.")
  }

  point <- vapply(
    score_cols,
    function(sc) pr_auc(d$Positive, d[[sc]]),
    numeric(1)
  )

  set.seed(seed)

  boot <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(score_cols),
    dimnames = list(NULL, score_cols)
  )

  for (bb in seq_len(B)) {
    idx0 <- sample(seq_len(nrow(null)), nrow(null), replace = TRUE)
    idx1 <- sample(seq_len(nrow(alt)), nrow(alt), replace = TRUE)

    db <- rbind(
      null[idx0],
      alt[idx1]
    )

    boot[bb, ] <- vapply(
      score_cols,
      function(sc) pr_auc(db$Positive, db[[sc]]),
      numeric(1)
    )
  }

  ans <- rbindlist(
    lapply(
      score_cols,
      function(sc) {
        data.table(
          Score = sc,
          AUPR = point[sc],
          AUPR_low = quantile(
            boot[, sc],
            0.025,
            na.rm = TRUE,
            names = FALSE
          ),
          AUPR_high = quantile(
            boot[, sc],
            0.975,
            na.rm = TRUE,
            names = FALSE
          ),
          N_null = nrow(null),
          N_alt = nrow(alt),
          Positive_prevalence =
            nrow(alt) / (nrow(null) + nrow(alt))
        )
      }
    )
  )

  # Paired difference against the complete score.
  complete_boot <- boot[, "score_complete"]

  ans[
    ,
    `:=`(
      Delta_vs_complete =
        AUPR - point["score_complete"],
      Delta_vs_complete_low =
        vapply(
          Score,
          function(sc) {
            quantile(
              boot[, sc] - complete_boot,
              0.025,
              na.rm = TRUE,
              names = FALSE
            )
          },
          numeric(1)
        ),
      Delta_vs_complete_high =
        vapply(
          Score,
          function(sc) {
            quantile(
              boot[, sc] - complete_boot,
              0.975,
              na.rm = TRUE,
              names = FALSE
            )
          },
          numeric(1)
        )
    )
  ]

  ans
}

config_cols <- c(
  "N_tsamples",
  "N_replicates",
  "Tsteps",
  "Platform",
  "Exprs_noise"
)

configs <- unique(dt[, ..config_cols])
setorder(
  configs,
  N_replicates,
  Tsteps,
  N_tsamples,
  Platform,
  Exprs_noise
)

run_one <- function(ii) {
  cfg <- configs[ii]

  d <- dt[
    N_tsamples == cfg$N_tsamples &
      N_replicates == cfg$N_replicates &
      Tsteps == cfg$Tsteps &
      Platform == cfg$Platform &
      Exprs_noise == cfg$Exprs_noise
  ]

  seed_i <- as.integer(
    (
      SEED +
        ii * 104729L
    ) %%
      2000000000L
  )

  cbind(
    cfg,
    analyze_config(
      d,
      B = N_BOOT,
      seed = seed_i
    )
  )
}

cat(
  "\n============================================================\n",
  "COMPOSITE SCORE ABLATION\n",
  "============================================================\n",
  "Configurations: ", nrow(configs), "\n",
  "Bootstrap resamples/configuration: ", N_BOOT, "\n",
  "Workers: ", min(N_WORKERS, nrow(configs)), "\n",
  sep = ""
)

if (.Platform$OS.type == "unix" && nrow(configs) > 1L) {
  result_list <- mclapply(
    seq_len(nrow(configs)),
    run_one,
    mc.cores = min(N_WORKERS, nrow(configs)),
    mc.preschedule = TRUE
  )
} else {
  result_list <- lapply(
    seq_len(nrow(configs)),
    run_one
  )
}

summary_dt <- rbindlist(result_list)

summary_dt[
  ,
  Score_label := factor(
    score_labels[Score],
    levels = unname(score_labels)
  )
]

# -----------------------------------------------------------------------------
# Overview across configurations.
# -----------------------------------------------------------------------------

overview <- summary_dt[
  ,
  .(
    N_configurations = .N,
    Median_AUPR = median(AUPR, na.rm = TRUE),
    Q25_AUPR = quantile(AUPR, 0.25, na.rm = TRUE, names = FALSE),
    Q75_AUPR = quantile(AUPR, 0.75, na.rm = TRUE, names = FALSE),
    Median_delta_vs_complete =
      median(Delta_vs_complete, na.rm = TRUE)
  ),
  by = .(
    Score,
    Score_label
  )
]

setorder(
  overview,
  Score_label
)

# -----------------------------------------------------------------------------
# Save tables.
# -----------------------------------------------------------------------------

raw_out <- dt[
  ,
  c(
    config_cols,
    "Gene",
    "Positive",
    "p_internal",
    "q_internal",
    "sigma_internal",
    "rss0_internal",
    "rss1_internal",
    "IR",
    "evidence",
    score_cols
  ),
  with = FALSE
]

fwrite(
  raw_out,
  OUT_RAW,
  sep = "\t"
)

fwrite(
  summary_dt,
  OUT_SUMMARY,
  sep = "\t"
)

fwrite(
  overview,
  OUT_OVERVIEW,
  sep = "\t"
)

# -----------------------------------------------------------------------------
# Figure: distribution of PR-AUC across corrected SHUTOFF configurations.
# -----------------------------------------------------------------------------

p <- ggplot(
  summary_dt,
  aes(
    x = Score_label,
    y = AUPR
  )
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.13,
    height = 0,
    alpha = 0.25,
    size = 0.9
  ) +
  geom_hline(
    yintercept = unique(summary_dt$Positive_prevalence)[1],
    linetype = "dotted",
    linewidth = 0.5
  ) +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = "PR-AUC"
  ) +
  theme_classic(
    base_size = 10.5
  ) +
  theme(
    panel.grid.major.y = element_line(
      colour = "grey92",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

ggsave(
  OUT_PDF,
  p,
  width = 9.5,
  height = 5.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  OUT_PNG,
  p,
  width = 9.5,
  height = 5.8,
  units = "in",
  dpi = 400,
  bg = "white"
)

score_ablation <- list(
  formula = list(
    evidence = "min(-log10(q), 6)",
    IR = "max(RSS0-RSS1,0)/RSS0",
    complete = "sigma_c_hat * IR * min(-log10(q),6)"
  ),
  summary = summary_dt,
  overview = overview
)

save(
  score_ablation,
  file = OUT_RDATA
)

cat(
  "\nOverview:\n"
)

print(overview)

cat(
  "\nGenerated in:\n  ",
  OUTPUT_DIR,
  "\n",
  sep = ""
)
