# =============================================================================
# FINAL Supplementary Figure S19
# Empirical real-data p-value QQ plots from the FINAL audited result table
#
# Input:
#   final_real_data_audit/final_all_realdata_standardized.tsv
#
# Outputs:
#   final_real_data_audit/FigS19_realdata_QQ_final.pdf
#   final_real_data_audit/FigS19_realdata_QQ_final.png
#
# IMPORTANT:
# This is a DESCRIPTIVE plot only. It is not used to establish calibration.
# =============================================================================

setwd("~/postexport-kinetics/real_datasets")

library(data.table)
library(ggplot2)

INFILE <- "final_real_data_audit/final_all_realdata_standardized.tsv"

OUT_PDF <- "final_real_data_audit/FigS19_realdata_QQ_final.pdf"
OUT_PNG <- "final_real_data_audit/FigS19_realdata_QQ_final.png"

dt <- fread(INFILE)

req <- c("dataset_final", "p_final")
miss <- setdiff(req, names(dt))
if (length(miss) > 0L) {
  stop("Missing columns: ", paste(miss, collapse = ", "))
}

dt <- dt[
  is.finite(p_final) &
    p_final >= 0 &
    p_final <= 1
]

dataset_levels <- c("Kc167", "K562", "NIH-3T3", "mESC")
dt[, dataset_final := factor(dataset_final, levels = dataset_levels)]

make_qq <- function(d) {

  p <- sort(d$p_final)
  n <- length(p)

  # Expected uniform order-statistic probabilities.
  expected_p <- ppoints(n)

  # Pointwise 95% envelope for Uniform(0,1) order statistics.
  i <- seq_len(n)
  low_p  <- qbeta(0.025, shape1 = i, shape2 = n + 1L - i)
  high_p <- qbeta(0.975, shape1 = i, shape2 = n + 1L - i)

  # Convert to -log10(p). Note reversal of lower/upper after transformation.
  data.table(
    expected = -log10(expected_p),
    observed = -log10(pmax(p, 1e-300)),
    band_low = -log10(pmax(high_p, 1e-300)),
    band_high = -log10(pmax(low_p, 1e-300))
  )
}

qq_list <- lapply(
  split(dt, dt$dataset_final, drop = TRUE),
  make_qq
)

qq <- rbindlist(
  lapply(
    names(qq_list),
    function(nm) {
      x <- qq_list[[nm]]
      x[, dataset_final := nm]
      x
    }
  )
)

counts <- dt[
  ,
  .(
    N = .N,
    N_p_lt_005 = sum(p_final < 0.05),
    N_p_lt_001 = sum(p_final < 0.01)
  ),
  by = dataset_final
]

counts[, label := sprintf(
  "n = %d\np < 0.05: %d\np < 0.01: %d",
  N,
  N_p_lt_005,
  N_p_lt_001
)]

qq[, dataset_final := factor(dataset_final, levels = dataset_levels)]
counts[, dataset_final := factor(dataset_final, levels = dataset_levels)]

max_xy <- max(c(
  qq$expected,
  qq$observed,
  qq$band_high
), na.rm = TRUE)

p <- ggplot(
  qq,
  aes(x = expected, y = observed)
) +
  geom_ribbon(
    aes(ymin = band_low, ymax = band_high),
    alpha = 0.14
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.55
  ) +
  geom_point(
    size = 0.9,
    alpha = 0.70
  ) +
  geom_text(
    data = counts,
    aes(
      x = Inf,
      y = -Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1.08,
    vjust = -0.35,
    size = 3.2
  ) +
  facet_wrap(
    ~ dataset_final,
    ncol = 2,
    scales = "free"
  ) +
  labs(
    x = expression(Expected~~-log[10](p)~~under~~Uniform(0,1)),
    y = expression(Observed~~-log[10](p))
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
    panel.grid.minor = element_blank()
  )

ggsave(
  OUT_PDF,
  p,
  width = 7.5,
  height = 6.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  OUT_PNG,
  p,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

cat("\n============================================================\n")
cat("FINAL REAL-DATA QQ FIGURE COMPLETE\n")
cat("============================================================\n\n")
print(counts)
cat("\nGenerated:\n  ", OUT_PDF, "\n  ", OUT_PNG, "\n", sep = "")
