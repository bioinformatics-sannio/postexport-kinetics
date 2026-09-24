# =============================================================================
# Supplementary figures for fraction-separation misspecification
#
# Purpose:
#   Generate ONLY the two figures needed for the Supplement:
#
#     1. Cytoplasmic/nuclear relative-scaling sensitivity
#     2. Cross-fraction contamination sensitivity
#
# Each figure is a single plot (no multi-panel assembly), showing the primary
# inferential endpoint: empirical Type-I error under sigma_c = 0.
#
# Input:
#   fraction_robustness_benchmark/fraction_robustness_summary.tsv
#
# Outputs:
#   fraction_robustness_benchmark/FigS_scaling_typeI.pdf
#   fraction_robustness_benchmark/FigS_scaling_typeI.png
#
#   fraction_robustness_benchmark/FigS_contamination_typeI.pdf
#   fraction_robustness_benchmark/FigS_contamination_typeI.png
#
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/synthetic_dataset")

library(data.table)
library(ggplot2)


# =============================================================================
# 1. Files
# =============================================================================

INPUT_FILE <-
  "fraction_robustness_benchmark/fraction_robustness_summary.tsv"

OUT_SCALING_PDF <-
  "fraction_robustness_benchmark/FigS_scaling_typeI.pdf"

OUT_SCALING_PNG <-
  "fraction_robustness_benchmark/FigS_scaling_typeI.png"

OUT_CONTAM_PDF <-
  "fraction_robustness_benchmark/FigS_contamination_typeI.pdf"

OUT_CONTAM_PNG <-
  "fraction_robustness_benchmark/FigS_contamination_typeI.png"


# =============================================================================
# 2. Settings
# =============================================================================

NOISE_LEVELS <- c(
  "Very low",
  "Low",
  "Medium",
  "High"
)

PLATFORM_LEVELS <- c(
  "RT-qPCR",
  "GAUSS",
  "RNA-seq"
)

DESIGN_LEVELS <- c(
  "5tp_3rep_dt10",
  "5tp_10rep_dt10"
)

DESIGN_LABELS <- c(
  "5tp_3rep_dt10" =
    "5 time points, 3 replicates",
  "5tp_10rep_dt10" =
    "5 time points, 10 replicates"
)

SCALING_LEVELS <- c(
  0.50,
  0.75,
  1.00,
  1.25,
  1.50,
  2.00
)

CONTAMINATION_PCT <- c(
  0,
  1,
  2.5,
  5,
  10,
  20
)

NOMINAL_ALPHA <- 0.05
PRACTICAL_LOW <- 0.025
PRACTICAL_HIGH <- 0.075


# =============================================================================
# 3. Load and validate
# =============================================================================

dt <- fread(
  INPUT_FILE
)

required <- c(
  "Analysis",
  "Design",
  "Platform",
  "Exprs_noise",
  "TypeI_005",
  "Cytoplasmic_scale",
  "Contamination_fraction"
)

missing <- setdiff(
  required,
  names(
    dt
  )
)

if (
  length(
    missing
  ) > 0L
) {
  stop(
    paste0(
      "Missing required columns:\n  ",
      paste(
        missing,
        collapse = "\n  "
      )
    )
  )
}


# =============================================================================
# 4. Ordered factors
# =============================================================================

dt[
  ,
  Exprs_noise :=
    factor(
      Exprs_noise,
      levels =
        NOISE_LEVELS
    )
]

dt[
  ,
  Platform :=
    factor(
      Platform,
      levels =
        PLATFORM_LEVELS
    )
]

dt[
  ,
  Design :=
    factor(
      Design,
      levels =
        DESIGN_LEVELS
    )
]


# =============================================================================
# 5. Publication theme
# =============================================================================

theme_robustness <- theme_classic(
  base_size =
    10.5
) +
  theme(

    panel.grid.major =
      element_line(
        colour =
          "grey92",
        linewidth =
          0.25
      ),

    panel.grid.minor =
      element_blank(),

    strip.background =
      element_blank(),

    strip.text =
      element_text(
        face =
          "bold",
        size =
          9.5
      ),

    axis.title =
      element_text(
        size =
          10.5
      ),

    axis.text =
      element_text(
        size =
          9
      ),

    legend.position =
      "bottom",

    legend.title =
      element_text(
        size =
          9.5
      ),

    legend.text =
      element_text(
        size =
          9
      ),

    plot.margin =
      margin(
        5,
        7,
        5,
        5
      )
  )


# =============================================================================
# 6. Relative-scaling figure
# =============================================================================

scaling_dt <- dt[
  Analysis ==
    "SCALING"
]


p_scaling <- ggplot(
  scaling_dt,
  aes(
    x =
      Cytoplasmic_scale,
    y =
      TypeI_005,
    group =
      Exprs_noise,
    linetype =
      Exprs_noise,
    shape =
      Exprs_noise
  )
) +

  annotate(
    "rect",
    xmin =
      -Inf,
    xmax =
      Inf,
    ymin =
      PRACTICAL_LOW,
    ymax =
      PRACTICAL_HIGH,
    alpha =
      0.045
  ) +

  geom_hline(
    yintercept =
      NOMINAL_ALPHA,
    linetype =
      "dashed",
    linewidth =
      0.5,
    colour =
      "grey35"
  ) +

  geom_vline(
    xintercept =
      1,
    linetype =
      "dotted",
    linewidth =
      0.45,
    colour =
      "grey50"
  ) +

  geom_line(
    linewidth =
      0.78
  ) +

  geom_point(
    size =
      1.9,
    stroke =
      0.4
  ) +

  facet_grid(
    Platform ~ Design,
    labeller =
      labeller(
        Design =
          DESIGN_LABELS
      )
  ) +

  scale_x_continuous(
    breaks =
      SCALING_LEVELS
  ) +

  labs(
    x =
      "Cytoplasmic relative scaling factor",
    y =
      "Empirical Type-I error",
    linetype =
      "Noise",
    shape =
      "Noise"
  ) +

  theme_robustness


ggsave(
  filename =
    OUT_SCALING_PDF,
  plot =
    p_scaling,
  width =
    10.5,
  height =
    7.6,
  units =
    "in",
  device =
    cairo_pdf
)

ggsave(
  filename =
    OUT_SCALING_PNG,
  plot =
    p_scaling,
  width =
    10.5,
  height =
    7.6,
  units =
    "in",
  dpi =
    400,
  bg =
    "white"
)


# =============================================================================
# 7. Cross-contamination figure
# =============================================================================

contam_dt <- dt[
  Analysis ==
    "CONTAMINATION"
]

contam_dt[
  ,
  Contamination_pct :=
    100 *
      Contamination_fraction
]


p_contam <- ggplot(
  contam_dt,
  aes(
    x =
      Contamination_pct,
    y =
      TypeI_005,
    group =
      Exprs_noise,
    linetype =
      Exprs_noise,
    shape =
      Exprs_noise
  )
) +

  annotate(
    "rect",
    xmin =
      -Inf,
    xmax =
      Inf,
    ymin =
      PRACTICAL_LOW,
    ymax =
      PRACTICAL_HIGH,
    alpha =
      0.045
  ) +

  geom_hline(
    yintercept =
      NOMINAL_ALPHA,
    linetype =
      "dashed",
    linewidth =
      0.5,
    colour =
      "grey35"
  ) +

  geom_vline(
    xintercept =
      0,
    linetype =
      "dotted",
    linewidth =
      0.45,
    colour =
      "grey50"
  ) +

  geom_line(
    linewidth =
      0.78
  ) +

  geom_point(
    size =
      1.9,
    stroke =
      0.4
  ) +

  facet_grid(
    Platform ~ Design,
    labeller =
      labeller(
        Design =
          DESIGN_LABELS
      )
  ) +

  scale_x_continuous(
    breaks =
      CONTAMINATION_PCT
  ) +

  labs(
    x =
      "Cross-fraction contamination (%)",
    y =
      "Empirical Type-I error",
    linetype =
      "Noise",
    shape =
      "Noise"
  ) +

  theme_robustness


ggsave(
  filename =
    OUT_CONTAM_PDF,
  plot =
    p_contam,
  width =
    10.5,
  height =
    7.6,
  units =
    "in",
  device =
    cairo_pdf
)

ggsave(
  filename =
    OUT_CONTAM_PNG,
  plot =
    p_contam,
  width =
    10.5,
  height =
    7.6,
  units =
    "in",
  dpi =
    400,
  bg =
    "white"
)


cat(
  "\nGenerated:\n",
  "  ",
  OUT_SCALING_PDF,
  "\n",
  "  ",
  OUT_SCALING_PNG,
  "\n",
  "  ",
  OUT_CONTAM_PDF,
  "\n",
  "  ",
  OUT_CONTAM_PNG,
  "\n",
  sep = ""
)
