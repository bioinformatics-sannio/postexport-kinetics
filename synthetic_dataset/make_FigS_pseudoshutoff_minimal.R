# =============================================================================
# Minimal supplementary figures for pseudo-shutoff misspecification
#
# Purpose:
#   Generate only the two figures needed to document sensitivity to incomplete
#   transcriptional shutoff:
#
#     1. Empirical Type-I error vs residual transcription
#     2. Mean estimated sigma_c under H0 vs residual transcription
#
# Power is intentionally omitted from the minimal Supplement because the
# principal issue is loss of null calibration, and apparent power gains under
# misspecification are not interpretable independently of Type-I inflation.
#
# Input:
#   pseudoshutoff_benchmark/pseudoshutoff_summary.tsv
#
# Outputs:
#   pseudoshutoff_benchmark/FigS_pseudoshutoff_typeI.pdf
#   pseudoshutoff_benchmark/FigS_pseudoshutoff_typeI.png
#
#   pseudoshutoff_benchmark/FigS_pseudoshutoff_sigma_null.pdf
#   pseudoshutoff_benchmark/FigS_pseudoshutoff_sigma_null.png
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
  "pseudoshutoff_benchmark/pseudoshutoff_summary.tsv"

OUT_TYPEI_PDF <-
  "pseudoshutoff_benchmark/FigS_pseudoshutoff_typeI.pdf"

OUT_TYPEI_PNG <-
  "pseudoshutoff_benchmark/FigS_pseudoshutoff_typeI.png"

OUT_SIGMA_PDF <-
  "pseudoshutoff_benchmark/FigS_pseudoshutoff_sigma_null.pdf"

OUT_SIGMA_PNG <-
  "pseudoshutoff_benchmark/FigS_pseudoshutoff_sigma_null.png"


# =============================================================================
# 2. Settings
# =============================================================================

RHO_BREAKS <- c(
  0,
  5,
  10,
  25,
  50,
  75,
  90,
  100
)

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
  "Design",
  "Platform",
  "Exprs_noise",
  "Residual_transcription_pct",
  "TypeI_005",
  "Null_sigma_hat_mean"
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
# 4. Ordering
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
# 5. Common publication theme
# =============================================================================

theme_pseudo <- theme_classic(
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
# 6. Figure 1: empirical Type-I error
# =============================================================================

p_typeI <- ggplot(
  dt,
  aes(
    x =
      Residual_transcription_pct,
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
      RHO_BREAKS,
    limits =
      c(
        0,
        100
      ),
    expand =
      expansion(
        mult =
          c(
            0.015,
            0.015
          )
      )
  ) +

  labs(
    x =
      "Residual transcription after nominal shutoff (%)",
    y =
      "Empirical Type-I error",
    linetype =
      "Noise",
    shape =
      "Noise"
  ) +

  theme_pseudo


ggsave(
  filename =
    OUT_TYPEI_PDF,
  plot =
    p_typeI,
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
    OUT_TYPEI_PNG,
  plot =
    p_typeI,
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
# 7. Figure 2: apparent sigma_c under H0
# =============================================================================

p_sigma <- ggplot(
  dt,
  aes(
    x =
      Residual_transcription_pct,
    y =
      Null_sigma_hat_mean,
    group =
      Exprs_noise,
    linetype =
      Exprs_noise,
    shape =
      Exprs_noise
  )
) +

  geom_hline(
    yintercept =
      0,
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
      ),
    scales =
      "free_y"
  ) +

  scale_x_continuous(
    breaks =
      RHO_BREAKS,
    limits =
      c(
        0,
        100
      ),
    expand =
      expansion(
        mult =
          c(
            0.015,
            0.015
          )
      )
  ) +

  labs(
    x =
      "Residual transcription after nominal shutoff (%)",
    y =
      expression(
        "Mean " *
          hat(
            sigma
          )[c] *
          " under " *
          H[0]
      ),
    linetype =
      "Noise",
    shape =
      "Noise"
  ) +

  theme_pseudo


ggsave(
  filename =
    OUT_SIGMA_PDF,
  plot =
    p_sigma,
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
    OUT_SIGMA_PNG,
  plot =
    p_sigma,
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
  OUT_TYPEI_PDF,
  "\n",
  "  ",
  OUT_TYPEI_PNG,
  "\n",
  "  ",
  OUT_SIGMA_PDF,
  "\n",
  "  ",
  OUT_SIGMA_PNG,
  "\n",
  sep = ""
)
