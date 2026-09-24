# =============================================================================
# Supplementary figure: null calibration of full vs cytoplasmic-only model
#
# Purpose
# -------
# Compare empirical Type-I error of:
#   1. the full compartment-resolved model;
#   2. the cytoplasmic-only kinetic baseline.
#
# The figure focuses ONLY on null calibration. Power is intentionally omitted
# because power comparisons are not interpretable when a method is strongly
# anti-conservative.
#
# Input
# -----
# cytoplasmic_only_baseline/cytoplasmic_only_summary.tsv
#
# Output
# ------
# cytoplasmic_only_baseline/FigS_full_vs_cytoplasmic_only_typeI.pdf
# cytoplasmic_only_baseline/FigS_full_vs_cytoplasmic_only_typeI.png
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

INPUT_FILE <- file.path(
  "cytoplasmic_only_baseline",
  "cytoplasmic_only_summary.tsv"
)

OUT_PDF <- file.path(
  "cytoplasmic_only_baseline",
  "FigS_full_vs_cytoplasmic_only_typeI.pdf"
)

OUT_PNG <- file.path(
  "cytoplasmic_only_baseline",
  "FigS_full_vs_cytoplasmic_only_typeI.png"
)


# =============================================================================
# 2. Load and validate
# =============================================================================

dt <- fread(
  INPUT_FILE
)

required <- c(
  "Design",
  "N_replicates",
  "Platform",
  "Exprs_noise",
  "Cyto_TypeI_005",
  "Cyto_Wilson_low",
  "Cyto_Wilson_high",
  "Cyto_CI_contains_005",
  "Full_TypeI_005",
  "Full_Wilson_low",
  "Full_Wilson_high",
  "Full_CI_contains_005"
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
# 3. Ordering and labels
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
# 4. Reshape to long format
# =============================================================================

full_dt <- dt[
  ,
  .(
    Design,
    N_replicates,
    Platform,
    Exprs_noise,
    Method =
      "Full compartment model",
    TypeI =
      Full_TypeI_005,
    Wilson_low =
      Full_Wilson_low,
    Wilson_high =
      Full_Wilson_high,
    CI_contains_005 =
      Full_CI_contains_005
  )
]

cyto_dt <- dt[
  ,
  .(
    Design,
    N_replicates,
    Platform,
    Exprs_noise,
    Method =
      "Cytoplasmic-only model",
    TypeI =
      Cyto_TypeI_005,
    Wilson_low =
      Cyto_Wilson_low,
    Wilson_high =
      Cyto_Wilson_high,
    CI_contains_005 =
      Cyto_CI_contains_005
  )
]

plot_dt <- rbind(
  full_dt,
  cyto_dt
)

plot_dt[
  ,
  Method :=
    factor(
      Method,
      levels = c(
        "Full compartment model",
        "Cytoplasmic-only model"
      )
    )
]


# =============================================================================
# 5. Diagnostics printed to console
# =============================================================================

calibration_counts <- plot_dt[
  ,
  .(
    N_scenarios =
      .N,
    N_calibrated =
      sum(
        CI_contains_005,
        na.rm =
          TRUE
      ),
    Fraction_calibrated =
      mean(
        CI_contains_005,
        na.rm =
          TRUE
      )
  ),
  by =
    Method
]

cat(
  "\nCalibration counts:\n"
)

print(
  calibration_counts
)


# =============================================================================
# 6. Figure
# =============================================================================

NOMINAL_ALPHA <- 0.05
PRACTICAL_LOW <- 0.025
PRACTICAL_HIGH <- 0.075


p <- ggplot(
  plot_dt,
  aes(
    x =
      Exprs_noise,
    y =
      TypeI,
    group =
      Method,
    linetype =
      Method,
    shape =
      Method
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
      0.05
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

  geom_errorbar(
    aes(
      ymin =
        Wilson_low,
      ymax =
        Wilson_high
    ),
    position =
      position_dodge(
        width =
          0.18
      ),
    width =
      0.08,
    linewidth =
      0.45
  ) +

  geom_line(
    linewidth =
      0.8,
    position =
      position_dodge(
        width =
          0.18
      )
  ) +

  geom_point(
    size =
      2.1,
    stroke =
      0.45,
    position =
      position_dodge(
        width =
          0.18
      )
  ) +

  facet_grid(
    Platform ~ Design,
    labeller =
      labeller(
        Design =
          DESIGN_LABELS
      )
  ) +

  scale_linetype_manual(
    values = c(
      "Full compartment model" =
        "solid",
      "Cytoplasmic-only model" =
        "longdash"
    )
  ) +

  scale_shape_manual(
    values = c(
      "Full compartment model" =
        16,
      "Cytoplasmic-only model" =
        17
    )
  ) +

  coord_cartesian(
    ylim =
      c(
        0,
        max(
          plot_dt$Wilson_high,
          na.rm =
            TRUE
        ) *
          1.05
      )
  ) +

  labs(
    x =
      "Measurement-noise regime",
    y =
      "Empirical Type-I error",
    linetype =
      NULL,
    shape =
      NULL
  ) +

  theme_classic(
    base_size =
      10.5
  ) +

  theme(
    panel.grid.major.y =
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

    axis.text.x =
      element_text(
        angle =
          30,
        hjust =
          1
      ),

    legend.position =
      "bottom",

    legend.key.width =
      grid::unit(
        1.5,
        "cm"
      )
  )


# =============================================================================
# 7. Save
# =============================================================================

ggsave(
  filename =
    OUT_PDF,
  plot =
    p,
  width =
    10.5,
  height =
    7.8,
  units =
    "in",
  device =
    cairo_pdf
)

ggsave(
  filename =
    OUT_PNG,
  plot =
    p,
  width =
    10.5,
  height =
    7.8,
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
  OUT_PDF,
  "\n",
  "  ",
  OUT_PNG,
  "\n",
  sep = ""
)
