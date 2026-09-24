# =============================================================================
# Title:
# Two-Gene ODE Example: NONE vs SHUTOFF vs PSEUDO-SHUTOFF
#
# Purpose:
#   Publication-style illustration of two genes with different transcriptional
#   onset times but a common pharmacological intervention time.
#
#   Conditions:
#     - NONE
#     - SHUTOFF (R_post = 0)
#     - PSEUDO-SHUTOFF (R_post = POST_R_FRACTION * R_pre)
#
# Figure requirements:
#   - same y-axis range across the three conditions for each gene;
#   - gene-specific kinetic parameters shown at left;
#   - only two items in the bottom key:
#       sampling time
#       common shutoff
#   - no global figure title/subtitle;
#   - common t_star for both genes;
#   - explicit validation that pseudo-shutoff retains residual R.
#
# Requires:
#   corrected ode.r containing simulate_scheduled_trajectory()
#
# Outputs:
#   ode_two_genes_three_conditions.pdf
#   ode_two_genes_three_conditions.png
#   ode_two_genes_three_conditions_metadata.tsv
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/ode_model")

source("ode.r")

library(data.table)
library(ggplot2)
library(patchwork)
library(grid)


# =============================================================================
# 1. Global settings
# =============================================================================

set.seed(20260911)

TIME_MIN <- 0
TIME_MAX <- 300

TIMES <- seq(
  TIME_MIN,
  TIME_MAX,
  by = 1
)

T_STAR <- 120

SAMPLE_TIMES <- c(
  60,
  90,
  120,
  180,
  240
)

# Make pseudo-shutoff visually distinguishable while remaining a partial shutoff.
POST_R_FRACTION <- 0.20


# =============================================================================
# 2. Two illustrative genes
#
# Both genes have the same experimental shutoff time.
# Their transcriptional onset times differ.
# =============================================================================

GENES <- data.table(
  gene = c(
    "Gene 1",
    "Gene 2"
  ),

  onset_time = c(
    20,
    60
  ),

  R = c(
    85,
    70
  ),

  tau = c(
    0.020,
    0.030
  ),

  tau_s = c(
    0.070,
    0.100
  ),

  sigma_n = c(
    0.040,
    0.030
  ),

  sigma_c = c(
    0.030,
    0.000
  ),

  alpha = c(
    0.040,
    0.020
  ),

  alpha_s = c(
    0.020,
    0.040
  )
)


# =============================================================================
# 3. Initial state
# =============================================================================

Y0 <- c(
  N = 0,
  N_s = 0,
  C = 0,
  C_s = 0
)


# =============================================================================
# 4. Conditions
# =============================================================================

CONDITIONS <- data.table(
  condition = c(
    "NONE",
    "SHUTOFF",
    "PSEUDO-SHUTOFF"
  ),

  post_R_fraction = c(
    1,
    0,
    POST_R_FRACTION
  )
)


# =============================================================================
# 5. Helpers
# =============================================================================

make_parameter_list <- function(
  gene_row
) {

  list(
    R = gene_row$R,
    tau = gene_row$tau,
    tau_s = gene_row$tau_s,
    sigma_n = gene_row$sigma_n,
    sigma_c = gene_row$sigma_c,
    alpha = gene_row$alpha,
    alpha_s = gene_row$alpha_s
  )
}


simulate_one <- function(
  gene_row,
  condition_row
) {

  p <- make_parameter_list(
    gene_row
  )

  condition_i <- as.character(
    condition_row$condition
  )

  t_star_i <- if (
    condition_i == "NONE"
  ) {
    NULL
  } else {
    T_STAR
  }

  post_fraction_i <- condition_row$post_R_fraction

  sim <- simulate_scheduled_trajectory(
    y0 = Y0,
    times = TIMES,
    params = p,
    onset_time = gene_row$onset_time,
    t_star = t_star_i,
    post_R_fraction = post_fraction_i,
    model_kinetics = rna_kinetics
  )

  sim <- as.data.table(
    sim
  )

  sim[
    ,
    `:=`(
      gene = gene_row$gene,
      condition = condition_i,
      onset_time = gene_row$onset_time,
      R_pre = gene_row$R,
      post_R_fraction = post_fraction_i
    )
  ]

  sim
}


# =============================================================================
# 6. Simulate all 2 x 3 combinations
# =============================================================================

sim_list <- list()

kk <- 1L

for (
  gi in seq_len(
    nrow(
      GENES
    )
  )
) {

  for (
    ci in seq_len(
      nrow(
        CONDITIONS
      )
    )
  ) {

    sim_list[[kk]] <- simulate_one(
      GENES[gi],
      CONDITIONS[ci]
    )

    kk <- kk + 1L
  }
}

sim_all <- rbindlist(
  sim_list,
  use.names = TRUE,
  fill = TRUE
)


# =============================================================================
# 7. Hard validation of R
# =============================================================================

# Complete shutoff must have R=0 from t_star onward.
bad_shutoff <- sim_all[
  condition == "SHUTOFF" &
    time >= T_STAR &
    abs(R) > 1e-12
]

if (
  nrow(
    bad_shutoff
  ) > 0L
) {
  stop(
    "SHUTOFF validation failed: R is non-zero at/after T_STAR."
  )
}


# Pseudo-shutoff must retain the requested fraction of gene-specific R.
pseudo_validation <- merge(
  sim_all[
    condition == "PSEUDO-SHUTOFF" &
      time >= T_STAR,
    .(
      gene,
      time,
      R_observed = R
    )
  ],
  GENES[
    ,
    .(
      gene,
      R_pre = R
    )
  ],
  by = "gene"
)

pseudo_validation[
  ,
  R_expected :=
    POST_R_FRACTION *
    R_pre
]

if (
  any(
    abs(
      pseudo_validation$R_observed -
        pseudo_validation$R_expected
    ) >
      1e-12
  )
) {
  stop(
    "PSEUDO-SHUTOFF validation failed: residual R is incorrect."
  )
}


# =============================================================================
# 8. Long format
# =============================================================================

STATE_LEVELS <- c(
  "N",
  "N_s",
  "C",
  "C_s"
)

STATE_COLORS <- c(
  N = "#00A6A6",
  N_s = "#8F46D3",
  C = "#72A800",
  C_s = "#F36C6C"
)

STATE_LABELS <- c(
  N = "N",
  N_s = "N\u209b",
  C = "C",
  C_s = "C\u209b"
)

sim_long <- melt(
  sim_all,
  id.vars = c(
    "time",
    "R",
    "gene",
    "condition",
    "onset_time",
    "R_pre",
    "post_R_fraction"
  ),
  measure.vars = STATE_LEVELS,
  variable.name = "State",
  value.name = "Value"
)

sim_long[
  ,
  State :=
    factor(
      State,
      levels = STATE_LEVELS
    )
]

sample_points <- sim_long[
  time %in% SAMPLE_TIMES
]


# =============================================================================
# 9. Shared y limits PER GENE
#
# This is the important alignment requested:
# all three conditions of a given gene use exactly the same y range.
# =============================================================================

gene_y_limits <- sim_long[
  ,
  .(
    ymin = 0,
    ymax =
      1.05 *
      max(
        Value,
        na.rm = TRUE
      )
  ),
  by = gene
]


# =============================================================================
# 10. Parameter labels
# =============================================================================

parameter_text <- function(
  gene_i
) {

  x <- GENES[
    gene == gene_i
  ][1]

  sprintf(
    paste0(
      "onset = %.0f min\n",
      "R = %.1f\n",
      "\u03c4 = %.3f min\u207b\u00b9\n",
      "\u03c4\u209b = %.3f min\u207b\u00b9\n",
      "\u03c3\u2099 = %.3f min\u207b\u00b9\n",
      "\u03c3c = %.3f min\u207b\u00b9\n",
      "\u03b1 = %.3f min\u207b\u00b9\n",
      "\u03b1\u209b = %.3f min\u207b\u00b9"
    ),
    x$onset_time,
    x$R,
    x$tau,
    x$tau_s,
    x$sigma_n,
    x$sigma_c,
    x$alpha,
    x$alpha_s
  )
}


# =============================================================================
# 11. One trajectory panel
# =============================================================================

make_panel <- function(
  gene_i,
  condition_i,
  show_y_title = FALSE,
  show_x_title = FALSE
) {

  d <- sim_long[
    gene == gene_i &
      condition == condition_i
  ]

  dp <- sample_points[
    gene == gene_i &
      condition == condition_i
  ]

  yy <- gene_y_limits[
    gene == gene_i
  ]

  g <- ggplot(
    d,
    aes(
      x = time,
      y = Value,
      colour = State
    )
  ) +

    # Sampling times.
    geom_vline(
      xintercept = SAMPLE_TIMES,
      linewidth = 0.35,
      linetype = "dashed",
      colour = "grey67"
    ) +

    # Common shutoff time only in intervention panels.
    {
      if (
        condition_i != "NONE"
      ) {
        geom_vline(
          xintercept = T_STAR,
          linewidth = 0.85,
          linetype = "solid",
          colour = "#B2182B"
        )
      }
    } +

    geom_line(
      linewidth = 1.10,
      lineend = "round"
    ) +

    geom_point(
      data = dp,
      aes(
        x = time,
        y = Value,
        colour = State
      ),
      shape = 21,
      fill = "white",
      size = 2.2,
      stroke = 0.45
    ) +

    scale_colour_manual(
      values = STATE_COLORS,
      breaks = STATE_LEVELS,
      labels = STATE_LABELS
    ) +

    scale_x_continuous(
      breaks = c(
        0,
        60,
        120,
        180,
        240,
        300
      ),
      limits = c(
        TIME_MIN,
        TIME_MAX
      ),
      expand = expansion(
        mult = c(
          0.01,
          0.02
        )
      )
    ) +

    coord_cartesian(
      ylim = c(
        yy$ymin,
        yy$ymax
      )
    ) +

    labs(
      x = if (
        show_x_title
      ) {
        "Time (min)"
      } else {
        NULL
      },

      y = if (
        show_y_title
      ) {
        "Relative abundance"
      } else {
        NULL
      },

      colour = NULL
    ) +

    theme_classic(
      base_size = 10.5
    ) +

    theme(
      legend.position = "none",

      axis.title.x = element_text(
        size = 10,
        margin = margin(
          t = 5
        )
      ),

      axis.title.y = element_text(
        size = 10,
        margin = margin(
          r = 5
        )
      ),

      axis.text.x = element_text(
        size = 8.5
      ),

      axis.text.y = element_text(
        size = 8.5
      ),

      panel.grid.major.x = element_line(
        colour = "grey94",
        linewidth = 0.25
      ),

      panel.grid.major.y = element_line(
        colour = "grey94",
        linewidth = 0.25
      ),

      panel.grid.minor = element_blank(),

      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )

  g
}


# =============================================================================
# 12. Parameter blocks at left
# =============================================================================

make_gene_label <- function(
  gene_i
) {

  ggplot() +

    annotate(
      "text",
      x = 0,
      y = 0.93,
      label = gene_i,
      hjust = 0,
      vjust = 1,
      size = 5,
      fontface = "bold"
    ) +

    annotate(
      "text",
      x = 0,
      y = 0.82,
      label =
        parameter_text(
          gene_i
        ),
      hjust = 0,
      vjust = 1,
      size = 3.2,
      lineheight = 1.18
    ) +

    xlim(
      0,
      1
    ) +

    ylim(
      0,
      1
    ) +

    theme_void()
}


# =============================================================================
# 13. Column headers
# =============================================================================

make_header <- function(
  txt
) {

  ggplot() +

    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = txt,
      size = 5,
      fontface = "bold"
    ) +

    xlim(
      0,
      1
    ) +

    ylim(
      0,
      1
    ) +

    theme_void()
}


header_none <- make_header(
  "No perturbation"
)

header_shutoff <- make_header(
  "Shutoff"
)

header_pseudo <- make_header(
  paste0(
    "Pseudo-shutoff (",
    round(
      100 *
        POST_R_FRACTION
    ),
    "% residual R)"
  )
)


# =============================================================================
# 14. Six trajectory panels
# =============================================================================

p11 <- make_panel(
  "Gene 1",
  "NONE",
  show_y_title = FALSE,
  show_x_title = FALSE
)

p12 <- make_panel(
  "Gene 1",
  "SHUTOFF",
  show_y_title = FALSE,
  show_x_title = FALSE
)

p13 <- make_panel(
  "Gene 1",
  "PSEUDO-SHUTOFF",
  show_y_title = FALSE,
  show_x_title = FALSE
)

p21 <- make_panel(
  "Gene 2",
  "NONE",
  show_y_title = FALSE,
  show_x_title = TRUE
)

p22 <- make_panel(
  "Gene 2",
  "SHUTOFF",
  show_y_title = FALSE,
  show_x_title = TRUE
)

p23 <- make_panel(
  "Gene 2",
  "PSEUDO-SHUTOFF",
  show_y_title = FALSE,
  show_x_title = TRUE
)


# =============================================================================
# 15. Bottom key
#
# ONLY:
#   sampling time
#   common shutoff
# =============================================================================

bottom_key <- ggplot() +

  annotate(
    "segment",
    x = 0.16,
    xend = 0.23,
    y = 0.5,
    yend = 0.5,
    colour = "grey67",
    linewidth = 0.55,
    linetype = "dashed"
  ) +

  annotate(
    "text",
    x = 0.24,
    y = 0.5,
    label = "sampling time",
    hjust = 0,
    size = 3.6
  ) +

  annotate(
    "segment",
    x = 0.55,
    xend = 0.62,
    y = 0.5,
    yend = 0.5,
    colour = "#B2182B",
    linewidth = 0.90
  ) +

  annotate(
    "text",
    x = 0.63,
    y = 0.5,
    label = paste0(
      "common shutoff (t* = ",
      T_STAR,
      " min)"
    ),
    hjust = 0,
    size = 3.6
  ) +

  xlim(
    0,
    1
  ) +

  ylim(
    0,
    1
  ) +

  theme_void()


# =============================================================================
# 16. Assemble figure
#
# No global title.
# No global subtitle.
# =============================================================================

header_row <- (
  plot_spacer() |
    header_none |
    header_shutoff |
    header_pseudo
) +
  plot_layout(
    widths = c(
      0.28,
      1,
      1,
      1
    )
  )

gene1_row <- (
  make_gene_label(
    "Gene 1"
  ) |
    p11 |
    p12 |
    p13
) +
  plot_layout(
    widths = c(
      0.28,
      1,
      1,
      1
    )
  )

gene2_row <- (
  make_gene_label(
    "Gene 2"
  ) |
    p21 |
    p22 |
    p23
) +
  plot_layout(
    widths = c(
      0.28,
      1,
      1,
      1
    )
  )

final_plot <- (
  header_row /
    gene1_row /
    gene2_row /
    bottom_key
) +
  plot_layout(
    heights = c(
      0.09,
      1,
      1,
      0.11
    )
  )


print(
  final_plot
)


# =============================================================================
# 17. Metadata
# =============================================================================

metadata <- merge(
  CJ(
    gene =
      GENES$gene,
    condition =
      CONDITIONS$condition,
    unique =
      TRUE
  ),
  GENES,
  by =
    "gene"
)

metadata <- merge(
  metadata,
  CONDITIONS,
  by =
    "condition"
)

metadata[
  ,
  `:=`(
    common_T_star =
      T_STAR,

    sampling_times =
      paste(
        SAMPLE_TIMES,
        collapse = ","
      )
  )
]

setcolorder(
  metadata,
  c(
    "gene",
    "condition",
    "onset_time",
    "common_T_star",
    "post_R_fraction",
    "sampling_times",
    "R",
    "tau",
    "tau_s",
    "sigma_n",
    "sigma_c",
    "alpha",
    "alpha_s"
  )
)

fwrite(
  metadata,
  "ode_two_genes_three_conditions_metadata.tsv",
  sep = "\t"
)


# =============================================================================
# 18. Console validation
# =============================================================================

cat(
  "\n============================================================\n",
  "TWO-GENE THREE-CONDITION FIGURE\n",
  "============================================================\n",
  sep = ""
)

cat(
  "\nPseudo-shutoff residual transcription:\n",
  "  requested fraction = ",
  POST_R_FRACTION,
  "\n",
  sep = ""
)

print(
  unique(
    pseudo_validation[
      ,
      .(
        gene,
        R_pre,
        R_expected,
        R_observed
      )
    ]
  )
)

cat(
  "\nR around common t_star:\n"
)

print(
  sim_all[
    time %in%
      c(
        T_STAR - 1,
        T_STAR,
        T_STAR + 1
      ),
    .(
      gene,
      condition,
      time,
      R
    )
  ][
    order(
      gene,
      condition,
      time
    )
  ]
)

cat(
  "\nShared y-axis limits by gene:\n"
)

print(
  gene_y_limits
)


# =============================================================================
# 19. Export
# =============================================================================

ggsave(
  filename =
    "ode_two_genes_three_conditions.pdf",
  plot =
    final_plot,
  width =
    16,
  height =
    9.3,
  units =
    "in",
  device =
    cairo_pdf
)

ggsave(
  filename =
    "ode_two_genes_three_conditions.png",
  plot =
    final_plot,
  width =
    16,
  height =
    9.3,
  units =
    "in",
  dpi =
    400,
  bg =
    "white"
)

cat(
  "\nSaved:\n",
  "  ode_two_genes_three_conditions.pdf\n",
  "  ode_two_genes_three_conditions.png\n",
  "  ode_two_genes_three_conditions_metadata.tsv\n",
  sep = ""
)
