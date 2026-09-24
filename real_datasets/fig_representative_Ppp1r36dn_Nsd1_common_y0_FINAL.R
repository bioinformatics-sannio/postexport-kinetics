# =============================================================================
# Title:
# Publication-Ready Representative mESC Dynamics
#
# Representative events:
#
#   Ppp1r36dn
#   Nsd1
#
# Description:
#   Generate a publication-quality main figure from the FINAL
#   20k-bootstrap mESC analysis.
#
# Input:
#
#   mesc_20k_results.rdata
#
# Objects used:
#
#   results_mesc_20k
#   observation-level trajectory table (auto-detected; historically rmats_all)
#
# Important:
#
#   - kinetic coefficients are NOT refitted;
#   - no bootstrap is rerun;
#   - continuous trajectories are ODE reconstructions from the FINAL
#     fitted kinetic coefficients;
#   - full and null trajectories use the SAME observed initial state;
#   - y0 = (N0, N_s0, C0, C_s0) is the replicate mean at the first
#     sampled time point;
#   - no initial-state optimization is performed;
#   - all kinetic parameters remain fixed;
#   - sigma_c is fixed to zero in the null model.
#
# Outputs:
#
#   main_representative_Ppp1r36dn_Nsd1.pdf
#   main_representative_Ppp1r36dn_Nsd1.png
#   main_representative_Ppp1r36dn_Nsd1.tiff
#
#   main_representative_Ppp1r36dn_Nsd1_parameters.tsv
#   main_representative_Ppp1r36dn_Nsd1_y0.tsv
#   main_representative_Ppp1r36dn_Nsd1_data.rdata
#
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/real_datasets")

library(data.table)
library(ggplot2)
library(deSolve)
library(patchwork)
library(grid)

data.table::setDTthreads(1)


# =============================================================================
# 1. Settings
# =============================================================================

DATASET_MAIN <- "mESC | GSE256335"

TARGET_GENES <- c(
  "Ppp1r36dn",
  "Nsd1"
)


# Transcriptional shutoff.
T_STAR <- 0


# Dense ODE grid for smooth trajectories.
ODE_STEP <- 0.20


# -----------------------------------------------------------------------------
# Initial state
# -----------------------------------------------------------------------------
# Full and null models are propagated from the same observed replicate-mean
# state at the first sampled time point.

# -----------------------------------------------------------------------------
# Colours
# -----------------------------------------------------------------------------

COL_FULL <- "#0072B2"

COL_NULL <- "#D55E00"

COL_REPLICATE <- "grey45"

COL_MEAN <- "black"

COL_ERROR <- "grey42"


# -----------------------------------------------------------------------------
# Figure dimensions
#
# Approximately full two-column width.
# -----------------------------------------------------------------------------

FIG_WIDTH_MM <- 183

FIG_HEIGHT_MM <- 118


# =============================================================================
# 2. Load FINAL mESC 20k analysis
# =============================================================================

load(
  "mesc_20k_results.rdata"
)


if (!exists("results_mesc_20k")) {

  stop(
    "Object 'results_mesc_20k' not found in mesc_20k_results.rdata."
  )
}


results_tested <- copy(
  as.data.table(
    results_mesc_20k
  )
)


# -----------------------------------------------------------------------------
# Observed trajectory data
# -----------------------------------------------------------------------------
# mesc_20k_results.rdata contains the final test results but, in the current
# analysis export, it does not contain the observation-level object rmats_all.
# Recover that object automatically from a separate RData/RDS/TSV/CSV file in
# the working directory.  No object is loaded into .GlobalEnv while searching.

required_observed_columns <- c(
  "dataset",
  "event",
  "time",
  "replicate",
  "N",
  "N_s",
  "C",
  "C_s"
)


has_observed_columns <- function(x) {

  is.data.frame(x) &&
    all(
      required_observed_columns %in% names(x)
    )
}


find_observed_data <- function() {

  # If the object already exists in the calling environment, use it.
  if (
    exists(
      "rmats_all",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  ) {

    x <- get(
      "rmats_all",
      envir = .GlobalEnv
    )

    if (has_observed_columns(x)) {

      cat(
        "Using existing object 'rmats_all' from the R session.\n"
      )

      return(
        copy(
          as.data.table(x)
        )
      )
    }
  }


  # Search serialized R files.  Files with informative names are inspected
  # first.  The final result file itself is excluded because it is known not to
  # contain the observation-level trajectories.
  r_files <- list.files(
    pattern = "\\.(rdata|RData|rds)$",
    full.names = TRUE
  )

  r_files <- r_files[
    basename(r_files) != "mesc_20k_results.rdata"
  ]

  if (length(r_files) > 0L) {

    priority <- grepl(
      "rmats|mesc|real|observ|trajectory|input|data",
      basename(r_files),
      ignore.case = TRUE
    )

    r_files <- c(
      r_files[priority],
      r_files[!priority]
    )

    r_files <- unique(r_files)
  }


  for (ff in r_files) {

    # Avoid accidentally opening extremely large unrelated workspaces.
    size_bytes <- file.info(ff)$size

    if (
      is.finite(size_bytes) &&
        size_bytes > 2 * 1024^3
    ) {
      next
    }

    ext <- tolower(
      tools::file_ext(ff)
    )

    if (ext == "rds") {

      obj <- tryCatch(
        readRDS(ff),
        error = function(e) NULL
      )

      if (has_observed_columns(obj)) {

        cat(
          "Observed trajectories recovered from ",
          ff,
          " (RDS).\n",
          sep = ""
        )

        return(
          copy(
            as.data.table(obj)
          )
        )
      }

    } else {

      ee <- new.env(
        parent = emptyenv()
      )

      object_names <- tryCatch(
        load(
          ff,
          envir = ee
        ),
        error = function(e) character(0)
      )

      if (length(object_names) == 0L) {
        next
      }

      # Prefer an object literally called rmats_all when present.
      object_names <- c(
        intersect(
          "rmats_all",
          object_names
        ),
        setdiff(
          object_names,
          "rmats_all"
        )
      )

      for (obj_name in object_names) {

        obj <- get(
          obj_name,
          envir = ee
        )

        if (has_observed_columns(obj)) {

          cat(
            "Observed trajectories recovered from ",
            ff,
            " :: ",
            obj_name,
            ".\n",
            sep = ""
          )

          return(
            copy(
              as.data.table(obj)
            )
          )
        }
      }
    }
  }


  # Search delimited tables by reading headers first.
  text_files <- list.files(
    pattern = "\\.(tsv|txt|csv)$",
    full.names = TRUE
  )

  if (length(text_files) > 0L) {

    priority <- grepl(
      "rmats|mesc|observ|trajectory|state|data",
      basename(text_files),
      ignore.case = TRUE
    )

    text_files <- unique(
      c(
        text_files[priority],
        text_files[!priority]
      )
    )
  }


  for (ff in text_files) {

    header <- tryCatch(
      fread(
        ff,
        nrows = 0L,
        showProgress = FALSE
      ),
      error = function(e) NULL
    )

    if (
      !is.null(header) &&
        all(
          required_observed_columns %in% names(header)
        )
    ) {

      cat(
        "Observed trajectories recovered from ",
        ff,
        ".\n",
        sep = ""
      )

      return(
        fread(
          ff,
          showProgress = TRUE
        )
      )
    }
  }


  stop(
    paste0(
      "Could not locate the observation-level trajectory table.\n",
      "The final file mesc_20k_results.rdata contains only: ",
      "results_mesc_20k, FDR05, FDR10, and mesc_summary.\n",
      "Place in ~/postexport-kinetics/real_datasets either an RData/RDS ",
      "object or a TSV/CSV table containing these columns:\n  ",
      paste(
        required_observed_columns,
        collapse = ", "
      ),
      "\nThe script will detect it automatically."
    )
  )
}


rmats_all <- find_observed_data()


cat(
  "\n============================================================\n",
  "FINAL mESC 20k ANALYSIS LOADED\n",
  "============================================================\n",
  sep = ""
)


cat(
  "Tested events: ",
  nrow(results_tested),
  "\n",
  sep = ""
)


cat(
  "Observed rows: ",
  nrow(rmats_all),
  "\n",
  sep = ""
)


# =============================================================================
# 3. Required columns
# =============================================================================

required_result_columns <- c(

  "dataset",
  "event",
  "gene_symbol",

  "p.value",
  "q.value",

  "Sigma",
  "IR",

  # Full model.
  "R_hat",
  "tau_hat",
  "tau_s_hat",
  "sigma_c_hat",
  "sigma_n_hat",
  "alpha_hat",
  "alpha_s_hat",

  # Null model.
  "R_hat_null",
  "tau_hat_null",
  "tau_s_hat_null",
  "sigma_n_hat_null",
  "alpha_hat_null",
  "alpha_s_hat_null"
)


missing_result_columns <- setdiff(
  required_result_columns,
  names(results_tested)
)


if (length(missing_result_columns) > 0L) {

  stop(
    paste0(
      "Missing result columns:\n  ",
      paste(
        missing_result_columns,
        collapse = "\n  "
      )
    )
  )
}


required_data_columns <- c(

  "dataset",
  "event",
  "time",
  "replicate",

  "N",
  "N_s",
  "C",
  "C_s"
)


missing_data_columns <- setdiff(
  required_data_columns,
  names(rmats_all)
)


if (length(missing_data_columns) > 0L) {

  stop(
    paste0(
      "Missing rmats_all columns:\n  ",
      paste(
        missing_data_columns,
        collapse = "\n  "
      )
    )
  )
}


# =============================================================================
# 4. Select representative events
# =============================================================================

selected_events <- results_tested[
  dataset == DATASET_MAIN &
    gene_symbol %in% TARGET_GENES &
    is.finite(q.value) &
    q.value < 0.05 &
    is.finite(Sigma) &
    Sigma > 0 &
    is.finite(IR)
]


if (
  !all(
    TARGET_GENES %in%
      selected_events$gene_symbol
  )
) {

  stop(
    "Not all requested target genes contain an FDR-supported event."
  )
}


# If more than one FDR-supported event is available for a gene,
# retain the one with the largest IR.

setorder(
  selected_events,
  gene_symbol,
  -IR,
  q.value
)


selected_events <- selected_events[
  ,
  .SD[1],
  by = gene_symbol
]


selected_events[
  ,
  gene_order :=
    match(
      gene_symbol,
      TARGET_GENES
    )
]


setorder(
  selected_events,
  gene_order
)


cat(
  "\nSelected representative events:\n"
)


print(
  selected_events[
    ,
    .(
      gene_symbol,
      event,
      q.value,
      Sigma,
      half_time_min =
        log(2) / Sigma,
      IR
    )
  ]
)


# =============================================================================
# 5. Four-state ODE model
# =============================================================================

ode_rhs <- function(
  t,
  y,
  p
) {

  N <-
    y[["N"]]

  Ns <-
    y[["N_s"]]

  C <-
    y[["C"]]

  Cs <-
    y[["C_s"]]


  # ---------------------------------------------------------------------------
  # Transcriptional shutoff
  # ---------------------------------------------------------------------------

  R <- if (

    !is.null(p$t_star) &&

      is.finite(p$t_star) &&

      t >= p$t_star

  ) {

    0

  } else {

    p$R
  }


  # ---------------------------------------------------------------------------
  # ODE system
  # ---------------------------------------------------------------------------

  dN <-
    R -
    p$sigma_n * N -
    p$tau * N


  dNs <-
    p$sigma_n * N -
    p$tau_s * Ns


  dC <-
    p$tau * N -
    p$sigma_c * C -
    p$alpha * C


  dCs <-
    p$tau_s * Ns +
    p$sigma_c * C -
    p$alpha_s * Cs


  list(
    c(
      dN,
      dNs,
      dC,
      dCs
    )
  )
}


# =============================================================================
# 6. ODE simulation
# =============================================================================

simulate_kinetic_model <- function(
  coefficients,
  y0,
  times,
  t_star = 0
) {

  required <- c(
    "R",
    "tau",
    "tau_s",
    "sigma_c",
    "sigma_n",
    "alpha",
    "alpha_s"
  )


  missing <- setdiff(
    required,
    names(coefficients)
  )


  if (length(missing) > 0L) {

    stop(
      paste(
        "Missing kinetic parameters:",
        paste(
          missing,
          collapse = ", "
        )
      )
    )
  }


  if (
    any(
      !is.finite(
        coefficients[
          required
        ]
      )
    )
  ) {

    stop(
      "Non-finite kinetic coefficient."
    )
  }


  if (
    any(!is.finite(y0)) ||
      any(y0 < 0)
  ) {

    stop(
      "Invalid initial state."
    )
  }


  p <- as.list(
    coefficients
  )


  p$t_star <-
    t_star


  out <- deSolve::ode(

    y =
      y0,

    times =
      times,

    func =
      ode_rhs,

    parms =
      p,

    method =
      "lsoda",

    rtol =
      1e-9,

    atol =
      1e-11
  )


  out <- as.data.table(
    out
  )


  setnames(
    out,
    c(
      "time",
      "N",
      "N_s",
      "C",
      "C_s"
    )
  )


  out
}


# =============================================================================
# 7. Extract final fitted coefficients
# =============================================================================

extract_full_coefficients <- function(
  result_row
) {

  c(

    R =
      result_row$R_hat,

    tau =
      result_row$tau_hat,

    tau_s =
      result_row$tau_s_hat,

    sigma_c =
      result_row$sigma_c_hat,

    sigma_n =
      result_row$sigma_n_hat,

    alpha =
      result_row$alpha_hat,

    alpha_s =
      result_row$alpha_s_hat
  )
}


extract_null_coefficients <- function(
  result_row
) {

  c(

    R =
      result_row$R_hat_null,

    tau =
      result_row$tau_hat_null,

    tau_s =
      result_row$tau_s_hat_null,

    # H0.
    sigma_c =
      0,

    sigma_n =
      result_row$sigma_n_hat_null,

    alpha =
      result_row$alpha_hat_null,

    alpha_s =
      result_row$alpha_s_hat_null
  )
}


# =============================================================================
# 8. Observed summaries
# =============================================================================

make_observed_summary <- function(
  ts_gene
) {

  states <- c(
    "N",
    "N_s",
    "C",
    "C_s"
  )


  long <- melt(

    copy(
      ts_gene
    ),

    id.vars =
      c(
        "time",
        "replicate"
      ),

    measure.vars =
      states,

    variable.name =
      "State",

    value.name =
      "Value"
  )


  summary <- long[
    ,
    .(

      Mean =
        mean(
          Value,
          na.rm = TRUE
        ),

      SD =
        sd(
          Value,
          na.rm = TRUE
        ),

      Nrep =
        sum(
          is.finite(Value)
        )
    ),

    by =
      .(
        time,
        State
      )
  ]


  list(

    long =
      long,

    summary =
      summary
  )
}


# =============================================================================
# 9. Starting point for y0 optimization
# =============================================================================

get_first_observed_y0 <- function(
  ts_gene
) {

  t0 <- min(
    ts_gene$time,
    na.rm = TRUE
  )


  yy <- ts_gene[
    time == t0,
    .(

      N =
        mean(
          N,
          na.rm = TRUE
        ),

      N_s =
        mean(
          N_s,
          na.rm = TRUE
        ),

      C =
        mean(
          C,
          na.rm = TRUE
        ),

      C_s =
        mean(
          C_s,
          na.rm = TRUE
        )
    )
  ]


  y0 <- c(

    N =
      yy$N[1],

    N_s =
      yy$N_s[1],

    C =
      yy$C[1],

    C_s =
      yy$C_s[1]
  )


  y0[
    !is.finite(y0) |
      y0 < 0
  ] <- 0


  y0
}


# =============================================================================
# 10. Human-readable genomic coordinates
# =============================================================================


format_event_label <- function(
  event_string
) {

  z <- strsplit(
    event_string,
    ":",
    fixed = TRUE
  )[[1]]


  if (length(z) == 7L) {

    start_position <- format(
      as.numeric(
        z[4]
      ),
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    )


    end_position <- format(
      as.numeric(
        z[7]
      ),
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    )


    return(
      paste0(
        z[2],
        ":",
        start_position,
        "\u2013",
        end_position,
        " (",
        z[3],
        ")"
      )
    )
  }


  event_string
}


# =============================================================================
# 11. Build one publication panel
# =============================================================================

make_publication_panel <- function(
  result_row
) {

  gene_i <- as.character(
    result_row$gene_symbol
  )


  event_i <- as.character(
    result_row$event
  )


  # ---------------------------------------------------------------------------
  # Retrieve observations
  # ---------------------------------------------------------------------------

  ts_gene <- rmats_all[
    dataset == DATASET_MAIN &
      event == event_i
  ]


  if (nrow(ts_gene) == 0L) {

    stop(
      paste(
        "No observations found for",
        gene_i
      )
    )
  }


  ts_gene <- copy(
    ts_gene
  )


  setorder(
    ts_gene,
    time,
    replicate
  )


  # ---------------------------------------------------------------------------
  # Final fitted kinetic parameters
  # ---------------------------------------------------------------------------

  coef_full <- extract_full_coefficients(
    result_row
  )


  coef_null <- extract_null_coefficients(
    result_row
  )


  # ---------------------------------------------------------------------------
  # Common observed initial state
  # ---------------------------------------------------------------------------

  # The two fitted models are propagated from exactly the same initial state,
  # defined as the replicate mean at the first observed time point.  This keeps
  # the visual comparison aligned with the inferential setup and avoids a
  # post-hoc model-specific optimization of y0.

  y0_common <- get_first_observed_y0(
    ts_gene
  )


  # ---------------------------------------------------------------------------
  # Time grids
  # ---------------------------------------------------------------------------

  observed_times <- sort(
    unique(
      ts_gene$time
    )
  )


  dense_times <- seq(

    min(
      observed_times
    ),

    max(
      observed_times
    ),

    by =
      ODE_STEP
  )


  dense_times <- sort(
    unique(
      c(
        dense_times,
        observed_times
      )
    )
  )


  # ---------------------------------------------------------------------------
  # ODE reconstruction: full
  # ---------------------------------------------------------------------------

  prediction_full <- simulate_kinetic_model(

    coefficients =
      coef_full,

    y0 =
      y0_common,

    times =
      dense_times,

    t_star =
      T_STAR
  )


  prediction_full[
    ,
    Model :=
      "Full"
  ]


  # ---------------------------------------------------------------------------
  # ODE reconstruction: null
  # ---------------------------------------------------------------------------

  prediction_null <- simulate_kinetic_model(

    coefficients =
      coef_null,

    y0 =
      y0_common,

    times =
      dense_times,

    t_star =
      T_STAR
  )


  prediction_null[
    ,
    Model :=
      "Null"
  ]


  predictions <- rbindlist(
    list(
      prediction_full,
      prediction_null
    )
  )


  # ---------------------------------------------------------------------------
  # Observed data
  # ---------------------------------------------------------------------------

  observed <- make_observed_summary(
    ts_gene
  )


  observation_long <-
    observed$long


  observation_summary <-
    observed$summary


  prediction_long <- melt(

    predictions,

    id.vars =
      c(
        "time",
        "Model"
      ),

    measure.vars =
      c(
        "N",
        "N_s",
        "C",
        "C_s"
      ),

    variable.name =
      "State",

    value.name =
      "Value"
  )


  # ---------------------------------------------------------------------------
  # Mathematical state labels
  # ---------------------------------------------------------------------------

  state_levels <- c(
    "N",
    "N_s",
    "C",
    "C_s"
  )


  state_labels <- c(
    "N",
    "N[s]",
    "C",
    "C[s]"
  )


  observation_long[
    ,
    State_plot :=
      factor(

        State,

        levels =
          state_levels,

        labels =
          state_labels
      )
  ]


  observation_summary[
    ,
    State_plot :=
      factor(

        State,

        levels =
          state_levels,

        labels =
          state_labels
      )
  ]


  prediction_long[
    ,
    State_plot :=
      factor(

        State,

        levels =
          state_levels,

        labels =
          state_labels
      )
  ]


  # ---------------------------------------------------------------------------
  # Statistics
  # ---------------------------------------------------------------------------

  half_time <-
    log(2) /
    result_row$Sigma


  coordinate_text <- format_event_label(
    event_i
  )


  q_text <- formatC(
    result_row$q.value,
    format = "g",
    digits = 3
  )


  sigma_text <- formatC(
    result_row$Sigma,
    format = "g",
    digits = 3
  )


  half_text <- sprintf(
    "%.1f",
    half_time
  )


  IR_text <- sprintf(
    "%.2f",
    result_row$IR
  )


  # ---------------------------------------------------------------------------
  # Metadata strings
  # ---------------------------------------------------------------------------

  metrics_line_1 <- paste0(
    "qBH = ", q_text,
    "    sigma_c = ", sigma_text, " min^-1"
  )

  metrics_line_2 <- paste0(
    "t1/2 = ", half_text, " min",
    "    IR = ", IR_text
  )

  # ---------------------------------------------------------------------------
  # SD error-bar width
  # ---------------------------------------------------------------------------

  error_width <- if (
    length(observed_times) >= 2L
  ) {

    0.08 *
      min(
        diff(
          observed_times
        )
      )

  } else {

    1
  }


  # =============================================================================
  # 12. Dynamics plot
  # =============================================================================

  dynamics_plot <- ggplot() +

    # -------------------------------------------------------------------------
    # Model trajectories
    # -------------------------------------------------------------------------

    geom_line(

      data =
        prediction_long,

      aes(
        x = time,
        y = Value,
        colour = Model,
        linetype = Model
      ),

      linewidth =
        0.90,

      lineend =
        "round"
    ) +

    # -------------------------------------------------------------------------
    # Standard deviation
    # -------------------------------------------------------------------------

    geom_errorbar(

      data =
        observation_summary,

      aes(
        x = time,

        ymin =
          pmax(
            Mean - SD,
            0
          ),

        ymax =
          Mean + SD
      ),

      width =
        error_width,

      linewidth =
        0.32,

      colour =
        COL_ERROR
    ) +

    # -------------------------------------------------------------------------
    # Individual biological replicates
    # -------------------------------------------------------------------------

    geom_point(

      data =
        observation_long,

      aes(
        x = time,
        y = Value
      ),

      size =
        1.25,

      colour =
        COL_REPLICATE,

      alpha =
        0.65
    ) +

    # -------------------------------------------------------------------------
    # Replicate mean
    # -------------------------------------------------------------------------

    geom_point(

      data =
        observation_summary,

      aes(
        x = time,
        y = Mean
      ),

      shape =
        21,

      size =
        2.70,

      stroke =
        0.50,

      colour =
        COL_MEAN,

      fill =
        "white"
    ) +

    # -------------------------------------------------------------------------
    # Four states
    # -------------------------------------------------------------------------

    facet_wrap(

      ~State_plot,

      ncol =
        2,

      scales =
        "free_y",

      labeller =
        label_parsed
    ) +

    # -------------------------------------------------------------------------
    # Model colours
    # -------------------------------------------------------------------------

    scale_colour_manual(

      values =
        c(
          Full =
            COL_FULL,

          Null =
            COL_NULL
        ),

      breaks =
        c(
          "Full",
          "Null"
        ),

      labels =
        expression(

          Full~(sigma[c] > 0),

          Null~(sigma[c] == 0)
        )
    ) +

    # -------------------------------------------------------------------------
    # Model line types
    # -------------------------------------------------------------------------

    scale_linetype_manual(

      values =
        c(
          Full =
            "solid",

          Null =
            "22"
        ),

      breaks =
        c(
          "Full",
          "Null"
        ),

      labels =
        expression(

          Full~(sigma[c] > 0),

          Null~(sigma[c] == 0)
        )
    ) +

    # -------------------------------------------------------------------------
    # X axis
    # -------------------------------------------------------------------------

    scale_x_continuous(

      breaks =
        observed_times,

      expand =
        expansion(
          mult =
            c(
              0.015,
              0.025
            )
        )
    ) +

    # -------------------------------------------------------------------------
    # Y axis
    # -------------------------------------------------------------------------

    scale_y_continuous(

      expand =
        expansion(
          mult =
            c(
              0.025,
              0.065
            )
        )
    ) +

    labs(

      x =
        "Time after transcriptional shutoff (min)",

      y =
        NULL,

      colour =
        NULL,

      linetype =
        NULL
    ) +

    # -------------------------------------------------------------------------
    # Publication styling
    # -------------------------------------------------------------------------

    theme_classic(
      base_size = 8.5
    ) +

    theme(

      # -----------------------------------------------------------------------
      # Very light HORIZONTAL grid only.
      #
      # No vertical grid to retain a cleaner publication appearance.
      # -----------------------------------------------------------------------

      panel.grid.major.y =
        element_line(
          colour = "grey92",
          linewidth = 0.25
        ),

      panel.grid.major.x =
        element_blank(),

      panel.grid.minor =
        element_blank(),


      # -----------------------------------------------------------------------
      # State names
      # -----------------------------------------------------------------------

      strip.background =
        element_blank(),

      strip.text =
        element_text(
          size = 10.7,
          face = "plain",
          colour = "black",
          margin =
            margin(
              b = 4
            )
        ),

      panel.spacing =
        unit(
          0.85,
          "lines"
        ),

      # -----------------------------------------------------------------------
      # Axes
      # -----------------------------------------------------------------------

      axis.line =
        element_line(
          linewidth = 0.38,
          colour = "black"
        ),

      axis.ticks =
        element_line(
          linewidth = 0.32,
          colour = "black"
        ),

      axis.ticks.length =
        unit(
          1.4,
          "mm"
        ),

      axis.text =
        element_text(
          size = 7.7,
          colour = "black"
        ),

      axis.title.x =
        element_text(
          size = 8.7,
          colour = "black",
          margin =
            margin(
              t = 5
            )
        ),

      # -----------------------------------------------------------------------
      # Legend
      # -----------------------------------------------------------------------

      legend.position =
        "bottom",

      legend.direction =
        "horizontal",

      legend.text =
        element_text(
          size = 8.2
        ),

      legend.key.width =
        unit(
          11,
          "mm"
        ),

      legend.key.height =
        unit(
          3,
          "mm"
        ),

      plot.margin =
        margin(
          t = 1,
          r = 5,
          b = 1,
          l = 4
        )
    )


  # =============================================================================
  # 13. Gene + genomic-coordinate header
  # =============================================================================

  # ---------------------------------------------------------------------------
  # Header with four explicitly separated text rows.
  #
  # Using a taller dedicated header prevents the genomic coordinate and
  # statistics from colliding with the gene name or the faceted dynamics.
  # ---------------------------------------------------------------------------

  header_plot <- ggplot() +

    annotate(
      "text",
      x = 0.00, y = 0.94,
      label = gene_i,
      hjust = 0, vjust = 1,
      fontface = "bold",
      size = 4.45,
      colour = "black"
    ) +

    annotate(
      "text",
      x = 0.00, y = 0.69,
      label = coordinate_text,
      hjust = 0, vjust = 1,
      size = 2.45,
      colour = "grey30"
    ) +

    annotate(
      "text",
      x = 0.00, y = 0.43,
      label = metrics_line_1,
      hjust = 0, vjust = 1,
      size = 2.35,
      colour = "grey10"
    ) +

    annotate(
      "text",
      x = 0.00, y = 0.18,
      label = metrics_line_2,
      hjust = 0, vjust = 1,
      size = 2.35,
      colour = "grey10"
    ) +

    xlim(0, 1) +
    ylim(0, 1) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(t = 1.5, r = 8, b = 1.5, l = 4)
    )


  # =============================================================================
  # 14. Combine header + dynamics
  # =============================================================================

  final_panel <- header_plot /

    dynamics_plot +

    plot_layout(

      heights =
        c(
          0.205,
          0.795
        )
    )


  list(

    plot =
      final_panel,

    data =
      ts_gene,

    prediction_full =
      prediction_full,

    prediction_null =
      prediction_null,

    coefficients_full =
      coef_full,

    coefficients_null =
      coef_null,

    y0_common =
      y0_common
  )
}


# =============================================================================
# 15. Generate representative panels
# =============================================================================

result_ppp <- selected_events[
  gene_symbol == "Ppp1r36dn"
]


result_nsd1 <- selected_events[
  gene_symbol == "Nsd1"
]


panel_ppp <- make_publication_panel(
  result_ppp
)


panel_nsd1 <- make_publication_panel(
  result_nsd1
)


# =============================================================================
# 16. Combine main figure
# =============================================================================

main_figure <- (

  panel_ppp$plot |

  panel_nsd1$plot

) +

  plot_layout(

    widths =
      c(
        1,
        1
      ),

    guides =
      "collect"
  ) &

  theme(

    legend.position =
      "bottom"
  )


print(
  main_figure
)


# =============================================================================
# 17. Export vector PDF
# =============================================================================

ggsave(

  filename =
    "main_representative_Ppp1r36dn_Nsd1.pdf",

  plot =
    main_figure,

  width =
    FIG_WIDTH_MM,

  height =
    FIG_HEIGHT_MM,

  units =
    "mm",

  device =
    cairo_pdf,

  bg =
    "white"
)


# =============================================================================
# 18. Export high-resolution PNG
# =============================================================================

ggsave(

  filename =
    "main_representative_Ppp1r36dn_Nsd1.png",

  plot =
    main_figure,

  width =
    FIG_WIDTH_MM,

  height =
    FIG_HEIGHT_MM,

  units =
    "mm",

  dpi =
    600,

  bg =
    "white"
)


# =============================================================================
# 19. Export high-resolution TIFF
# =============================================================================

ggsave(

  filename =
    "main_representative_Ppp1r36dn_Nsd1.tiff",

  plot =
    main_figure,

  width =
    FIG_WIDTH_MM,

  height =
    FIG_HEIGHT_MM,

  units =
    "mm",

  dpi =
    600,

  compression =
    "lzw",

  bg =
    "white"
)


# =============================================================================
# 20. Export exact kinetic parameters
# =============================================================================

parameter_table <- rbindlist(

  list(

    data.table(

      gene =
        "Ppp1r36dn",

      model =
        "Full",

      parameter =
        names(
          panel_ppp$coefficients_full
        ),

      estimate =
        as.numeric(
          panel_ppp$coefficients_full
        )
    ),

    data.table(

      gene =
        "Ppp1r36dn",

      model =
        "Null",

      parameter =
        names(
          panel_ppp$coefficients_null
        ),

      estimate =
        as.numeric(
          panel_ppp$coefficients_null
        )
    ),

    data.table(

      gene =
        "Nsd1",

      model =
        "Full",

      parameter =
        names(
          panel_nsd1$coefficients_full
        ),

      estimate =
        as.numeric(
          panel_nsd1$coefficients_full
        )
    ),

    data.table(

      gene =
        "Nsd1",

      model =
        "Null",

      parameter =
        names(
          panel_nsd1$coefficients_null
        ),

      estimate =
        as.numeric(
          panel_nsd1$coefficients_null
        )
    )
  )
)


fwrite(

  parameter_table,

  "main_representative_Ppp1r36dn_Nsd1_parameters.tsv",

  sep = "\t"
)


# =============================================================================
# 21. Export common observed initial states
# =============================================================================

y0_table <- rbindlist(
  list(
    data.table(
      gene = "Ppp1r36dn",
      state = names(panel_ppp$y0_common),
      common_observed_y0 = as.numeric(panel_ppp$y0_common)
    ),
    data.table(
      gene = "Nsd1",
      state = names(panel_nsd1$y0_common),
      common_observed_y0 = as.numeric(panel_nsd1$y0_common)
    )
  )
)


fwrite(
  y0_table,
  "main_representative_Ppp1r36dn_Nsd1_y0.tsv",
  sep = "\t"
)


# =============================================================================
# 22. Save exact figure objects/data
# =============================================================================

main_representative_figure <- list(

  source =
    "mesc_20k_results.rdata",

  results_object =
    "results_mesc_20k",

  selected_events =
    selected_events,

  Ppp1r36dn =
    panel_ppp,

  Nsd1 =
    panel_nsd1,

  figure =
    main_figure
)


save(

  main_representative_figure,

  file =
    "main_representative_Ppp1r36dn_Nsd1_data.rdata"
)


# =============================================================================
# 23. Final report
# =============================================================================

cat(
  "\n============================================================\n",
  "PUBLICATION-READY REPRESENTATIVE FIGURE COMPLETE\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nSelected events:\n"
)


print(
  selected_events[
    ,
    .(
      gene_symbol,
      event,
      q.value,
      Sigma,
      half_time_min =
        log(2) / Sigma,
      IR
    )
  ]
)


cat(
  "\nCommon observed initial states used for BOTH full and null models:\n"
)

print(y0_table)


cat(
  "\nGenerated:\n",

  "  main_representative_Ppp1r36dn_Nsd1.pdf\n",

  "  main_representative_Ppp1r36dn_Nsd1.png\n",

  "  main_representative_Ppp1r36dn_Nsd1.tiff\n",

  "  main_representative_Ppp1r36dn_Nsd1_parameters.tsv\n",

  "  main_representative_Ppp1r36dn_Nsd1_y0.tsv\n",

  "  main_representative_Ppp1r36dn_Nsd1_data.rdata\n",

  sep = ""
)


cat(
  "\nFigure conventions:\n",

  "  grey dots    = biological replicates\n",

  "  open circles = replicate means\n",

  "  error bars   = +/- 1 SD\n",

  "  blue solid   = full model, sigma_c > 0\n",

  "  orange dash  = null model, sigma_c = 0\n",

  "  light grid   = horizontal major grid only\n",

  "\n",

  "Continuous curves are ODE reconstructions using FINAL 20k kinetic\n",

  "estimates. Full and null trajectories use the SAME replicate-mean\n",

  "initial state at the first sampled time; kinetic parameters are not refitted.\n",

  "============================================================\n",

  sep = ""
)