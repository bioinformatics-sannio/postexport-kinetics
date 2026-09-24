# =============================================================================
# FINAL REAL-DATA AUDIT
#
# Purpose
# -------
# Build a single source of truth for the final revision using:
#
#   1. results_realdatasets_revision.rdata
#        -> Kc167, K562, NIH-3T3
#
#   2. mesc_20k_results.rdata
#        -> final mESC 20k analysis
#
# The script:
#   - identifies the result tables automatically;
#   - standardizes dataset labels;
#   - recomputes BH q-values if needed;
#   - reports:
#         N tested
#         N p < 0.05
#         N q < 0.10
#         N q < 0.05
#         sigma_c boundary fraction
#         number of unique genes
#   - exports final mESC FDR05 / FDR10 event tables;
#   - extracts Ppp1r36dn and Nsd1 from the same final mESC source;
#   - saves all audit objects in one RData file.
#
# IMPORTANT
# ---------
# This script DOES NOT modify any source result file.
# It only reads, standardizes, audits, and exports final tables.
#
# Run from:
#   ~/postexport-kinetics/real_datasets
#
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/real_datasets")

library(data.table)


# =============================================================================
# 1. Input files
# =============================================================================

OLD_REALDATA_FILE <- "results_realdatasets_revision.rdata"
MESC_FILE <- "mesc_20k_results.rdata"

OUT_DIR <- "final_real_data_audit"

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}


# =============================================================================
# 2. Output files
# =============================================================================

OUT_COUNTS <- file.path(
  OUT_DIR,
  "final_dataset_counts.tsv"
)

OUT_BOUNDARY <- file.path(
  OUT_DIR,
  "final_boundary_summary.tsv"
)

OUT_MESC_FDR05 <- file.path(
  OUT_DIR,
  "final_mesc_FDR05_events.tsv"
)

OUT_MESC_FDR10 <- file.path(
  OUT_DIR,
  "final_mesc_FDR10_events.tsv"
)

OUT_REP <- file.path(
  OUT_DIR,
  "final_representative_events.tsv"
)

OUT_ALL_STANDARDIZED <- file.path(
  OUT_DIR,
  "final_all_realdata_standardized.tsv"
)

OUT_RDATA <- file.path(
  OUT_DIR,
  "final_real_data_audit.rdata"
)


# =============================================================================
# 3. Helpers
# =============================================================================

pick_col <- function(
  candidates,
  nm,
  required = TRUE,
  label = NULL
) {

  hit <- candidates[
    candidates %in% nm
  ]

  if (length(hit) > 0L) {
    return(hit[1])
  }

  if (required) {
    stop(
      paste0(
        "Could not find ",
        ifelse(
          is.null(label),
          "required column",
          label
        ),
        ". Tried: ",
        paste(
          candidates,
          collapse = ", "
        )
      )
    )
  }

  NULL
}


safe_unique_n <- function(x) {

  x <- x[
    !is.na(x) &
      nzchar(
        as.character(x)
      )
  ]

  uniqueN(x)
}


normalize_dataset_label <- function(x) {

  y <- as.character(x)
  z <- tolower(y)

  out <- y

  out[
    grepl(
      "kc167|gse83620|drosophila",
      z
    )
  ] <- "Kc167"

  out[
    grepl(
      "k562",
      z
    )
  ] <- "K562"

  out[
    grepl(
      "3t3|nih-?3t3|nih3t3",
      z
    )
  ] <- "NIH-3T3"

  out[
    grepl(
      "mesc|gse256335|esc",
      z
    )
  ] <- "mESC"

  out
}


extract_table_from_env <- function(
  env,
  preferred = NULL,
  must_have_any = c(
    "p.value",
    "p_value",
    "pvalue"
  )
) {

  objs <- ls(env)

  if (!is.null(preferred)) {
    for (nm in preferred) {
      if (nm %in% objs) {
        x <- env[[nm]]

        if (
          is.data.frame(x) ||
            is.data.table(x)
        ) {
          return(
            list(
              name = nm,
              data = as.data.table(x)
            )
          )
        }
      }
    }
  }

  candidates <- list()

  for (nm in objs) {

    x <- env[[nm]]

    if (
      is.data.frame(x) ||
        is.data.table(x)
    ) {

      xx <- as.data.table(x)

      if (
        any(
          must_have_any %in%
            names(xx)
        )
      ) {
        candidates[[nm]] <- xx
      }
    }
  }

  if (length(candidates) == 0L) {
    stop(
      "No plausible result table found in loaded environment."
    )
  }

  sizes <- vapply(
    candidates,
    nrow,
    numeric(1)
  )

  best <- names(
    sizes
  )[
    which.max(
      sizes
    )
  ]

  list(
    name = best,
    data = candidates[[best]]
  )
}


standardize_result_table <- function(
  x,
  dataset_source = NULL,
  force_dataset = NULL
) {

  x <- copy(
    as.data.table(x)
  )


  # ---------------------------------------------------------------------------
  # Core columns
  # ---------------------------------------------------------------------------

  p_col <- pick_col(
    c(
      "p.value",
      "p_value",
      "pvalue",
      "P.value",
      "PValue",
      "pval"
    ),
    names(x),
    required = TRUE,
    label = "p-value"
  )

  q_col <- pick_col(
    c(
      "q.value",
      "q_value",
      "qvalue",
      "BH_q",
      "BH.q",
      "padj",
      "FDR",
      "fdr"
    ),
    names(x),
    required = FALSE
  )

  sigma_col <- pick_col(
    c(
      "Sigma",
      "sigma_c_hat",
      "sigma_hat",
      "sigma_c",
      "Sigma_full",
      "sigma_c_est"
    ),
    names(x),
    required = FALSE
  )

  gene_col <- pick_col(
    c(
      "gene_symbol",
      "Gene",
      "gene",
      "gene_name",
      "symbol",
      "GeneSymbol"
    ),
    names(x),
    required = FALSE
  )

  event_col <- pick_col(
    c(
      "event",
      "Event",
      "event_id",
      "EventID",
      "ID",
      "ri_event"
    ),
    names(x),
    required = FALSE
  )

  dataset_col <- pick_col(
    c(
      "dataset",
      "Dataset",
      "dataset_name",
      "source",
      "study"
    ),
    names(x),
    required = FALSE
  )

  rss0_col <- pick_col(
    c(
      "RSS0",
      "rss0",
      "RSS_null",
      "rss_null"
    ),
    names(x),
    required = FALSE
  )

  rss1_col <- pick_col(
    c(
      "RSS1",
      "rss1",
      "RSS_full",
      "rss_full"
    ),
    names(x),
    required = FALSE
  )


  # ---------------------------------------------------------------------------
  # Standardized fields
  # ---------------------------------------------------------------------------

  x[
    ,
    p_final :=
      as.numeric(
        get(
          p_col
        )
      )
  ]


  if (!is.null(q_col)) {

    x[
      ,
      q_final :=
        as.numeric(
          get(
            q_col
          )
        )
    ]

  } else {

    # q-values will be recomputed later within dataset
    x[
      ,
      q_final :=
        NA_real_
    ]
  }


  if (!is.null(sigma_col)) {

    x[
      ,
      sigma_c_final :=
        as.numeric(
          get(
            sigma_col
          )
        )
    ]

  } else {

    x[
      ,
      sigma_c_final :=
        NA_real_
    ]
  }


  if (!is.null(gene_col)) {

    x[
      ,
      gene_final :=
        as.character(
          get(
            gene_col
          )
        )
    ]

  } else {

    x[
      ,
      gene_final :=
        NA_character_
    ]
  }


  if (!is.null(event_col)) {

    x[
      ,
      event_final :=
        as.character(
          get(
            event_col
          )
        )
    ]

  } else {

    x[
      ,
      event_final :=
        NA_character_
    ]
  }


  if (!is.null(force_dataset)) {

    x[
      ,
      dataset_final :=
        force_dataset
    ]

  } else if (!is.null(dataset_col)) {

    x[
      ,
      dataset_final :=
        normalize_dataset_label(
          get(
            dataset_col
          )
        )
    ]

  } else if (!is.null(dataset_source)) {

    x[
      ,
      dataset_final :=
        dataset_source
    ]

  } else {

    stop(
      "Could not determine dataset label."
    )
  }


  if (
    !is.null(rss0_col) &&
      !is.null(rss1_col)
  ) {

    x[
      ,
      IR_final :=
        fifelse(
          as.numeric(
            get(
              rss0_col
            )
          ) >
            0,
          pmax(
            as.numeric(
              get(
                rss0_col
              )
            ) -
              as.numeric(
                get(
                  rss1_col
                )
              ),
            0
          ) /
            as.numeric(
              get(
                rss0_col
              )
            ),
          NA_real_
        )
    ]

  } else {

    ir_col <- pick_col(
      c(
        "IR",
        "Relative_RSS",
        "relative_rss",
        "RSS_improvement",
        "RelativeRSS"
      ),
      names(x),
      required = FALSE
    )

    if (!is.null(ir_col)) {
      x[
        ,
        IR_final :=
          as.numeric(
            get(
              ir_col
            )
          )
      ]
    } else {
      x[
        ,
        IR_final :=
          NA_real_
      ]
    }
  }


  x
}


# =============================================================================
# 4. Load revision real-data results
# =============================================================================

env_old <- new.env()

load(
  OLD_REALDATA_FILE,
  envir = env_old
)

cat(
  "\n============================================================\n",
  "LOADED: ",
  OLD_REALDATA_FILE,
  "\n",
  "Objects:\n  ",
  paste(
    ls(env_old),
    collapse = "\n  "
  ),
  "\n",
  sep = ""
)


old_obj <- extract_table_from_env(
  env_old,
  preferred = c(
    "results_realdatasets_revision",
    "results_revision",
    "results",
    "res"
  )
)

cat(
  "\nUsing object for Kc167/K562/NIH-3T3:\n  ",
  old_obj$name,
  "\nRows: ",
  nrow(
    old_obj$data
  ),
  "\n",
  sep = ""
)


old_dt <- standardize_result_table(
  old_obj$data
)


# Remove mESC from the revision table if present.
old_dt <- old_dt[
  dataset_final %chin%
    c(
      "Kc167",
      "K562",
      "NIH-3T3"
    )
]


# =============================================================================
# 5. Load FINAL mESC 20k results
# =============================================================================

env_mesc <- new.env()

load(
  MESC_FILE,
  envir = env_mesc
)

if (!"results_mesc_20k" %in% ls(env_mesc)) {
  stop(
    "Object 'results_mesc_20k' not found in mesc_20k_results.rdata."
  )
}

mesc_dt <- standardize_result_table(
  env_mesc$results_mesc_20k,
  force_dataset = "mESC"
)


cat(
  "\n============================================================\n",
  "LOADED FINAL mESC 20k\n",
  "Rows: ",
  nrow(
    mesc_dt
  ),
  "\n",
  sep = ""
)


# =============================================================================
# 6. Combine final sources
# =============================================================================

keep_cols <- c(
  "dataset_final",
  "gene_final",
  "event_final",
  "p_final",
  "q_final",
  "sigma_c_final",
  "IR_final"
)


all_dt <- rbindlist(
  list(
    old_dt[
      ,
      ..keep_cols
    ],
    mesc_dt[
      ,
      ..keep_cols
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)


# =============================================================================
# 7. Sanity checks
# =============================================================================

all_dt <- all_dt[
  is.finite(
    p_final
  )
]


if (
  any(
    all_dt$p_final <
      0 |
      all_dt$p_final >
        1,
    na.rm = TRUE
  )
) {
  stop(
    "Found p-values outside [0,1]."
  )
}


# Recompute BH q-values WITHIN DATASET if:
#   - missing;
#   - non-finite;
#   - or obviously inconsistent in length.
#
# We retain source q-values only when complete.
all_dt[
  ,
  q_source_complete :=
    all(
      is.finite(
        q_final
      )
    ),
  by = dataset_final
]


all_dt[
  ,
  q_recomputed :=
    p.adjust(
      p_final,
      method = "BH"
    ),
  by = dataset_final
]


all_dt[
  !is.finite(
    q_final
  ),
  q_final :=
    q_recomputed
]


# =============================================================================
# 8. Audit source-vs-recomputed q-values
# =============================================================================

q_audit <- all_dt[
  ,
  .(
    N =
      .N,

    Source_q_complete =
      all(
        q_source_complete
      ),

    Max_abs_q_difference =
      if (
        all(
          is.finite(
            q_final
          )
        )
      ) {
        max(
          abs(
            q_final -
              q_recomputed
          ),
          na.rm = TRUE
        )
      } else {
        NA_real_
      }
  ),
  by = dataset_final
]


cat(
  "\n============================================================\n",
  "Q-VALUE AUDIT\n",
  "============================================================\n",
  sep = ""
)

print(
  q_audit
)


# =============================================================================
# 9. Final dataset counts
# =============================================================================

dataset_counts <- all_dt[
  ,
  .(
    N_tested =
      .N,

    N_p_lt_005 =
      sum(
        p_final <
          0.05,
        na.rm = TRUE
      ),

    N_q_lt_010 =
      sum(
        q_final <
          0.10,
        na.rm = TRUE
      ),

    N_q_lt_005 =
      sum(
        q_final <
          0.05,
        na.rm = TRUE
      ),

    N_unique_genes =
      safe_unique_n(
        gene_final
      ),

    N_unique_events =
      safe_unique_n(
        event_final
      )
  ),
  by = dataset_final
]


dataset_order <- c(
  "Kc167",
  "K562",
  "NIH-3T3",
  "mESC"
)

dataset_counts[
  ,
  dataset_final :=
    factor(
      dataset_final,
      levels = dataset_order
    )
]

setorder(
  dataset_counts,
  dataset_final
)

dataset_counts[
  ,
  dataset_final :=
    as.character(
      dataset_final
    )
]


# =============================================================================
# 10. Boundary summary
# =============================================================================

boundary_summary <- all_dt[
  ,
  .(
    N_with_sigma =
      sum(
        is.finite(
          sigma_c_final
        )
      ),

    N_sigma_boundary =
      sum(
        is.finite(
          sigma_c_final
        ) &
          abs(
            sigma_c_final
          ) <
            1e-12
      ),

    Boundary_fraction =
      mean(
        abs(
          sigma_c_final[
            is.finite(
              sigma_c_final
            )
          ]
        ) <
          1e-12
      )
  ),
  by = dataset_final
]


boundary_summary[
  ,
  dataset_final :=
    factor(
      dataset_final,
      levels = dataset_order
    )
]

setorder(
  boundary_summary,
  dataset_final
)

boundary_summary[
  ,
  dataset_final :=
    as.character(
      dataset_final
    )
]


# =============================================================================
# 11. mESC FDR tables
# =============================================================================

mesc_final <- all_dt[
  dataset_final ==
    "mESC"
]


mesc_fdr05 <- mesc_final[
  q_final <
    0.05
][
  order(
    q_final,
    p_final,
    -IR_final
  )
]


mesc_fdr10 <- mesc_final[
  q_final <
    0.10
][
  order(
    q_final,
    p_final,
    -IR_final
  )
]


# =============================================================================
# 12. Representative events from FINAL mESC source
# =============================================================================

REP_GENES <- c(
  "Ppp1r36dn",
  "Nsd1"
)


representative_events <- mesc_final[
  gene_final %chin%
    REP_GENES
][
  order(
    match(
      gene_final,
      REP_GENES
    ),
    q_final,
    p_final
  )
]


if (
  !all(
    REP_GENES %chin%
      representative_events$gene_final
  )
) {
  warning(
    paste0(
      "One or more representative genes were not found: ",
      paste(
        setdiff(
          REP_GENES,
          representative_events$gene_final
        ),
        collapse = ", "
      )
    )
  )
}


# =============================================================================
# 13. Print final audit
# =============================================================================

cat(
  "\n============================================================\n",
  "FINAL REAL-DATA AUDIT\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nDataset-level discovery counts:\n"
)

print(
  dataset_counts
)


cat(
  "\nBoundary summary:\n"
)

print(
  boundary_summary
)


cat(
  "\nFINAL mESC FDR < 0.05 events: ",
  nrow(
    mesc_fdr05
  ),
  "\n",
  sep = ""
)


cat(
  "FINAL mESC FDR < 0.10 events: ",
  nrow(
    mesc_fdr10
  ),
  "\n",
  sep = ""
)


cat(
  "\nRepresentative events from the same FINAL mESC source:\n"
)

print(
  representative_events[
    ,
    .(
      gene_final,
      event_final,
      p_final,
      q_final,
      sigma_c_final,
      IR_final
    )
  ]
)


# =============================================================================
# 14. Internal consistency checks
# =============================================================================

cat(
  "\n============================================================\n",
  "CONSISTENCY CHECKS\n",
  "============================================================\n",
  sep = ""
)


expected_rep_q <- data.table(
  gene_final = c(
    "Ppp1r36dn",
    "Nsd1"
  ),
  expected_q = c(
    0.019720,
    0.012325
  )
)


rep_check <- merge(
  expected_rep_q,
  representative_events[
    ,
    .(
      gene_final,
      q_final
    )
  ],
  by = "gene_final",
  all.x = TRUE
)


rep_check[
  ,
  abs_difference :=
    abs(
      q_final -
        expected_q
    )
]


print(
  rep_check
)


if (
  any(
    rep_check$abs_difference >
      1e-5,
    na.rm = TRUE
  )
) {

  warning(
    paste0(
      "Representative-event q-values differ from the figure values. ",
      "Do NOT freeze Table 5/Figure 5 until resolved."
    )
  )

} else {

  cat(
    "\nRepresentative-event q-values match the finalized figure.\n"
  )
}


# =============================================================================
# 15. Save
# =============================================================================

fwrite(
  dataset_counts,
  OUT_COUNTS,
  sep = "\t"
)

fwrite(
  boundary_summary,
  OUT_BOUNDARY,
  sep = "\t"
)

fwrite(
  mesc_fdr05,
  OUT_MESC_FDR05,
  sep = "\t"
)

fwrite(
  mesc_fdr10,
  OUT_MESC_FDR10,
  sep = "\t"
)

fwrite(
  representative_events,
  OUT_REP,
  sep = "\t"
)

fwrite(
  all_dt,
  OUT_ALL_STANDARDIZED,
  sep = "\t"
)


final_real_data_audit <- list(
  source_files = list(
    pseudo_shutoff =
      OLD_REALDATA_FILE,
    mesc =
      MESC_FILE
  ),

  q_audit =
    q_audit,

  dataset_counts =
    dataset_counts,

  boundary_summary =
    boundary_summary,

  mesc_FDR05 =
    mesc_fdr05,

  mesc_FDR10 =
    mesc_fdr10,

  representative_events =
    representative_events,

  all_standardized =
    all_dt
)


save(
  final_real_data_audit,
  file = OUT_RDATA
)


cat(
  "\n============================================================\n",
  "SAVED\n",
  "============================================================\n",
  "  ",
  OUT_COUNTS,
  "\n  ",
  OUT_BOUNDARY,
  "\n  ",
  OUT_MESC_FDR05,
  "\n  ",
  OUT_MESC_FDR10,
  "\n  ",
  OUT_REP,
  "\n  ",
  OUT_ALL_STANDARDIZED,
  "\n  ",
  OUT_RDATA,
  "\n",
  sep = ""
)



library(data.table)

x <- fread(
  "final_real_data_audit/final_mesc_FDR05_events.tsv"
)

x <- x[
  order(q_final, p_final, -IR_final)
]

x[
  ,
  half_time_min :=
    fifelse(
      sigma_c_final > 0,
      log(2) / sigma_c_final,
      Inf
    )
]

cat(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\small\n",
  "\\caption{\\textbf{FDR-supported retained-intron events in the pharmacological-shutoff mESC dataset.} ",
  "The table reports the final event-level analysis after Benjamini--Hochberg correction. ",
  "$\\widehat{\\sigma}_c$ is the estimated post-export conversion rate, ",
  "$t_{1/2,c}=\\log(2)/\\widehat{\\sigma}_c$ is the corresponding characteristic half-time, ",
  "and $I_R$ is the relative RSS improvement of the full model over the constrained null model.}\n",
  "\\label{tab:mesc-fdr}\n",
  "\\begin{tabular}{lrrrrr}\n",
  "\\hline\n",
  "Gene & $\\widehat{\\sigma}_c$ (min$^{-1}$) & $t_{1/2,c}$ (min) & $I_R$ & $p$ & BH $q$ \\\\\n",
  "\\hline\n",
  sep = ""
)

for (i in seq_len(nrow(x))) {

  cat(
    sprintf(
      "%s & %.5f & %.1f & %.3f & %.5g & %.5g \\\\\n",
      x$gene_final[i],
      x$sigma_c_final[i],
      x$half_time_min[i],
      x$IR_final[i],
      x$p_final[i],
      x$q_final[i]
    )
  )
}

cat(
  "\\hline\n",
  "\\end{tabular}\n",
  "\\end{table}\n",
  sep = ""
)