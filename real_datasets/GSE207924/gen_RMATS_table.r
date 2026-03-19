# =============================================================================
# Title: Preprocessing of the GSE207924 Human and Mouse Datasets
# Description:
#   This script prepares the GSE207924 datasets (K562 and 3T3) for downstream
#   nested-test analysis.
#
#   Main steps for each cell type:
#     - load compartment-specific fraction files,
#     - reshape time-course fraction data into long format,
#     - load normalized rMATS-based count data,
#     - compute sample-level library-size normalization factors,
#     - join the fraction estimates to the event-level count table,
#     - construct event identifiers,
#     - parse time and replicate metadata,
#     - map compartment/version combinations to the four model states,
#     - aggregate counts into event-by-time-by-replicate measurements,
#     - attach Ensembl-based gene annotations,
#     - save the processed event-level dataset.
#
# Intended use:
#   Dataset-specific preprocessing for the K562 and 3T3 real-data analyses.
#
# Copyright (c) 2026 Luigi Cerulo
#
# Permission is hereby granted to use, copy, modify, and distribute this
# software for academic and research purposes, provided that this notice is
# retained in all copies.
#
# This software is provided "as is", without warranty of any kind, express
# or implied, including but not limited to the warranties of merchantability,
# fitness for a particular purpose, and noninfringement. In no event shall
# the authors be liable for any claim, damages, or other liability arising
# from, out of, or in connection with the software or its use.
# =============================================================================


# -----------------------------------------------------------------------------
# Set working directory for the GSE207924 datasets.
# -----------------------------------------------------------------------------
setwd("~/postexport-kinetics/real_datasets/GSE207924")

library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(biomaRt)


# -----------------------------------------------------------------------------
# Define the cell types to process.
#
# K562:
#   human dataset
#
# 3T3:
#   mouse dataset
# -----------------------------------------------------------------------------
cells = c("K562", "3T3")

# Species-specific annotation databases.
orgdb = c(
  "3T3" = org.Mm.eg.db,
  "K562" = org.Hs.eg.db
)

# Species-specific Ensembl dataset names.
ensdb = c(
  "3T3" = "mmusculus_gene_ensembl",
  "K562" = "hsapiens_gene_ensembl"
)


# -----------------------------------------------------------------------------
# Process each cell type independently.
# -----------------------------------------------------------------------------
for (cc in cells) {

  # ---------------------------------------------------------------------------
  # Load compartment-specific time-course fraction files.
  #
  # Expected file naming pattern:
  #   <cell>_(nuclear|cytoplasm)_WT_rep<rep>.tsv
  # ---------------------------------------------------------------------------
  pattern <- paste0("^", cc, "_(nuclear|cytoplasm)_WT_rep\\d+\\.tsv$")
  files <- list.files(path = ".", pattern = pattern, full.names = TRUE)

  # ---------------------------------------------------------------------------
  # Read all matching files and attach metadata inferred from file names.
  #
  # For each file:
  #   - keep Gene / Symbol plus Mean columns,
  #   - rename time columns to t0, t15, t30, t60, t120,
  #   - add compartment and replicate metadata.
  # ---------------------------------------------------------------------------
  dt_all <- rbindlist(lapply(files, function(f) {
    d <- fread(f)

    mean_cols <- grep("Mean", names(d), value = TRUE)

    # Parse metadata from the filename.
    bn <- basename(f)

    fraction <- sub(
      paste0("^", cc, "_(nuclear|cytoplasm)_WT_rep\\d+\\.tsv$"),
      "\\1",
      bn
    )

    replicate <- as.integer(sub(
      paste0("^", cc, "_(?:nuclear|cytoplasm)_WT_rep(\\d+)\\.tsv$"),
      "\\1",
      bn
    ))

    # Keep only identifier columns and time-wise mean columns.
    base_cols <- intersect(c("Gene", "Symbol"), names(d))
    d <- d[, c(base_cols, mean_cols), with = FALSE]

    # Rename mean columns to explicit timepoint labels.
    tp_levels <- c("t0", "t15", "t30", "t60", "t120")
    setnames(
      d,
      old = mean_cols,
      new = tp_levels[seq_along(mean_cols)]
    )

    # Add compartment and replicate metadata.
    d[, `:=`(compartment = fraction, replicate = replicate)]

    d
  }), use.names = TRUE, fill = TRUE)


  # ---------------------------------------------------------------------------
  # Convert the fraction table into long format.
  #
  # Each row represents one:
  #   - gene
  #   - timepoint
  #   - compartment
  #   - replicate
  #   - new_fraction
  # ---------------------------------------------------------------------------
  dt_long <- melt(
    dt_all,
    id.vars = c("Gene", "Symbol", "compartment", "replicate"),
    measure.vars = c("t0", "t15", "t30", "t60", "t120"),
    variable.name = "timepoint",
    value.name = "new_fraction"
  )

  # Harmonize gene column name.
  setnames(dt_long, "Gene", "gene")

  # Map long compartment names to the shorter convention used elsewhere.
  map <- c(
    nuclear = "nuc",
    cytoplasm = "cyt"
  )
  dt_long[, compartment := map[compartment]]

  # Convert replicate to the same string format used in the rMATS input.
  dt_long[, replicate := paste0("r", replicate)]


  # ---------------------------------------------------------------------------
  # Load normalized rMATS-based event table.
  #
  # An alternative raw CSV is left commented for reference.
  # ---------------------------------------------------------------------------
  # d_rmats = fread(paste0("rmats_gse207924_", cc, ".csv"))
  d_rmats = fread(paste0("rmats_gse207924_", cc, "_normalized.csv"))


  # ---------------------------------------------------------------------------
  # Build a sample identifier and compute library-size normalization factors.
  # ---------------------------------------------------------------------------
  d_rmats[, sample_id := paste(timepoint, compartment, replicate, sep = "_")]

  libDT <- d_rmats[, .(libsize = sum(count, na.rm = TRUE)), by = sample_id]
  libDT[, lib_factor := libsize / median(libsize)]

  # Join normalization factors back and compute normalized counts.
  d_rmats <- libDT[d_rmats, on = "sample_id"]
  d_rmats[, count_norm := count / lib_factor]


  # ---------------------------------------------------------------------------
  # Join the time-course fraction estimates to the rMATS event table.
  #
  # The intended key is:
  #   gene + timepoint + compartment + replicate
  # ---------------------------------------------------------------------------
  dt_long2 <- unique(dt_long[, .(gene, timepoint, compartment, replicate, new_fraction)])

  # Check whether the join key is unique in the fraction table.
  dt_long2[duplicated(dt_long2, by = c("gene", "timepoint", "compartment", "replicate"))]

  setkey(d_rmats, gene, timepoint, compartment, replicate)
  setkey(dt_long2, gene, timepoint, compartment, replicate)

  # Left join: preserve all rows from the rMATS table.
  d_rmats <- dt_long2[d_rmats]

  # Retain only rows for which a fraction estimate is available.
  d_rmats <- d_rmats[!is.na(new_fraction)]


  # ---------------------------------------------------------------------------
  # Build a unique event identifier from genomic coordinates and gene ID.
  # ---------------------------------------------------------------------------
  d_rmats[, event := paste(
    gene, chr, strand,
    upstreamES, upstreamEE, downstreamES, downstreamEE,
    sep = ":"
  )]

  # Parse numeric time and replicate variables.
  d_rmats[, time := as.numeric(sub("^t", "", timepoint))]
  d_rmats[, replicate := as.integer(sub("^r", "", replicate))]


  # ---------------------------------------------------------------------------
  # Map each observation to one of the four model states:
  #   - N   : nuclear inclusion
  #   - N_s : nuclear skipping
  #   - C   : cytoplasmic inclusion
  #   - C_s : cytoplasmic skipping
  # ---------------------------------------------------------------------------
  d_rmats[, State := fifelse(
    compartment == "nuc" & version == "inclusion", "N",
    fifelse(
      compartment == "nuc" & version == "skipping", "N_s",
      fifelse(
        compartment == "cyt" & version == "inclusion", "C",
        fifelse(
          compartment == "cyt" & version == "skipping", "C_s",
          NA_character_
        )
      )
    )
  )]


  # ---------------------------------------------------------------------------
  # Aggregate counts into event-level measurements.
  #
  # The current script weights counts by (1 - new_fraction), effectively
  # down-weighting observations according to the joined fraction estimate.
  # ---------------------------------------------------------------------------
  agg = d_rmats[, .(
    Value = (1 - new_fraction) * sum(count, na.rm = TRUE)
  ), by = .(event, time, replicate, State)]


  # ---------------------------------------------------------------------------
  # Reshape into wide format with one column per state:
  #   N, N_s, C, C_s
  # ---------------------------------------------------------------------------
  rmats_dt = dcast(
    agg,
    event + time + replicate ~ State,
    value.var = "Value",
    fill = 0
  )

  # Extract Ensembl gene ID from the event string and remove version suffixes.
  rmats_dt[, ensembl := sub(":.*", "", event)]
  rmats_dt[, ensembl := sub("\\..*", "", ensembl)]


  # ===========================================================================
  # Retrieve gene metadata from Ensembl
  # ===========================================================================

  gene_ids <- unique(rmats_dt$ensembl)

  mart <- useEnsembl(
    biomart = "genes",
    dataset = ensdb[[cc]],
    mirror = "useast"
  )

  gene_info <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "description",
      "gene_biotype"
    ),
    filters = "ensembl_gene_id",
    values  = gene_ids,
    mart    = mart
  )

  rownames(gene_info) <- gene_info$ensembl_gene_id

  # Join Ensembl annotations back to the event-level table.
  rmats_dt[gene_info, on = .(ensembl = ensembl_gene_id),
    `:=`(
      gene_symbol = i.external_gene_name,
      description = i.description,
      gene_biotype = i.gene_biotype
    )]


  # ---------------------------------------------------------------------------
  # Save the processed dataset for this cell type.
  # ---------------------------------------------------------------------------
  save(rmats_dt, file = paste0("rmats_dt_", cc, ".rdata"))
}