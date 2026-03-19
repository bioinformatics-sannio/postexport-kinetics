# =============================================================================
# Title: Preprocessing of the GSE83620 Drosophila Dataset
# Description:
#   This script prepares the Drosophila GSE83620 dataset for downstream
#   nested-test analysis.
#
#   Main steps:
#     - load normalized/count-level event data,
#     - compute sample-level library-size normalization factors,
#     - normalize counts across replicates,
#     - construct event identifiers,
#     - parse time and replicate metadata,
#     - map compartment/version combinations to the four model states,
#     - aggregate counts into event-by-time-by-replicate measurements,
#     - reshape to wide format,
#     - attach Ensembl-based gene annotations,
#     - save the processed event-level table.
#
# Intended use:
#   Dataset-specific preprocessing for the Kc167 / GSE83620 real-data analysis.
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
# Set working directory for the Drosophila real dataset.
# -----------------------------------------------------------------------------
setwd("~/postexport-kinetics/real_datasets/GSE83620")

# Dataset: Drosophila
require(data.table)
library(org.Dm.eg.db)
library(biomaRt)


# -----------------------------------------------------------------------------
# Load the preprocessed/normalized input table.
#
# An alternative raw rMATS-based file is left commented for reference.
# -----------------------------------------------------------------------------
# d_rmats = fread("rmats_gse83620_drosophila.csv")
d_rmats = fread("gse83620_drosophila_normalized.csv")


# -----------------------------------------------------------------------------
# Construct a unique sample identifier from:
#   - timepoint
#   - compartment
#   - replicate
#
# This is used to compute library-size normalization factors.
# -----------------------------------------------------------------------------
d_rmats[, sample_id := paste(timepoint, compartment, replicate, sep = "_")]


# -----------------------------------------------------------------------------
# Compute sample-level library sizes and corresponding normalization factors.
#
# Each sample is scaled relative to the median library size.
# -----------------------------------------------------------------------------
libDT <- d_rmats[, .(libsize = sum(count, na.rm = TRUE)), by = sample_id]
libDT[, lib_factor := libsize / median(libsize)]


# -----------------------------------------------------------------------------
# Join the normalization factor back to the main table and normalize counts.
#
# Normalization is applied at the replicate/sample level.
# -----------------------------------------------------------------------------
d_rmats <- libDT[d_rmats, on = "sample_id"]
d_rmats[, count_norm := count / lib_factor]


# -----------------------------------------------------------------------------
# Build a unique event identifier from genomic coordinates and gene ID.
#
# Event format:
#   gene:chr:strand:upstreamES:upstreamEE:downstreamES:downstreamEE
# -----------------------------------------------------------------------------
d_rmats[, event := paste(
  gene, chr, strand,
  upstreamES, upstreamEE, downstreamES, downstreamEE,
  sep = ":"
)]


# -----------------------------------------------------------------------------
# Parse time and replicate variables into numeric form.
#
# Example transformations:
#   - timepoint "t1_5" -> 1.5
#   - replicate "r2"   -> 2
# -----------------------------------------------------------------------------
d_rmats[, time := as.numeric(sub("_", ".", sub("^t", "", timepoint)))]
d_rmats[, replicate := as.integer(sub("^r", "", replicate))]


# -----------------------------------------------------------------------------
# Map each observation to one of the four model states:
#   - N   : nuclear inclusion
#   - N_s : nuclear skipping
#   - C   : cytoplasmic inclusion
#   - C_s : cytoplasmic skipping
# -----------------------------------------------------------------------------
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


# -----------------------------------------------------------------------------
# Aggregate counts at the level of:
#   - event
#   - time
#   - replicate
#   - model state
#
# The script currently sums raw counts rather than count_norm.
# -----------------------------------------------------------------------------
agg = d_rmats[, .(Value = sum(count, na.rm = TRUE)),
              by = .(event, time, replicate, State)]


# -----------------------------------------------------------------------------
# Reshape into wide format with one column per state:
#   N, N_s, C, C_s
# -----------------------------------------------------------------------------
rmats_dt = dcast(
  agg,
  event + time + replicate ~ State,
  value.var = "Value",
  fill = 0
)


# -----------------------------------------------------------------------------
# Extract Ensembl gene ID from the event identifier.
#
# Also remove Ensembl version suffixes, if present.
# -----------------------------------------------------------------------------
rmats_dt[, ensembl := sub(":.*", "", event)]
rmats_dt[, ensembl := sub("\\..*", "", ensembl)]


# -----------------------------------------------------------------------------
# Convert time from hours to minutes for consistency with downstream modeling.
# -----------------------------------------------------------------------------
rmats_dt[, time := time * 60]


# =============================================================================
# Retrieve gene metadata from Ensembl
# =============================================================================

# -----------------------------------------------------------------------------
# Query Ensembl for:
#   - external gene symbol
#   - description
#   - gene biotype
# -----------------------------------------------------------------------------
gene_ids <- unique(rmats_dt$ensembl)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "dmelanogaster_gene_ensembl",
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


# -----------------------------------------------------------------------------
# Join Ensembl annotations back into the event-level table.
# -----------------------------------------------------------------------------
rmats_dt[gene_info, on = .(ensembl = ensembl_gene_id),
  `:=`(
    gene_symbol = i.external_gene_name,
    description = i.description,
    gene_biotype = i.gene_biotype
  )]


# -----------------------------------------------------------------------------
# Save the processed dataset for downstream real-data analysis.
# -----------------------------------------------------------------------------
save(rmats_dt, file = paste0("rmats_dt_Kc167.rdata"))