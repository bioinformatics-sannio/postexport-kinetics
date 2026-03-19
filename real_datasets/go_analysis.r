# =============================================================================
# Title: Functional Enrichment and Cross-Dataset Comparison on Real Data Results
# Description:
#   This script performs downstream biological interpretation of the nested-test
#   results obtained on real datasets.
#
#   Main analyses:
#     - assign each dataset to a species,
#     - perform GO over-representation analysis (ORA) on selected events,
#     - compare enriched GO terms across datasets,
#     - compare gene sets across datasets through ortholog mapping,
#     - run GSEA-style enrichment analyses (GO, KEGG, Reactome),
#     - inspect enriched pathways and recover contributing genes/events,
#     - generate illustrative plots for selected positive examples.
#
# Intended use:
#   Downstream enrichment analysis and biological interpretation of the
#   post-export kinetics results on real datasets.
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


setwd("~/postexport-kinetics/real_datasets")

require(data.table)
require(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(org.Dm.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(ReactomePA)


# -----------------------------------------------------------------------------
# Load the real-dataset nested-test results.
# -----------------------------------------------------------------------------
load("results_realdatasets.rdata", verbose = TRUE)


# -----------------------------------------------------------------------------
# Assign each dataset to a species code.
#
# Species labels used here:
#   - hs    : human
#   - dm    : drosophila
#   - mm3t3 : mouse 3T3
#   - mmesc : mouse ESC
# -----------------------------------------------------------------------------
results[, species := fcase(
  grepl("^K562", dataset),  "hs",
  grepl("^Kc167", dataset), "dm",
  grepl("^3T3", dataset),   "mm3t3",
  grepl("^mESC", dataset),  "mmesc",
  default = "unknown"
)]


# =============================================================================
# GO over-representation analysis on selected event-associated genes
# =============================================================================

# -----------------------------------------------------------------------------
# Run GO Biological Process ORA separately for each dataset/species.
#
# For each dataset:
#   - select the organism-specific annotation database,
#   - load the expressed-gene background,
#   - restrict candidate genes to that expressed background,
#   - perform enrichGO(),
#   - return the enrichment table.
#
# Note:
#   The input gene universe here is defined by the expressed genes from each
#   dataset-specific background file, rather than by the full genome.
# -----------------------------------------------------------------------------
ora_ir <- results[, {

  org_db = switch(
    species,
    "mm3t3" = org.Mm.eg.db,
    "mmesc" = org.Mm.eg.db,
    "hs"    = org.Hs.eg.db,
    "dm"    = org.Dm.eg.db
  )

  bg_file = switch(
    species,
    "mm3t3" = "GSE207924/gse207924_mouse_genes.csv",
    "mmesc" = "GSE256335/gse256335_mouse_genes.csv",
    "hs"    = "GSE207924/gse207924_human_genes.csv",
    "dm"    = "GSE83620/gse83620_drosophila_genes.csv"
  )

  # Dataset-specific expressed-gene background.
  bg_d = fread(bg_file)
  bg_genome = unique(bg_d$gene)

  # Candidate genes restricted to the background.
  genes <- intersect(ensembl, bg_genome)

  ego <- enrichGO(
    gene          = genes,
    universe      = bg_genome,
    OrgDb         = org_db,
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.01,
    minGSSize     = 10,
    readable      = TRUE
  )

  df <- as.data.table(ego@result)

  df[, .(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count)]
}, by = .(dataset, species)]


# Example inspection: significant drosophila GO terms.
ora_ir[species == "dm" & p.adjust < 0.01][order(Description)]


# -----------------------------------------------------------------------------
# Collect significant GO term sets per dataset for overlap analysis.
# -----------------------------------------------------------------------------
fdr_thr <- 0.01
go_sets <- ora_ir[p.adjust < fdr_thr, .(go_ids = list(unique(ID))), by = dataset]


# -----------------------------------------------------------------------------
# Helper to compute pairwise overlap matrices between set collections.
#
# Returns:
#   - n_intersect: number of shared terms
#   - jaccard: Jaccard similarity
# -----------------------------------------------------------------------------
overlap_mats <- function(setsDT, col = "go_ids") {
  ds <- setsDT$dataset
  n  <- length(ds)

  mat_n <- matrix(0L, n, n, dimnames = list(ds, ds))
  mat_j <- matrix(0,  n, n, dimnames = list(ds, ds))

  for (i in seq_len(n)) {
    A <- setsDT[[col]][[i]]
    for (j in seq_len(n)) {
      B <- setsDT[[col]][[j]]
      inter <- length(intersect(A, B))
      uni   <- length(union(A, B))
      mat_n[i, j] <- inter
      mat_j[i, j] <- if (uni == 0) NA_real_ else inter / uni
    }
  }

  list(n_intersect = mat_n, jaccard = mat_j)
}

go_mat <- overlap_mats(go_sets, "go_ids")

# Pairwise counts of shared GO terms.
go_mat$n_intersect

# Pairwise Jaccard similarities.
go_mat$jaccard


# -----------------------------------------------------------------------------
# Export overlap count matrix as LaTeX table if needed.
# -----------------------------------------------------------------------------
kbl <- knitr::kable(
  go_mat$n_intersect,
  format = "latex",
  booktabs = TRUE,
  longtable = FALSE,
  escape = TRUE,
  digits = 3
)

kbl


# -----------------------------------------------------------------------------
# Example: inspect GO terms shared by two selected datasets.
# -----------------------------------------------------------------------------
go_sets <- ora_ir[p.adjust < fdr_thr, .(go_ids = list(unique(ID))), by = dataset]

get_set <- function(ds) go_sets[dataset == ds, go_ids][[1]]

A <- get_set("Kc167 | GSE83620")
B <- get_set("mESC | GSE256335")

inter_ids <- intersect(A, B)

ora_ir[
  ID %in% inter_ids &
  dataset %chin% c("Kc167 | GSE83620", "mESC | GSE256335"),
  .(dataset, ID, Description, p.adjust)
][order(ID, dataset)]


# -----------------------------------------------------------------------------
# Build a dataset-by-term summary table for enriched GO terms.
#
# For each GO term, the table stores:
#   - number of datasets in which the term is enriched,
#   - gene ratios for each dataset.
# -----------------------------------------------------------------------------
ds_sets <- ora_ir[p.adjust < fdr_thr, {
  ds_gratio = rep("", 4)
  names(ds_gratio) = c("Kc167", "K562", "3T3", "mESC")
  ds_gratio[gsub(" \\|.*$", "", dataset)] = GeneRatio

  .(
    term   = unique(Description),
    N      = .N,
    G_Kc167 = ds_gratio["Kc167"],
    G_K562  = ds_gratio["K562"],
    G_3T3   = ds_gratio["3T3"],
    G_mESC  = ds_gratio["mESC"]
  )
}, by = ID]

require(openxlsx)
write.xlsx(ds_sets[order(-N)], file = "ora_go_terms.xlsx")


# =============================================================================
# Cross-dataset gene overlap via ortholog mapping
# =============================================================================

library(data.table)
library(biomaRt)

# -----------------------------------------------------------------------------
# Gene sets per dataset.
#
# Here all tested genes are used; alternative commented code shows how one could
# restrict to only significant candidates.
# -----------------------------------------------------------------------------
ir_sets <- results[, .(genes = list(unique(ensembl))), by = .(dataset, species)]

# Example alternative:
# ir_sets <- results[deltaAIC > 0 & p.value < 0.05 & Sigma > 0.001,
#                    .(genes = list(unique(ensembl))), by = .(dataset, species)]

ir_sets[, .(l = lengths(genes)), by = .(dataset)]


# -----------------------------------------------------------------------------
# Initialize Ensembl BioMart connections.
#
# A fixed mirror is used to reduce host mismatch issues.
# -----------------------------------------------------------------------------
mir <- "useast"
mart_hs <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl", mirror = mir)
mart_mm <- useEnsembl("genes", dataset = "mmusculus_gene_ensembl", mirror = mir)
mart_dm <- useEnsembl("genes", dataset = "dmelanogaster_gene_ensembl", mirror = mir)


# -----------------------------------------------------------------------------
# Retry wrapper for getBM() to make BioMart access more robust.
# -----------------------------------------------------------------------------
safe_getBM <- function(..., tries = 5, sleep = 2) {
  for (k in seq_len(tries)) {
    out <- try(getBM(...), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    Sys.sleep(sleep * k)
  }
  stop("getBM failed after retries")
}


# -----------------------------------------------------------------------------
# Map dataset-specific gene identifiers to a common human ENSG space.
#
# Strategy:
#   - human: retain ENSG directly
#   - mouse: use hsapiens_homolog_ensembl_gene
#   - drosophila: map FlyBase -> Ensembl -> human homolog ENSG
# -----------------------------------------------------------------------------
map_to_hs_ensg <- function(sp, ids) {
  ids <- unique(ids)
  ids <- ids[!is.na(ids) & ids != ""]
  if (!length(ids)) return(character())

  if (sp == "hs") {
    hs <- ids[grepl("^ENSG", ids)]
    return(unique(hs))
  }

  if (sp == "mm3t3" | sp == "mmesc") {
    out <- safe_getBM(
      attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"),
      filters    = "ensembl_gene_id",
      values     = ids,
      mart       = mart_mm
    )
    out <- as.data.table(out)
    hs <- unique(out$hsapiens_homolog_ensembl_gene)
    return(hs[!is.na(hs) & hs != ""])
  }

  if (sp == "dm") {
    # Two-step drosophila -> Ensembl -> human homolog mapping.
    dm1 <- safe_getBM(
      attributes = c("flybase_gene_id", "ensembl_gene_id"),
      filters    = "flybase_gene_id",
      values     = ids,
      mart       = mart_dm
    )
    dm1 <- as.data.table(dm1)

    dm_ens <- unique(dm1$ensembl_gene_id)
    dm_ens <- dm_ens[!is.na(dm_ens) & dm_ens != ""]
    if (!length(dm_ens)) return(character())

    dm2 <- safe_getBM(
      attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"),
      filters    = "ensembl_gene_id",
      values     = dm_ens,
      mart       = mart_dm
    )
    dm2 <- as.data.table(dm2)

    hs <- unique(dm2$hsapiens_homolog_ensembl_gene)
    return(hs[!is.na(hs) & hs != ""])
  }

  character()
}


# -----------------------------------------------------------------------------
# Build human-ortholog sets for each dataset.
# -----------------------------------------------------------------------------
ir_hs_sets <- ir_sets[, .(
  hs_orth = list(unique(map_to_hs_ensg(species, genes[[1]]))),
  genes   = list(genes[[1]]),
  n_orig  = length(genes[[1]])
), by = .(dataset, species)]


# -----------------------------------------------------------------------------
# Compute a symmetric overlap matrix:
#   - diagonal: original gene overlap within dataset
#   - off-diagonal: overlap between human ortholog sets
# -----------------------------------------------------------------------------
ds <- ir_hs_sets$dataset
n  <- length(ds)

mat_orth <- matrix(0L, n, n, dimnames = list(ds, ds))

for (i in seq_len(n)) {
  A <- ir_hs_sets$hs_orth[[i]]
  A_g = unique(ir_hs_sets$genes[[i]])

  for (j in seq_len(n)) {
    B <- ir_hs_sets$hs_orth[[j]]
    B_g = unique(ir_hs_sets$genes[[j]])

    if (i == j) {
      mat_orth[i, j] <- length(intersect(A_g, B_g))
    } else {
      mat_orth[i, j] <- length(intersect(A, B))
    }
  }
}

mat_orth


# Export ortholog-overlap matrix as LaTeX table if needed.
kbl <- knitr::kable(
  mat_orth,
  format = "latex",
  booktabs = TRUE,
  longtable = FALSE,
  escape = TRUE,
  digits = 3
)

kbl


# =============================================================================
# GSEA-like enrichment analysis on ranked nested-test results
# =============================================================================

# -----------------------------------------------------------------------------
# Build ranked gene lists per dataset.
#
# For each gene, retain the maximum score across events.
# -----------------------------------------------------------------------------
gene_list = results[, {
  .(score = max(score_sig_IR_logq, na.rm = TRUE))
}, by = .(dataset, ensembl, species)]


# -----------------------------------------------------------------------------
# Run GSEA separately for each dataset/species.
#
# Enrichment databases:
#   - GO Biological Process
#   - KEGG (except drosophila skipped here)
#   - Reactome (except drosophila skipped here)
# -----------------------------------------------------------------------------
gsea_results <- gene_list[, {

  org_db <- switch(
    species,
    "mm3t3" = org.Mm.eg.db,
    "mmesc" = org.Mm.eg.db,
    "hs"    = org.Hs.eg.db,
    "dm"    = org.Dm.eg.db
  )

  key_type <- if (species == "dm") "FLYBASE" else "ENSEMBL"
  kegg_org <- switch(species, "mmesc" = "mmu", "mm3t3" = "mmu", "hs" = "hsa", "dm" = "dme")
  react_org <- switch(species, "mm3t3" = "mouse", "mmesc" = "mouse", "hs" = "human", "dm" = "fly")

  gl <- score
  names(gl) <- ensembl
  gl <- gl[is.finite(gl) & !is.na(names(gl)) & names(gl) != ""]
  gl <- sort(gl, decreasing = TRUE)

  # Keep only genes mappable in the corresponding OrgDb.
  gl <- gl[names(gl) %in% keys(org_db, keytype = key_type)]

  out_list <- list()

  # ---------------------------------------------------------------------------
  # GO Biological Process GSEA
  # ---------------------------------------------------------------------------
  gsea_bp <- gseGO(
    geneList      = gl,
    OrgDb         = org_db,
    keyType       = key_type,
    ont           = "BP",
    pAdjustMethod = "BH",
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = 0.2,
    eps           = 1e-10,
    verbose       = FALSE
  )
  dt_bp <- as.data.table(gsea_bp)[, .(DB = "GO-BP", ID, Description, NES, p.adjust)]
  out_list[["GO_BP"]] <- dt_bp

  # ---------------------------------------------------------------------------
  # Map ranked genes to ENTREZ IDs for KEGG / Reactome.
  # ---------------------------------------------------------------------------
  entrez_map <- mapIds(
    org_db,
    keys      = names(gl),
    column    = "ENTREZID",
    keytype   = key_type,
    multiVals = "first"
  )

  gl_entrez <- gl
  names(gl_entrez) <- unname(entrez_map[names(gl)])
  gl_entrez <- gl_entrez[!is.na(names(gl_entrez)) & names(gl_entrez) != ""]

  # Aggregate duplicate ENTREZ IDs using the maximum score.
  gl_entrez <- tapply(gl_entrez, names(gl_entrez), max, na.rm = TRUE)
  gl_entrez <- sort(gl_entrez, decreasing = TRUE)

  gl_entrez_n = names(gl_entrez)
  gl_entrez = as.numeric(gl_entrez)
  names(gl_entrez) = gl_entrez_n

  # ---------------------------------------------------------------------------
  # KEGG and Reactome GSEA
  #
  # Skipped for drosophila in this script.
  # ---------------------------------------------------------------------------
  if (species != "dm") {

    gsea_kegg <- gseKEGG(
      geneList      = gl_entrez,
      organism      = kegg_org,
      pAdjustMethod = "BH",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.2,
      eps           = 1e-10,
      verbose       = TRUE
    )
    dt_kegg <- as.data.table(gsea_kegg)[, .(DB = "KEGG", ID, Description, NES, p.adjust)]
    out_list[["KEGG"]] <- dt_kegg

    gsea_react <- gsePathway(
      geneList      = gl_entrez,
      organism      = react_org,
      pAdjustMethod = "BH",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.2,
      eps           = 1e-10,
      verbose       = FALSE
    )
    dt_react <- as.data.table(gsea_react)[, .(DB = "Reactome", ID, Description, NES, p.adjust)]
    out_list[["Reactome"]] <- dt_react
  }

  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}, by = .(dataset, species)]


# Inspect positively enriched significant terms.
gsea_results[p.adjust < 0.1 & NES > 0]

