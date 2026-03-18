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



load("results_realdatasets.rdata",verbose=T)

results[, species := fcase(
  grepl("^K562", dataset),  "hs",
  grepl("^Kc167", dataset), "dm",
  grepl("^3T3", dataset),   "mm3t3",
  grepl("^mESC", dataset),  "mmesc",
  default = "unknown"
)]

# =========================
# GO analysis on IR events
# =========================

ora_ir <- results[, {
  org_db = switch(species,
         "mm3t3" = org.Mm.eg.db,
         "mmesc" = org.Mm.eg.db,
         "hs" = org.Hs.eg.db,
         "dm" = org.Dm.eg.db)
  
  bg_file = switch(species,
         "mm3t3" = "GSE207924/gse207924_mouse_genes.csv",
         "mmesc" = "GSE256335/gse256335_mouse_genes.csv",
         "hs" = "GSE207924/gse207924_human_genes.csv",
         "dm" = "GSE83620/gse83620_drosophila_genes.csv")

  #bg_genome <- keys(org_db, keytype = "ENSEMBL") # full genome

  bg_d = fread(bg_file)
  bg_genome = unique(bg_d$gene)
  # only expressed genes

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


ora_ir[species=="dm" & p.adjust<0.01][order(Description)]


fdr_thr <- 0.01  # cambia soglia
go_sets <- ora_ir[p.adjust < fdr_thr, .(go_ids = list(unique(ID))), by = dataset]

# helper per matrix overlap
overlap_mats <- function(setsDT, col="go_ids"){
  ds <- setsDT$dataset
  n  <- length(ds)
  mat_n <- matrix(0L, n, n, dimnames=list(ds, ds))
  mat_j <- matrix(0,  n, n, dimnames=list(ds, ds))

  for(i in seq_len(n)){
    A <- setsDT[[col]][[i]]
    for(j in seq_len(n)){
      B <- setsDT[[col]][[j]]
      inter <- length(intersect(A,B))
      uni   <- length(union(A,B))
      mat_n[i,j] <- inter
      mat_j[i,j] <- if (uni==0) NA_real_ else inter/uni
    }
  }
  list(n_intersect=mat_n, jaccard=mat_j)
}

go_mat <- overlap_mats(go_sets, "go_ids")

go_mat$n_intersect  # numero GO comuni
go_mat$jaccard      # Jaccard GO

kbl <- knitr::kable(
  go_mat$n_intersect,
  format = "latex",
  booktabs = TRUE,
  longtable = FALSE,
  escape = TRUE,
  digits = 3
)

kbl

# set di GO per ciascun dataset
go_sets <- ora_ir[p.adjust < fdr_thr, .(go_ids=list(unique(ID))), by=dataset]

get_set <- function(ds) go_sets[dataset==ds, go_ids][[1]]
A <- get_set("Kc167 | GSE83620")
B <- get_set("mESC | GSE256335")
inter_ids <- intersect(A,B)
ora_ir[ID %in% inter_ids & dataset %chin% c("Kc167 | GSE83620","mESC | GSE256335"),
      .(dataset, ID, Description, p.adjust)][order(ID, dataset)]


ds_sets <- ora_ir[p.adjust < fdr_thr, {
  ds_gratio = rep("",4)
  names(ds_gratio) = c("Kc167","K562","3T3","mESC")
  ds_gratio[gsub(" \\|.*$","",dataset)] = GeneRatio
  
  .(term=unique(Description),
      N = .N,
      G_Kc167 = ds_gratio["Kc167"],
      G_K562 = ds_gratio["K562"],
      G_3T3 = ds_gratio["3T3"],
      G_mESC = ds_gratio["mESC"])
}, by=ID]

require(openxlsx)
write.xlsx(ds_sets[order(-N)], file="ora_go_terms.xlsx")

# =======================
# Matrice geni comuni dataset vs dataset
# =======================

library(data.table)
library(biomaRt)

# IR gene set per dataset (originale)
ir_sets <- results[, .(genes = list(unique(ensembl))), by = .(dataset,species)]
# only significant candidates
#ir_sets <- results[deltaAIC > 0 & p.value < 0.05 & Sigma > 0.001, 
#              .(genes = list(unique(ensembl))), by = .(dataset,species)]
ir_sets[, .(l=lengths(genes)),by=.(dataset)]

# --- biomaRt marts (stesso mirror per evitare host mismatch) ---
mir <- "useast"  # prova "www" o "asia" se 502
mart_hs <- useEnsembl("genes", dataset="hsapiens_gene_ensembl", mirror=mir)
mart_mm <- useEnsembl("genes", dataset="mmusculus_gene_ensembl", mirror=mir)
mart_dm <- useEnsembl("genes", dataset="dmelanogaster_gene_ensembl", mirror=mir)

# retry helper (riduce 502)
safe_getBM <- function(..., tries=5, sleep=2) {
  for (k in seq_len(tries)) {
    out <- try(getBM(...), silent=TRUE)
    if (!inherits(out, "try-error")) return(out)
    Sys.sleep(sleep * k)
  }
  stop("getBM failed after retries")
}

# map -> ENSG (spazio comune umano)
map_to_hs_ensg <- function(sp, ids){
  ids <- unique(ids)
  ids <- ids[!is.na(ids) & ids!=""]
  if (!length(ids)) return(character())

  if (sp == "hs") {
    # già ENSG
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
    return(hs[!is.na(hs) & hs!=""])
  }

  if (sp == "dm") {
    # FBgn -> hs homolog direttamente (se disponibile) è spesso ok,
    # altrimenti fai 2-step (qui metto 2-step robusto).
    dm1 <- safe_getBM(
      attributes = c("flybase_gene_id", "ensembl_gene_id"),
      filters    = "flybase_gene_id",
      values     = ids,
      mart       = mart_dm
    )
    dm1 <- as.data.table(dm1)
    dm_ens <- unique(dm1$ensembl_gene_id)
    dm_ens <- dm_ens[!is.na(dm_ens) & dm_ens!=""]
    if (!length(dm_ens)) return(character())

    dm2 <- safe_getBM(
      attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"),
      filters    = "ensembl_gene_id",
      values     = dm_ens,
      mart       = mart_dm
    )
    dm2 <- as.data.table(dm2)
    hs <- unique(dm2$hsapiens_homolog_ensembl_gene)
    return(hs[!is.na(hs) & hs!=""])
  }

  character()
}

# costruisci ortologhi umani per dataset
ir_hs_sets <- ir_sets[, .(
  hs_orth = list(unique(map_to_hs_ensg(species, genes[[1]]))),
  genes = list(genes[[1]]),
  n_orig  = length(genes[[1]])
), by = .(dataset, species)]

# --- matrice simmetrica: #ortologhi umani condivisi ---
ds <- ir_hs_sets$dataset
n  <- length(ds)
mat_orth <- matrix(0L, n, n, dimnames=list(ds, ds))

for (i in seq_len(n)) {
  A <- ir_hs_sets$hs_orth[[i]]
  A_g = unique(ir_hs_sets$genes[[i]])
  for (j in seq_len(n)) {
    B <- ir_hs_sets$hs_orth[[j]]
    B_g = unique(ir_hs_sets$genes[[j]])
    if (i==j) {
      mat_orth[i, j] <- length(intersect(A_g, B_g))
    } else {
      mat_orth[i, j] <- length(intersect(A, B))
    }
  }
}

mat_orth 

kbl <- knitr::kable(
  mat_orth,
  format = "latex",
  booktabs = TRUE,
  longtable = FALSE,
  escape = TRUE,
  digits = 3
)

kbl


# =========================
# GO enrichment for nested test results
# =========================
# mESC K562 Kc167 3T3 score_siglogp score_aic score_aiclogp
gene_list = results[, {
    .(score=max(score_sig_IR_logq, na.rm = TRUE))
}, by = .(dataset,ensembl,species)]


gsea_results <- gene_list[, {

  org_db <- switch(species,
                   "mm3t3" = org.Mm.eg.db,
                   "mmesc" = org.Mm.eg.db,
                   "hs" = org.Hs.eg.db,
                   "dm" = org.Dm.eg.db)

  key_type <- if (species == "dm") "FLYBASE" else "ENSEMBL"
  kegg_org <- switch(species, "mmesc"="mmu", "mm3t3"="mmu", "hs"="hsa", "dm"="dme")
  react_org <- switch(species, "mm3t3"="mouse","mmesc"="mouse", "hs"="human", "dm"="fly")

  gl <- score
  
  names(gl) <- ensembl
  gl <- gl[is.finite(gl) & !is.na(names(gl)) & names(gl) != ""]
  gl <- sort(gl, decreasing = TRUE)

  # tieni solo geni mappabili nel relativo OrgDb
  gl <- gl[names(gl) %in% keys(org_db, keytype = key_type)]

  out_list <- list()

  # ---------------- GO BP ----------------
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
  dt_bp <- as.data.table(gsea_bp)[,.(DB="GO-BP",ID, Description, NES,p.adjust)]
  out_list[["GO_BP"]] <- dt_bp

  # ---------------- map to ENTREZ ----------------
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

  # aggrega duplicati ENTREZ: max score
  gl_entrez <- tapply(gl_entrez, names(gl_entrez), max, na.rm=TRUE)
  gl_entrez <- sort(gl_entrez, decreasing = TRUE)
  gl_entrez_n = names(gl_entrez)
  gl_entrez = as.numeric(gl_entrez)
  names(gl_entrez) = gl_entrez_n
  # ---------------- KEGG + Reactome (skip dm) ----------------
  if (species != "dm") {

    gsea_kegg <- gseKEGG(
      geneList      = gl_entrez,
      organism      = kegg_org,
      #keyType       = key_type,
      pAdjustMethod = "BH",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.2,
      eps           = 1e-10,
      verbose       = T
    )
    dt_kegg <- as.data.table(gsea_kegg)[,.(DB="KEGG",ID, Description, NES,p.adjust)]
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
    dt_react <- as.data.table(gsea_react)[,.(DB="Reactome",ID, Description, NES,p.adjust)]
    out_list[["Reactome"]] <- dt_react
  }

  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}, by = .(dataset, species)]

# risultati significativi
gsea_results[p.adjust < 0.1 & NES>0]




results[, score_scaled:=scale(score_siglogq), by=.(dataset)]

library(ReactomePA)
gene_list = results[, {
    .(score=max(score_sig_IR_logq, na.rm = TRUE))
}, by = .(dataset,ensembl,species)]

gl = gene_list[species=="mmesc",.(score,ensembl)]
gl_ens = gl$score
names(gl_ens) = gl$ensembl
gl_ens <- sort(gl_ens, decreasing = TRUE)

entrez_map <- mapIds(
    org.Mm.eg.db,
    keys      = gl$ensembl,
    column    = "ENTREZID",
    keytype   = "ENSEMBL",
    multiVals = "first"
)

gl_entrez <- gl$score
names(gl_entrez) <- unname(entrez_map[gl$ensembl])
gl_entrez <- gl_entrez[!is.na(names(gl_entrez)) & names(gl_entrez) != ""]
gl_entrez <- tapply(gl_entrez, names(gl_entrez), max, na.rm=TRUE)
gl_entrez <- sort(gl_entrez, decreasing = TRUE)
gl_entrez_n = names(gl_entrez)
gl_entrez = as.numeric(gl_entrez)
names(gl_entrez) = gl_entrez_n

gsea_bp <- gseGO(
    geneList      = gl_ens,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = 0.2,
    eps           = 1e-10,
    verbose       = FALSE
)
dotplot(gsea_bp, showCategory=15, split=".sign") +
  facet_grid(. ~ .sign) +
  theme_bw(base_size=11)
gsea_react <- gsePathway(
      geneList      = gl_entrez,
      organism      = "mouse",
      pAdjustMethod = "BH",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.5,
      eps           = 1e-10,
      verbose       = FALSE
)

library(enrichplot)
p =dotplot(gsea_react, showCategory=15, split=".sign") +
  facet_grid(. ~ .sign) +
  theme_bw(base_size=11)
ggsave(p, file="go_plot.pdf",width = 10,height = 8)
gseaplot(gsea_react)

gseaplot2(gsea_react, geneSetID = "R-MMU-68877")


dt_react <- as.data.table(gsea_react)
dt_react[NES>0]



res_esc = results[grepl("ESC", dataset)]
genes_long <- dt_react[NES > 0, {
  entrez <- unique(unlist(strsplit(core_enrichment, "/", fixed = TRUE)))

  sym <- mapIds(org.Mm.eg.db,
                keys = entrez,
                column = "SYMBOL",
                keytype = "ENTREZID",
                multiVals = "first")
  sym <- unique(na.omit(unname(sym)))

  # prendi TUTTI gli eventi per quei geni e tieni anche event
  dA <- res_esc[gene_symbol %chin% sym,
                .(event, ensembl, T.obs, gene_symbol, p.value, Sigma)]

  # attacca metadati del termine (ripetuti su tutte le righe)
  dA[, `:=`(ID = ID,
            Description = Description,
            p.adjust = p.adjust,
            NES = NES)]
  dA
}, by = .(ID, Description, p.adjust, NES)]

# esempio: solo termini sig.
genes_long[order(-deltaAIC)][p.value<0.07]

#esempi positivi

mygene="ENSMUSG00000007564.16:chr17:+:21179202:21179272:21179651:21179785"
mygene="ENSMUSG00000025474.10:chr7:-:139584712:139584892:139585183:139585363"
mygene="ENSMUSG00000036672.6:chr8:-:106571892:106572068:106572948:106573023"
mygene="ENSMUSG00000007564.16:chr17:+:21179202:21179272:21179651:21179785"
mygene="ENSMUSG00000015120.17:chr17:-:25484093:25484151:25487515:25487587"
load("rmats_all.rdata")
source("../commons/nested_test.r")
source("../commons/plot.r")
ts_gene <- rmats_all[event==mygene]
ss_gene <- rmats_all[1:3,]
ss_gene$replicate <- 1:3
res_nested = test_sigma_nested(ts_gene, ss_gene, scaling_A=T, usa_W=T,bootstrap_type="wild",
                          wild_dist="rademacher", #rademacher normal mammen
                               add_SS=FALSE, t_star=0, B_n=5000)


evento = unlist(strsplit(ts_gene$event[1],":"))
evento = paste0(evento[2],":",evento[3],":",evento[5],"-",evento[6])
titolo = paste0("(",ts_gene$gene_symbol[1],") ",evento)
sottotitolo = paste0("(nested test p=",signif(res_nested$p.value, 2),", DeltaAIC=",signif(res_nested$deltaAIC, 2),")")

p_fit = plot_fit_two_models_paper(ts_gene, res_nested, t_star=0, show_shutoff = T,
          step=1) + ggtitle(titolo, subtitle = sottotitolo) + 
          theme(plot.title = element_text(size=14, hjust = 0.5),
          plot.subtitle = element_text(size=10, hjust = 0.5))
ggsave("Mis12.pdf", p_fit, width=8, height=5)

