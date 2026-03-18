
setwd("~/postexport-kinetics/real_datasets/GSE256335")
# dataset Mouse con shutoff


require(data.table)
library(org.Mm.eg.db)
library(biomaRt)

#d_rmats = fread("rmats_counts.csv")
d_rmats = fread("gse256335_mouse_normalized.csv")

d_rmats[, sample_id := paste(timepoint, compartment, replicate, sep = "_")]

libDT <- d_rmats[, .(libsize = sum(count, na.rm = TRUE)), by = sample_id]
libDT[, lib_factor := libsize / median(libsize)]

    ## 3) join del fattore e normalizzazione counts (a livello replicate)
d_rmats <- libDT[d_rmats, on = "sample_id"]
d_rmats[, count_norm := count / lib_factor]

d_rmats[, event := paste(gene, chr, strand,
                        upstreamES, upstreamEE, downstreamES, downstreamEE,
                        sep=":")]
# fix time and replicate variables
d_rmats[, time := as.numeric(sub("_",".",sub("^t", "", timepoint)))]
d_rmats[, replicate := as.integer(sub("^r", "", replicate))]

d_rmats[, State := fifelse(compartment=="nuc" & version=="inclusion", "N",
                 fifelse(compartment=="nuc" & version=="skipping", "N_s",
                 fifelse(compartment=="cyt" & version=="inclusion", "C",
                 fifelse(compartment=="cyt" & version=="skipping", "C_s", NA_character_))))]

agg = d_rmats[, .(Value = sum(count, na.rm=TRUE)),
            by = .(event, time, replicate, State)]


rmats_dt = dcast(agg, event + time + replicate ~ State, value.var="Value", fill=0)

rmats_dt[,ensembl := sub(":.*", "", event)]
rmats_dt[,ensembl := sub("\\..*", "", ensembl)]


# ======================
# get metadata from ensembl
# ======================
gene_ids <- unique(rmats_dt$ensembl)
mart <- useEnsembl(biomart="genes",dataset="mmusculus_gene_ensembl", mirror = "useast")
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

rmats_dt[gene_info, on = .(ensembl = ensembl_gene_id),
        `:=`(
          gene_symbol = i.external_gene_name,
          description = i.description,
          gene_biotype = i.gene_biotype
        )]
    
save(rmats_dt,file=paste0("rmats_dt_ECS.rdata"))
