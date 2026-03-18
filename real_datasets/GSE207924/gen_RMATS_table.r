
setwd("~/postexport-kinetics/real_datasets/GSE207924")

library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(biomaRt)

cells = c("K562", "3T3")

orgdb = c("3T3" = org.Mm.eg.db, "K562" = org.Hs.eg.db)
ensdb = c("3T3" = "mmusculus_gene_ensembl", "K562" = "hsapiens_gene_ensembl")


for (cc in cells) {
    pattern <- paste0("^", cc, "_(nuclear|cytoplasm)_WT_rep\\d+\\.tsv$")
    files <- list.files(path = ".", pattern=pattern, full.names = TRUE)

    dt_all <- rbindlist(lapply(files, function(f) {
        d <- fread(f)
        mean_cols <- grep("Mean", names(d), value = TRUE)
        # metadati dal nome file
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

        base_cols <- intersect(c("Gene", "Symbol"), names(d))
        d <- d[, c(base_cols, mean_cols), with = FALSE]
        tp_levels <- c("t0","t15","t30","t60","t120")
        setnames(
            d,
            old = mean_cols,
            new = tp_levels[seq_along(mean_cols)]
        )

        d[, `:=`(compartment = fraction, replicate = replicate)]

        d
    }), use.names = TRUE, fill = TRUE)

    dt_long <- melt(
        dt_all,
        id.vars = c("Gene","Symbol","compartment","replicate"),
        measure.vars = c("t0","t15","t30","t60","t120"),
        variable.name = "timepoint",
        value.name = "new_fraction"
    )
    setnames(dt_long, "Gene", "gene")
    map <- c(
        nuclear = "nuc",
        cytoplasm = "cyt"
    )
    dt_long[, compartment := map[compartment]]
    dt_long[, replicate := paste0("r",replicate)]


    #d_rmats = fread(paste0("rmats_gse207924_",cc,".csv"))
    d_rmats = fread(paste0("rmats_gse207924_",cc,"_normalized.csv"))

    d_rmats[, sample_id := paste(timepoint, compartment, replicate, sep = "_")]

    libDT <- d_rmats[, .(libsize = sum(count, na.rm = TRUE)), by = sample_id]
    libDT[, lib_factor := libsize / median(libsize)]

    ## 3) join del fattore e normalizzazione counts (a livello replicate)
    d_rmats <- libDT[d_rmats, on = "sample_id"]
    d_rmats[, count_norm := count / lib_factor]

    ## 5) Join new_fraction (deve essere 1:1 su gene+timepoint+compartment+replicate)
    # tieni solo colonne necessarie e rendi univoca la chiave
    dt_long2 <- unique(dt_long[, .(gene, timepoint, compartment, replicate, new_fraction)])

    # check duplicati (se >0, c'è un problema nei dati di dt_long)
    dt_long2[duplicated(dt_long2, by = c("gene","timepoint","compartment","replicate"))]

    setkey(d_rmats, gene, timepoint, compartment, replicate)
    setkey(dt_long2, gene, timepoint, compartment, replicate)

    # aggiungi new_fraction a d_rmats (left join: mantieni tutte le righe rMATS)
    d_rmats <- dt_long2[d_rmats]

    # filtra NA se vuoi solo righe con new_fraction
    d_rmats <- d_rmats[!is.na(new_fraction)]


    d_rmats[, event := paste(gene, chr, strand,
                        upstreamES, upstreamEE, downstreamES, downstreamEE,
                        sep=":")]
    # fix time and replicate variables
    d_rmats[, time := as.numeric(sub("^t", "", timepoint))]
    d_rmats[, replicate := as.integer(sub("^r", "", replicate))]

    d_rmats[, State := fifelse(compartment=="nuc" & version=="inclusion", "N",
                 fifelse(compartment=="nuc" & version=="skipping", "N_s",
                 fifelse(compartment=="cyt" & version=="inclusion", "C",
                 fifelse(compartment=="cyt" & version=="skipping", "C_s", NA_character_))))]

    agg = d_rmats[, .(Value = (1-new_fraction)*sum(count, na.rm=TRUE)), #(1-new_fraction)*
            by = .(event, time, replicate, State)]


    rmats_dt = dcast(agg, event + time + replicate ~ State, value.var="Value", fill=0)
    rmats_dt[,ensembl := sub(":.*", "", event)]
    rmats_dt[,ensembl := sub("\\..*", "", ensembl)]

# ======================
# get metadata from ensembl
# ======================

    gene_ids <- unique(rmats_dt$ensembl)
    mart <- useEnsembl(biomart="genes",dataset=ensdb[[cc]], mirror = "useast")
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
    
    save(rmats_dt,file=paste0("rmats_dt_",cc,".rdata"))
}
