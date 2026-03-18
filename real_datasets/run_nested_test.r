setwd("~/postexport-kinetics/real_datasets")
source("../commons/nested_test.r")


require(data.table)

load("GSE83620/rmats_dt_Kc167.rdata")
dt_kc167 <- copy(rmats_dt)
dt_kc167[, dataset := "Kc167 | GSE83620"]

load("GSE207924/rmats_dt_K562.rdata")
dt_k562 <- copy(rmats_dt)
dt_k562[, dataset := "K562 | GSE207924"]

load("GSE207924/rmats_dt_3T3.rdata")
dt_3t3 <- copy(rmats_dt)
dt_3t3[, dataset := "3T3 | GSE207924"]

load("GSE256335/rmats_dt_ECS.rdata")
dt_esc <- copy(rmats_dt)
dt_esc[, dataset := "mESC | GSE256335"]  

# 2) combina tutto (fill=TRUE gestisce colonne diverse)
rmats_all <- rbindlist(list(dt_kc167, dt_k562, dt_3t3, dt_esc),
                       use.names = TRUE, fill = TRUE)

save(rmats_all, file="rmats_all.rdata")

# metriche per evento (prima del nested test)
evt_qc <- rmats_all[, .(
  cov_total = sum(N + N_s + C + C_s, na.rm=TRUE),
  cov_cyt   = sum(C + C_s, na.rm=TRUE),
  rng_C     = max(C, na.rm=TRUE) - min(C, na.rm=TRUE),
  rng_Cs    = max(C_s, na.rm=TRUE) - min(C_s, na.rm=TRUE),
  nz_C      = sum(C  > 0, na.rm=TRUE),
  nz_Cs     = sum(C_s> 0, na.rm=TRUE),
  frac_nzC  = mean(C  > 0, na.rm=TRUE),
  frac_nzCs = mean(C_s> 0, na.rm=TRUE)
), by=.(dataset, event)]

# soglie relative per dataset (esempio: 20° percentile)
thr <- evt_qc[, .(
  cov_total_min = quantile(cov_total, 0.01, na.rm=TRUE),
  cov_cyt_min   = quantile(cov_cyt,   0.01, na.rm=TRUE),
  rng_C_min     = quantile(rng_C,     0.01, na.rm=TRUE),
  rng_Cs_min    = quantile(rng_Cs,    0.01, na.rm=TRUE)
), by=dataset]

# lambda_var
# R=2  0.7
# R=3  0.5
# R>4  0.3
lambda = c("2"=0.7,"3"=0.5,"4"=0.3,"5"=0.3)

results <- rmats_all[, {

    d <- .SD

    min_pos      <- 3      # già tuo (cyt)
    min_frac     <- 0.3    # già tuo
    min_pos_nuc  <- 3
    max_zero_frac_C  <- 0.7
    max_zero_frac_Cs <- 0.8
    th <- thr[dataset==dataset[1]]
    n_rep = length(unique(d$replicate))

    bad <- (all(d$N == 0) || all(d$N_s == 0) || all(d$C == 0) || all(d$C_s == 0)) ||
       # informatività citoplasma
        (sum(d$C   > 0, na.rm=TRUE) < min_pos) ||
        (sum(d$C_s > 0, na.rm=TRUE) < min_pos) ||
        (mean(d$C  > 0, na.rm=TRUE) < min_frac) ||
        # informatività nucleo
        (sum(d$N   > 0, na.rm=TRUE) < min_pos_nuc) ||
        (sum(d$N_s > 0, na.rm=TRUE) < min_pos_nuc) ||        
        (sum(d$N + d$N_s + d$C + d$C_s, na.rm=TRUE) < th$cov_total_min) ||
        (sum(d$C + d$C_s, na.rm=TRUE) < th$cov_cyt_min) ||
        ((max(d$C,  na.rm=TRUE) - min(d$C,  na.rm=TRUE)) < th$rng_C_min) ||
        ((max(d$C_s,na.rm=TRUE) - min(d$C_s,na.rm=TRUE)) < th$rng_Cs_min) ||
       # zero-inflation esplicita
        (mean(d$C  == 0, na.rm=TRUE) > max_zero_frac_C) ||
        (mean(d$C_s== 0, na.rm=TRUE) > max_zero_frac_Cs)


    if (bad) {} else {
        WWW <- test_sigma_nested(
            tsampled_data = d,
            #ss_data       = d[1:3,],
            scaling_A     = TRUE,
            boot_mode = "parametric_whitened", #"wild_whitened", #"parametric_whitened"
            #bootstrap_type = "wild",
            #wild_dist     = "rademacher", #"rademacher", #mammen
            #usa_W         = TRUE,
            #add_SS        = FALSE,
            B_n           = 5000,
            lambda_var = lambda[n_rep],
            t_star        = 0
         )

        .(p.value = WWW$p.value,
        Sigma   = WWW$Sigma,
        Alpha   = WWW$Alpha,
        T.obs = WWW$T.obs,
        RSS0 = WWW$RSS0,
        RSS1 = WWW$RSS1, 
        IR = WWW$IR,
        #n_obs = WWW$N_obs, 
        #deltaAIC = WWW$deltaAIC,
        gene_symbol = gene_symbol[1],
        description = description[1],
        gene_biotype = gene_biotype[1])
    }
}, by = .(dataset, event, ensembl)]

results[, .(m=.N, ng=length(unique(ensembl)), pmin=min(p.value, na.rm=TRUE)), by=dataset]

a=results[, .(N=.N),by=dataset]
b=rmats_all[, .(N=length(unique(event))),by=dataset]
1-a$N/b$N

eps = 1e-10
results[, q.value := p.adjust(p.value, method="fdr"), 
    by = .(dataset)]
results[, neglogq := -log10(pmax(q.value, eps))]
results[, neglogp := -log10(pmax(p.value, eps))]   # optional diagnostic
# cap evidence to avoid extreme domination
cap <- 6
results[, neglogq_cap := pmin(neglogq, cap)]
results[, neglogp_cap := pmin(neglogp, cap)]

# model-improvement metrics
results[, DeltaRSS := pmax(RSS0 - RSS1, 0)]
results[, IR := fifelse(RSS0 > 0, DeltaRSS / RSS0, 0)]      # improvement ratio in [0,1]
results[, logRSSratio := fifelse(RSS1 > 0, log(RSS0 / RSS1), 0)]

# 3) robust variant (less sensitive to Sigma outliers)
# optional robust scaling for Sigma (prevents outliers dominating)
results[, s0 := median(abs(Sigma), na.rm = TRUE),
 by = .(dataset)] 
results[s0==0,s0:=1]

results[, sig_log := log1p(abs(Sigma+eps) / s0)]

# 1) old-school: Sigma * evidence
results[, score_siglogq := Sigma * neglogp_cap]

# 2) recommended: Sigma * improvement * evidence
results[, score_sig_IR_logq := Sigma * IR * neglogp_cap]

# 3) robust variant (less sensitive to Sigma outliers)
results[, score_siglog_IR_logq := sig_log * IR * neglogp_cap]


summary(results)

save(results,file="results_realdatasets.rdata")
load("results_realdatasets.rdata")
results[order(-score_sig_IR_logq)][1:10]

library(openxlsx)

# 1) tabella completa: solo eventi testati (p.value non NA)
supp_all <- results[!is.na(p.value),
  .(dataset, event, ensembl,
    gene_symbol, gene_biotype,
    T.obs, RSS0, RSS1,score_sig_IR_logq,
    p.value,q.value)
]
# 3) ordina (dataset, p)
setorder(supp_all, dataset, -score_sig_IR_logq)




# top n score_siglogq score_sig_IR_logq
results[IR>0.1 & p.value<0.01 & Sigma>0.01]
results[order(-score_sig_IR_logq)][1:10]

results[order(-score_siglogq)][1:10]

supp_all = results[,{
  .SD[order(-score_sig_IR_logq),.(gene_symbol,gene_biotype,event,description,IR,p.value,Sigma)][IR>0.1 & p.value<0.05]
},by=dataset]

require(openxlsx)
write.xlsx(supp_all, file = "nested_test_results.xlsx", overwrite = TRUE)
