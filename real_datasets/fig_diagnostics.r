setwd("~/postexport-kinetics/real_datasets")
require(data.table)
require(ggplot2)

load("results_realdatasets.rdata",verbose=T)


# =====================
# hist di T.obs
# =====================
# prepara logT
eps <- 1e-8
results[!is.na(T.obs) & T.obs >= 0, logT := log10(T.obs + eps)]

# statistiche per dataset (median/95/99, Pr(T>0), n)
T_stats <- results[!is.na(T.obs) & T.obs >= 0,
                   .(n = .N,
                     pr_pos = mean(T.obs > 0),
                     q50 = quantile(T.obs, 0.50),
                     q95 = quantile(T.obs, 0.95),
                     q99 = quantile(T.obs, 0.99)),
                   by = dataset]

T_stats[, `:=`(
  v_med = log10(q50 + eps),
  v_95  = log10(q95 + eps),
  v_99  = log10(q99 + eps),
  ann_txt = sprintf("n=%d\nPr(T>0)=%.3f\nmedian=%.3g\n95%%=%.3g\n99%%=%.3g",
                    n, pr_pos, q50, q95, q99)
)]

# posizione label per dataset (in coordinate del pannello)
lab_pos_T <- results[!is.na(logT),
                     .(x_lab = quantile(logT, 0.70),
                       y_lab = Inf),
                     by = dataset]

T_stats <- lab_pos_T[T_stats, on="dataset"]

p_T_facet <- ggplot(results[!is.na(logT)],
                    aes(x = logT)) +
  geom_histogram(bins = 60, linewidth = 0.2, fill = "grey80", color = "grey30") +
  geom_vline(data = T_stats, aes(xintercept = v_med), linetype = 2, linewidth = 0.5) +
  geom_vline(data = T_stats, aes(xintercept = v_95),  linetype = 3, linewidth = 0.5) +
  geom_vline(data = T_stats, aes(xintercept = v_99),  linetype = 1, linewidth = 0.5) +
  geom_label(data = T_stats,
             aes(x = x_lab, y = y_lab, label = ann_txt),
             hjust = 0, vjust = 1.1,
             size = 3.1, linewidth = 0.25) +
  facet_wrap(~ dataset, scales = "free_y") +
  labs(
    title = "Distribution of fit improvement statistic by dataset",
    subtitle = expression(T[obs] == RSS[0] - RSS[1] ~ " (full vs null)"),
    x = expression(log[10](T[obs] + epsilon)),
    y = "Count"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_T_facet

# =====================
# plot deltaAIC
# =====================

#results[!is.na(RSS0) & !is.na(RSS1) & RSS0 > 0 & RSS1 > 0,
#        deltaAIC := n_obs * log(RSS0 / RSS1) - 2]

p_simple_aic <- ggplot(results,
                       aes(x = Sigma, y = deltaAIC)) +
  geom_hline(yintercept = 0,  linetype = 2, linewidth = 0.4) +
  geom_hline(yintercept = 2,  linetype = 3, linewidth = 0.4) +
  geom_hline(yintercept = 10, linetype = 3, linewidth = 0.4) +
  geom_point(size = 0.7, alpha = 0.6) +
  facet_wrap(~ dataset, ncol = 2, scales = "free_x") +
  labs(
    x = expression(hat(sigma)[c]),
    y = expression(Delta*AIC),
    #title = expression(paste(Delta, "AIC vs ", hat(sigma)[c], " across datasets"))
  ) +
  theme_bw(base_size = 11) +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

p_simple_aic

ggsave(p_simple_aic, filename = "aic_by_dataset.pdf",
  width    = 10,  
  height   = 8
)

# =======================
# q-q plot
# =======================

# ---------- 1) prepara tabella QQ per dataset ----------
qq_dt <- results[
  !is.na(p.value) & p.value > 0 & p.value <= 1,
  .(p = p.value),
  by = dataset
]

setorder(qq_dt, dataset, p)

# expected p per rank (per dataset)
qq_dt[, exp_p := (seq_len(.N) - 0.5) / .N, by = dataset]
qq_dt[, obs_logp := -log10(p)]
qq_dt[, exp_logp := -log10(exp_p)]

# banda 95% per dataset (Beta order-statistics)
qq_dt[, i := seq_len(.N), by = dataset]
qq_dt[, m := .N, by = dataset]
alpha <- 0.05
qq_dt[, lo := qbeta(alpha/2, i, m - i + 1)]
qq_dt[, hi := qbeta(1 - alpha/2, i, m - i + 1)]
qq_dt[, lo_logp := -log10(lo)]
qq_dt[, hi_logp := -log10(hi)]

# ---------- 2) label per facet: n, n(p<0.05), n(p<0.01) ----------
lab_dt <- qq_dt[, .(
  n = .N,
  n_p05 = sum(p < 0.05),
  n_p01 = sum(p < 0.01),
  x_lab = quantile(exp_logp, 0.05),   # posizione label (in coordinate pannello)
  y_lab = quantile(obs_logp, 0.95)
), by = dataset]

lab_dt[, label := sprintf("n=%d\np<0.05: %d\np<0.01: %d", n, n_p05, n_p01)]

# (Opzionale) aggiungi anche n(ΔAIC>10) se deltaAIC è già in results:
# aic_lab <- results[!is.na(deltaAIC), .(n_dAIC10 = sum(deltaAIC > 10)), by=dataset]
# lab_dt <- aic_lab[lab_dt, on="dataset"]
# lab_dt[, label := sprintf("%s\nΔAIC>10: %d", label, n_dAIC10)]

# ---------- 3) QQ plot facet "paper-ready" ----------
p_qq_facet <- ggplot(qq_dt, aes(x = exp_logp, y = obs_logp)) +
  geom_ribbon(aes(ymin = lo_logp, ymax = hi_logp), alpha = 0.2) +
  geom_point(size = 0.6, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
  geom_label(
    data = lab_dt,
    aes(x = x_lab, y = y_lab, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = -2,
    size = 3.1,
    linewidth = 0.25
  ) +
  facet_wrap(~ dataset, ncol = 2, scales = "fixed") +
  labs(
    #title = "QQ-plots of nested-test p-values by dataset",
    x = "Expected -log10(p) under Uniform(0,1)",
    y = "Observed -log10(p)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )

p_qq_facet
ggsave(p_qq_facet, filename = "qq_by_dataset.pdf",
  width    = 10,   # inches: 2x2 facet leggibile
  height   = 8
)


# =======================
# Bias audit (ΔAIC vs coverage + C-zero fraction)
# =======================

# if not exists run_nested_test.r
load("rmats_all.rdata")

# 1) frazioni di zeri e coverage proxy per evento
cov_dt <- rmats_all[, .(
  fracC0  = mean(C == 0, na.rm=TRUE),
  fracCs0 = mean(C_s == 0, na.rm=TRUE),
  meanC   = mean(C, na.rm=TRUE),
  meanCs  = mean(C_s, na.rm=TRUE),
  meanN   = mean(N, na.rm=TRUE),
  meanNs  = mean(N_s, na.rm=TRUE),
  # coverage proxy: somma media delle 4 specie (puoi cambiarla)
  mean_total = mean(N + N_s + C + C_s, na.rm=TRUE)
), by = .(dataset, event)]

# 2) unisci a results
setkeyv(results, c("dataset", "event"))
setkeyv(cov_dt,  c("dataset", "event"))

res_cov <- cov_dt[results]   # mantiene tutte le righe di results


res_cov[, log_cov := log10(mean_total + 1)]

p_cov <- ggplot(res_cov[!is.na(deltaAIC) & !is.na(log_cov)],
                aes(x = log_cov, y = deltaAIC)) +
  geom_hline(yintercept = 0, linetype=2, linewidth=0.35) +
  geom_hline(yintercept = 2, linetype=3, linewidth=0.35) +
  geom_hline(yintercept = 10, linetype=3, linewidth=0.35) +
  geom_point(size=0.6, alpha=0.5) +
  geom_smooth(method="loess", se=TRUE, linewidth=0.6) +
  facet_wrap(~dataset, ncol=2, scales="free_x") +
  theme_bw(base_size = 11) +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  ) +
  labs(#title = expression(paste(Delta,"AIC vs coverage (near-null audit)")),
       x = "log10(mean total abundance + 1)",
       y = expression(Delta*AIC))

p_cov
ggsave(p_cov, filename = "deltaAIC_abdundance.pdf",
  width    = 10,   # inches: 2x2 facet leggibile
  height   = 8
)

# ===============
# ΔAIC vs zero-inflation in C/Cs
# ===============

p_zero <- ggplot(res_cov[!is.na(deltaAIC)],
                 aes(x = fracC0, y = deltaAIC)) +
  geom_hline(yintercept = 0, linetype=2, linewidth=0.35) +
  geom_hline(yintercept = 2, linetype=3, linewidth=0.35) +
  geom_hline(yintercept = 10, linetype=3, linewidth=0.35) +
  geom_point(size=0.6, alpha=0.5) +
  geom_smooth(method="loess", se=TRUE, linewidth=0.6) +
  facet_wrap(~dataset, ncol=2, scales="fixed") +
  theme_bw(base_size = 11) +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )+
  labs(#title = expression(paste(Delta,"AIC vs zero fraction in cytoplasmic unspliced (C)")),
       x = "Fraction of timepoints with C=0",
       y = expression(Delta*AIC))

p_zero
ggsave(p_zero, filename = "deltaAIC_inflactionC.pdf",
  width    = 10,   # inches: 2x2 facet leggibile
  height   = 8
)


topN_gene <- 10

tab_top_genes <- results[!is.na(deltaAIC)][
  order(p.value)
][, .SD[1], by = .(dataset, ensembl)]  # prende il primo dopo l'order: max ΔAIC per gene

tab_top_genes <- tab_top_genes[
  order(p.value)
][, head(.SD, topN_gene), by = dataset,
  .SDcols = c("ensembl","gene_symbol","gene_biotype","deltaAIC","Sigma","p.value","event")]

tab_top_genes_print <- copy(tab_top_genes)
tab_top_genes_print[, `:=`(
  deltaAIC = round(deltaAIC, 2),
  Sigma    = signif(Sigma, 3),
  p.value  = signif(p.value, 3)
)]

kable(tab_top_genes_print, format="latex", booktabs=TRUE,
      caption=paste0("Top ", topN_gene, " genes per dataset ranked by maximum event-level $\\Delta$AIC."))
results[deltaAIC > 2 & p.value < 0.05]
results[deltaAIC > 10]


