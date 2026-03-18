setwd("~/postexport-kinetics/synthetic_dataset")
source("../commons/nested_test.r")
source("../commons/psi_test.r")
source("../commons/platforms.r")


require(data.table)
library(parallel)

load("ode_states_2k_20p.rdata")

range_platform <- c("RT-qPCR", "GAUSS", "RNA-seq")
range_perturbation <- c("NONE", "SHUTOFF")
range_technical_noise <- c("Very low", "Low", "Medium", "High") # technical noise simulation
range_lambda <-  c(0,0.3,0.5,0.7)
scaleA <- T #c(TRUE, FALSE)
#test_method <- c("Nested-parametric", "Nested-wild-rademacher", "Nested-wild-mammen", "Nested-wild-normal", "Baseline", "OneSample")
#test_method <- c("Nested-wild-rademacher", "PSIBaseline", "Nested-wild-mammen")
#test_method <- c("PSIBaseline")
test_method <- c("Nested-parametric_whitened", "Nested-wild_whitened", "PSIBaseline")


range_gauss_noise <- c(0.02, 0.05, 0.1, 0.2) # technical noise simulation (for GAUSS only)
names(range_gauss_noise) <- range_technical_noise


range_n_time_samples <- unique(ode_states$N_time_samples)
range_n_replicates <- unique(ode_states$N_replicates)
range_tsteps <- 10 #unique(ode_states$T_step)

genes <- unique(ode_states$Gene)

run_for_gene <- function(gene_i,N_boot=5000) {
  list_results_gene <- list()
  ki <- 1

  for (time_samples_i in range_n_time_samples) {
    for (perturbation_i in range_perturbation) {
      for (tstep_i in range_tsteps) {
        for (n_replicates_i in range_n_replicates) {
          for (platform_i in range_platform) {
            for (noise_i in range_technical_noise) {
              for (scaling_A_i in scaleA) {
                for (lambda_i in range_lambda) {
                    for (test_method_i in test_method) {
                      s_test_i <- unlist(strsplit(test_method_i, "-"))
                      test_i <- s_test_i[1]
                      b_type_i <- if (length(s_test_i) > 1) s_test_i[2]
                      #wild_type_i <- if (length(s_test_i) > 1) s_test_i[3]

                      # sottoinsieme per questo gene e scenario
                      dt_results_sub <- ode_states[
                        Gene == gene_i &
                          N_replicates == n_replicates_i &
                          N_time_samples == time_samples_i &
                          T_step == tstep_i
                      ]

                      # se per qualche motivo non esiste questa combinazione, salta
                      if (nrow(dt_results_sub) == 0L) next

                      dt_results_sub_timecourse <- dt_results_sub[Perturbation == perturbation_i]
                      dt_results_sub_steadystate <- dt_results_sub[Perturbation == "SteadyState"]

                      if (nrow(dt_results_sub_timecourse) == 0L ||
                        nrow(dt_results_sub_steadystate) == 0L) {
                        next
                      }

                      cols_to_noise <- c("N", "C", "C_s", "N_s")

                      #--------------------------------
                      # Simulazione rumore per piattaforma
                      #--------------------------------
                      if (platform_i == "GAUSS") {
                        noise_sd_i <- range_gauss_noise[noise_i]
                        dt_results_sub_timecourse <- add_gaussian_noise(dt_results_sub_timecourse, cols_to_noise, noise_sd_i)
                        dt_results_sub_steadystate <- add_gaussian_noise(dt_results_sub_steadystate, cols_to_noise, noise_sd_i)
                      } else if (platform_i == "RT-qPCR") {
                        cd_sd_i <- switch(noise_i,
                          "Very low" = 0.01,
                          "Low"      = 0.05,
                          "Medium"   = 0.1,
                          "High"     = 0.25
                        )
                        scale_copies_i <- switch(noise_i,
                          "Very low" = 20,
                          "Low"      = 15,
                          "Medium"   = 10,
                          "High"     = 5
                        )
                        dt_results_sub_timecourse <- simulate_rt_qpcr(
                          dt_results_sub_timecourse,
                          targets      = cols_to_noise,
                          ct_sd        = cd_sd_i,
                          scale_copies = scale_copies_i
                        )
                        dt_results_sub_steadystate <- simulate_rt_qpcr(
                          dt_results_sub_steadystate,
                          targets      = cols_to_noise,
                          ct_sd        = cd_sd_i,
                          scale_copies = scale_copies_i
                        )
                      } else if (platform_i == "RNA-seq") {
                        scale_counts_i <- switch(noise_i,
                          "Very low" = 10000,
                          "Low" = 5000,
                          "Medium" = 1000,
                          "High" = 200
                        )
                        mean_disp_i <- switch(noise_i,
                          "Very low" = 0.01,
                          "Low" = 0.05,
                          "Medium" = 0.10,
                          "High" = 0.25
                        ) * 0.25
                        cv_disp_i <- switch(noise_i,
                          "Very low" = 0.5,
                          "Low" = 0.7,
                          "Medium" = 0.8,
                          "High" = 1.0
                        )
                        dt_results_sub_timecourse <- simulate_rnaseq(
                          dt_results_sub_timecourse,
                          targets = cols_to_noise,
                          scale_counts = scale_counts_i,
                          mean_disp = mean_disp_i, cv_disp = cv_disp_i
                        )
                        dt_results_sub_steadystate <- simulate_rnaseq(
                          dt_results_sub_steadystate,
                          targets = cols_to_noise,
                          scale_counts = scale_counts_i,
                          mean_disp = mean_disp_i, cv_disp = cv_disp_i
                        )
                      }

                      #--------------------------------
                      # Test statistico
                      #--------------------------------
                      if (test_i == "Nested") {
                        WW <- test_sigma_nested(
                          tsampled_data = dt_results_sub_timecourse,
                          #ss_data = dt_results_sub_steadystate,
                          scaling_A = scaling_A_i,
                          #usa_W = usa_W_i,
                          #add_SS = usa_SS_i,
                          boot_mode = b_type_i, #"parametric_whitened",
                          lambda_var = lambda_i,
                          B_n = N_boot,
                          #bootstrap_type = b_type_i,
                          #wild_dist = wild_type_i,
                          t_star = if (perturbation_i == "NONE") NULL else unique(dt_results_sub_timecourse$T_star)
                        )
                      } else if (test_i == "PSIBaseline") {
                        if (lambda_i>0) next
                        dt_results_sub_timecourse[, C := C+1]
                        dt_results_sub_timecourse[, C_s := C_s+1]
                        WW <- baseline_psi_logit_ttest_pre_vs_post(
                          tsampled_data = dt_results_sub_timecourse,
                          t_star = NULL #if (perturbation_i == "NONE") NULL else unique(dt_results_sub_timecourse$T_star)
                        )
                      } else {
                        next
                      }

                      r_i <- data.frame(
                        Gene         = gene_i,
                        N_replicates = n_replicates_i,
                        N_tsamples   = time_samples_i,
                        Lambda = lambda_i,
                        Tsteps       = tstep_i,
                        Perturbation = perturbation_i,
                        Exprs_noise  = noise_i,
                        Platform     = platform_i,
                        Scaling_A    = scaling_A_i,
                        Test_method  = test_method_i,
                        Positive     = unique(dt_results_sub$truth_pos),
                        p.value      = WW$p.value,
                        T.obs        = WW$T.obs,
                        Sigma        = WW$Sigma,
                        Alpha        = WW$Alpha,
                        RSS0 = WW$RSS0,
                        RSS1 = WW$RSS1, 
                        IR = WW$IR,
                        #deltaAIC     = WW$deltaAIC,
                        sigma_true_mean = mean(dt_results_sub_steadystate$Base_sigma_c),
                        sigma_true_sd = sd(dt_results_sub_steadystate$Base_sigma_c),
                        sigma_true_median = median(dt_results_sub_steadystate$Base_sigma_c),
                        sigma_true_min = min(dt_results_sub_steadystate$Base_sigma_c),
                        sigma_true_max = max(dt_results_sub_steadystate$Base_sigma_c)
                      )

                      list_results_gene[[ki]] <- r_i
                      ki <- ki + 1
                    } # test_i
                } # usa_W_i
              } # scaling_A_i
            } # noise_i
          } # platform_i
        } # n_replicates_i
      } # tstep_i
    } # perturbation_i
  } # time_samples_i

  if (length(list_results_gene) == 0L) {
    return(NULL)
  } else {
    return(list_results_gene)
  }
}

#--------------------------------------------------
# Esecuzione parallela sui geni (Linux/macOS)
#--------------------------------------------------
n_cores <- max(1L, detectCores() - 90)

currentTs <- Sys.time()
run_for_gene(1)
elapsed <- Sys.time() - currentTs
elapsed*length(genes)/n_cores


res_par <- mclapply(genes, run_for_gene, mc.cores = n_cores)



# Unisci tutti i risultati
res_par <- Filter(Negate(is.null), res_par)
res_flat <- unlist(res_par, recursive = FALSE)
length(res_flat) # numero totale di test eseguiti
dt_test_results <- data.table::rbindlist(res_flat, use.names = TRUE, fill = TRUE)

dt_test_results$Exprs_noise <- factor(dt_test_results$Exprs_noise, levels = c("Very low", "Low", "Medium", "High"), ordered=TRUE)
dt_test_results$Perturbation <- factor(dt_test_results$Perturbation, levels = c("NONE", "SHUTOFF"))
dt_test_results$Platform <- factor(dt_test_results$Platform, levels = c("RT-qPCR", "GAUSS", "RNA-seq"))
dt_test_results$Test_method <- factor(dt_test_results$Test_method)


eps = 1e-10
dt_test_results[, q.value := p.adjust(p.value, method="fdr"), 
    by = .(Lambda,N_replicates, Tsteps, Test_method, N_tsamples, Platform, Exprs_noise, Perturbation)]
dt_test_results[, neglogq := -log10(pmax(q.value, eps))]
dt_test_results[, neglogp := -log10(pmax(p.value, eps))]   # optional diagnostic
# cap evidence to avoid extreme domination
cap <- 6
dt_test_results[, neglogq_cap := pmin(neglogq, cap)]

# model-improvement metrics
dt_test_results[, DeltaRSS := pmax(RSS0 - RSS1, 0)]
dt_test_results[, IR := fifelse(RSS0 > 0, DeltaRSS / RSS0, 0)]      # improvement ratio in [0,1]
dt_test_results[, logRSSratio := fifelse(RSS1 > 0, log(RSS0 / RSS1), 0)]

# 3) robust variant (less sensitive to Sigma outliers)
# optional robust scaling for Sigma (prevents outliers dominating)
dt_test_results[, s0 := median(abs(Sigma), na.rm = TRUE),
 by = .(Lambda,N_replicates, Tsteps, Test_method, N_tsamples, Platform, Exprs_noise, Perturbation)] 
dt_test_results[s0==0,s0:=1]

dt_test_results[, sig_log := log1p(abs(Sigma+eps) / s0)]

# 1) old-school: Sigma * evidence
dt_test_results[, score_siglogq := Sigma * neglogq_cap]

# 2) recommended: Sigma * improvement * evidence
dt_test_results[, score_sig_IR_logq := Sigma * IR * neglogq_cap]

# 3) robust variant (less sensitive to Sigma outliers)
dt_test_results[, score_siglog_IR_logq := sig_log * IR * neglogq_cap]

summary(dt_test_results)


save(dt_test_results, file = "test_sigma_2k_new.rdata")






require(plotROC)
gg <- ggplot(dt_test_results, aes(d = Positive, m = Score, color = Test_method)) +
  geom_roc(n.cuts = 0, size = 0.5) +
  facet_grid(Exprs_noise ~ N_replicates, labeller = label_both)
calc_auc(gg)
