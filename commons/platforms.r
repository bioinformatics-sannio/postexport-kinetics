
library(MASS)


add_gaussian_noise <- function(df, cols, noise_sd) {
  for (col in cols) {
    df[[col]] <- df[[col]] + stats::rnorm(
      nrow(df),
      mean = 0,
      sd   = pmax(0, df[[col]]) * noise_sd
    )
  }
  return(df)
}


sample_dispersion_gamma <- function(n_genes, mean_disp, cv_disp) {
  shape <- 1 / (cv_disp^2)
  scale <- mean_disp * cv_disp^2
  rgamma(n_genes, shape = shape, scale = scale)
}

simulate_rnaseq <- function(
  dt,
  targets       = c("N","N_s","C","C_s"),
  scale_counts  = 200,    # come 10^3, 10^2, 50, ecc.
  mean_disp=0.1, cv_disp=0.5
){
  dt <- copy(as.data.table(dt))
  stopifnot(all(targets %in% names(dt)))
  
  for (tg in targets) {
    latent <- pmax(dt[[tg]], 0)
    latent[is.na(latent)] <- 0
    
    # media di conteggio attesa
    mu_counts <- latent * scale_counts
    
    disp = sample_dispersion_gamma(1, mean_disp, cv_disp) 
    # simulazione
    theta <- if (disp > 0) 1/disp else 1e6
    counts <- ifelse(
        mu_counts > 0,
        MASS::rnegbin(length(mu_counts), mu = mu_counts, theta = theta),
        0L
    )
    
    # normalizziamo i valori per riportarli in scala lineare ~latent
    dt[[tg]] <- counts / scale_counts
  }
  
  dt[]
}

simulate_rt_qpcr <- function(dt,
  targets      = c("N", "C", "C_s", "N_s"),
  scale_copies = 10,   # copie per unità di stato latente
  ct_intercept = 35,   # Ct per ~1 copia
  ct_sd        = 0.25, # rumore tecnico qPCR
  lod_ct       = 40    # limite detection
){
  dt <- copy(as.data.table(dt))
  stopifnot(all(targets %in% names(dt)))
  
  for (tg in targets) {
    
    ## 1) lambda = numero medio di copie generate dalle molecole simulate dall’ODE
    lambda <- pmax(dt[[tg]], 0) * scale_copies
    lambda[is.na(lambda)] <- 0
    
    ## 2) rumore di campionamento (Poisson)
    copies <- rpois(nrow(dt), lambda)
    copies[copies < 1] <- 1L   # evita log2(0)
    
    ## 3) Ct ideale + rumore tecnico qPCR
    Ct <- ct_intercept - log2(copies)
    Ct <- Ct + rnorm(nrow(dt), 0, ct_sd)
    
    ## 4) limite di detection
    Ct[Ct > lod_ct] <- lod_ct + 0.5
    #Ct[Ct > lod_ct] <- NA_real_
    
    ## 5) sovrascrivo direttamente la colonna latente
    dt[[tg]] <- 2^(-Ct)
  }
  
  dt[]
}
