# =============================================================================
# Title: Measurement Noise and Assay Simulation Utilities
# Description:
#   This script provides helper functions to perturb latent continuous signals
#   and to simulate simplified RNA-seq and RT-qPCR measurements from model-based
#   latent quantities.
#
# Included utilities:
#   - additive Gaussian perturbation with signal-dependent scale,
#   - gamma sampling for dispersion parameters,
#   - negative-binomial count simulation for RNA-seq-like data,
#   - Poisson + Ct transformation for RT-qPCR-like data.
#
# Intended use:
#   Research code accompanying the paper and shared for transparency/review.
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


library(MASS)


# -----------------------------------------------------------------------------
# Add multiplicative-scale Gaussian noise to selected columns.
#
# For each selected column:
#   - noise is sampled from a Gaussian distribution with mean 0,
#   - the standard deviation is proportional to the non-negative signal level,
#   - the perturbed value is written back into the data frame.
#
# This provides a simple way to mimic measurement variability that increases
# with signal magnitude.
# -----------------------------------------------------------------------------
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


# -----------------------------------------------------------------------------
# Sample dispersion values from a Gamma distribution.
#
# The Gamma parameters are chosen so that the sampled dispersions have:
#   - approximate mean = mean_disp
#   - approximate coefficient of variation = cv_disp
#
# This is useful for generating heterogeneous dispersion values across genes
# or simulated features.
# -----------------------------------------------------------------------------
sample_dispersion_gamma <- function(n_genes, mean_disp, cv_disp) {
  shape <- 1 / (cv_disp^2)
  scale <- mean_disp * cv_disp^2
  rgamma(n_genes, shape = shape, scale = scale)
}


# -----------------------------------------------------------------------------
# Simulate RNA-seq-like observations from latent continuous abundances.
#
# For each target variable:
#   1. latent values are truncated below at zero,
#   2. expected counts are obtained by multiplying by scale_counts,
#   3. one dispersion value is sampled from a Gamma distribution,
#   4. counts are drawn from a Negative Binomial distribution,
#   5. counts are rescaled back to the original continuous scale.
#
# Notes:
#   - scale_counts controls sequencing depth / count resolution.
#   - mean_disp and cv_disp control overdispersion.
#   - a single dispersion value is sampled per target column in each call.
# -----------------------------------------------------------------------------
simulate_rnaseq <- function(
  dt,
  targets       = c("N","N_s","C","C_s"),
  scale_counts  = 200,
  mean_disp=0.1, cv_disp=0.5
){
  dt <- copy(as.data.table(dt))
  stopifnot(all(targets %in% names(dt)))
  
  for (tg in targets) {
    # Enforce non-negative latent abundances.
    latent <- pmax(dt[[tg]], 0)
    latent[is.na(latent)] <- 0
    
    # Expected count scale for the target.
    mu_counts <- latent * scale_counts
    
    # Sample one dispersion value for the current target.
    disp = sample_dispersion_gamma(1, mean_disp, cv_disp)
    
    # Convert dispersion to theta parameterization used by rnegbin().
    # Large theta corresponds to weaker overdispersion.
    theta <- if (disp > 0) 1/disp else 1e6
    
    # Draw counts; zero expected mean gives zero observed counts.
    counts <- ifelse(
        mu_counts > 0,
        MASS::rnegbin(length(mu_counts), mu = mu_counts, theta = theta),
        0L
    )
    
    # Rescale counts back to a continuous scale comparable to the latent input.
    dt[[tg]] <- counts / scale_counts
  }
  
  dt[]
}


# -----------------------------------------------------------------------------
# Simulate RT-qPCR-like observations from latent continuous abundances.
#
# For each target variable:
#   1. latent abundance is converted into an expected copy number,
#   2. observed copies are sampled with Poisson noise,
#   3. copy number is transformed into Ct space,
#   4. Gaussian technical noise is added in Ct space,
#   5. a limit of detection is applied,
#   6. Ct is mapped back to a positive signal scale via 2^(-Ct).
#
# Notes:
#   - scale_copies controls the expected number of molecules/copies.
#   - ct_intercept defines the Ct corresponding to approximately one copy.
#   - ct_sd controls technical variability.
#   - lod_ct is the detection threshold in Ct units.
# -----------------------------------------------------------------------------
simulate_rt_qpcr <- function(dt,
  targets      = c("N", "C", "C_s", "N_s"),
  scale_copies = 10,
  ct_intercept = 35,
  ct_sd        = 0.25,
  lod_ct       = 40
){
  dt <- copy(as.data.table(dt))
  stopifnot(all(targets %in% names(dt)))
  
  for (tg in targets) {
    
    # Convert latent abundance into expected copy number.
    lambda <- pmax(dt[[tg]], 0) * scale_copies
    lambda[is.na(lambda)] <- 0
    
    # Poisson sampling noise on molecular copies.
    copies <- rpois(nrow(dt), lambda)
    
    # Enforce at least one copy to avoid log2(0) in Ct conversion.
    copies[copies < 1] <- 1L
    
    # Ideal Ct transformation plus Gaussian technical noise.
    Ct <- ct_intercept - log2(copies)
    Ct <- Ct + rnorm(nrow(dt), 0, ct_sd)
    
    # Apply a simple limit of detection.
    # Values above the threshold are censored just beyond the LOD.
    Ct[Ct > lod_ct] <- lod_ct + 0.5
    # Alternative option:
    # Ct[Ct > lod_ct] <- NA_real_
    
    # Map Ct back to a positive signal-like quantity.
    dt[[tg]] <- 2^(-Ct)
  }
  
  dt[]
}