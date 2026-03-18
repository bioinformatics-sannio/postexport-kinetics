# Post-export RNA kinetic  
### Predict post-export RNA processing with a nested ODE kinetic test

Compartment-resolved kinetic inference and **nested hypothesis testing** for detecting **post-export RNA processing** from time-resolved nuclear/cytoplasmic transcriptomic data.

This repository contains the R implementation of the modeling and inference framework described in:

> Faretra L., Napolitano F., Pancione M., Cerulo L. (2026). *A kinetic modeling framework to detect post-export RNA processing from time-resolved transcriptomic data.* Manuscript.

---

## Overview

Many transcriptomic analyses assume that RNA processing (e.g., splicing) occurs exclusively in the nucleus. However, increasing evidence suggests that **RNA processing and remodeling may continue after nuclear export**.

This repository provides:

- A **mechanistic ODE-based model** of RNA kinetics  
- A **nested hypothesis testing framework**  
- Tools to assess whether **cytoplasmic processing (σ_c)** is supported by the data  
- Scripts to **fully reproduce all results** presented in the paper  

---

## Model inputs

The method operates on time-resolved measurements of four RNA pools:

- **N(t)**: nuclear *unprocessed* RNA  
- **N_s(t)**: nuclear *processed* RNA  
- **C(t)**: cytoplasmic *unprocessed* RNA  
- **C_s(t)**: cytoplasmic *processed* RNA  

The model estimates kinetic parameters and compares:

- **Null model**: no cytoplasmic processing  
- **Alternative model**: includes cytoplasmic processing (σ_c)

---

## Repository structure

The repository is organized into four main directories:

### `commons/`
Utility functions shared across the project:

- Data handling  
- Plotting utilities  
- General helper functions  

### `ode_model/`
Core implementation of the kinetic model:

- ODE system definition  
- Parameter estimation  
- Nested hypothesis testing  
- Model fitting routines  

### `real_datasets/`
Scripts to reproduce results on **real datasets**.

- Each dataset (GSE) has its own subdirectory  
- Inside each GSE directory:
  - Scripts to generate **time-course RNA states** starting from **rMATS outputs**

Main entry point:

- `run_nested_test.r` → runs the full analysis across all 4 datasets

### `synthetic_dataset/`
Scripts for **simulation and benchmarking** on synthetic data.

- `gen_synthetic_ODE_states.r`  (run first)
  → generates synthetic RNA time-course data from the ODE model  

- `run_tests.r`  
  → runs the nested tests on synthetic data  

- `fig_pr_roc.r`  
  → generates Precision-Recall and ROC curves  

- `fig_calibration_power.r`  
  → generates calibration and statistical power plots  

### Settings
All scripts include explicit documentation in their headers, describing:
-	Required input files
-	Output formats
-	Model parameters
-	Optional arguments

### Notes
-	Synthetic data must be generated before running synthetic tests
-	Real dataset preprocessing depends on rMATS outputs
-	Plotting scripts reproduce the figures from the paper