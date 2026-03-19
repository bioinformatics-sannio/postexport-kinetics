# Post-export RNA kinetic  
### Predict post-export RNA processing with a nested ODE kinetic test

Compartment-resolved kinetic inference and **nested hypothesis testing** for detecting **post-export RNA processing** from time-resolved nuclear/cytoplasmic transcriptomic data.

This repository contains the R implementation of the modeling and inference framework described in:

> Faretra L., Napolitano F., Pancione M., Cerulo L. (2026). *A kinetic modeling framework to detect post-export RNA processing from time-resolved transcriptomic data.* Submitted manuscript.

---

## Overview

Many transcriptomic analyses assume that RNA processing (e.g., splicing) occurs exclusively in the nucleus. However, increasing evidence suggests that **RNA processing and remodeling may continue after nuclear export**.

This repository basically provides scripts to reproduce results reported in the manuscript. Specifically it provides:

- A **mechanistic ODE-based model** of RNA kinetics  
- A **nested hypothesis testing framework**  
- Tools to assess whether **cytoplasmic processing (σ_c)** is supported by the data  in real and synthetic datasets

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

## Required packages

Required packages can be installe as follows

```r
install.packages(c(
  "data.table","deSolve","nnls","MASS","parallel",
  "ggplot2","patchwork","scales","plotROC","pROC","PRROC",
  "biomaRt","AnnotationDbi","clusterProfiler","ReactomePA","enrichplot",
  "openxlsx","knitr","tidyr"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "org.Hs.eg.db","org.Mm.eg.db","org.Dm.eg.db"
))
```

### Repository organization

The repository is organized into four main directories:

- `commons/` utility functions shared across the project for: data handling ,plotting utilities, and general helper functions.

- `ode_model/` core implementation of the kinetic model: ODE system definition ,parameter estimation, nested hypothesis testing, model fitting routines.

- `real_datasets/` scripts to reproduce results on **real datasets**. Each dataset (GSE) has its own subdirectory. Inside each GSE directory a script `gen_RMATS_table.r` generates **time-course RNA states** starting from **rMATS outputs** (.csv provided in each directory).

- `synthetic_dataset/` scripts for **simulation and benchmarking** on synthetic data.

### Synthetic data

Synthetic data must be generated first with `gen_synthetic_ODE_states.r`. The script generates synthetic RNA time-course data from the ODE model. Generation parameters can be set in the header of the script.

The test (nested and PSI) can then be executed with `run_tests.r`. The script runs the nested tests on synthetic data. The execution can be performed in parallel using the parallel R package.

Figures and tables reported in the paper can be obtained with `fig_pr_roc.r`, which generates Precision-Recall and ROC curves, and `fig_calibration_power.r`, which generates calibration and statistical power plots.

### Real data

Real datasets are derived from RNA-seq experiments processed with rMATS, starting from raw FASTQ files. Raw reads are aligned to the reference genome using tools such as STAR, producing BAM files. The aligned reads are analyzed with rMATS to obtain event-level counts for inclusion/skipping isoforms across compartments and time points. The rMATS outputs are then processed using `rmats_preprocessing.r`. This step converts rMATS outputs into time-resolved abundance estimates of the four RNA species required by the ODE model. For convenience, in each GSE directory the processed rMATS outputs are already provided as .csv files.

The nested test can be applied to all datasets using `run_nested_test.r`. This script merges all datasets, applies quality filtering, and runs the nested kinetic test on each event.

The main figures and tables reported in the manuscript can be reproduced using: `fig_diagnostics.r`, `fig_representative_dynamics.r`, and `go_analysis.r`.

### Settings

All scripts assume that the repository has been cloned in the home directory (~) and include explicit documentation in their headers, describing:

-	Required input files
-	Output formats
-	Model parameters
-	Optional arguments

### Notes

The implementation includes:
- variance shrinkage for stability
- weighted regression (heteroskedastic noise)
- non-negativity constraints (NNLS)
- bootstrap-based inference due to boundary parameters

The full pipeline is fully reproducible starting from: synthetic simulations or rMATS-derived real datasets

-	Synthetic data must be generated before running synthetic tests
-	Plotting scripts to reproduce the figures from the paper are provided