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
- A **cytoplasmic processing ($\sigma_c$) nested hypothesis testing framework**  with variance shrinkage for stability, weighted regression (heteroskedastic noise), non-negativity constraints (NNLS), and bootstrap-based inference due to boundary parameters.
- **Synthetic** and **real** datasets used to test the framework.

---

## Model inputs

The method operates on time-resolved measurements of four RNA pools:

- $N(t)$: nuclear *unprocessed* RNA  
- $N_s(t)$: nuclear *processed* RNA  
- $C(t)$: cytoplasmic *unprocessed* RNA  
- $C_s(t)$: cytoplasmic *processed* RNA  

The model estimates kinetic parameters and compares:

- **Null model**: no cytoplasmic processing  
- **Alternative model**: includes cytoplasmic processing ($\sigma_c$)

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

Real datasets are derived from RNA-seq experiments, deposited by the respective studies. Raw FASTQ reads are aligned to the reference genomes using STAR (v2.7) and processed according to authors directives. The aligned reads then analyzed with rMATS to obtain event-level counts for inclusion/skipping isoforms across compartments and time points. The rMATS outputs are then converted into time-resolved abundance estimates of the four RNA species required by the ODE model. For convenience, in each GSE directory the processed rMATS outputs are already provided as .csv files.

The nested test can be applied to all datasets using `run_nested_test.r`. This script merges all datasets, applies quality filtering, and runs the nested kinetic test on each event.

The main figures and tables reported in the manuscript can be reproduced using: `fig_diagnostics.r`, `fig_representative_dynamics.r`, and `go_analysis.r`.

#### Details on raw data processing

To reproduce the processing from raw FASTQ files we applied the following procedure:

- For the GSE207924 (human K562 and mouse NIH3T3 cells), we followed the workflow described by Chen et al. Adapter sequences and low-quality bases were removed using `cutadapt`, reads were aligned to a chimeric reference genome, retrieved from Zenodo, using STAR (v2.7), and aligned reads were then filtered using `samtools` to remove secondary alignments, PCR duplicates, and unmapped mates, and subsequently partitioned into host-specific, mitochondrial, and spike-in subsets based on chromosome identifiers. The shell script `GSE207924.sh` implements the above described FASTQ processing pipeline. 

- For the GSE83620 (D. melanogaster), initial quality control and adapter trimming were performed using `fastp`, including poly-G tail removal. Reads were aligned to the D. melanogaster (BDGP6) and S. cerevisiae (R64-1) genomes, the latter serving as a normalization spike-in. The shell script `GSE83620.sh` implements the above described FASTQ processing pipeline.

Intron retention (RI) events were quantified using rMATS (v4.1), configured for paired-end reads in GSE207924 and single-end reads in GSE83620, both with support for variable read lengths. This analysis provided event-level counts for inclusion isoforms (reads supporting intron retention) and skipping isoforms (reads supporting exon-exon junctions) across nuclear and cytoplasmic compartments at all time points.

Downstream processing was performed in the R environment using the script rmats_preprocessing.r, which converts rMATS outputs into time-resolved abundance estimates of the four RNA species required by the ODE model. To account for technical variability in sequencing depth and sample recovery, raw counts were normalized using exogenous spike-ins: ERCC counts for human and mouse datasets, and yeast-mapped reads for Drosophila. Counts were first adjusted for effective isoform length—accounting for differences between Inclusion and Skipping forms—and subsequently scaled by the total number of uniquely mapped spike-in reads.

### Settings

All scripts assume that the repository has been cloned in the home directory (~) and include explicit documentation in their headers, describing required input files, output formats, model parameters, and any other parameter for the reproduction of results.



