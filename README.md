# Post export RNA kinetic
### Predict post-export RNA processing with a nested ODE kinetic test

Compartment-resolved kinetic inference and **nested hypothesis testing** for detecting **post-export RNA processing** from time-resolved nuclear/cytoplasmic transcriptomic data.

This repository contains the R implementation of the modeling and inference framework described in:

> Faretra L., Napolitano F., Pancione M., Cerulo L. (2026). *A kinetic modeling framework to detect post-export RNA processing from time-resolved transcriptomic data.* Manuscript.

---

## What problem does it solve?

Many analyses treat RNA processing (e.g., splicing) as purely nuclear. However, in several biological contexts, **RNA processing and remodeling may continue after nuclear export**.

Given time-resolved measurements of four RNA pools for each gene/event:

- **N(t)**: nuclear *unprocessed* RNA  
- **N\_s(t)**: nuclear *processed* RNA  
- **C(t)**: cytoplasmic *unprocessed* RNA  
- **C\_s(t)**: cytoplasmic *processed* RNA  

`nestedRNA` fits a mechanistic kinetic model and tests whether a **cytoplasmic conversion term** (σ\_c) is supported by the data.

---

## Model (mechanistic core)

The framework models the four pools with a first-order linear ODE system:

$$
\begin{aligned}
\frac{dN(t)}{dt} &= R - \sigma_n N(t) - \tau N(t) \\
\frac{dN_s(t)}{dt} &= \sigma_n N(t) - \tau_s N_s(t) \\
\frac{dC(t)}{dt} &= \tau N(t) - \sigma_c C(t) - \alpha C(t) \\
\frac{dC_s(t)}{dt} &= \tau_s N_s(t) + \sigma_c C(t) - \alpha_s C_s(t)
\end{aligned}
$$

- The parameter of interest is **σ\_c ≥ 0** (post-export conversion).
- The **nested null model** sets **σ\_c = 0**.

---

## Inference & statistical test (high-level)

1) **Integral balance discretization (Crank–Nicolson)** over each interval \([t_k, t_{k+1}]\)  
   → avoids numerical differentiation of noisy trajectories.

2) **Weighted non-negative least squares (NNLS)** for parameter estimation  
   → enforces non-negativity of kinetic rates and stabilizes estimation.

3) **Nested model comparison** (Full vs Null) using RSS difference  
   $$T = \mathrm{RSS}_0 - \mathrm{RSS}_1$$

4) **Rademacher wild bootstrap** under the null  
   → provides calibrated one-sided p-values in the presence of boundary constraints (σ\_c ≥ 0).

---

## Requirements

- R (>= 4.1 recommended)
- Suggested R packages (depending on your implementation choices):
  - `nnls`, `Matrix`, `stats`
  - `data.table` / `dplyr` for data handling
  - `ggplot2` for plotting
  - `BiocManager` if integrating with Bioconductor infrastructure

> If you are packaging this as an R package, keep dependencies minimal and list them in `DESCRIPTION`.

---

## Installation

### From GitHub (development)
```r
install.packages("remotes")
remotes::install_github("<GITHUB_ORG_OR_USER>/nested-rna")
