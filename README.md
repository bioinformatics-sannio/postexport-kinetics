# Post-export RNA kinetics

Compartment-resolved kinetic modeling and constrained nested-model comparison for identifying RNA trajectories consistent with an additional post-export conversion component.

This repository contains the analysis and reproducibility code associated with:

> Faretra L., Napolitano F., Pancione M., Cerulo L. (2026). *Kinetic model comparison identifies RNA trajectories consistent with post-export processing.* Revised manuscript submitted to **Bioinformatics** (Manuscript BIOINF-2026-0699).

## Scientific scope

The framework models four compartment-resolved RNA states: `N(t)` (nuclear unprocessed), `N_s(t)` (nuclear processed), `C(t)` (cytoplasmic unprocessed), and `C_s(t)` (cytoplasmic processed).

The parameter `sigma_c` is a phenomenological post-export conversion rate from `C` to `C_s`. A positive or statistically supported `sigma_c` does **not** by itself identify a biochemical mechanism and should not be interpreted as direct evidence of cytoplasmic splicing.

Inference compares a null model (`sigma_c = 0`) with a full constrained model (`sigma_c >= 0`), with scientific one-sided alternative `sigma_c > 0`. Significance is assessed by a replicate-level generative bootstrap with covariance propagation, reconstruction of observation-derived quantities, non-negativity constraints, explicit boundary handling, and an add-one p-value correction.

## Reproducibility note

This repository contains historical development scripts as well as the finalized revised-manuscript analysis. Only the scripts listed under **Final manuscript pipeline** define the authoritative analysis path.

Frozen manuscript tag:

`manuscript-revision-v1.0`

A reusable R package is developed separately as `postexportKinetics`.

## Final manuscript pipeline

### Core inference
- `commons/nested_test2.r`
- `ode_model/ode.r`

### Corrected synthetic benchmark
1. `synthetic_dataset/gen_synthetic_ODE_states_corrected_onset.R`
2. `synthetic_dataset/run_benchmark_main_corrected_onset_revision.R`
3. `synthetic_dataset/analyze_benchmark_corrected_onset_final.R`

The corrected factorial benchmark contains 1,152 experimental configurations and uses 1,999 generative-bootstrap replicates per test.

### Pseudo-shutoff and model misspecification
- `synthetic_dataset/run_benchmark_pseudoshutoff_revision.R`
- `synthetic_dataset/run_pseudoshutoff_misspecification_benchmark.R`
- `synthetic_dataset/run_fraction_separation_robustness_benchmark.R`
- `synthetic_dataset/make_FigS_pseudoshutoff_minimal.R`
- `synthetic_dataset/make_FigS_fraction_misspecification.R`

### Practical identifiability and numerical validation
- `synthetic_dataset/analyze_practical_identifiability.R`
- `synthetic_dataset/make_FigS_two_representative_synthetic_refits_FINAL.R`
- `synthetic_dataset/analyze_CN_vs_exact_transition.R`

### Comparator and ranking analyses
- `synthetic_dataset/run_cytoplasmic_only_baseline.R`
- `synthetic_dataset/make_FigS_full_vs_cytoplasmic_only_typeI.R`
- `synthetic_dataset/analyze_three_way_AUPR.R`
- `synthetic_dataset/analyze_composite_score_ablation.R`
- `synthetic_dataset/regenerate_corrected_discrimination_figures.R`

The composite score is an exploratory prioritization summary, not an inferential statistic or an optimized ranking rule.

### mESC matched-design sensitivity
- `synthetic_dataset/run_GSE256335_matched_corrected_onset_benchmark.R`
- `synthetic_dataset/analyze_GSE256335_effect_size_power.R`

### Real-data analysis
- `real_datasets/run_real_datasets_revision.R`
- `real_datasets/run_mESC_20k.R`
- `real_datasets/audit_final_real_data.R`
- `real_datasets/fig_representative_Ppp1r36dn_Nsd1_common_y0_FINAL.R`
- `real_datasets/make_FigS19_realdata_QQ_final.R`
- `real_datasets/make_S5testresults_final.py`

Final Supplementary Table S5:
- `real_datasets/S5testresults_final.xlsx`

## Final real-data audit

Benjamini-Hochberg correction is performed separately within each dataset.

| Dataset | Tested RI events | Unique genes | p < 0.05 | q < 0.10 | q < 0.05 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Kc167 | 337 | 286 | 16 | 0 | 0 |
| K562 | 2696 | 1635 | 121 | 0 | 0 |
| NIH-3T3 | 1746 | 1288 | 79 | 0 | 0 |
| mESC | 1972 | 1302 | 179 | 28 | 14 |

The mESC pharmacological-shutoff analysis contains FDR-supported events. The three pseudo-shutoff datasets are interpreted as exploratory rankings rather than FDR-controlled discoveries.

Machine-readable audit summaries are under `real_datasets/final_real_data_audit/`.

## Real datasets

| Dataset | System | Design |
| --- | --- | --- |
| GSE83620 | Drosophila Kc167 | pseudo-shutoff |
| GSE207924 | human K562 | pseudo-shutoff |
| GSE207924 | mouse NIH-3T3 | pseudo-shutoff |
| GSE256335 | mouse embryonic stem cells | pharmacological shutoff |

GSE256335 uses pharmacological transcriptional inhibition. GSE83620 and GSE207924 use metabolic-labeling-based pseudo-shutoff constructions and are not treated as equivalent to direct inhibition.

## RNA-seq preprocessing

Retained-intron events are quantified with rMATS and nuclear/cytoplasmic inclusion and skipping measurements are converted to the four kinetic states. Dataset-specific preprocessing and processed event-level inputs are retained where required for manuscript reproduction. Raw public sequencing data are not duplicated here.

Raw-read workflows use standard tools including SRA Toolkit, cutadapt, STAR, samtools, rMATS, fastp, and RSEM. Large genomes, transcriptomes, FASTQ/BAM files, third-party installations, checkpoints, and temporary intermediates are intentionally excluded.

## Software dependencies

Analyses use CRAN packages including `data.table`, `deSolve`, `nnls`, `MASS`, `parallel`, `ggplot2`, `patchwork`, `scales`, `pROC`, `PRROC`, `openxlsx`, `knitr`, `tidyr`, and `readr`, plus Bioconductor annotation/enrichment packages including `biomaRt`, `AnnotationDbi`, `clusterProfiler`, `ReactomePA`, `enrichplot`, `org.Hs.eg.db`, `org.Mm.eg.db`, and `org.Dm.eg.db`.

Not every dependency is required for every analysis. Final computational runs record session information where available.

## Repository organization

- `commons/`: statistical inference and shared utilities.
- `ode_model/`: ODE definitions and simulation utilities.
- `synthetic_dataset/`: synthetic generation, benchmarking, robustness, comparators, identifiability, and figures.
- `real_datasets/`: preprocessing, final real-data inference, audits, and supplementary outputs.

## Detailed real-data preprocessing

The public RNA-seq datasets are processed from raw sequencing data to event-level retained-intron measurements and then converted into the four kinetic states required by the model. Raw FASTQ files are not redistributed in this repository.

### GSE256335 — mouse embryonic stem cells

Reads were aligned to the mouse reference genome with STAR (v2.7) in two-pass mode. Alignment allowed a maximum of three mismatches per read, end-to-end alignment, and a maximum intron length of 299,999 bp. Retained-intron events were quantified with rMATS. Event-level counts were corrected for effective isoform length; sample-specific normalization used total-expression estimates from RSEM, with nuclear/cytoplasmic recovery scaling as described in the manuscript and Supplementary Methods.

This dataset uses direct pharmacological transcriptional inhibition and is analyzed separately from the pseudo-shutoff datasets.

### GSE207924 — human K562 and mouse NIH-3T3

Adapter and low-quality sequence removal used cutadapt. Reads were aligned with STAR to the appropriate host/spike-in reference and filtered with samtools. Retained-intron events were quantified with rMATS. Event-level measurements were combined with gene-level metabolic-labeling information to reconstruct the pre-existing RNA component used as a pseudo-shutoff approximation.

Because labeling fractions are estimated at gene level whereas retained-intron measurements are event-specific, this reconstruction is treated as an approximation rather than as equivalent to pharmacological shutoff.

### GSE83620 — Drosophila Kc167

Initial quality control and adapter trimming used fastp, including poly-G tail removal. Reads were aligned to the Drosophila reference genome and to the Saccharomyces cerevisiae spike-in reference used for normalization. Retained-intron events were quantified with rMATS. The experimentally isolated unlabeled/pre-existing RNA fraction is used as the pseudo-shutoff measurement.

### rMATS quantification and kinetic states

rMATS v4.1 was used for retained-intron quantification, with paired-end configuration for GSE207924 and GSE256335 and single-end configuration for GSE83620, with variable read lengths supported as appropriate.

Nuclear/cytoplasmic inclusion and skipping measurements are transformed into:

- nuclear inclusion -> `N`
- nuclear skipping -> `N_s`
- cytoplasmic inclusion -> `C`
- cytoplasmic skipping -> `C_s`

Raw event counts are adjusted for effective isoform length. GSE207924 and GSE83620 use exogenous spike-in information for normalization. GSE256335 uses sample-specific total-expression information from RSEM together with compartment-recovery scaling. The manuscript and Supplementary Methods remain the authoritative description of preprocessing assumptions.

## Reproducibility workflow

The repository supports two levels of reproduction.

### Level 1 — inference from processed event-level inputs

This is the recommended route for reproducing the statistical results without repeating raw-read alignment and rMATS quantification.

```text
processed event-level measurements
        |
        v
construction of N, N_s, C, C_s
        |
        +--> real_datasets/run_real_datasets_revision.R
        |       |
        |       +--> Kc167 / K562 / NIH-3T3 final results
        |
        +--> real_datasets/run_mESC_20k.R
                |
                +--> final mESC results
        |
        v
real_datasets/audit_final_real_data.R
        |
        +--> final dataset counts
        +--> FDR-supported mESC events
        +--> boundary summaries
        +--> representative-event audit
        |
        +--> real_datasets/make_FigS19_realdata_QQ_final.R
        +--> real_datasets/fig_representative_Ppp1r36dn_Nsd1_common_y0_FINAL.R
        +--> real_datasets/make_S5testresults_final.py
```

Final machine-readable audit summaries provide lightweight reference outputs against which a reproduced analysis can be checked.

### Level 2 — processed inputs from public sequencing data

Raw sequencing data can be retrieved from the original public accessions. Dataset-specific preprocessing workflows document alignment, retained-intron quantification, normalization, and construction of the model states.

This route additionally requires the appropriate genome/transcriptome references and external command-line tools. Large references, indexes, FASTQ/BAM files, and third-party installations are intentionally not version-controlled.

## Synthetic benchmark reproducibility

The corrected benchmark follows:

```text
ode_model/ode.r
        |
        v
synthetic_dataset/gen_synthetic_ODE_states_corrected_onset.R
        |
        v
synthetic_dataset/run_benchmark_main_corrected_onset_revision.R
        |
        v
synthetic_dataset/analyze_benchmark_corrected_onset_final.R
```

The complete benchmark evaluates 1,152 configurations with 1,999 generative-bootstrap replicates per test and is intended for parallel Linux/HPC execution. Small final summaries, run metadata, and session information are included so that key manuscript results can be checked without distributing all large intermediate benchmark objects.

## Reference outputs for verification

A successful reproduction of the final real-data audit should recover:

| Dataset | Tested RI events | Unique genes | p < 0.05 | q < 0.10 | q < 0.05 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Kc167 | 337 | 286 | 16 | 0 | 0 |
| K562 | 2696 | 1635 | 121 | 0 | 0 |
| NIH-3T3 | 1746 | 1288 | 79 | 0 | 0 |
| mESC | 1972 | 1302 | 179 | 28 | 14 |

Final representative mESC reference events include:

- `Ppp1r36dn`: BH q = 0.019720; estimated `sigma_c` = 0.01899230 min^-1
- `Nsd1`: BH q = 0.012325; estimated `sigma_c` = 0.04742467 min^-1

These are compact verification targets for the frozen manuscript analysis.

## Execution assumptions

Most manuscript analyses were developed for Linux/HPC execution and document expected inputs and outputs in script headers. Some scripts retain project-relative working-directory assumptions because this repository preserves the exact research workflow.

For exact reproduction, preserve the repository directory structure. The separate `postexportKinetics` package is intended to provide a portable user-facing interface without manuscript-specific path assumptions.

## Random-number generation and parallel execution

Simulation and bootstrap scripts explicitly control random-number generation. Parallelization and seed handling are documented in relevant script headers. When validating the frozen analysis, preserve the script-level seed strategy unless numerical equivalence of an alternative strategy has been established.

## Session information

Session information is retained for key final computational runs where available, including the corrected factorial benchmark, to document the R environment and package versions used.

## Intentionally excluded artifacts

The Git repository intentionally excludes large or regenerable artifacts such as:

- raw FASTQ and BAM files;
- genome/transcriptome indexes;
- large intermediate `.rdata` objects;
- benchmark checkpoints and progress files;
- local third-party software installations;
- diagnostic figures not used in the manuscript;
- wet-lab primer-design intermediates.

The frozen Git/Zenodo release is intended to archive the scientific code, processed inputs required for documented workflows, lightweight verification summaries, and manuscript-facing supplementary outputs.

## External software

Raw-read preprocessing requires recent versions of:

- SRA Toolkit
- cutadapt
- STAR
- samtools
- rMATS
- fastp
- RSEM

Dataset-specific assumptions are documented in the preprocessing workflows and in the manuscript Supplementary Methods.
