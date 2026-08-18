# sars-cov2-upper-airway-age-signature

Cross-cohort transcriptomic analysis of age-associated gene expression in the
SARS-CoV-2-positive upper airway.

Repository: https://github.com/joseovandor/sars-cov2-upper-airway-age-signature

This repository contains the complete analysis code and result tables
supporting the manuscript:

> **Cross-dataset Identification of Age-Associated Upper-Airway Transcriptomic
> Signatures in SARS-CoV-2 Infection**
> Submitted to *In Silico Research in Biomedicine* (manuscript
> INSILI-D-26-00226R1; currently under revision).

This is a reorganized, cleaned-up version of the original analysis
repository; the underlying data, methods, and results are unchanged.

## Overview

Pediatric and adult patients infected with SARS-CoV-2 differ markedly in
disease severity, but reproducible, cross-cohort evidence for age-associated
transcriptional differences in the infected upper airway has been limited.
This project integrates three independent, publicly available bulk RNA-seq
cohorts of SARS-CoV-2-positive upper airway samples and asks which
age-associated genes are reproducible across cohorts, rather than relying on
a single dataset or a small, arbitrarily thresholded gene list.

The pipeline:

1. Runs an independent DESeq2 differential expression analysis (Pediatric vs.
   Adult) in each of the three cohorts.
2. Harmonizes the cohort-level effect sizes by gene symbol and combines them
   with a random-effects (REML) meta-analysis, correcting for multiple
   testing (Benjamini-Hochberg FDR) across the complete eligible gene
   universe, without preselecting genes by effect size or direction.
3. Assesses the robustness of the resulting signature with a
   leave-one-cohort-out sensitivity analysis and a descriptive
   robustness classification (Robust / Supported / Heterogeneous).
4. Performs pre-ranked GSEA independently within each cohort on the complete,
   unthresholded ranked gene list (MSigDB Hallmark and Reactome/GO
   collections), and evaluates which pathways are enriched with a consistent
   direction across cohorts.
5. Validates the cellular origin of the signature genes in an independent,
   publicly available single-cell RNA-seq dataset of the SARS-CoV-2-infected
   nasal mucosa (Yoshida et al., *Nature*, 2022).

## Data sources

| Cohort | Accession | Samples used | Notes |
|---|---|---|---|
| GSE172274 | [GSE172274](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172274) | 13 adults / 5 pediatric | Age groups defined directly from sample metadata (Pediatric: age < 18) |
| GSE179277 | [GSE179277](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179277) | 45 adults / 38 pediatric | Age groups as labeled in the original study ("Peds"/"Adult") |
| GSE231409 | [GSE231409](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231409) | 20 adults / 104 pediatric | Sequencing batch included as a covariate when estimable |
| Yoshida et al. (single-cell) | [Nature 602:321-327 (2022)](https://doi.org/10.1038/s41586-021-04345-x) | 10 adult / 18 pediatric donors | CELLxGENE H5AD object converted to a Seurat object for orthogonal validation |

Total bulk RNA-seq samples across the three cohorts: 225 (78 adults, 147
pediatric).

Raw sequencing data and processed count matrices are **not** redistributed in
this repository; they are downloaded directly from GEO / CELLxGENE by the
scripts in `scripts/01_deseq2/` and `scripts/04_single_cell_rna_seq/`.

## Repository structure

```
sars-cov2-upper-airway-age-signature/
├── scripts/
│   ├── 01_deseq2/                     # Per-cohort differential expression (DESeq2)
│   │   ├── 01_GSE172274.R
│   │   ├── 02_GSE179277.R
│   │   └── 03_GSE231409.R
│   ├── 02_meta_analysis/              # Cross-cohort random-effects meta-analysis
│   │   ├── 01_harmonize_meta_inputs.R
│   │   ├── 02_random_effects_meta_analysis.R
│   │   ├── 03_meta_analysis_robustness.R
│   │   ├── 04_meta_analysis_figures.R
│   │   └── 05_previous_signature_concordance.R
│   ├── 03_gsea/                       # Pre-ranked GSEA, by cohort and cross-cohort
│   │   ├── 01_GSEA_by_cohort.R
│   │   ├── 02_GSEA_cross_cohort_concordance.R
│   │   └── 03_GSEA_metabolism_and_figures.R
│   ├── 04_single_cell_rna_seq/        # Orthogonal single-cell validation
│   │   ├── 00_convert_Yoshida_h5ad_to_Seurat.R
│   │   └── 01_Yoshida_airway_single_cell_validation.R
│   └── 05_figures/                    # Main-figure assembly scripts
│       ├── 02_Figure 2 Forest plot of cohort-specific and pooled random-effects estimates.R
│       ├── 03_Figure 3 Leave-one-cohort-out sensitivity.R
│       └── 04_Figure 4 Cross-dataset GSEA.R
└── results/
    ├── deseq2/<GSE...>/{primary,sensitivity,diagnostics,metadata}/
    ├── meta_analysis/{harmonization,primary,sensitivity,figures,signature_concordance,tables}/
    ├── gsea/{by_cohort,cross_cohort,figures_final,figures_final_v2,diagnostics}/
    ├── scRNAseq_COVID_nasal_age/{all_tissues,cellular_composition,tables}/
    └── scRNAseq_validation_reviewer/{MAIN_compact,TABLES}/
```

Each `results/` subfolder contains the numeric outputs and a
`diagnostics/sessionInfo_*.txt` file recording the exact R and package
versions used to generate it.

This structure mirrors the original working repository at the time of the
round-3 manuscript resubmission. If you reorganize folders while cleaning up
this repo, please update this section (and the paths listed under
"Reproducing the analysis" below) to match, so the two never drift apart.

Note: the script that generated the PCA/forest-plot panels used in main
Figure 1 was not yet present under `scripts/05_figures/` at the time this
README was written; that figure was assembled from
`results/deseq2/*/diagnostics/*_PCA_scores.csv` and the pooled forest-plot
panel. Add the corresponding script here when it is finalized, so that all
four main figures are fully reproducible from this repository.

## Reproducing the analysis

The project is organized as an RStudio Project and scripts resolve all paths
with `here::here()` — clone the repository and open the `.Rproj` file (or
otherwise set the working directory to the repository root) rather than
hard-coding machine-specific paths.

Recommended run order:

```
scripts/01_deseq2/01_GSE172274.R
scripts/01_deseq2/02_GSE179277.R
scripts/01_deseq2/03_GSE231409.R
scripts/02_meta_analysis/01_harmonize_meta_inputs.R
scripts/02_meta_analysis/02_random_effects_meta_analysis.R
scripts/02_meta_analysis/03_meta_analysis_robustness.R
scripts/02_meta_analysis/04_meta_analysis_figures.R
scripts/02_meta_analysis/05_previous_signature_concordance.R
scripts/03_gsea/01_GSEA_by_cohort.R
scripts/03_gsea/02_GSEA_cross_cohort_concordance.R
scripts/03_gsea/03_GSEA_metabolism_and_figures.R
scripts/04_single_cell_rna_seq/00_convert_Yoshida_h5ad_to_Seurat.R
scripts/04_single_cell_rna_seq/01_Yoshida_airway_single_cell_validation.R
scripts/05_figures/02_Figure 2 Forest plot of cohort-specific and pooled random-effects estimates.R
scripts/05_figures/03_Figure 3 Leave-one-cohort-out sensitivity.R
scripts/05_figures/04_Figure 4 Cross-dataset GSEA.R
```

Each script is intentionally linear (no custom function definitions) and
writes its outputs to the corresponding `results/` subfolder shown above.

### Software environment

All analyses were run under:

- R 4.3.2 (2023-10-31)
- DESeq2 1.42.1
- metafor 4.8-0
- fgsea 1.28.0
- Seurat 5.5.0 / SeuratObject 5.4.0
- GEOquery, tidyverse, org.Hs.eg.db, AnnotationDbi, apeglm, edgeR, msigdbr, here

Complete `sessionInfo()` output for every analysis stage is stored alongside
its results (see `results/**/diagnostics/sessionInfo_*.txt`), and is also
distributed as part of the manuscript's supplementary materials package
(`reproducibility_environment/`) so the exact environment can be checked
without re-running the pipeline.

## Key methodological conventions

These conventions are used consistently across scripts and are important for
interpreting any intermediate file in `results/`:

- Contrast is always **Pediatric vs. Adult**: a positive log2 fold change (or
  positive GSEA NES) means higher expression / enrichment in pediatric
  samples; negative means higher in adults.
- Primary low-count filtering: raw count ≥ 10 in at least as many samples as
  the smaller age group in that cohort, per cohort. A CPM ≥ 1
  sensitivity filter is also run in parallel for comparison
  (`results/deseq2/*/sensitivity/`).
- Cross-cohort meta-analysis uses **unshrunken** DESeq2 log2FC/SE, combined
  by gene symbol, with a random-effects (REML) model as primary and a
  fixed-effect (inverse-variance) model as a sensitivity check.
- Multiple-testing correction (Benjamini-Hochberg) is applied across the
  complete set of genes eligible for meta-analysis (present in ≥ 2 of 3
  cohorts), with no preselection of genes before FDR correction.
- Pathway enrichment uses pre-ranked GSEA on the complete, unthresholded
  per-cohort ranking (DESeq2 Wald statistic), not over-representation
  analysis of a small significant-gene list.

## Supplementary materials mapping

The tables and figures referenced from the manuscript and from the point-by
point response to reviewers are direct, unmodified exports from this
pipeline (only file names were standardized to the `Supplementary_Table_S_*`
labels used in the manuscript text). See the `README_index.md` distributed
with the manuscript's Supplementary Materials package for the complete
file-by-file mapping to specific manuscript sections, figures, and reviewer
comments.

## License

No license file is currently included in this repository. Until one is
added, please treat the code and result tables as "all rights reserved" by
the authors and contact them before reuse or redistribution; the authors are
encouraged to add an explicit license (e.g., MIT for code, CC-BY for data
tables) appropriate to their institution's policy.

## Citation

If you use this code or these results, please cite the manuscript above
(citation details to be finalized upon acceptance) and the original data
sources: GSE172274, GSE179277, GSE231409 (NCBI GEO), and Yoshida et al.,
*Nature* 602:321-327 (2022) for the single-cell validation dataset.

## Contact

For questions about this repository or the associated manuscript, please
contact the corresponding author listed on the manuscript title page, or
open an issue in this repository.
