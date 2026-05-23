# Computational Analysis of tRNA-Derived Fragments in Papillary Renal Cell Carcinoma

## Overview

This repository contains all analysis scripts for a computational study of tRNA-derived fragment (tRF) dysregulation in papillary renal cell carcinoma (KIRP) using publicly available TCGA data. The pipeline integrates MINTbase v2.0 tRF expression profiles with GDC RNA-sequencing gene expression and clinical data to identify candidate tRF biomarkers and tRF–mRNA regulatory interactions.

**Key findings:**
- 1,585 differentially expressed tRFs identified (670 up, 915 down in tumor vs. normal)
- 104 significant anticorrelated tRF–mRNA pairs; top 5 all involve *MET*
- tRF-23-7SIR3DR2DV achieves AUC = 0.949 for internal tumor-normal discrimination
- tRF-23-79MP9P9MDD has predicted RNA22 binding to the *MET* 3′UTR (ΔG = −19.0 kcal/mol, p = 0.036)
- tRF-30-81PV6RRNLNK8 shows exploratory association with overall survival (HR = 1.16, FDR = 0.033)

## Repository Structure

```
KIRP_tRF_project/
├── scripts/                  # All analysis scripts (see Pipeline below)
│   ├── config.R              # Set ROOT path here (R scripts)
│   ├── config.py             # Set ROOT path here (Python scripts)
│   ├── install_r_packages.R  # Install core R dependencies
│   ├── install_phase8_packages.R
│   ├── phase3_sample_metadata.py
│   ├── phase4a_build_trf_matrix.py
│   ├── phase4b_build_mrna_matrix.py
│   ├── phase4b_gdc_uuid_map.py
│   ├── phase4_trf_DE_edgeR.R
│   ├── phase5_gene_DE_edgeR.R
│   ├── phase6_anticorrelation.R
│   ├── phase7a_get_sequences.py
│   ├── phase7b_parse_rna22.py
│   ├── phase8_strengthening.R
│   ├── generate_table1.R
│   └── build_publication_figures.R
├── data_raw/                 # Raw downloaded data (not tracked — see Data below)
│   ├── mintbase_trf/KIRP/    # MINTbase v2.0 per-sample expression files
│   ├── gdc_rnaseq/           # GDC STAR-count TSV files (one folder per UUID)
│   ├── gdc_clinical/         # GDC clinical.tsv, follow_up.tsv, etc.
│   ├── gdc_biospecimen/      # GDC sample.tsv, aliquot.tsv
│   └── data_raw_rna22_*_results - Sheet1.tsv   # RNA22 v2 output (manual step)
├── data_processed/           # Intermediate matrices (not tracked — large files)
│   ├── trf_count_matrix.tsv          # 30,076 tRFs × 325 samples
│   ├── mrna_count_matrix.tsv         # 60,660 genes × 322 samples
│   ├── trf_sample_metadata.tsv
│   ├── mrna_sample_metadata.tsv
│   ├── trf_sequences.tsv
│   ├── mrna_3utr_sequences.tsv
│   └── rna22_binding_results.tsv
└── results/
    ├── tables/               # Output TSV tables (tracked)
    │   ├── table1_cohort.tsv
    │   ├── table2_candidates.tsv
    │   ├── trf_DE_full.tsv
    │   ├── trf_DE_significant.tsv
    │   ├── gene_DE_full.tsv
    │   ├── gene_DE_significant.tsv
    │   ├── anticorrelation_significant.tsv
    │   ├── rna22_binding_significant.tsv
    │   ├── roc_auc_top_trfs.tsv
    │   ├── stage_kruskal_results.tsv
    │   └── cox_univariate_results.tsv
    └── figures/              # Output PDFs/PNGs (not tracked)
```

## Data Availability

All input data are publicly available:

| Dataset | Source | Access |
|---|---|---|
| TCGA-KIRP tRF expression | MINTbase v2.0 | https://cm.jefferson.edu/tcga-mintmap-profiles/ |
| TCGA-KIRP RNA-seq counts | GDC Data Portal | https://portal.gdc.cancer.gov/projects/TCGA-KIRP (phs000178) |
| TCGA-KIRP clinical data | GDC Data Portal | Same as above |

## Setup

### 1. Configure project root

Edit `scripts/config.R` and `scripts/config.py` to point to wherever you cloned this repository:

```r
# scripts/config.R
ROOT <- "/path/to/your/KIRP_tRF_project"
```

```python
# scripts/config.py
ROOT = "/path/to/your/KIRP_tRF_project"
```

### 2. Install dependencies

**R** (≥ 4.4; tested on R 4.5.1):
```r
Rscript scripts/install_r_packages.R
Rscript scripts/install_phase8_packages.R
```

Core R packages: `edgeR`, `ggplot2`, `ggrepel`, `pheatmap`, `RColorBrewer`, `pROC`, `survival`, `survminer`, `patchwork`, `dplyr`, `tidyr`, `readr`

**Python** (≥ 3.9; tested on 3.12):
```bash
pip install pandas requests
```

### 3. Download raw data

**MINTbase tRF profiles:**  
Download TCGA-KIRP per-sample expression files from https://cm.jefferson.edu/tcga-mintmap-profiles/ and place them in `data_raw/mintbase_trf/KIRP/`. Each sample is a folder named by TCGA barcode containing a `*-exclusive-tRFs.expression.txt` file.

**GDC RNA-seq and clinical data:**  
Using the GDC Data Transfer Tool (`gdc-client v2.3`), download TCGA-KIRP STAR-count files (Workflow Type: STAR - Counts) and place in `data_raw/gdc_rnaseq/`. Download clinical and biospecimen TSV bundles from the GDC portal and place in `data_raw/gdc_clinical/` and `data_raw/gdc_biospecimen/`.

## Pipeline

Run the scripts in the following order:

| Step | Script | Description |
|---|---|---|
| 1 | `phase3_sample_metadata.py` | Parse TCGA barcodes; build sample metadata tables |
| 2 | `phase4a_build_trf_matrix.py` | Build tRF raw count matrix from MINTbase files |
| 3 | `phase4b_gdc_uuid_map.py` | Map GDC UUIDs to TCGA barcodes via GDC REST API |
| 4 | `phase4b_build_mrna_matrix.py` | Build mRNA raw count matrix from GDC STAR-count files |
| 5 | `phase4_trf_DE_edgeR.R` | Differential expression of tRFs (edgeR, TMM+QL) |
| 6 | `phase5_gene_DE_edgeR.R` | Differential expression of 16 KIRP candidate genes |
| 7 | `phase6_anticorrelation.R` | Spearman anticorrelation: DE tRFs × up-regulated genes |
| 8 | `phase7a_get_sequences.py` | Fetch tRF sequences and 3′UTR sequences via Ensembl API |
| 9 | RNA22 v2 (manual) | Submit sequences to RNA22 v2 (https://cm.jefferson.edu/rna22/); save results to `data_raw/data_raw_rna22_{GENE}_results - Sheet1.tsv` |
| 10 | `phase7b_parse_rna22.py` | Parse RNA22 output; produce binding table |
| 11 | `phase8_strengthening.R` | ROC/AUC, Kruskal–Wallis stage association, Cox survival |
| 12 | `generate_table1.R` | Generate Table 1 (cohort) and Table 2 (candidate tRFs) |
| 13 | `build_publication_figures.R` | Assemble publication-ready figures |

## Methods Summary

- **tRF differential expression:** edgeR with TMM normalization and quasi-likelihood F-tests; filter CPM > 1 in ≥ 34 samples; FDR < 0.05, |log₂FC| ≥ 1
- **Gene differential expression:** Same edgeR workflow on 16 KIRP-relevant genes from GDC STAR-counts
- **Anticorrelation analysis:** Spearman correlation (log₂CPM+1) across 290 matched tumor samples; Benjamini–Hochberg correction; 4,755 tRF–gene pairs tested
- **RNA22 binding prediction:** Sensitivity 63%, specificity 61%, seed size 7, max 1 unpaired in seed, min 12 paired bases, max folding energy −12 kcal/mol, unrestricted G:U wobble
- **ROC/AUC:** pROC package; tumor vs. solid normal as binary outcome
- **Stage association:** Kruskal–Wallis test across AJCC stages I–IV; Benjamini–Hochberg correction
- **Survival:** Univariable Cox regression (survival package); log₂CPM+1 as continuous predictor; 251 tumor samples with valid follow-up, 44 death events; interpreted as exploratory

## Citation

If you use these scripts, please cite the manuscript or github repo.

## License

Scripts are released under the MIT License.
