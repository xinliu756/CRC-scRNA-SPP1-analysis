# SPP1+ Macrophage Enrichment in Colorectal Cancer Tumor Microenvironment

A single-cell RNA-seq analysis exploring SPP1 expression across immune cell types in colorectal cancer (CRC), using publicly available data from GEO (GSE146771).

---

## Background

SPP1 (Osteopontin) has emerged as a key marker of tumor-associated macrophages (TAMs) in colorectal cancer. SPP1+ macrophages are associated with immunosuppression, tumor progression, and poor prognosis. This analysis investigates the distribution and enrichment of SPP1-expressing immune cells across tumor and normal tissues in CRC patients.

---

## Dataset

| Item | Detail |
|------|--------|
| GEO Accession | [GSE146771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146771) |
| Platform | 10x Chromium |
| Samples | 10 CRC patients |
| Total Cells | 43,817 leukocytes |
| Tissue Types | Tumor (T), Normal (N), Peripheral Blood (P) |
| Cell Types | B cell, CD4 T cell, CD8 T cell, Myeloid cell, ILC |

---

## Key Finding

> **SPP1+ myeloid cells are dramatically enriched in tumor tissue (45.0%) compared to normal tissue (1.46%) — a ~30-fold increase.**

| Tissue | Cell Type | Total Cells | SPP1+ Cells | SPP1+ Ratio |
|--------|-----------|-------------|-------------|-------------|
| Normal | Myeloid cell | 1,025 | 15 | 1.46% |
| Tumor | Myeloid cell | 2,559 | 1,152 | **45.0%** |
| Normal | B cell | 6,016 | 21 | 0.35% |
| Tumor | B cell | 1,956 | 177 | 9.05% |
| Normal | CD4 T cell | 2,808 | 2 | 0.07% |
| Tumor | CD4 T cell | 3,827 | 109 | 2.85% |

---

## Visualization

![SPP1 Analysis](results/figures/SPP1_analysis.png)

**Left:** SPP1+ cell ratio across immune cell types in tumor vs normal tissue.  
**Right:** SPP1 expression distribution in myeloid cells (tumor vs normal).

---

## Analysis Pipeline

```
Raw Data (GEO GSE146771)
        ↓
00_setup.R          — Environment setup & package installation
        ↓
01_load_data.R      — Load and inspect metadata
        ↓
02_load_expression.R — Load expression matrix, verify target genes
        ↓
03_seurat_object.R  — Create Seurat object
        ↓
04_SPP1_analysis.R  — SPP1 expression analysis across cell types & tissues
        ↓
05_visualization.R  — Generate publication-quality figures
```

---

## Requirements

```r
# R version 4.5+
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(harmony)
```

Install all dependencies:
```r
install.packages(c("Seurat", "tidyverse", "ggplot2", "patchwork", "harmony"))
```

---

## How to Reproduce

```r
# 1. Clone this repository
# 2. Download raw data from GEO (GSE146771) into data/ folder
# 3. Run scripts in order:
source("scripts/00_setup.R")
source("scripts/01_load_data.R")
source("scripts/02_load_expression.R")
source("scripts/03_seurat_object.R")
source("scripts/04_SPP1_analysis.R")
source("scripts/05_visualization.R")
```

---

## Biological Interpretation

The striking enrichment of SPP1+ macrophages in the tumor microenvironment suggests:

- **Immunosuppressive remodeling**: SPP1+ TAMs are known to suppress cytotoxic T cell activity
- **Tumor-promoting functions**: SPP1 promotes angiogenesis, ECM remodeling, and tumor invasion
- **Potential therapeutic target**: Targeting SPP1+ macrophages may improve immunotherapy response in CRC

These findings are consistent with recent literature reporting SPP1+ TAMs as a dominant immunosuppressive population in CRC tumors.

---

## Author

**Xin Liu**  
Molecular Biology Researcher | Bioinformatics | Medical AI  
GitHub: [xinliu756](https://github.com/xinliu756)

---

## Reference

Che LH, Liu JW, Huo JP, et al. A single-cell atlas of liver metastases of colorectal cancer reveals reprogramming of the tumor microenvironment in response to preoperative chemotherapy. *Cell Discovery*, 2021.

Data source: GSE146771
