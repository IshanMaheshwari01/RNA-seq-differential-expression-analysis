# RNA-Seq Differential Expression Analysis: Tumor vs Normal

[![R](https://img.shields.io/badge/R-4.3+-276DC3?style=flat&logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.17+-1881C2?style=flat)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **A comprehensive RNA-seq analysis pipeline demonstrating differential gene expression analysis, quality control, and pathway enrichment in cancer research.**

---

## Project Overview

This project implements a complete **end-to-end RNA-seq differential expression analysis pipeline** comparing tumor tissue versus normal tissue. The analysis identifies differentially expressed genes, performs quality control, and reveals enriched biological pathways characteristic of cancer.

### Objectives

- Identify genes differentially expressed between tumor and normal tissue
- Perform comprehensive quality control on RNA-seq data
- Discover enriched biological pathways and processes
- Create publication-quality visualizations
- Demonstrate reproducible bioinformatics research

---

## Biological Context

**Research Question:** Which genes and pathways are dysregulated in tumor tissue compared to normal tissue?

**Findings:**
- **3,296 differentially expressed genes** identified (16.7% of transcriptome)
- **1,654 genes upregulated** in tumor (proliferation, cell cycle)
- **1,642 genes downregulated** in tumor (immune response, defense)
- Fold changes ranging from **0.04x to 28x**

**Biological Interpretation:**
- Tumor shows strong **proliferation signature** (cell cycle genes activated)
- Evidence of **immune evasion** (immune pathways suppressed)
- Consistent with known cancer hallmarks

---

## Project Structure

```r
RNA-seq-differential-expression-analysis/
│
├── README.md # Project documentation
├── environment.R # Package installation script
├── analysis_report.html # Interactive HTML report
│
├── app/
│ └── app.R # Shiny dashboard application
│
├── data/
│ ├── raw/ # Original count data
│ │ ├── count_matrix.rds
│ │ └── sample_metadata.rds
│ └── processed/ # Filtered & normalized data
│ ├── counts_filtered.rds
│ └── vsd_matrix.rds
│
├── scripts/
│ ├── 01_data_preparation.R # Data simulation & setup
│ ├── 02_quality_control.R # QC analysis & visualization
│ ├── 03_differential_expression.R # DESeq2 analysis
│ └── 04_pathway_enrichment.R # GO enrichment analysis
│
└── results/
├── figures/ # 16 publication-quality plots
│ ├── 01_library_sizes.png
│ ├── 02_gene_detection.png
│ ├── 03_count_distribution.png
│ ├── 04_sample_correlation.png
│ ├── 05_pca_plot.png
│ ├── 06_sample_distances.png
│ ├── 07_volcano_plot.png
│ ├── 08_MA_plot.png
│ ├── 09_heatmap_top_genes.png
│ ├── 10_fold_change_distribution.png
│ ├── 11_pvalue_distribution.png
│ ├── 12_GO_upregulated_barplot.png
│ ├── 13_GO_upregulated_dotplot.png
│ ├── 14_GO_downregulated_barplot.png
│ ├── 15_GO_downregulated_dotplot.png
│ └── 16_GO_combined_comparison.png
│
└── tables/ # 9 result tables (CSV)
├── qc_summary.csv
├── DE_results_full.csv
├── DE_results_significant.csv
├── DE_results_top100.csv
├── DE_genes_upregulated.csv
├── DE_genes_downregulated.csv
├── GO_enrichment_upregulated.csv
└── GO_enrichment_downregulated.csv
```
---

## Technologies & Methods

### **Core Technologies**
- **R 4.3+** - Statistical computing
- **Bioconductor 3.17+** - Bioinformatics packages
- **DESeq2** - Differential expression analysis
- **clusterProfiler** - Pathway enrichment

### **R Packages**
```r
# Differential Expression
library(DESeq2)           # DE analysis (Love et al., 2014)
library(edgeR)            # Alternative DE method

# Visualization  
library(EnhancedVolcano)  # Volcano plots
library(pheatmap)         # Heatmaps
library(ggplot2)          # Publication-quality graphics

# Pathway Analysis
library(clusterProfiler)  # GO/KEGG enrichment
library(org.Hs.eg.db)     # Human gene annotations

# Data Manipulation
library(tidyverse)        # Data wrangling

```

## Statistical Methods
- Normalization: DESeq2 median-of-ratios
- Variance stabilization: VST transformation
- Batch correction: Linear model accounting for batch effects
- Multiple testing: Benjamini-Hochberg FDR correction
- Significance thresholds: Adjusted p-value < 0.05, |log2FC| ≥ 1

## Analysis Workflow

```r
1. Data Preparation

source("scripts/01_data_preparation.R")

- Simulated realistic RNA-seq count data
- 20,000 genes × 12 samples (6 tumor, 6 normal)
- Included batch effects and biological variability

2. Quality Control

source("scripts/02_quality_control.R")

Metrics Evaluated:

- Library size distribution (2.3M reads average)
- Gene detection rate (91.4% of transcriptome)
- Sample correlation (within > between groups)
- PCA analysis (PC1 = 24.2% variance)
- Batch effect assessment (PC2 = 7.9%)

3.  Differential Expression

source("scripts/03_differential_expression.R")

DESeq2 Analysis:

- Design: ~ batch + condition
- Identified 3,296 significant DE genes
- Top upregulated: 21.8-fold increase
- Top downregulated: 27-fold decrease

4. Pathway Enrichment

source("scripts/04_pathway_enrichment.R")

GO Enrichment Results:

- Upregulated: Cell cycle, DNA replication, mitosis
- Downregulated: Immune response, defense mechanisms
- 16 significantly enriched pathways (p < 0.05)

```

## Results

Differential Expression Summary

- Total genes analyzed - 19,779
- Significant DE genes - 3,296 (16.7%)
- Upregulated in tumor - 1,654
- Downregulated in tumor - 1,642
- Max fold change (up) - 28.1x
- Max fold change (down) - 0.04x (27x decrease)
- Min adjusted p-value - 1.36e-14

Top Enriched Pathways

Upregulated (Tumor > Normal):
- Cell cycle (95 genes, p = 2.4e-12)
- Cell division (78 genes, p = 6.8e-11)
- Mitotic cell cycle (72 genes, p = 1.1e-09)

Downregulated (Tumor < Normal):
- Immune response (89 genes, p = 1.8e-10)
- Immune system process (82 genes, p = 2.4e-10)
- Innate immune response (71 genes, p = 6.8e-09)

## Visualizations

### Quality Control

**PCA Analysis**

![PCA Plot](results/figures/05_pca_plot.png)

*PCA showing clear separation between tumor (red) and normal (blue) samples. PC1 captures 24.2% of variance.*

**Sample Correlation**

![Correlation Heatmap](results/figures/04_sample_correlation.png)

*Sample-to-sample correlation heatmap showing samples cluster by condition.*

---

### Differential Expression

**Volcano Plot**

![Volcano Plot](results/figures/07_volcano_plot.png)

*Volcano plot highlighting 3,296 differentially expressed genes. Red = upregulated, Green = downregulated.*

**Heatmap of Top DE Genes**

![Heatmap](results/figures/09_heatmap_top_genes.png)

*Heatmap of top 50 DE genes showing clear separation between tumor and normal samples.*

---

### Pathway Enrichment

![GO Enrichment](results/figures/16_GO_combined_comparison.png)

*GO enrichment analysis showing upregulated proliferation pathways (red) and downregulated immune pathways (blue).*

## How to Reproduce:

```r
- Prerequisites

# R version 4.3 or higher
# Bioconductor 3.17 or higher

# Install required packages
source("environment.R")

- Run Complete Pipeline

# Step 1: Data preparation
source("scripts/01_data_preparation.R")

# Step 2: Quality control
source("scripts/02_quality_control.R")

# Step 3: Differential expression
source("scripts/03_differential_expression.R")

# Step 4: Pathway enrichment
source("scripts/04_pathway_enrichment.R")

```

## Output Files:

Figures (16 total)
- QC plots: Library sizes, gene detection, PCA, correlation
- DE plots: Volcano, MA, heatmaps, distributions
- Pathway plots: GO enrichment bar plots, dot plots, comparisons

Tables (9 total)
- DE_results_full.csv - All genes with statistics
- DE_results_significant.csv - 3,296 significant genes
- DE_genes_upregulated.csv - 1,654 upregulated genes
- DE_genes_downregulated.csv - 1,642 downregulated genes
- GO_enrichment_upregulated.csv - Enriched pathways (up)
- GO_enrichment_downregulated.csv - Enriched pathways (down)

Additional QC and summary tables

## Skills Demonstrated:

Bioinformatics
- RNA-seq data analysis
- Differential expression analysis (DESeq2)
- Quality control & normalization
- Batch effect correction
- Pathway enrichment analysis
- Functional annotation

Statistical Analysis
- Negative binomial modeling
- Multiple testing correction (FDR)
- Variance stabilization
- Dimensionality reduction (PCA)
- Hierarchical clustering

Data Science
- R programming
- Data wrangling (tidyverse)
- Data visualization (ggplot2)
- Reproducible research
- Version control (Git/GitHub)

Domain Knowledge
- Cancer biology
- Molecular biology
- Genomics
- Gene regulation
- Biological pathway analysis

## References:

Methods & Packages
- Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550.
- Yu, G., Wang, L.G., Han, Y., He, Q.Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS, 16(5):284-287.
- Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.

Biological Context
- Hanahan, D., Weinberg, R.A. (2011). Hallmarks of cancer: the next generation. Cell, 144(5):646-674.

## Author:
Ishan Maheshwari
- MSc in Genomics Data Science | University of Galway
- LinkedIn: www.linkedin.com/in/ishanmaheshwari2001
- Email: ishanmaheshwari02@gmail.com

## License:
This project is licensed under the MIT License.

## Acknowledgements:
- Bioconductor community for excellent bioinformatics tools
- DESeq2 developers for robust statistical methods
- R community for comprehensive data science ecosystem
