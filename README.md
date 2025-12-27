# Single-Cell Melanoma Pathway Analysis (SCMPA)

[![Reproducibility](https://img.shields.io/badge/Reproducibility-Guaranteed-brightgreen.svg)](REPRODUCIBILITY_GUIDE.md)
[![R Version](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains R scripts for analyzing single-cell RNA sequencing (scRNA-seq) data with a focus on pathway analysis and drug treatment effects. The analysis pipeline processes data from the GSE164897 dataset, which examines melanoma cells under different drug treatments.

## 🎯 Overview

The analysis consists of four main components that work together to provide comprehensive insights into cellular responses to drug treatments:

1. **Data Processing and Integration** (`GSE164897.R`)
2. **Single-Cell Pathway Analysis** (`SCPA.R`)
3. **Advanced Visualizations** (`visualizations.R`)
4. **Trajectory and Pseudotime Analysis** (`Psudotime.R`)

## 📋 Table of Contents

- [File Descriptions](#file-descriptions)
- [Reproducibility](#reproducibility)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Workflow](#workflow)
- [Documentation](#documentation)
- [Data Source](#data-source)

## 📁 File Descriptions

### 1. GSE164897.R - Data Processing and Integration

This script handles the complete preprocessing pipeline for the GSE164897 dataset:

**Key Functions:**
- **Data Loading**: Reads multiple gzipped count files from the dataset directory
- **Seurat Object Creation**: Converts raw count matrices into Seurat objects for each sample
- **Quality Control**: Applies standard filtering (min.cells = 3, min.features = 200)
- **Sample Balancing**: Downsamples treatment groups to equal sizes (reproducible with seed = 42)
- **Normalization**: Performs log-normalization and scaling
- **Dimensionality Reduction**: Runs PCA and UMAP
- **Batch Correction**: Uses CCA (Canonical Correlation Analysis) integration
- **Clustering**: Graph-based clustering with Louvain algorithm
- **Cell Type Annotation**: Automatic melanoma cell state annotation using module scoring

**Treatment Groups:**
- `untreated`: Control cells
- `Vemurafenib`: BRAF inhibitor monotherapy
- `vem_cob`: Vemurafenib + Cobimetinib combination
- `vem_tram`: Vemurafenib + Trametinib combination

**Cell Types Annotated:**
- `Differentiated`: Melanocytic/pigmented state (MLANA, TYR, DCT)
- `Undifferentiated`: Neural crest-like state (SOX10, MITF, PAX3)
- `Invasive`: Mesenchymal/resistant state (AXL, NGFR, VIM)
- `Proliferative`: Actively dividing cells (MKI67, TOP2A)
- `Resistant`: Adaptive resistance state (EGFR, PDGFRB)
- `Hypoxic`: Stress response state (HIF1A, VEGFA)
- `ImmuneEvasion`: Immune checkpoint expression (CD274, PDCD1LG2)

**Output**: 
- Integrated Seurat object (`GSE164897`) with cell type annotations
- `GSE164897$celltype` - Dominant melanoma cell state per cell
- `GSE164897$celltype_score` - Confidence score for cell type assignment
- `GSE164897$<State>_score` - Individual module scores for each state
- `session_info_GSE164897.txt`
- UMAP visualizations including cell type plots

### 2. SCPA.R - Single-Cell Pathway Analysis

This script performs pathway-level differential analysis using SCPA:

**Key Functions:**
- **Gene ID Conversion**: Automatic detection and conversion of Ensembl IDs to gene symbols
- **Sample Extraction**: Separates cells by treatment condition
- **Pathway Database**: Loads MSigDB Hallmark gene sets (50 pathways)
- **SCPA Analysis**: Rank-based pathway enrichment with FDR correction
- **Visualization**: Creates bar plots, volcano plots, and dot plots

**Key Features:**
- Parallel processing (4 cores)
- Multiple testing correction (q-values)
- Publication-ready visualizations

**Output**: 
- `SCPA_untreated_vs_vem.csv` - Complete results
- `SCPA_BarPlot_*.pdf/png` - Top pathways
- `SCPA_VolcanoPlot_*.pdf/png` - FC vs significance
- `SCPA_DotPlot_Significant.pdf/png` - Significant pathways
- `session_info_SCPA.txt`

### 3. visualizations.R - Advanced Pathway Visualization

This script creates comprehensive cross-treatment, cross-cluster pathway heatmaps:

**Key Functions:**
- **Pathway Database**: MSigDB Hallmark gene sets (50 pathways)
- **Cluster-Specific Analysis**: Analyzes each cluster separately for each treatment
- **Signed Pathway Scores**: Calculates sign(FC) × -log10(qval) for directionality
- **Complex Heatmaps**: Treatment × Cluster pathway activity matrices
- **Variance Analysis**: Identifies most variable pathways across conditions

**Key Features:**
- Cluster-specific volcano plots
- Pathway variance ranking
- Custom color schemes for treatments and clusters

**Output**: 
- Complex heatmaps (displayed in R)
- Pathway variance plots
- Cluster-specific volcano plots
- `session_info_visualizations.txt`

### 4. Psudotime.R - Trajectory and Pseudotime Analysis

This script infers cellular trajectories and analyzes pathway dynamics:

**Key Functions:**
- **Trajectory Inference**: Uses Slingshot algorithm with CLARA clustering
- **Pseudotime Calculation**: Orders cells along developmental/response trajectory
- **Milestone Assignment**: Discretizes trajectory into states
- **Pathway Dynamics**: Analyzes pathway activity along pseudotime
- **Per-Cell Scoring**: Calculates module scores for all pathways
- **UMAP Visualization**: Spatial distribution of pathway activity

**Key Features:**
- Reproducible trajectory inference (seed controlled)
- Pathway activity heatmaps across pseudotime
- Top variable pathway identification
- Viridis color scales for accessibility

**Output**: 
- Trajectory plots (pseudotime, treatment, milestones)
- Pseudotime × pathway heatmaps
- Pathway ranking plots
- UMAP feature plots for all variable pathways
- `session_info_Psudotime.txt`

## 🔬 Reproducibility

**This analysis is fully reproducible.** All random operations use a fixed seed (42) to ensure identical results across runs.

### Quick Start

```r
# All scripts automatically load reproducibility configuration
source("GSE164897.R")      # Preprocessing
source("SCPA.R")           # Pathway analysis
source("visualizations.R") # Visualizations
source("Psudotime.R")      # Trajectory analysis
```

### Key Features

- **Global Seed**: All random operations use seed = 42
- **Session Tracking**: Every script saves package versions and system info
- **Parallel Reproducibility**: L'Ecuyer-CMRG RNG for multi-core operations
- **Automated Verification**: Built-in package version checks

### Documentation

- 📖 **[REPRODUCIBILITY_GUIDE.md](REPRODUCIBILITY_GUIDE.md)** - Complete reproducibility documentation
- 📋 **[REPRODUCIBILITY_QUICKREF.md](REPRODUCIBILITY_QUICKREF.md)** - Quick reference card
- 📊 **[Technical_Documentation.md](Technical_Documentation.md)** - Detailed technical documentation

### Verification

After running analyses, check:
```r
# List session info files
list.files(pattern = "session_info_.*\\.txt")

# List reproducibility reports
list.files(pattern = "reproducibility_.*\\.txt")
```

## 📦 Dependencies

### Core Packages

| Package | Version | Purpose |
|---------|---------|---------|
| **Seurat** | ≥5.0.0 | Single-cell analysis framework |
| **SeuratObject** | Latest | Seurat data structures |
| **SCPA** | Any | Pathway analysis for scRNA-seq |
| **dyno** | Any | Trajectory inference |
| **msigdbr** | Any | MSigDB pathway databases |
| **tidyverse** | Any | Data manipulation |
| **ComplexHeatmap** | Any | Advanced heatmaps |
| **org.Hs.eg.db** | Any | Human gene annotation |
| **ggrepel** | Any | Non-overlapping text labels |

### Installation

```r
# CRAN packages
install.packages(c("Seurat", "tidyverse", "ggrepel", "BiocManager"))

# Bioconductor packages
BiocManager::install(c("org.Hs.eg.db", "ComplexHeatmap"))

# GitHub packages (if needed)
remotes::install_github("jackbibby1/SCPA")
remotes::install_github("dynverse/dyno")
```

## 🚀 Usage

### 1. Setup

```r
# Clone repository
git clone https://github.com/SreeSatyaGit/Single-Cell-Melanoma-Pathway-Analysis.git
cd Single-Cell-Melanoma-Pathway-Analysis

# Update data path in GSE164897.R
base_dir <- "/path/to/your/data"  # Line 7 in GSE164897.R
```

### 2. Run Analysis

```r
# Run in order (each script sources reproducibility_config.R automatically)
source("GSE164897.R")      # ~10-15 min
source("SCPA.R")           # ~5-10 min
source("visualizations.R") # ~15-20 min
source("Psudotime.R")      # ~10-15 min
```

### 3. Verify Reproducibility

```r
# Check that session info files were created
list.files(pattern = "session_info")
# Should show: session_info_GSE164897.txt, session_info_SCPA.txt, etc.

# Verify package versions
source("reproducibility_config.R")
verify_package_versions()
```

## 🔄 Workflow

```
Raw Data (GSE164897)
    ↓
[GSE164897.R] Preprocessing & Integration
    ↓
Integrated Seurat Object
    ↓
    ├─→ [SCPA.R] Treatment-level pathway analysis
    │       ↓
    │   Pathway fold changes & significance
    │
    ├─→ [visualizations.R] Cluster-specific pathway analysis
    │       ↓
    │   Treatment × Cluster pathway matrix
    │
    └─→ [Psudotime.R] Trajectory & temporal analysis
            ↓
        Pathway dynamics along pseudotime
```

## 📚 Documentation

- **[Technical_Documentation.md](Technical_Documentation.md)** - Complete technical documentation of all analyses
- **[REPRODUCIBILITY_GUIDE.md](REPRODUCIBILITY_GUIDE.md)** - Comprehensive reproducibility guide
- **[REPRODUCIBILITY_QUICKREF.md](REPRODUCIBILITY_QUICKREF.md)** - Quick reference for reproducibility

## 📊 Data Source

This analysis uses the **GSE164897** dataset from GEO (Gene Expression Omnibus):
- **Organism**: *Homo sapiens* (Human)
- **Cell Type**: Melanoma cells
- **Platform**: Single-cell RNA sequencing
- **Treatments**: BRAF/MEK inhibitor combinations

## 🔑 Key Insights

The analysis pipeline enables:
- **Pathway-level understanding** of drug responses
- **Heterogeneity detection** across cell populations
- **Temporal dynamics** of pathway activation
- **Cross-treatment comparison** of therapeutic effects
- **Cluster-specific responses** to different drug combinations

## 📝 Citation

If you use this code, please cite:
```
[Your Citation Information]
```

## 📧 Contact

**Maintainer**: Bharadwaja Nandivada  
**Repository**: https://github.com/SreeSatyaGit/Single-Cell-Melanoma-Pathway-Analysis

## 📄 License

[Specify License]

---

**Last Updated**: December 1, 2025  
**Version**: 1.0