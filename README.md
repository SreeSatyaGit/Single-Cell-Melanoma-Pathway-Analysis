# Single-Cell Multi-Pathway Analysis (SCMPA)

This repository contains R scripts for analyzing single-cell RNA sequencing (scRNA-seq) data with a focus on pathway analysis and drug treatment effects. The analysis pipeline processes data from the GSE164897 dataset, which examines melanoma cells under different drug treatments.

## Overview

The analysis consists of three main components that work together to provide comprehensive insights into cellular responses to drug treatments:

1. **Data Processing and Integration** (`GSE164897.R`)
2. **Pathway Analysis** (`PathwayAnalysis.R`) 
3. **AUCell Scoring** (`AUcell.R`)

## File Descriptions

### 1. GSE164897.R - Data Processing and Integration

This script handles the complete preprocessing pipeline for the GSE164897 dataset:

**Key Functions:**
- **Data Loading**: Reads multiple gzipped count files from the dataset directory
- **Seurat Object Creation**: Converts raw count matrices into Seurat objects for each sample
- **Quality Control**: Applies standard filtering (min.cells = 3, min.features = 200)
- **Normalization**: Performs data normalization and scaling
- **Dimensionality Reduction**: Runs PCA analysis
- **Clustering**: Identifies cell clusters using unintegrated data
- **Data Integration**: Uses CCA (Canonical Correlation Analysis) integration to merge samples
- **UMAP Visualization**: Creates 2D embeddings for visualization

**Treatment Groups:**
- `untreated`: Control cells
- `Vemurafenib`: BRAF inhibitor treatment
- `vem_cob`: Vemurafenib + Cobimetinib combination
- `vem_tram`: Vemurafenib + Trametinib combination

**Output**: Integrated Seurat object (`GSE164897`) ready for downstream analysis

### 2. PathwayAnalysis.R - Gene Set Enrichment Analysis

This script performs comprehensive pathway analysis using multiple approaches:

**Key Functions:**
- **Gene ID Conversion**: Maps Ensembl gene IDs to gene symbols using `org.Hs.eg.db`
- **SCTransform**: Applies variance stabilizing transformation for improved normalization
- **PCA Analysis**: Runs principal component analysis with reversed PCA for pathway scoring
- **GESECA Analysis**: Performs Gene Set Enrichment Analysis for Cellular Components and Activities using `fgsea`
- **Pathway Filtering**: Focuses on oncogenic pathways (AKT, RAF, KRAS, MTOR signaling)
- **Visualization**: Creates coregulation profile plots showing pathway activity across UMAP space

**Key Libraries Used:**
- `msigdbr`: Gene set collections (C6 oncogenic signatures)
- `fgsea`: Fast gene set enrichment analysis
- `sctransform`: Single-cell transformation
- `AnnotationDbi`: Gene annotation database

**Output**: 
- GESECA results table with pathway enrichment scores
- Coregulation profile plots showing pathway activity patterns

### 3. AUcell.R - AUCell Scoring Analysis

This script implements AUCell (Area Under the Curve) analysis to quantify pathway activity:

**Key Functions:**
- **Ranking Generation**: Builds gene expression rankings for each cell using `AUCell_buildRankings`
- **AUC Calculation**: Computes AUC scores for each pathway in each cell
- **Global Pathway Activity**: Identifies the top 10 most active pathways across all cells
- **Visualization**: Creates horizontal bar plots showing mean AUCell scores

**Key Features:**
- Uses multi-core processing (4 cores) for efficiency
- Focuses on globally active pathways rather than cell-specific patterns
- Provides quantitative scoring of pathway activity levels

**Output**: 
- AUC matrix (pathways × cells)
- Bar plot visualization of top 10 highly expressed pathways

## Workflow Integration

The three scripts work together in a sequential pipeline:

1. **GSE164897.R** → Processes raw data and creates integrated Seurat object
2. **PathwayAnalysis.R** → Performs gene set enrichment analysis on the processed data
3. **AUcell.R** → Quantifies pathway activity using AUCell scoring

## Dependencies

### Required R Packages:
- `Seurat` - Single-cell analysis framework
- `AUCell` - Area under the curve analysis
- `msigdbr` - Gene set collections
- `fgsea` - Fast gene set enrichment analysis
- `AnnotationDbi` - Annotation database interface
- `org.Hs.eg.db` - Human gene annotation database
- `sctransform` - Single-cell transformation
- `glmGamPoi` - Generalized linear model for single-cell data
- `data.table` - Fast data manipulation
- `Matrix` - Sparse matrix operations
- `ggplot2` - Data visualization
- `cowplot` - Plot arrangement

## Usage

1. Ensure all required packages are installed
2. Update the `base_dir` path in `GSE164897.R` to point to your data directory
3. Run scripts in order: `GSE164897.R` → `PathwayAnalysis.R` → `AUcell.R`
4. Results will include integrated Seurat objects, pathway enrichment scores, and visualization plots

## Data Source

This analysis uses the GSE164897 dataset, which contains single-cell RNA-seq data from melanoma cells treated with various drug combinations targeting the MAPK pathway.

## Key Insights

The analysis pipeline enables:
- Identification of treatment-specific cellular responses
- Quantification of oncogenic pathway activity
- Visualization of pathway co-regulation patterns
- Comparison of drug treatment effects at the single-cell level