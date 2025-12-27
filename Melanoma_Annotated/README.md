# GSE115978 - Melanoma scRNA-seq Dataset

## 📋 Overview

This directory contains the **GSE115978** dataset - a pre-annotated melanoma single-cell RNA-seq dataset with existing cell type labels.

## 📊 Dataset Information

**Source:** GEO (Gene Expression Omnibus) - GSE115978  
**Organism:** *Homo sapiens* (Human)  
**Cell Type:** Melanoma cells with annotations  
**Platform:** Single-cell RNA sequencing

## 📁 Files

| File | Size | Description |
|------|------|-------------|
| `GSE115978_counts.csv.gz` | ~51 MB | Raw count matrix (genes × cells) |
| `GSE115978_tpm.csv.gz` | ~238 MB | TPM normalized expression matrix |
| `GSE115978_cell.annotations.csv.gz` | ~101 KB | Cell metadata and annotations |
| `GES115978.R` | - | Processing script |

## 🚀 Quick Start

### Run the Processing Pipeline

```r
# From the SCMPA directory
source("Melanoma_Annotated/GES115978.R")

# Or from within Melanoma_Annotated directory
setwd("Melanoma_Annotated")
source("GES115978.R")
```

### What the Script Does

The processing pipeline includes:

1. **Data Loading**
   - Reads count matrix (~23,000 genes)
   - Loads cell annotations
   - Creates Seurat object

2. **Quality Control**
   - Calculates QC metrics (genes/cell, UMIs/cell, % mitochondrial)
   - Visualizes QC distributions
   - Optional filtering

3. **Normalization & Scaling**
   - Log-normalization
   - Variable feature selection
   - Data scaling

4. **Dimensionality Reduction**
   - PCA
   - UMAP

5. **Clustering**
   - Graph-based clustering
   - Visualization

6. **Existing Annotations**
   - Visualizes pre-existing cell type labels
   - Shows annotation distributions

7. **Melanoma State Annotation**
   - Adds melanoma-specific cell states
   - Module scoring for 7 states
   - State assignment

## 📈 Expected Output

After running the script, you'll have:

### Seurat Object: `GSE115978`

**Metadata columns:**
- Original annotations from the dataset
- `seurat_clusters` - Seurat clustering results
- `melanoma_state` - Melanoma cell state (Differentiated, Invasive, etc.)
- `melanoma_state_score` - Confidence score
- Individual state scores (`Differentiated_score`, `Invasive_score`, etc.)

### Visualizations

- QC metrics (violin plots, scatter plots)
- PCA plots
- Elbow plot
- UMAP by Seurat clusters
- UMAP by existing annotations
- UMAP by melanoma states

## 🔬 Melanoma Cell States

The script annotates cells with 7 melanoma-specific states:

| State | Description | Key Markers |
|-------|-------------|-------------|
| **Differentiated** | Melanocytic/pigmented | MLANA, TYR, DCT |
| **Undifferentiated** | Neural crest-like | SOX10, MITF, PAX3 |
| **Invasive** | Mesenchymal/resistant | AXL, NGFR, VIM |
| **Proliferative** | Actively dividing | MKI67, TOP2A |
| **Resistant** | Adaptive resistance | EGFR, PDGFRB |
| **Hypoxic** | Stress response | HIF1A, VEGFA |
| **ImmuneEvasion** | Immune checkpoints | CD274, PDCD1LG2 |

## 💡 Usage Examples

### Basic Exploration

```r
# After running GES115978.R

# Check cell counts
ncol(GSE115978)  # Number of cells
nrow(GSE115978)  # Number of genes

# View metadata
head(GSE115978@meta.data)
colnames(GSE115978@meta.data)

# Cell state distribution
table(GSE115978$melanoma_state)
```

### Visualizations

```r
# UMAP by melanoma state
DimPlot(GSE115978, group.by = "melanoma_state", label = TRUE)

# UMAP by original annotations (if available)
# Replace 'celltype' with actual column name
DimPlot(GSE115978, group.by = "celltype", label = TRUE)

# Feature plots for key markers
FeaturePlot(GSE115978, features = c("MLANA", "AXL", "MKI67", "SOX10"))

# Violin plots
VlnPlot(GSE115978, features = "Invasive_score", group.by = "melanoma_state")
```

### Subset Analysis

```r
# Extract specific cell states
invasive_cells <- subset(GSE115978, melanoma_state == "Invasive")
differentiated_cells <- subset(GSE115978, melanoma_state == "Differentiated")

# Find markers for a specific state
invasive_markers <- FindMarkers(GSE115978, 
                                ident.1 = "Invasive",
                                group.by = "melanoma_state")
```

### Compare with Original Annotations

```r
# Cross-tabulation (replace 'original_celltype' with actual column)
table(GSE115978$melanoma_state, GSE115978$original_celltype)

# Visualize both
p1 <- DimPlot(GSE115978, group.by = "original_celltype", label = TRUE)
p2 <- DimPlot(GSE115978, group.by = "melanoma_state", label = TRUE)

library(patchwork)
p1 | p2
```

## 🔄 Integration with Other Datasets

### Compare with GSE164897

```r
# After processing both datasets
source("GSE164897.R")
source("Melanoma_Annotated/GES115978.R")

# Compare cell state distributions
table(GSE164897$celltype)
table(GSE115978$melanoma_state)

# Merge datasets (if needed)
# Note: Requires careful batch correction
merged <- merge(GSE164897, GSE115978, 
                add.cell.ids = c("GSE164897", "GSE115978"))
```

## 📊 Downstream Analysis

### Pathway Analysis

```r
library(SCPA)

# Extract cells by state
invasive <- seurat_extract(GSE115978,
                          meta1 = "melanoma_state",
                          value_meta1 = "Invasive")

differentiated <- seurat_extract(GSE115978,
                                 meta1 = "melanoma_state",
                                 value_meta1 = "Differentiated")

# Compare pathways
pathways <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- format_pathways(pathways)

scpa_out <- compare_pathways(
  samples = list(invasive, differentiated),
  pathways = pathways
)
```

### Differential Expression

```r
# Find markers for each melanoma state
Idents(GSE115978) <- "melanoma_state"
all_markers <- FindAllMarkers(GSE115978, only.pos = TRUE)

# Top markers per state
top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Heatmap
DoHeatmap(GSE115978, features = top_markers$gene)
```

## 🐛 Troubleshooting

### Memory Issues

If you encounter memory issues with the large files:

```r
# Use data.table for faster reading
library(data.table)
counts <- fread("GSE115978_counts.csv.gz", nrows = 1000)  # Test with subset

# Or use sparse matrix from the start
library(Matrix)
```

### Cell Name Mismatches

If cell names don't match between counts and annotations:

```r
# Check cell name formats
head(colnames(counts))
head(rownames(annotations))

# Adjust as needed
colnames(counts) <- gsub("\\.", "-", colnames(counts))
```

## 📚 References

- **GEO Accession:** GSE115978
- **Publication:** [Add publication details if available]

## 📧 Questions?

See main repository README or contact the maintainer.

---

**Last Updated:** December 6, 2025
