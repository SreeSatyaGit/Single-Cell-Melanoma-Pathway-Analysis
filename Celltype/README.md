# Cell Type Annotation for GSE164897 Melanoma Dataset

This directory contains scripts for annotating cell types and melanoma cell states in the GSE164897 dataset.

## 📋 Overview

Two complementary approaches are provided:

1. **SingleR.R** - Automated reference-based cell type annotation
2. **MelanomaStates.R** - Melanoma-specific cell state profiling

## 🚀 Quick Start

### Prerequisites

Make sure you have run `GSE164897.R` first to create the Seurat object.

### Running the Analysis

```r
# Option 1: Run SingleR for broad cell type annotation
source("Celltype/SingleR.R")

# Option 2: Run melanoma state annotation (recommended for drug resistance analysis)
source("Celltype/MelanomaStates.R")

# Option 3: Run both for comprehensive annotation
source("Celltype/SingleR.R")
source("Celltype/MelanomaStates.R")
```

## 📁 File Descriptions

### 1. SingleR.R - Reference-Based Cell Type Annotation

**Purpose:** Identifies broad cell types using reference datasets

**Method:** SingleR algorithm with Human Primary Cell Atlas reference

**Output:**
- `GSE164897$celltype_hpca` - Cell type labels
- `GSE164897$celltype_score` - Annotation confidence scores
- UMAP plots with cell type labels
- Score heatmaps showing annotation quality

**Best for:**
- Identifying contaminating cell types (immune cells, fibroblasts, etc.)
- Quality control
- Mixed cell type datasets

**Runtime:** ~5-10 minutes (first run downloads reference data)

### 2. MelanomaStates.R - Melanoma Cell State Profiling

**Purpose:** Identifies melanoma-specific functional states relevant for drug resistance

**Method:** Module scoring with melanoma state-specific marker genes

**Cell States Identified:**
- **Differentiated** - Melanocytic/pigmented state (MLANA, TYR, DCT)
- **Undifferentiated** - Neural crest-like state (SOX10, MITF, PAX3)
- **Invasive** - Mesenchymal/resistant state (AXL, NGFR, VIM)
- **Proliferative** - Actively dividing cells (MKI67, TOP2A)
- **Resistant** - Adaptive resistance state (EGFR, PDGFRB, AXL)
- **Hypoxic** - Stress response state (HIF1A, VEGFA)
- **ImmuneEvasion** - Immune checkpoint expression (CD274, PDCD1LG2)

**Output:**
- `GSE164897$dominant_state` - Dominant cell state per cell
- `GSE164897$<State>_score` - Individual state scores
- `melanoma_state_annotations.csv` - Full annotation table
- `melanoma_state_summary.rds` - Summary statistics
- Multiple visualizations (UMAPs, heatmaps, violin plots)

**Best for:**
- Drug resistance analysis
- Understanding treatment effects
- Melanoma-specific biology
- **Recommended for this dataset**

**Runtime:** ~3-5 minutes

## 🔬 Which Method Should I Use?

### For GSE164897 Melanoma Dataset:

**Primary Recommendation: MelanomaStates.R**

This dataset contains melanoma cells treated with BRAF/MEK inhibitors. The most biologically relevant annotation is **melanoma cell states** rather than broad cell types.

**Why MelanomaStates.R is better for this analysis:**
- ✅ Identifies drug resistance mechanisms (invasive/resistant states)
- ✅ Tracks treatment-induced state transitions
- ✅ Directly relevant to BRAF/MEK inhibitor biology
- ✅ Integrates with pathway analysis (SCPA)
- ✅ Provides actionable biological insights

**When to also use SingleR.R:**
- If you suspect contaminating cell types
- For quality control
- To validate dataset purity
- For publication (shows comprehensive annotation)

### Workflow Recommendation:

```r
# Step 1: Quick check with SingleR (optional but recommended)
source("Celltype/SingleR.R")
table(GSE164897$celltype_hpca)  # Check for non-melanoma cells

# Step 2: Melanoma state annotation (essential)
source("Celltype/MelanomaStates.R")

# Step 3: Integrate with pathway analysis
source("SCPA.R")  # Now you can analyze pathways by cell state
```

## 📊 Integration with Existing Workflow

### With SCPA Analysis

```r
# After running MelanomaStates.R, you can subset by state
library(SCPA)

# Extract cells by state
invasive_cells <- seurat_extract(GSE164897,
                                 meta1 = "dominant_state", 
                                 value_meta1 = "Invasive")

differentiated_cells <- seurat_extract(GSE164897,
                                       meta1 = "dominant_state", 
                                       value_meta1 = "Differentiated")

# Compare pathways between states
scpa_out <- compare_pathways(samples = list(invasive_cells, differentiated_cells),
                             pathways = pathways)
```

### With Visualizations

```r
# After running MelanomaStates.R
source("visualizations.R")

# Create state-specific pathway heatmaps
# Modify visualizations.R to use dominant_state instead of seurat_clusters
```

### With Pseudotime Analysis

```r
# After running MelanomaStates.R
source("Psudotime.R")

# Analyze state transitions along pseudotime
# Color trajectory by dominant_state
```

## 📈 Expected Results

### SingleR Results

For a pure melanoma dataset, you should see:
- Majority: "Melanocytes" or similar melanoma-related cell type
- Possible minor populations: Immune cells, fibroblasts (if present)
- High annotation scores (>0.7) indicate confident assignments

### MelanomaStates Results

For GSE164897 with drug treatments, expect:
- **Untreated**: Higher differentiated state
- **Vemurafenib**: Shift toward resistant/invasive states
- **Combination therapies**: Variable state distributions
- **Key insight**: Invasive/resistant states correlate with drug resistance

## 🔍 Quality Control

### Check Annotation Quality

```r
# After SingleR
summary(apply(pred_hpca$scores, 1, max))  # Should be >0.5
table(GSE164897$celltype_hpca)  # Check distribution

# After MelanomaStates
summary(GSE164897$dominant_state_score)  # Check confidence
table(GSE164897$dominant_state, GSE164897$treatment)  # Check balance
```

### Validate with Markers

```r
# Check key melanoma markers
FeaturePlot(GSE164897, features = c("MLANA", "SOX10", "AXL", "MKI67"))

# Check marker expression by state
VlnPlot(GSE164897, features = "AXL", group.by = "dominant_state")
```

## 📦 Dependencies

### Required Packages

```r
# For SingleR.R
BiocManager::install(c("SingleR", "celldex", "SingleCellExperiment"))

# For MelanomaStates.R
install.packages(c("Seurat", "ggplot2", "dplyr", "tidyr"))
BiocManager::install("ComplexHeatmap")
```

### Installation

All dependencies are automatically installed when running the scripts.

## 🐛 Troubleshooting

### Issue: BiocNeighbors error

**Solution:** The updated SingleR.R script now automatically installs all dependencies.

### Issue: "GSE164897 object not found"

**Solution:** Run `GSE164897.R` first to create the Seurat object.

### Issue: Missing marker genes

**Solution:** MelanomaStates.R automatically checks for available markers and adapts.

### Issue: Low annotation scores

**Solution:** 
- Check data normalization
- Try different reference datasets in SingleR
- For melanoma states, this is expected (states are continuous, not discrete)

## 📚 References

### SingleR
- Aran et al. (2019) "Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage" *Nature Immunology*

### Melanoma Cell States
- Rambow et al. (2018) "Toward Minimal Residual Disease-Directed Therapy in Melanoma" *Cell*
- Tirosh et al. (2016) "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq" *Science*
- Tsoi et al. (2018) "Multi-stage Differentiation Defines Melanoma Subtypes with Differential Vulnerability to Drug-Induced Iron-Dependent Oxidative Stress" *Cancer Cell*

## 💡 Tips for Analysis

1. **Start with MelanomaStates.R** - Most relevant for drug resistance
2. **Use SingleR for validation** - Confirms dataset purity
3. **Integrate with SCPA** - Analyze pathways by cell state
4. **Check treatment effects** - Compare state distributions across treatments
5. **Focus on invasive/resistant states** - Key for understanding drug resistance

## 📧 Questions?

See main repository README or contact the maintainer.

---

**Last Updated:** December 6, 2025
