
library(SCPA)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(statmod)
library(dplyr)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)
source("/projects/vanaja_lab/satya/SCPA/Reproduce.R")
# Check and convert gene IDs if needed
print("=== Converting Ensembl IDs to Gene Symbols ===")
current_genes <- rownames(GSE164897)
print(paste("Current genes in Seurat object:", length(current_genes)))
print("Sample genes:")
print(head(current_genes, 10))

# Check if Ensembl IDs
is_ensembl <- grepl("^ENSG", head(current_genes, 100))
ensembl_pct <- sum(is_ensembl) / length(is_ensembl) * 100

if (ensembl_pct > 50) {
  print(paste("Detected", round(ensembl_pct, 1), "% Ensembl IDs - converting to gene symbols..."))
  
  # Join layers first if needed
  tryCatch({
    GSE164897[["RNA"]] <- JoinLayers(GSE164897[["RNA"]])
    print("Joined RNA layers")
  }, error = function(e) {
    print("Layers already joined or not needed")
  })
  
  # Map Ensembl to symbols
  ensembl_ids <- rownames(GSE164897)
  gene_symbols <- mapIds(x = org.Hs.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  
  # Make unique and remove NAs
  gene_symbols <- make.unique(as.character(gene_symbols))
  
  # Update rownames
  rownames(GSE164897) <- gene_symbols
  
  # Remove NA genes
  GSE164897 <- GSE164897[!is.na(rownames(GSE164897)), ]
  
  print(paste("Conversion complete. Genes after conversion:", nrow(GSE164897)))
  print("New gene names:")
  print(head(rownames(GSE164897), 10))
} else {
  print("Gene symbols already present - no conversion needed")
}

DimPlot(GSE164897, split.by = "treatment") +
  theme(aspect.ratio = 1)



untreated <- seurat_extract(GSE164897,
                            meta1 = "treatment", value_meta1 = "untreated")

vem_cob <- seurat_extract(GSE164897,
                          meta1 = "treatment", value_meta1 = "vem_cob")

vem_tram <- seurat_extract(GSE164897,
                           meta1 = "treatment", value_meta1 = "vem_tram")

vem <- seurat_extract(GSE164897,
                      meta1 = "treatment", value_meta1 = "Vemurafenib")

# Debug: Check structure of extracted samples
print("Checking extracted sample structure:")
print(paste("untreated class:", class(untreated)))
print(paste("untreated dimensions:", paste(dim(untreated), collapse = " x ")))
print(paste("vem class:", class(vem)))
print(paste("vem dimensions:", paste(dim(vem), collapse = " x ")))

# Debug: Check gene names in Seurat object
print("\n=== Checking gene names in GSE164897 ===")
seurat_genes <- rownames(GSE164897)
print(paste("Total genes in Seurat object:", length(seurat_genes)))
print("First 10 genes in Seurat object:")
print(head(seurat_genes, 10))

# Check if they're Ensembl IDs or gene symbols
is_ensembl <- grepl("^ENSG", head(seurat_genes, 100))
print(paste("Percentage Ensembl IDs:", round(sum(is_ensembl) / length(is_ensembl) * 100, 1), "%"))

# Load Hallmark pathways from msigdbr  
hallmark_df <- msigdbr(species = "Homo sapiens", 
                       category = "H")

print(paste("\nTotal Hallmark gene sets retrieved:", length(unique(hallmark_df$gs_name))))

# Use SCPA's format_pathways helper function to convert msigdbr output
pathways <- format_pathways(hallmark_df)

print(paste("Formatted", length(pathways), "Hallmark pathways for SCPA"))
print("First 10 genes in first pathway:")
print(head(pathways[[1]], 10))

# Check overlap between Seurat genes and pathway genes
all_pathway_genes <- unique(unlist(pathways))
overlap <- intersect(seurat_genes, all_pathway_genes)
print(paste("\nGene overlap check:"))
print(paste("  Genes in pathways:", length(all_pathway_genes)))
print(paste("  Genes in Seurat object:", length(seurat_genes)))
print(paste("  Overlapping genes:", length(overlap)))
print(paste("  Overlap percentage:", round(length(overlap) / length(seurat_genes) * 100, 2), "%"))




