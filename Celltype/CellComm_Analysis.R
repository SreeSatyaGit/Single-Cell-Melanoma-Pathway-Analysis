# =============================================================================
# CELL-CELL COMMUNICATION ANALYSIS
# =============================================================================
# Analyzes intercellular communication networks in melanoma resistance
# Uses CellChat to identify ligand-receptor interactions
# Compares communication patterns across treatments
# =============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(CellChat)
library(ComplexHeatmap)
library(circlize)

# Load reproducibility configuration
source("/projects/vanaja_lab/satya/SCPA/Reproduce.R")
set_analysis_seed("cellcomm_analysis")

cat("\n=============================================================================\n")
cat("CELL-CELL COMMUNICATION ANALYSIS\n")
cat("=============================================================================\n\n")

# Create output directory
dir.create("CellComm_Results", showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA AND PREPARE FOR CELLCHAT
# =============================================================================

cat("=== 1. LOADING DATA ===\n")

# Load GSE164897
if (!exists("GSE164897")) {
  cat("  Loading GSE164897...\n")
  source("GSE164897.R")
}

cat(paste("  Loaded:", ncol(GSE164897), "cells\n"))
cat(paste("  Cell types:", length(unique(GSE164897$celltype)), "\n"))
cat(paste("  Treatments:", length(unique(GSE164897$treatment)), "\n"))

# =============================================================================
# 2. LOAD CELLCHAT DATABASE
# =============================================================================

cat("\n=== 2. LOADING CELLCHAT DATABASE ===\n")

# Load human ligand-receptor database
CellChatDB <- CellChatDB.human

# Show database info
cat(paste("  Total interactions:", nrow(CellChatDB$interaction), "\n"))
cat(paste("  Interaction types:\n"))
print(table(CellChatDB$interaction$annotation))

# Use all interactions (can be filtered if needed)
CellChatDB.use <- CellChatDB

cat("  Using full database for analysis\n")

# =============================================================================
# 3. CREATE CELLCHAT OBJECTS FOR EACH TREATMENT
# =============================================================================

cat("\n=== 3. CREATING CELLCHAT OBJECTS ===\n")

treatments <- c("untreated", "Vemurafenib", "vem_cob", "vem_tram")
cellchat_objects <- list()

for (trt in treatments) {
  cat(paste("\n  Processing:", trt, "\n"))
  
  # Subset to treatment
  seurat_subset <- subset(GSE164897, treatment == trt)
  
  # Ensure normalized data exists (fixes "data layer not found" warnings)
  seurat_subset <- NormalizeData(seurat_subset, verbose = FALSE)
  
  cat(paste("    Cells:", ncol(seurat_subset), "\n"))
  
  # Extract data for CellChat
  data.input <- GetAssayData(seurat_subset, layer = "data", assay = "RNA")
  
  # Create metadata
  meta <- data.frame(
    labels = seurat_subset$celltype,
    row.names = colnames(seurat_subset)
  )
  
  # Create CellChat object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  
  # Set database
  cellchat@DB <- CellChatDB.use
  
  # Preprocess - subset data to signaling genes
  cellchat <- subsetData(cellchat)
  
  # Identify overexpressed genes and interactions
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Compute communication probability
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
  
  # Filter out low-quality communications
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Infer signaling pathways
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Calculate aggregated network
  cellchat <- aggregateNet(cellchat)
  
  # Compute network centrality (required for computeNetSimilarity)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  cellchat_objects[[trt]] <- cellchat
  
  cat(paste("    Identified", nrow(cellchat@net$count), "interactions\n"))
}

cat("\n  CellChat objects created for all treatments\n")

# =============================================================================
# 4. ANALYZE OVERALL COMMUNICATION PATTERNS
# =============================================================================

cat("\n=== 4. OVERALL COMMUNICATION ANALYSIS ===\n")

# For each treatment, analyze communication
for (trt in treatments) {
  cat(paste("\n  Analyzing:", trt, "\n"))
  
  cellchat <- cellchat_objects[[trt]]
  
  # Number of interactions
  n_interactions <- sum(cellchat@net$count)
  n_pathways <- length(cellchat@netP$pathways)
  
  cat(paste("    Total interactions:", n_interactions, "\n"))
  cat(paste("    Active pathways:", n_pathways, "\n"))
  
  # Interaction strength
  total_strength <- sum(cellchat@net$weight)
  cat(paste("    Total interaction strength:", round(total_strength, 2), "\n"))
}

# =============================================================================
# 5. VISUALIZE COMMUNICATION NETWORKS
# =============================================================================

cat("\n=== 5. CREATING NETWORK VISUALIZATIONS ===\n")

# --- 5.1 Circle plots for each treatment ---
cat("  Creating circle plots...\n")

for (trt in treatments) {
  cellchat <- cellchat_objects[[trt]]
  
  # Create circle plot
  pdf(file.path("CellComm_Results", paste0("CirclePlot_", trt, ".pdf")),
      width = 10, height = 10)
  
  # Get cell counts per group
  group.size <- as.numeric(table(cellchat@idents))
  
  netVisual_circle(cellchat@net$count, 
                   vertex.weight = group.size, 
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = paste("Number of Interactions -", trt))
  
  dev.off()
  
  cat(paste("    Saved circle plot for", trt, "\n"))
}

# --- 5.2 Heatmaps of interaction counts ---
cat("  Creating interaction heatmaps...\n")

for (trt in treatments) {
  cellchat <- cellchat_objects[[trt]]
  
  pdf(file.path("CellComm_Results", paste0("Heatmap_Interactions_", trt, ".pdf")),
      width = 8, height = 7)
  
  netVisual_heatmap(cellchat, color.heatmap = "Reds")
  
  dev.off()
}

# =============================================================================
# 6. IDENTIFY KEY SIGNALING PATHWAYS
# =============================================================================

cat("\n=== 6. IDENTIFYING KEY SIGNALING PATHWAYS ===\n")

pathway_analysis <- list()

for (trt in treatments) {
  cat(paste("\n  Analyzing pathways in:", trt, "\n"))
  
  cellchat <- cellchat_objects[[trt]]
  
  # Get pathway information
  pathways <- cellchat@netP$pathways
  
  if (length(pathways) > 0) {
    # Calculate pathway strength
    pathway_strength <- sapply(pathways, function(pathway) {
      sum(cellchat@netP$prob[,,pathway], na.rm = TRUE)
    })
    
    # Sort by strength
    pathway_strength <- sort(pathway_strength, decreasing = TRUE)
    
    # Store results
    pathway_df <- data.frame(
      pathway = names(pathway_strength),
      strength = as.numeric(pathway_strength),
      treatment = trt,
      stringsAsFactors = FALSE
    )
    
    pathway_analysis[[trt]] <- pathway_df
    
    # Print top 10
    cat("    Top 10 pathways:\n")
    for (i in 1:min(10, length(pathway_strength))) {
      cat(sprintf("      %2d. %-20s (strength: %.3f)\n", 
                  i, names(pathway_strength)[i], pathway_strength[i]))
    }
  }
}

# Combine pathway results
all_pathways <- bind_rows(pathway_analysis)
write.csv(all_pathways,
          file.path("CellComm_Results", "Pathway_Strengths_All_Treatments.csv"),
          row.names = FALSE)

# =============================================================================
# 7. COMPARE TREATMENTS
# =============================================================================

cat("\n=== 7. COMPARING TREATMENTS ===\n")

# Merge CellChat objects for comparison
cat("  Merging CellChat objects...\n")

cellchat_merged <- mergeCellChat(cellchat_objects, add.names = treatments)

cat("  Comparison analysis...\n")

# --- 7.1 Compare number of interactions ---
pdf(file.path("CellComm_Results", "Comparison_Interaction_Numbers.pdf"),
    width = 10, height = 6)

compareInteractions(cellchat_merged, show.legend = TRUE, group = 1:4)

dev.off()

# --- 7.2 Compare interaction strength ---
pdf(file.path("CellComm_Results", "Comparison_Interaction_Strength.pdf"),
    width = 10, height = 6)

compareInteractions(cellchat_merged, show.legend = TRUE, group = 1:4, measure = "weight")

dev.off()

# --- 7.3 Differential interaction analysis ---
cat("  Identifying differential signaling pathways...\n")

# Compare treated vs untreated using rankNet
for (trt in c("Vemurafenib", "vem_cob", "vem_tram")) {
  
  # Get position of treatments in merged object
  pos_untreated <- which(treatments == "untreated")
  pos_treated <- which(treatments == trt)
  
  cat(paste("    Comparing", trt, "vs untreated...\n"))
  
  tryCatch({
    # 1. Compare signaling pathway information flow
    # Save plot
    pdf(file.path("CellComm_Results", paste0("Differential_Pathways_", trt, "_vs_untreated.pdf")), width = 8, height = 6)
    p_path <- rankNet(cellchat_merged, 
                     mode = "comparison", 
                     stacked = FALSE,
                     comparison = c(pos_untreated, pos_treated),
                     do.stat = TRUE)
    print(p_path)
    dev.off()
    
    # Save data
    df_path <- rankNet(cellchat_merged, 
                      mode = "comparison", 
                      stacked = FALSE,
                      comparison = c(pos_untreated, pos_treated),
                      do.stat = TRUE,
                      return.data = TRUE)
                      
    if (is.data.frame(df_path)) {
       write.csv(df_path,
                 file.path("CellComm_Results", 
                           paste0("Differential_Pathways_", trt, "_vs_untreated.csv")),
                 row.names = FALSE)
    }

    # 2. Compare interaction strength (L-R pairs)
    # Save plot
    pdf(file.path("CellComm_Results", paste0("Differential_LR_Pairs_", trt, "_vs_untreated.pdf")), width = 8, height = 10)
    p_lr <- rankNet(cellchat_merged,
                   mode = "comparison",
                   stacked = FALSE,
                   comparison = c(pos_untreated, pos_treated),
                   slot.name = "net",
                   do.stat = TRUE)
    print(p_lr)
    dev.off()
                   
    # Save data
    df_lr <- rankNet(cellchat_merged,
                    mode = "comparison",
                    stacked = FALSE,
                    comparison = c(pos_untreated, pos_treated),
                    slot.name = "net",
                    do.stat = TRUE,
                    return.data = TRUE)
    
    if (is.data.frame(df_lr)) {
      write.csv(df_lr,
                file.path("CellComm_Results",
                          paste0("Differential_LR_Pairs_", trt, "_vs_untreated.csv")),
                row.names = FALSE)
    }
    
  }, error = function(e) {
    cat(paste("      Error in differential analysis for", trt, ":", e$message, "\n"))
  })
}

# =============================================================================
# 8. CELL TYPE-SPECIFIC COMMUNICATION
# =============================================================================

cat("\n=== 8. CELL TYPE-SPECIFIC COMMUNICATION ===\n")

# Focus on resistance-related cell types
resistance_celltypes <- c("Invasive", "Resistant", "Hypoxic", "ImmuneEvasion")
resistance_celltypes <- intersect(resistance_celltypes, unique(GSE164897$celltype))

cat(paste("  Analyzing communication for:", 
          paste(resistance_celltypes, collapse = ", "), "\n"))

# For each treatment, analyze communication involving resistance cell types
celltype_comm_summary <- list()

for (trt in treatments) {
  cellchat <- cellchat_objects[[trt]]
  
  for (ct in resistance_celltypes) {
    # Incoming signals to this cell type
    incoming <- subsetCommunication(cellchat, targets.use = ct)
    
    # Outgoing signals from this cell type
    outgoing <- subsetCommunication(cellchat, sources.use = ct)
    
    if (nrow(incoming) > 0 || nrow(outgoing) > 0) {
      summary_df <- data.frame(
        treatment = trt,
        celltype = ct,
        n_incoming = nrow(incoming),
        n_outgoing = nrow(outgoing),
        incoming_strength = sum(incoming$prob, na.rm = TRUE),
        outgoing_strength = sum(outgoing$prob, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      
      celltype_comm_summary[[paste(trt, ct, sep = "_")]] <- summary_df
    }
  }
}

# Combine results
celltype_comm_df <- bind_rows(celltype_comm_summary)
write.csv(celltype_comm_df,
          file.path("CellComm_Results", "CellType_Communication_Summary.csv"),
          row.names = FALSE)

cat("\n  Cell type communication summary:\n")
print(celltype_comm_df)

# =============================================================================
# 9. LIGAND-RECEPTOR PAIR ANALYSIS
# =============================================================================

cat("\n=== 9. LIGAND-RECEPTOR PAIR ANALYSIS ===\n")

# Identify top L-R pairs for each treatment
top_lr_pairs <- list()

for (trt in treatments) {
  cat(paste("\n  Top L-R pairs in:", trt, "\n"))
  
  cellchat <- cellchat_objects[[trt]]
  
  # Get all L-R pairs
  lr_pairs <- subsetCommunication(cellchat)
  
  if (nrow(lr_pairs) > 0) {
    # Sort by probability
    lr_pairs <- lr_pairs %>%
      arrange(desc(prob)) %>%
      mutate(treatment = trt)
    
    top_lr_pairs[[trt]] <- head(lr_pairs, 20)
    
    # Print top 10
    cat("    Top 10 L-R pairs:\n")
    for (i in 1:min(10, nrow(lr_pairs))) {
      cat(sprintf("      %2d. %s -> %s: %s (prob: %.3f)\n",
                  i,
                  lr_pairs$source[i],
                  lr_pairs$target[i],
                  lr_pairs$interaction_name[i],
                  lr_pairs$prob[i]))
    }
  }
}

# Combine and save
all_lr_pairs <- bind_rows(top_lr_pairs)
write.csv(all_lr_pairs,
          file.path("CellComm_Results", "Top_LR_Pairs_All_Treatments.csv"),
          row.names = FALSE)

# =============================================================================
# 10. PATHWAY-SPECIFIC VISUALIZATIONS
# =============================================================================

cat("\n=== 10. PATHWAY-SPECIFIC VISUALIZATIONS ===\n")

# Select key pathways to visualize
key_pathways <- c("VEGF", "TGFb", "NOTCH", "WNT", "EGF", "PDGF", "FGF", "IGF")

for (pathway in key_pathways) {
  cat(paste("  Visualizing pathway:", pathway, "\n"))
  
  # Check which treatments have this pathway
  for (trt in treatments) {
    cellchat <- cellchat_objects[[trt]]
    
    if (pathway %in% cellchat@netP$pathways) {
      
      # Create pathway visualization
      pdf(file.path("CellComm_Results", 
                    paste0("Pathway_", pathway, "_", trt, ".pdf")),
          width = 8, height = 6)
      
      tryCatch({
        netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
      }, error = function(e) {
        cat(paste("    Error visualizing", pathway, "in", trt, "\n"))
      })
      
      dev.off()
    }
  }
}

# =============================================================================
# 10b. ADVANCED VISUALIZATIONS
# =============================================================================

cat("\n=== 10b. CREATING ADVANCED VISUALIZATIONS ===\n")

for (trt in treatments) {
  cellchat <- cellchat_objects[[trt]]
  
  # 1. Bubble plots of significant interactions
  cat(paste("  Creating bubble plots for", trt, "...\n"))
  pdf(file.path("CellComm_Results", paste0("BubblePlot_All_", trt, ".pdf")), width = 12, height = 15)
  tryCatch({
      print(netVisual_bubble(cellchat, remove.isolate = FALSE))
  }, error = function(e) { cat(paste("    Error (Bubble):", e$message, "\n")) })
  dev.off()
  
  # 2. Signaling Role Scatter Plot
  cat(paste("  Creating signaling role scatter plots for", trt, "...\n"))
  pdf(file.path("CellComm_Results", paste0("SignalingRole_Scatter_", trt, ".pdf")), width = 10, height = 8)
  tryCatch({
       p <- netAnalysis_signalingRole_scatter(cellchat)
       print(p)
  }, error = function(e) { cat(paste("    Error (Scatter):", e$message, "\n")) })
  dev.off()
  
  # 3. Violin plot of signaling gene expression
  if (length(cellchat@netP$pathways) > 0) {
      cat(paste("  Creating gene expression violin plots for", trt, "...\n"))
      top_paths <- cellchat@netP$pathways[1:min(5, length(cellchat@netP$pathways))]
      
      for (path in top_paths) {
        pdf(file.path("CellComm_Results", paste0("ViolinPlot_GeneExpr_", path, "_", trt, ".pdf")), width = 8, height = 6)
        tryCatch({
          print(plotGeneExpression(cellchat, signaling = path))
        }, error = function(e) { cat(paste("    Error (Violin for", path, "):", e$message, "\n")) })
        dev.off()
      }
  }
}

# =============================================================================
# 11. RESISTANCE-SPECIFIC COMMUNICATION PATTERNS
# =============================================================================

cat("\n=== 11. RESISTANCE-SPECIFIC COMMUNICATION ===\n")

# Compare communication in resistant vs sensitive cells
if ("pseudotime_clinical" %in% colnames(GSE164897@meta.data)) {
  
  cat("  Analyzing communication by resistance stage...\n")
  
  # Create CellChat objects for resistant vs sensitive
  resistant_cells <- subset(GSE164897, pseudotime_clinical == "Resistant")
  sensitive_cells <- subset(GSE164897, pseudotime_clinical == "Sensitive")
  
  # Process resistant cells
  if (ncol(resistant_cells) > 100) {
    cat("    Processing resistant cells...\n")
    
    data.input <- GetAssayData(resistant_cells, layer = "data", assay = "RNA")
    meta <- data.frame(labels = resistant_cells$celltype,
                       row.names = colnames(resistant_cells))
    
    cellchat_resistant <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    cellchat_resistant@DB <- CellChatDB.use
    cellchat_resistant <- subsetData(cellchat_resistant)
    cellchat_resistant <- identifyOverExpressedGenes(cellchat_resistant)
    cellchat_resistant <- identifyOverExpressedInteractions(cellchat_resistant)
    cellchat_resistant <- computeCommunProb(cellchat_resistant)
    cellchat_resistant <- filterCommunication(cellchat_resistant, min.cells = 10)
    cellchat_resistant <- computeCommunProbPathway(cellchat_resistant)
    cellchat_resistant <- aggregateNet(cellchat_resistant)
    cellchat_resistant <- netAnalysis_computeCentrality(cellchat_resistant, slot.name = "netP")
    
    # Save resistant communication network
    pdf(file.path("CellComm_Results", "CirclePlot_Resistant_Cells.pdf"),
        width = 10, height = 10)
    
    # Get cell counts per group
    group.size <- as.numeric(table(cellchat_resistant@idents))
    
    netVisual_circle(cellchat_resistant@net$count,
                     vertex.weight = group.size,
                     weight.scale = TRUE,
                     title.name = "Resistant Cells Communication")
    dev.off()
  }
  
  # Process sensitive cells
  if (ncol(sensitive_cells) > 100) {
    cat("    Processing sensitive cells...\n")
    
    data.input <- GetAssayData(sensitive_cells, layer = "data", assay = "RNA")
    meta <- data.frame(labels = sensitive_cells$celltype,
                       row.names = colnames(sensitive_cells))
    
    cellchat_sensitive <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    cellchat_sensitive@DB <- CellChatDB.use
    cellchat_sensitive <- subsetData(cellchat_sensitive)
    cellchat_sensitive <- identifyOverExpressedGenes(cellchat_sensitive)
    cellchat_sensitive <- identifyOverExpressedInteractions(cellchat_sensitive)
    cellchat_sensitive <- computeCommunProb(cellchat_sensitive)
    cellchat_sensitive <- filterCommunication(cellchat_sensitive, min.cells = 10)
    cellchat_sensitive <- computeCommunProbPathway(cellchat_sensitive)
    cellchat_sensitive <- aggregateNet(cellchat_sensitive)
    cellchat_sensitive <- netAnalysis_computeCentrality(cellchat_sensitive, slot.name = "netP")
    
    # Save sensitive communication network
    pdf(file.path("CellComm_Results", "CirclePlot_Sensitive_Cells.pdf"),
        width = 10, height = 10)
    
    # Get cell counts per group
    group.size <- as.numeric(table(cellchat_sensitive@idents))
    
    netVisual_circle(cellchat_sensitive@net$count,
                     vertex.weight = group.size,
                     weight.scale = TRUE,
                     title.name = "Sensitive Cells Communication")
    dev.off()
  }
}

# =============================================================================
# 12. SUMMARY STATISTICS
# =============================================================================

cat("\n=== 12. SUMMARY STATISTICS ===\n")

# Create comprehensive summary
summary_stats <- data.frame(
  Treatment = character(),
  Total_Interactions = integer(),
  Total_Strength = numeric(),
  Active_Pathways = integer(),
  Avg_Interactions_Per_CellType = numeric(),
  stringsAsFactors = FALSE
)

for (trt in treatments) {
  cellchat <- cellchat_objects[[trt]]
  
  n_celltypes <- length(unique(cellchat@meta$labels))
  
  summary_stats <- rbind(summary_stats, data.frame(
    Treatment = trt,
    Total_Interactions = sum(cellchat@net$count),
    Total_Strength = sum(cellchat@net$weight),
    Active_Pathways = length(cellchat@netP$pathways),
    Avg_Interactions_Per_CellType = sum(cellchat@net$count) / n_celltypes
  ))
}

cat("\n  Communication Summary:\n")
print(summary_stats)

write.csv(summary_stats,
          file.path("CellComm_Results", "Communication_Summary_Statistics.csv"),
          row.names = FALSE)

# =============================================================================
# 13. SAVE CELLCHAT OBJECTS
# =============================================================================

cat("\n=== 13. SAVING CELLCHAT OBJECTS ===\n")

# Save all CellChat objects
saveRDS(cellchat_objects, 
        file.path("CellComm_Results", "CellChat_Objects_All_Treatments.rds"))

cat("  CellChat objects saved\n")

# Save session info
save_session_info(file.path("CellComm_Results", "session_info_cellcomm.txt"))

# =============================================================================
# 14. FINAL SUMMARY
# =============================================================================

cat("\n=============================================================================\n")
cat("✓ CELL-CELL COMMUNICATION ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
cat("\nSummary:\n")
cat(paste("  Treatments analyzed:", length(treatments), "\n"))
cat(paste("  Total CellChat objects created:", length(cellchat_objects), "\n"))
cat("\nKey findings:\n")
for (i in 1:nrow(summary_stats)) {
  cat(sprintf("  %s: %d interactions, %d pathways\n",
              summary_stats$Treatment[i],
              summary_stats$Total_Interactions[i],
              summary_stats$Active_Pathways[i]))
}
cat("\nKey outputs:\n")
cat("  - Circle plots for each treatment\n")
cat("  - Interaction heatmaps\n")
cat("  - Pathway strength analysis\n")
cat("  - Treatment comparisons\n")
cat("  - Cell type-specific communication\n")
cat("  - Top L-R pairs\n")
cat("  - Pathway-specific visualizations\n")
cat("  - Resistance stage communication\n")
cat(paste("\nResults saved in:", file.path(getwd(), "CellComm_Results"), "\n"))
cat("=============================================================================\n\n")
