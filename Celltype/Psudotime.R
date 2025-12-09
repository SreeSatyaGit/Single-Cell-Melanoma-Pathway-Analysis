library(SCPA)
library(Seurat)
library(tidyverse)
library(magrittr)
library(dyno)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)
library(ggplot2)
library(patchwork)

source("/projects/vanaja_lab/satya/SCPA/Reproduce.R")

# Check if cell types are available
if (!"seurat_clusters" %in% colnames(GSE164897@meta.data)) {
  stop("seurat_clusters not found! Please ensure clustering has been performed.")
}

message("\n=== Pseudotime Analysis with Cell Type Dynamics ===")

var_genes <- VariableFeatures(GSE164897)

counts <- GetAssayData(GSE164897, layer = "counts", assay = "RNA")
expression <- GetAssayData(GSE164897, layer = "data", assay = "RNA")

counts_var <- t(as.matrix(counts[var_genes, ]))
expression_var <- t(as.matrix(expression[var_genes, ]))

print("Wrapping data for dyno...")
dataset <- wrap_expression(
  counts = counts_var,
  expression = expression_var
)


print("Inferring trajectory using Slingshot...")

model <- infer_trajectory(dataset, method = "slingshot", parameters = list(cluster_method = "clara"))

# Calculate Pseudotime
print("Calculating pseudotime...")
pseudotime <- calculate_pseudotime(model)


GSE164897$pseudotime <- pseudotime[colnames(GSE164897)]

# ============================================================================
# CELL TYPE VISUALIZATIONS IN PSEUDOTIME
# ============================================================================

message("\n=== Visualizing Cell Types in Pseudotime ===")

# 1. Pseudotime colored by cell type
print("Creating pseudotime plots with cell types...")

# Get cells that are in the model
model_cells <- rownames(model$dimred)
seurat_cells <- colnames(GSE164897)

# Find common cells
common_cells <- intersect(model_cells, seurat_cells)
message(paste("Model has", length(model_cells), "cells"))
message(paste("Seurat has", length(seurat_cells), "cells"))
message(paste("Common cells:", length(common_cells)))

# Create metadata for plotting that matches model cells
celltype_for_plot <- GSE164897$seurat_clusters[match(model_cells, seurat_cells)]
treatment_for_plot <- GSE164897$treatment[match(model_cells, seurat_cells)]

p_celltype <- plot_dimred(
  model, 
  grouping = celltype_for_plot,
  plot_trajectory = TRUE, 
  size_cells = 1.5, 
  alpha_cells = 0.7
) + 
  ggtitle("Cell Types in Pseudotime") +
  theme(aspect.ratio = 1, legend.position = "right")

# Print cell type plot separately
print("Cell Types in Pseudotime:")
print(p_celltype)

# 2. Pseudotime gradient
p_pseudotime <- plot_dimred(
  model, 
  "pseudotime", 
  pseudotime = pseudotime, 
  hex_cells = FALSE,
  plot_trajectory = TRUE, 
  size_cells = 1.5, 
  alpha_cells = 0.7
) + 
  ggtitle("Pseudotime Progression") +
  theme(aspect.ratio = 1)

# 3. Treatment groups
p_treatment <- plot_dimred(
  model, 
  grouping = treatment_for_plot,
  plot_trajectory = TRUE, 
  size_cells = 1.5, 
  alpha_cells = 0.7
) + 
  ggtitle("Treatment Groups") +
  theme(aspect.ratio = 1)

# Combine only pseudotime and treatment plots
print("Pseudotime and Treatment overview:")
print(p_pseudotime + p_treatment)

# 4. Cell types split by treatment - use ggplot instead of facet_wrap on dyno plot
print("Creating cell type plots split by treatment...")

# Extract dimred coordinates from model
dimred_df <- as.data.frame(model$dimred)
dimred_df$cell <- rownames(dimred_df)

# Match with Seurat metadata
dimred_df$celltype <- celltype_for_plot
dimred_df$treatment <- treatment_for_plot

# Remove any NA values
dimred_df <- dimred_df[!is.na(dimred_df$celltype) & !is.na(dimred_df$treatment), ]

# Create faceted plot with ggplot
p_celltype_split <- ggplot(dimred_df, aes(x = comp_1, y = comp_2, color = celltype)) +
  geom_point(size = 1, alpha = 0.6) +
  facet_wrap(~treatment, ncol = 2) +
  labs(title = "Cell Types by Treatment",
       x = "Component 1", 
       y = "Component 2",
       color = "Cell Type") +
  theme_minimal() +
  theme(aspect.ratio = 1, 
        legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

print(p_celltype_split)

# ============================================================================
# CELL TYPE COMPOSITION ACROSS PSEUDOTIME
# ============================================================================

message("\n=== Analyzing Cell Type Changes Across Pseudotime ===")

# Create pseudotime bins with biologically meaningful names
# Based on melanoma drug resistance progression stages
message("\n=== Creating Biologically Meaningful Pseudotime Stages ===")

# Option 1: Detailed 10-stage progression
GSE164897$pseudotime_stage <- cut(GSE164897$pseudotime, 
                                  breaks = 10, 
                                  labels = c(
                                    "1_Naive",           # 0.0-0.1: Treatment-naive
                                    "2_Early_Response",  # 0.1-0.2: Initial drug response
                                    "3_Adaptation",      # 0.2-0.3: Early adaptation
                                    "4_Stress",          # 0.3-0.4: Stress response
                                    "5_Transition",      # 0.4-0.5: Phenotype switching begins
                                    "6_Switching",       # 0.5-0.6: Active phenotype switching
                                    "7_Pre_Resistant",   # 0.6-0.7: Pre-resistant state
                                    "8_Resistant",       # 0.7-0.8: Resistant phenotype
                                    "9_Stable_Resistant",# 0.8-0.9: Stable resistance
                                    "10_Fully_Resistant" # 0.9-1.0: Fully resistant
                                  ))

# Option 2: Simplified 5-stage progression (also create this for broader analysis)
GSE164897$pseudotime_phase <- cut(GSE164897$pseudotime, 
                                  breaks = 5, 
                                  labels = c(
                                    "Naive",              # 0.0-0.2: Treatment-naive
                                    "Early_Adaptation",   # 0.2-0.4: Early adaptation
                                    "Transition",         # 0.4-0.6: Phenotype transition
                                    "Emerging_Resistance",# 0.6-0.8: Resistance emerging
                                    "Resistant"           # 0.8-1.0: Fully resistant
                                  ))

# Option 3: Clinical 3-stage progression
GSE164897$pseudotime_clinical <- cut(GSE164897$pseudotime, 
                                     breaks = 3, 
                                     labels = c(
                                       "Sensitive",       # 0.0-0.33: Drug-sensitive
                                       "Adaptive",        # 0.33-0.67: Adaptive/transitional
                                       "Resistant"        # 0.67-1.0: Drug-resistant
                                     ))

# Display stage distributions
message("\n=== Stage Distributions ===")
message("\n10-Stage Detailed Progression:")
print(table(GSE164897$pseudotime_stage))

message("\n5-Stage Simplified Progression:")
print(table(GSE164897$pseudotime_phase))

message("\n3-Stage Clinical Classification:")
print(table(GSE164897$pseudotime_clinical))

# Use the detailed 10-stage for main analysis
message("\nUsing 10-stage detailed progression for trajectory analysis...")

# Calculate cell type proportions across pseudotime stages
celltype_pseudotime <- GSE164897@meta.data %>%
  group_by(pseudotime_stage, seurat_clusters) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(pseudotime_stage) %>%
  mutate(proportion = count / sum(count) * 100)

# Plot cell type composition across pseudotime stages
p_composition <- ggplot(celltype_pseudotime, 
                        aes(x = pseudotime_stage, y = proportion, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell Type Composition Across Resistance Progression",
       x = "Resistance Stage", 
       y = "Proportion (%)",
       fill = "Seurat Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"))

print(p_composition)

# ============================================================================
# CELL TYPE CHANGES BY TREATMENT
# ============================================================================

message("\n=== Cell Type Dynamics by Treatment ===")

# Calculate cell type proportions across pseudotime stages for each treatment
celltype_treatment_pseudotime <- GSE164897@meta.data %>%
  group_by(treatment, pseudotime_stage, seurat_clusters) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(treatment, pseudotime_stage) %>%
  mutate(proportion = count / sum(count) * 100)

# Plot cell type composition by treatment
p_composition_treatment <- ggplot(celltype_treatment_pseudotime, 
                                  aes(x = pseudotime_stage, y = proportion, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~treatment, ncol = 2) +
  labs(title = "Cell Type Dynamics Across Resistance Stages by Treatment",
       x = "Resistance Stage", 
       y = "Proportion (%)",
       fill = "Seurat Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom",
        strip.text = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

print(p_composition_treatment)

# ============================================================================
# INDIVIDUAL CELL TYPE TRAJECTORIES
# ============================================================================

message("\n=== Individual Cell Type Trajectories ===")

# Line plot showing each cell type's proportion across pseudotime
# Convert stage names to numeric for plotting
stage_order <- c("1_Naive", "2_Early_Response", "3_Adaptation", "4_Stress", 
                 "5_Transition", "6_Switching", "7_Pre_Resistant", "8_Resistant",
                 "9_Stable_Resistant", "10_Fully_Resistant")

celltype_trajectory <- celltype_pseudotime %>%
  mutate(stage_numeric = match(pseudotime_stage, stage_order))

p_trajectory <- ggplot(celltype_trajectory, 
                       aes(x = stage_numeric, y = proportion, color = seurat_clusters, group = seurat_clusters)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:10, 
                     labels = c("Naive", "Early\nResponse", "Adaptation", "Stress",
                                "Transition", "Switching", "Pre-\nResistant", "Resistant",
                                "Stable\nResistant", "Fully\nResistant")) +
  labs(title = "Cell Type Trajectories Across Resistance Progression",
       x = "Resistance Stage", 
       y = "Proportion (%)",
       color = "Seurat Cluster") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(size = 14, face = "bold"))

print(p_trajectory)

# By treatment
celltype_treatment_trajectory <- celltype_treatment_pseudotime %>%
  mutate(stage_numeric = match(pseudotime_stage, stage_order))

p_trajectory_treatment <- ggplot(celltype_treatment_trajectory, 
                                 aes(x = stage_numeric, y = proportion, color = seurat_clusters, group = seurat_clusters)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~treatment, ncol = 2) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 10), 
                     labels = c("Naive", "Adapt", "Trans", "Pre-Res", "Stable", "Full")) +
  labs(title = "Cell Type Trajectories by Treatment",
       x = "Resistance Stage", 
       y = "Proportion (%)",
       color = "Seurat Cluster") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

print(p_trajectory_treatment)

# ============================================================================
# HEATMAP OF CELL TYPE SCORES ACROSS PSEUDOTIME
# ============================================================================

message("\n=== Cell Type Score Heatmap Across Resistance Stages ===")

# Get cell type score columns
score_cols <- grep("_score$", colnames(GSE164897@meta.data), value = TRUE)
score_cols <- score_cols[!grepl("celltype_score|prediction", score_cols)]

if (length(score_cols) > 0) {
  # Calculate average scores per pseudotime stage
  score_data <- GSE164897@meta.data %>%
    select(pseudotime_stage, all_of(score_cols)) %>%
    group_by(pseudotime_stage) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    column_to_rownames("pseudotime_stage")
  
  # Transpose for heatmap
  score_matrix <- t(score_data)
  rownames(score_matrix) <- gsub("_score", "", rownames(score_matrix))
  
  # Reorder columns by stage progression
  stage_order <- c("1_Naive", "2_Early_Response", "3_Adaptation", "4_Stress", 
                   "5_Transition", "6_Switching", "7_Pre_Resistant", "8_Resistant",
                   "9_Stable_Resistant", "10_Fully_Resistant")
  score_matrix <- score_matrix[, stage_order]
  
  # Simplify column names for display
  colnames(score_matrix) <- c("Naive", "Early\nResp", "Adapt", "Stress",
                              "Trans", "Switch", "Pre-Res", "Resist",
                              "Stable\nRes", "Fully\nRes")
  
  # Create heatmap
  library(pheatmap)
  pheatmap(score_matrix,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           main = "Cell Type Scores Across Resistance Progression",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           fontsize = 11,
           fontsize_row = 10,
           fontsize_col = 8,
           angle_col = "45",
           border_color = "grey60")
  
  message("Cell type score heatmap created")
}

# ============================================================================
# STATISTICAL ANALYSIS
# ============================================================================

message("\n=== Statistical Analysis of Cell Type Changes ===")

# Test for differences in cell type proportions across treatments
for (ct in unique(GSE164897$seurat_clusters)) {
  # Create binary variable for this cell type
  GSE164897@meta.data[[paste0(ct, "_binary")]] <- ifelse(GSE164897$seurat_clusters == ct, 1, 0)
  
  # Test correlation with pseudotime by treatment
  cor_results <- GSE164897@meta.data %>%
    group_by(treatment) %>%
    summarise(
      correlation = cor(pseudotime, get(paste0(ct, "_binary")), use = "complete.obs"),
      .groups = "drop"
    )
  
  message(paste("\n", ct, "correlation with pseudotime:"))
  print(cor_results)
}

# ============================================================================
# ORIGINAL MILESTONE ANALYSIS
# ============================================================================

message("\n=== Milestone Analysis ===")


mile_group <- data.frame(group_onto_nearest_milestones(model)) %>%
  set_colnames("milestone") %>%
  rownames_to_column("cell")


GSE164897$milestone <- mile_group$milestone



cd4_pseudo <- list()
for (i in 1:max(mile_group$milestone)) {
  cd4_pseudo[[i]] <- seurat_extract(GSE164897, meta1 = "milestone", value_meta1 = i)
}


# Load Hallmark Pathways
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- format_pathways(hallmark_df)

# Run SCPA Analysis across Pseudotime Milestones
print("Running SCPA across milestones...")
pathActivation <- compare_pathways(samples = cd4_pseudo, 
                                   pathways = pathways, 
                                   parallel = TRUE, 
                                   cores = 4) # Reduced cores for safety


pathActivation <- as.data.frame(pathActivation)

pathActivation <- pathActivation %>%
  data.frame() %>%
  select(Pathway, qval) %>%
  column_to_rownames("Pathway")

# Handle Inf values for color scale (replace Inf with max finite value)
qvals <- as.matrix(pathActivation)
max_val <- max(qvals[is.finite(qvals)])
qvals[is.infinite(qvals)] <- max_val

col_hm <- colorRamp2(colors = c("white", "red"), breaks = c(0, max_val))

Heatmap(t(qvals),
        name = "Qvalue",
        col = col_hm,
        border = TRUE,
        rect_gp = gpar(col = "white", lwd = 0.1),
        show_column_dend = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 8))



plot_rank(scpa_out = pathActivation, 
          pathway = "HALLMARK_FATTY_ACID_METABOLISM",
          base_point_size = 2.5,
          highlight_point_size = 3)


# === Visualize Pathway Activity Across Clusters ===

print("Calculating per-cell pathway scores...")

# Filter pathways to only include genes present in the Seurat object
seurat_features <- rownames(GSE164897)

# Check if pathways is a data frame (SCPA format) or list
if (is.data.frame(pathways)) {
  pathways_list <- split(pathways$Genes, pathways$Pathway)
} else {
  pathways_list <- pathways
}

# Filter to genes present in dataset
pathways_filtered <- lapply(pathways_list, function(genes) {
  g_vec <- if(is.data.frame(genes) || is.list(genes)) unlist(genes) else genes
  base::intersect(as.character(g_vec), as.character(seurat_features))
})

# Remove empty pathways
pathways_filtered <- pathways_filtered[sapply(pathways_filtered, length) > 0]
print(paste("Retained", length(pathways_filtered), "pathways with overlapping genes"))

# Calculate pathway scores
GSE164897 <- AddModuleScore(GSE164897, features = pathways_filtered, name = "PathwayScore")

# Rename columns to pathway names
meta <- GSE164897@meta.data
renamed_count <- 0
for(i in seq_along(pathways_filtered)) {
  old_name <- paste0("PathwayScore", i)
  new_name <- names(pathways_filtered)[i]
  
  # Check both old name exists and new name is valid
  if(old_name %in% colnames(meta) && !is.null(new_name) && nchar(new_name) > 0) {
    colnames(meta)[colnames(meta) == old_name] <- new_name
    renamed_count <- renamed_count + 1
  }
}
GSE164897@meta.data <- meta
print(paste("Successfully renamed", renamed_count, "pathway columns"))

# Find top 5 most variable pathways across clusters
print("Identifying top variable pathways across clusters...")
pathway_names <- names(pathways_filtered)

# Verify pathway names exist in metadata
pathway_names_in_meta <- pathway_names[pathway_names %in% colnames(GSE164897@meta.data)]
print(paste("Found", length(pathway_names_in_meta), "pathway columns in metadata"))

if(length(pathway_names_in_meta) == 0) {
  print("ERROR: No pathway columns found in metadata. Using PathwayScore columns instead.")
  pathway_names_in_meta <- grep("^PathwayScore", colnames(GSE164897@meta.data), value = TRUE)
}

cluster_means <- GSE164897@meta.data %>%
  select(seurat_clusters, all_of(pathway_names_in_meta)) %>%
  group_by(seurat_clusters) %>%
  summarise(across(everything(), mean))

pathway_vars <- cluster_means %>%
  select(-seurat_clusters) %>%
  summarise(across(everything(), var)) %>%
  t() %>%
  as.data.frame() %>%
  setNames("variance") %>%
  arrange(desc(variance))

top_pathways <- rownames(pathway_vars)
# Remove any NA values
top_pathways <- top_pathways[!is.na(top_pathways)]
print(paste("Plotting", length(top_pathways), "pathways in groups of 5..."))

# Create FeaturePlots for all pathways, 5 per combined plot
if(length(top_pathways) > 0) {
  print("Creating FeaturePlots...")
  
  # Split into groups of 5
  pathway_groups <- split(top_pathways, ceiling(seq_along(top_pathways) / 5))
  
  # Create a plot for each group
  for(group_idx in seq_along(pathway_groups)) {
    group_pathways <- pathway_groups[[group_idx]]
    
    # Create a mapping from PathwayScore names to actual pathway names
    pathway_names_for_plot <- character(length(group_pathways))
    
    for(i in seq_along(group_pathways)) {
      score_col <- group_pathways[i]
      
      if(grepl("^PathwayScore", score_col)) {
        # Extract the index number
        idx <- as.numeric(gsub("PathwayScore", "", score_col))
        # Get the actual pathway name from the pathways list
        if(!is.na(idx) && idx <= length(pathways)) {
          # pathways is a list of tibbles, extract the Pathway name from the first row
          pathway_name <- pathways[[idx]]$Pathway[1]
          # Check if pathway name is valid (not NULL, not empty)
          if(!is.null(pathway_name) && length(pathway_name) > 0 && !is.na(pathway_name) && nchar(pathway_name) > 0) {
            pathway_names_for_plot[i] <- pathway_name
          } else {
            pathway_names_for_plot[i] <- score_col  # Fallback to score name
          }
        } else {
          pathway_names_for_plot[i] <- score_col
        }
      } else {
        pathway_names_for_plot[i] <- score_col
      }
    }
    
    plot_list <- lapply(1:length(group_pathways), function(i) {
      score_col <- group_pathways[i]
      pathway_name <- pathway_names_for_plot[i]
      
      FeaturePlot(GSE164897, features = score_col, reduction = "umap") + 
        scale_color_viridis_c(option = "C") + 
        ggtitle(gsub("HALLMARK_", "", pathway_name)) +
        theme(aspect.ratio = 1, 
              plot.title = element_text(size = 9))
    })
    
    # Combine plots for this group
    library(patchwork)
    combined_plot <- wrap_plots(plot_list, ncol = 3)
    print(combined_plot)
    print(paste("Plotted group", group_idx, "of", length(pathway_groups)))
  }
} else {
  print("No valid pathways found for plotting")
}
