# Melanoma Cell State Annotation
# This script identifies melanoma cell states relevant for drug resistance analysis

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# Install pheatmap if needed (doesn't have the GlobalOptions corruption issue)
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)

# Check if GSE164897 object exists
if (!exists("GSE164897")) {
  stop("ERROR: GSE164897 Seurat object not found. Please run GSE164897.R first.")
}

message("\n=== Melanoma Cell State Annotation ===")
message("Analyzing melanoma-specific states for drug resistance profiling\n")

# Define melanoma-specific marker gene sets
melanoma_states <- list(
  # Differentiated melanocytic state
  Differentiated = c("MLANA", "TYR", "DCT", "PMEL", "TYRP1", "GPR143", "OCA2"),
  
  # Undifferentiated/Neural crest-like state
  Undifferentiated = c("SOX10", "MITF", "PAX3", "SOX9", "FOXD3"),
  
  # Invasive/Mesenchymal state (associated with resistance)
  Invasive = c("AXL", "NGFR", "VIM", "TWIST1", "ZEB1", "SNAI2", "FN1", "CDH2"),
  
  # Proliferative state
  Proliferative = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "CCNA2"),
  
  # Resistant/Adaptive state
  Resistant = c("EGFR", "PDGFRB", "AXL", "NGFR", "WNT5A", "JUN", "FOS"),
  
  # Hypoxic/Stress response
  Hypoxic = c("HIF1A", "VEGFA", "SLC2A1", "LDHA", "PDK1", "BNIP3"),
  
  # Immune evasion
  ImmuneEvasion = c("CD274", "PDCD1LG2", "IDO1", "HAVCR2", "LAG3")
)

message("Defined", length(melanoma_states), "melanoma cell states:")
for (state in names(melanoma_states)) {
  message(paste("  -", state, ":", length(melanoma_states[[state]]), "markers"))
}

# Check which markers are present in the dataset
message("\n=== Checking Marker Gene Availability ===")
all_markers <- unique(unlist(melanoma_states))
available_markers <- intersect(all_markers, rownames(GSE164897))
missing_markers <- setdiff(all_markers, rownames(GSE164897))

message(paste("Total markers defined:", length(all_markers)))
message(paste("Available in dataset:", length(available_markers), 
              paste0("(", round(100*length(available_markers)/length(all_markers), 1), "%)")))

if (length(missing_markers) > 0) {
  message(paste("\nMissing markers:", paste(head(missing_markers, 10), collapse = ", ")))
  if (length(missing_markers) > 10) {
    message(paste("  ... and", length(missing_markers) - 10, "more"))
  }
}

# Filter marker lists to only include available genes
melanoma_states_filtered <- lapply(melanoma_states, function(markers) {
  intersect(markers, rownames(GSE164897))
})

# Remove empty gene sets
melanoma_states_filtered <- melanoma_states_filtered[sapply(melanoma_states_filtered, length) > 0]

message("\n=== Calculating Module Scores ===")

# First, remove any existing score columns to prevent duplicates
existing_score_cols <- grep("_score1?$", colnames(GSE164897@meta.data), value = TRUE)
if (length(existing_score_cols) > 0) {
  message("Removing existing score columns to prevent duplicates...")
  cols_to_remove <- c()
  for (state in names(melanoma_states_filtered)) {
    # Check for both _score and _score1 versions
    pattern <- paste0("^", state, "_score1?$")
    matching_cols <- grep(pattern, colnames(GSE164897@meta.data), value = TRUE)
    cols_to_remove <- c(cols_to_remove, matching_cols)
  }
  if (length(cols_to_remove) > 0) {
    GSE164897@meta.data <- GSE164897@meta.data[, !colnames(GSE164897@meta.data) %in% cols_to_remove, drop = FALSE]
    message(paste("Removed", length(unique(cols_to_remove)), "existing score columns"))
  }
}

# Calculate module scores for each state
for (state in names(melanoma_states_filtered)) {
  if (length(melanoma_states_filtered[[state]]) > 0) {
    message(paste("Scoring", state, "state..."))
    GSE164897 <- AddModuleScore(
      GSE164897,
      features = list(melanoma_states_filtered[[state]]),
      name = paste0(state, "_score"),
      seed = 42
    )
  }
}

message("Module scoring complete!")

# Get all score column names (AddModuleScore adds "1" suffix)
score_cols <- grep("_score1$", colnames(GSE164897@meta.data), value = TRUE)

# Rename score columns for clarity (remove the "1" suffix)
for (i in seq_along(names(melanoma_states_filtered))) {
  old_name <- paste0(names(melanoma_states_filtered)[i], "_score1")
  new_name <- paste0(names(melanoma_states_filtered)[i], "_score")
  if (old_name %in% colnames(GSE164897@meta.data)) {
    colnames(GSE164897@meta.data)[colnames(GSE164897@meta.data) == old_name] <- new_name
  }
}

# Update score_cols with new names
score_cols <- paste0(names(melanoma_states_filtered), "_score")

message("\n=== Assigning Dominant Cell States ===")
# Assign dominant state to each cell
state_scores <- GSE164897@meta.data[, score_cols, drop = FALSE]
GSE164897$dominant_state <- colnames(state_scores)[apply(state_scores, 1, which.max)]
GSE164897$dominant_state <- gsub("_score", "", GSE164897$dominant_state)
GSE164897$dominant_state_score <- apply(state_scores, 1, max)

# Cell state distribution
message("\n=== Cell State Distribution ===")
state_table <- table(GSE164897$dominant_state)
print(state_table)

# Distribution by treatment
message("\n=== Cell States by Treatment ===")
state_by_treatment <- table(GSE164897$dominant_state, GSE164897$treatment)
print(state_by_treatment)

# Calculate proportions
state_proportions <- prop.table(state_by_treatment, margin = 2) * 100
message("\n=== Cell State Proportions (%) by Treatment ===")
print(round(state_proportions, 2))

message("\n=== Creating Visualizations ===")

# 1. UMAP colored by dominant state
p1 <- DimPlot(GSE164897, group.by = "dominant_state", label = TRUE, repel = TRUE) +
  ggtitle("Melanoma Cell States") +
  theme_minimal() +
  theme(legend.position = "right")

print(p1)

# 2. UMAP split by treatment
p2 <- DimPlot(GSE164897, group.by = "dominant_state", split.by = "treatment", 
              label = TRUE, repel = TRUE, ncol = 2) +
  ggtitle("Melanoma Cell States by Treatment") +
  theme_minimal()

print(p2)

# 3. Feature plots for key markers
message("Creating feature plots for key markers...")
key_markers <- c("MLANA", "AXL", "NGFR", "MKI67", "SOX10", "VIM")
available_key_markers <- intersect(key_markers, rownames(GSE164897))

if (length(available_key_markers) > 0) {
  p3 <- FeaturePlot(GSE164897, features = available_key_markers, ncol = 3) &
    theme_minimal() &
    theme(legend.position = "right")
  print(p3)
}

# 4. Violin plots of state scores by treatment
message("Creating violin plots...")
for (state in names(melanoma_states_filtered)) {
  score_col <- paste0(state, "_score")
  if (score_col %in% colnames(GSE164897@meta.data)) {
    p <- VlnPlot(GSE164897, features = score_col, group.by = "treatment", pt.size = 0) +
      ggtitle(paste(state, "Score by Treatment")) +
      theme_minimal()
    print(p)
  }
}

# 5. Heatmap of average state scores by treatment
message("Creating state score heatmap...")
avg_scores <- GSE164897@meta.data %>%
  group_by(treatment) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE)) %>%
  column_to_rownames("treatment")

# Transpose for better visualization
avg_scores_t <- t(avg_scores)
rownames(avg_scores_t) <- gsub("_score", "", rownames(avg_scores_t))

# Create heatmap using pheatmap (no GlobalOptions dependency)
pheatmap(avg_scores_t,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Average Melanoma State Scores by Treatment",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(min(avg_scores_t), max(avg_scores_t), length.out = 101),
         fontsize = 11,
         fontsize_row = 10,
         fontsize_col = 10,
         angle_col = 45,
         border_color = "grey60",
         cellwidth = 40,
         cellheight = 20)

message("Heatmap created successfully with pheatmap")

# 6. Stacked bar plot of state proportions
message("Creating stacked bar plot...")
state_prop_df <- as.data.frame(state_proportions)
colnames(state_prop_df) <- c("State", "Treatment", "Percentage")

p_bar <- ggplot(state_prop_df, aes(x = Treatment, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Melanoma Cell State Composition by Treatment",
       y = "Percentage (%)", x = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

print(p_bar)

# 7. Dot plot of marker expression by cluster
message("Creating marker dot plot...")
if (length(available_markers) > 0) {
  p_dot <- DotPlot(GSE164897, 
                   features = head(available_markers, 20),  # Top 20 markers
                   group.by = "seurat_clusters") +
    RotatedAxis() +
    ggtitle("Melanoma Marker Expression by Cluster") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  print(p_dot)
}

message("\n=== Statistical Analysis ===")
# Test for significant differences in state scores between treatments
message("Testing for treatment effects on cell states...")

for (state in names(melanoma_states_filtered)) {
  score_col <- paste0(state, "_score")
  if (score_col %in% colnames(GSE164897@meta.data)) {
    # Kruskal-Wallis test (non-parametric ANOVA)
    kw_test <- kruskal.test(GSE164897@meta.data[[score_col]] ~ GSE164897$treatment)
    
    if (kw_test$p.value < 0.05) {
      message(paste("\n", state, "state:"))
      message(paste("  Kruskal-Wallis p-value:", format(kw_test$p.value, scientific = TRUE)))
      message("  ** Significant difference between treatments **")
      
      # Calculate mean scores by treatment
      mean_scores <- tapply(GSE164897@meta.data[[score_col]], 
                            GSE164897$treatment, 
                            mean, na.rm = TRUE)
      message("  Mean scores by treatment:")
      for (trt in names(mean_scores)) {
        message(paste("    ", trt, ":", round(mean_scores[trt], 4)))
      }
    }
  }
}

message("\n=== Saving Results ===")


message("\n=== Melanoma State Annotation Complete! ===")
message("Annotations stored in:")
message("  - GSE164897$dominant_state (dominant cell state)")
message("  - GSE164897$dominant_state_score (confidence score)")
message("  - GSE164897$<State>_score (individual state scores)")
message("\nKey findings:")
message(paste("  - Identified", length(unique(GSE164897$dominant_state)), "distinct melanoma states"))
message(paste("  - Analyzed", ncol(GSE164897), "cells across", length(unique(GSE164897$treatment)), "treatments"))
message("  - Results saved to CSV and RDS files")




