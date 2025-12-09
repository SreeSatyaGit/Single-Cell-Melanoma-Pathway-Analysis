# =============================================================================
# Advanced Pathway Visualizations
# =============================================================================
# Load reproducibility configuration
source("/projects/vanaja_lab/satya/SCPA/Reproduce.R")

# Load required libraries
library(SCPA)
library(Seurat)
library(msigdbr)
library(magrittr)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(ggrepel)

# Set seed for parallel processing
set_parallel_seeds(n_cores = 20)


# Load Hallmark pathways only
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")

# Format pathways for SCPA
pathways <- format_pathways(hallmark_df)


split_treatment <- SplitObject(GSE164897, split.by = "treatment")

# Initialize lists
Pvem_tram <- list(); Pvem_cob <- list(); Pvem <- list()

# =============================================================================
# STEP 1: RUN SCPA ANALYSIS ACROSS CELL STATES (CLUSTERS)
# =============================================================================
if (!"celltype" %in% colnames(GSE164897@meta.data)) {
  # Fallback if annotation script wasn't run
  GSE164897$celltype <- GSE164897$seurat_clusters
}

cell_types <- unique(GSE164897$celltype)
print(paste("Processing", length(cell_types), "cell types (clusters)..."))

for (i in cell_types) {
  print(paste("Analyzing cluster:", i))
  
  cluster_obj <- subset(GSE164897, celltype == i)
  
  if (ncol(cluster_obj) < 5) {
    print(paste("Skipping cluster", i, "- too few cells (total < 5)"))
    next
  }
  
  tryCatch({
    untreated_sub <- seurat_extract(cluster_obj, meta1 = "treatment", value_meta1 = "untreated")
    vem_sub <- seurat_extract(cluster_obj, meta1 = "treatment", value_meta1 = "Vemurafenib")
    vem_cob_sub <- seurat_extract(cluster_obj, meta1 = "treatment", value_meta1 = "vem_cob")
    vem_tram_sub <- seurat_extract(cluster_obj, meta1 = "treatment", value_meta1 = "vem_tram")
    
    u_n <- if(is.null(untreated_sub)) 0 else ncol(untreated_sub)
    vt_n <- if(is.null(vem_tram_sub)) 0 else ncol(vem_tram_sub)
    vc_n <- if(is.null(vem_cob_sub)) 0 else ncol(vem_cob_sub)
    v_n <- if(is.null(vem_sub)) 0 else ncol(vem_sub)
    
    if (u_n < 3) {
      print(paste("  Skipping", i, "- insufficient Untreated cells for comparison"))
    } else {
      idx <- as.character(i)
      if (vt_n >= 3) try({ Pvem_tram[[idx]] <- compare_pathways(list(vem_tram_sub, untreated_sub), pathways, parallel = TRUE, cores = 20) }, silent = TRUE)
      if (vc_n >= 3) try({ Pvem_cob[[idx]] <- compare_pathways(list(vem_cob_sub, untreated_sub), pathways, parallel = TRUE, cores = 20) }, silent = TRUE)
      if (v_n >= 3) try({ Pvem[[idx]] <- compare_pathways(list(vem_sub, untreated_sub), pathways, parallel = TRUE, cores = 20) }, silent = TRUE)
    }
  }, error = function(e) {
    print(paste("Error in cluster", i, ":", e$message))
  })
}

# =============================================================================
# FIGURE 1: MAPK-PI3K CROSSTALK ANALYSIS
# =============================================================================

cat("\n=== Generating Figure 1 Key Plots: MAPK vs PI3K Crosstalk ===\n")

# Function to extract specific pathway scores
extract_pathways <- function(result_list, treatment_label, pathways_of_interest) {
  do.call(rbind, lapply(names(result_list), function(ct) {
    res <- result_list[[ct]]
    if (is.null(res)) return(NULL)
    
    # Filter for pathways of interest
    matches <- res[res$Pathway %in% pathways_of_interest, ]
    
    if (nrow(matches) > 0) {
      data.frame(
        CellType = ct,
        Treatment = treatment_label,
        Pathway = matches$Pathway,
        FC = matches$FC,
        Qval = matches$qval,
        Score = sign(matches$FC) * -log10(matches$qval + 1e-50)
      )
    } else {
      NULL
    }
  }))
}

# Check if data was generated
if (length(Pvem) == 0 && length(Pvem_cob) == 0 && length(Pvem_tram) == 0) {
  warning("WARNING: No pathway analysis results generated. Figures will differ or fail.")
}

# Define the crosstalk pathways
target_pathways <- c("HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_KRAS_SIGNALING_UP")

# Extract data
crosstalk_data <- rbind(
  extract_pathways(Pvem, "Vem", target_pathways),
  extract_pathways(Pvem_cob, "Vem+Cob", target_pathways),
  extract_pathways(Pvem_tram, "Vem+Tram", target_pathways)
)

# Clean pathway names for plotting
crosstalk_data$Pathway <- gsub("HALLMARK_", "", crosstalk_data$Pathway)
crosstalk_data$Pathway <- gsub("_SIGNALING_UP", " (MAPK)", crosstalk_data$Pathway)
crosstalk_data$Pathway <- gsub("_SIGNALING", "", crosstalk_data$Pathway)
crosstalk_data$Pathway <- gsub("_", " ", crosstalk_data$Pathway)

# Plot 1: Side-by-Side Comparison (Bar Plot)
p_mech_bar <- ggplot(crosstalk_data, aes(x = CellType, y = FC, fill = Pathway)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(~Treatment) +
  scale_fill_manual(values = c("KRAS (MAPK)" = "#3498DB", "PI3K AKT MTOR" = "#E74C3C")) +
  labs(title = "Differential Response: MAPK Suppression vs PI3K Persistence",
       subtitle = "Standard dual inhibition (Vem+Cob/Tram) fails to suppress PI3K as effectively as MAPK",
       y = "Pathway Enrichment (Fold Change vs Untreated)",
       x = NULL,
       fill = "Signaling Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        panel.border = element_rect(fill = NA, color = "grey60"))


# Plot 2: Crosstalk Scatter (The "Escape" Plot)
# Reshape data to wide format: One row per CellType+Treatment, columns for MAPK and PI3K scores
library(tidyr)
crosstalk_wide <- crosstalk_data %>%
  dplyr::select(CellType, Treatment, Pathway, FC) %>%
  pivot_wider(names_from = Pathway, values_from = FC)

names(crosstalk_wide) <- make.names(names(crosstalk_wide)) # Handle spaces
# Assuming names became something like PI3K.AKT.MTOR and KRAS..MAPK.

# Identify column names dynamically to be safe
pi3k_col <- grep("PI3K", names(crosstalk_wide), value = TRUE)
mapk_col <- grep("KRAS", names(crosstalk_wide), value = TRUE)

p_crosstalk_scatter <- ggplot(crosstalk_wide, aes_string(x = mapk_col, y = pi3k_col, color = "Treatment")) +
  geom_point(size = 4, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Vem" = "#E74C3C", "Vem+Cob" = "#3498DB", "Vem+Tram" = "#2ECC71")) +
  labs(title = "Mechanism of Resistance: The PI3K Escape",
       subtitle = "Points above the diagonal indicate PI3K dominance over MAPK suppression",
       x = "MAPK Inhibition Level (Fold Change)",
       y = "PI3K Activity Level (Fold Change)") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "bottom") +
  geom_text_repel(aes(label = CellType), size = 3, max.overlaps = 10)


print(p_mech_bar)
print(p_crosstalk_scatter)

# Assemble Figure 1 with these narrative elements
cat("Assembling Final Manuscript Figure 1...\n")
library(patchwork)
library(ggplotify)

# Convert ComplexHeatmap to grob
p_heatmap_grob <- grid.grabExpr(draw(ht))

# Layout: 
# A | B
# C C
# D D
# E E
layout_design <- "
AAABBB
CCCCCC
CCCCCC
DDDDDD
DDDDDD
EEEEEE
EEEEEE
"

# Note: Using p_mech_bar (Panel D) and p_crosstalk_scatter (Panel E)
fig1 <- (p_umap_tx | p_umap_ct) / 
  as.ggplot(p_heatmap_grob) / 
  p_mech_bar /
  p_crosstalk_scatter + 
  plot_layout(heights = c(1, 2, 1.2, 1.5)) +
  plot_annotation(tag_levels = 'A', 
                  title = "Figure 1: MAPK-PI3K Crosstalk Drives Resistance to BRAF+MEK Inhibition",
                  subtitle = "Dual inhibition of BRAF and PI3K works better than dual BRAF and MEK inhibition")

ggsave("Figure1_MAPK_PI3K_Crosstalk_Complete.pdf", fig1, width = 14, height = 24)
cat("✓ Use 'Figure1_MAPK_PI3K_Crosstalk_Complete.pdf' for your manuscript.\n")

# =============================================================================
# FIGURE 1 ASSEMBLY PREP
# =============================================================================

# 1. Generate Context Map (UMAP)
cat("Generating UMAP context plots...\n")
p_umap_tx <- DimPlot(GSE164897, group.by = "treatment", cols = c("untreated" = "grey70", "Vemurafenib" = "#E74C3C", "vem_cob" = "#3498DB", "vem_tram" = "#2ECC71"), pt.size = 0.5) +
  ggtitle("Treatment Groups") + theme_void() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_umap_ct <- DimPlot(GSE164897, group.by = "celltype", label = TRUE, label.size = 3, repel = TRUE) +
  ggtitle("Cell Types") + theme_void() + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")



# (Analysis loop moved to Step 1 at top of file)


# Helper function to get Signed Q-values (sign(FC) * -log10(qval))
get_signed_scores <- function(scpa_out, name) {
  df <- list()
  for (i in names(scpa_out)) {
    res <- scpa_out[[i]]
    # Add small epsilon to qval to avoid Inf
    res$score <- sign(res$FC) * -log10(res$qval + 1e-50)
    
    df[[i]] <- res %>% dplyr::select(Pathway, score)
  }
  
  col_names <- names(df)
  for (i in 1:length(df)) {
    df[[i]] <- setNames(df[[i]], c("pathway", paste(name, col_names[[i]], sep = "_")))
  }
  return(df)
}


scpa_results <- Reduce(full_join, c(get_signed_scores(Pvem_tram, "Vemtram"),
                                    get_signed_scores(Pvem_cob, "Vemcob"),
                                    get_signed_scores(Pvem, "Vem")))


all_data <- scpa_results %>%
  column_to_rownames("pathway")


all_data <- all_data[grep("HALLMARK_", rownames(all_data)), ]


all_data[is.na(all_data)] <- 0
mat <- as.matrix(all_data)


rownames(mat) <- gsub("HALLMARK_", "", rownames(mat))
rownames(mat) <- gsub("REACTOME_", "", rownames(mat))
rownames(mat) <- gsub("KEGG_", "", rownames(mat))
rownames(mat) <- gsub("BIOCARTA_", "", rownames(mat))
rownames(mat) <- gsub("WP_", "", rownames(mat))
rownames(mat) <- gsub("_", " ", rownames(mat)) 


row_var <- apply(mat, 1, var)
top_pathways <- names(sort(row_var, decreasing = TRUE))[1:20]
mat_subset <- mat[top_pathways, , drop = FALSE]


col_names <- colnames(mat)

treatment <- sapply(strsplit(col_names, "_"), `[`, 1)
cell_type <- sapply(strsplit(col_names, "_"), `[`, 2)

n_types <- length(unique(cell_types))
type_colors <- setNames(rainbow(n_types), unique(cell_types))
treatment_colors <- c("Vem" = "#E74C3C", "Vemcob" = "#3498DB", "Vemtram" = "#2ECC71")

col_an <- HeatmapAnnotation(
  Treatment = treatment,
  CellType = cell_type,
  col = list(Treatment = treatment_colors, CellType = type_colors),
  gp = gpar(col = "white", lwd = 0.05),
  simple_anno_size = unit(3, "mm")
)

print(paste("Plotting signed heatmap for top", nrow(mat_subset), "pathways"))


ht <- Heatmap(mat_subset,
              name = "Signed -log10(qval)",
              top_annotation = col_an,
              col = colorRamp2(c(-10, 0, 10), c("blue", "white", "red")), 
              show_row_names = TRUE,
              show_column_names = TRUE,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              row_names_gp = gpar(fontsize = 9),
              column_names_gp = gpar(fontsize = 8),
              rect_gp = gpar(col = "white", lwd = 0.5))
draw(ht)

p_var <- apply(mat, 1, var) %>%
  as.data.frame() %>% 
  setNames("variance") %>%
  arrange(desc(variance)) %>% 
  rownames_to_column("pathway") %>%
  ggplot(aes(reorder(pathway, variance), variance)) +
  geom_point(shape = 21, cex = 3, fill = "royalblue2", color = 'black', stroke = 0.2) +
  scale_x_discrete(expand = c(0.04, 0.04)) +
  labs(title = "Pathway Variance Across Conditions",
       x = "Pathway", 
       y = "Variance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"))
print(p_var)




