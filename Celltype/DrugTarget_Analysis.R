# =============================================================================
# DRUG TARGET IDENTIFICATION ANALYSIS
# =============================================================================
# Identifies druggable genes and pathways for melanoma resistance
# Prioritizes therapeutic targets based on expression and druggability
# Suggests combination therapy strategies
# =============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# Load reproducibility configuration
source("/projects/vanaja_lab/satya/SCPA/Reproduce.R")
set_analysis_seed("drug_target_analysis")

cat("\n=============================================================================\n")
cat("DRUG TARGET IDENTIFICATION ANALYSIS\n")
cat("=============================================================================\n\n")

# Create output directory
dir.create("DrugTarget_Results", showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA AND DEG RESULTS
# =============================================================================

cat("=== 1. LOADING DATA ===\n")

# Load GSE164897
if (!exists("GSE164897")) {
  cat("  Loading GSE164897...\n")
  source("GSE164897.R")
}

# Load DEG results if available
deg_results_file <- "DEG_Results/All_DEG_Results.RData"
if (file.exists(deg_results_file)) {
  cat("  Loading DEG results...\n")
  load(deg_results_file)
} else {
  cat("  Warning: DEG results not found. Run DEG_Analysis.R first.\n")
  cat("  Proceeding with limited analysis...\n")
  treatment_degs <- NULL
  stage_degs <- NULL
}

# =============================================================================
# 2. KNOWN DRUGGABLE GENE FAMILIES
# =============================================================================

cat("\n=== 2. DEFINING DRUGGABLE GENE FAMILIES ===\n")

# Druggable gene families
druggable_genes <- list(
  
  # Receptor Tyrosine Kinases (RTKs)
  RTKs = c("EGFR", "ERBB2", "ERBB3", "ERBB4", "MET", "ALK", "ROS1", 
           "AXL", "PDGFRA", "PDGFRB", "FGFR1", "FGFR2", "FGFR3", "FGFR4",
           "IGF1R", "INSR", "KIT", "FLT3", "RET", "NTRK1", "NTRK2", "NTRK3"),
  
  # Kinases
  Kinases = c("BRAF", "NRAS", "KRAS", "MEK1", "MEK2", "ERK1", "ERK2",
              "AKT1", "AKT2", "AKT3", "MTOR", "PI3K", "CDK4", "CDK6",
              "JAK1", "JAK2", "JAK3", "SRC", "ABL1", "BTK"),
  
  # Immune Checkpoints
  Immune_Checkpoints = c("PD1", "PDL1", "PDL2", "CTLA4", "LAG3", "TIM3",
                         "TIGIT", "VISTA", "B7H3", "ICOS", "OX40", "CD40"),
  
  # Transcription Factors
  Transcription_Factors = c("SOX10", "SOX9", "MITF", "FOXD3", "TFAP2A",
                            "JUN", "FOS", "STAT3", "NFKB1", "NFKB2", "REL"),
  
  # Epigenetic Regulators
  Epigenetic = c("EZH2", "HDAC1", "HDAC2", "HDAC3", "HDAC6", "DNMT1",
                 "DNMT3A", "DNMT3B", "KDM1A", "KDM4A", "BRD4"),
  
  # Metabolic Enzymes
  Metabolic = c("LDHA", "PKM", "HK2", "PFKFB3", "GLS", "GLUD1",
                "FASN", "ACLY", "IDH1", "IDH2"),
  
  # Apoptosis Regulators
  Apoptosis = c("BCL2", "BCL2L1", "BCL2L2", "MCL1", "BAX", "BAK1",
                "BID", "BIM", "PUMA", "NOXA"),
  
  # Cell Cycle Regulators
  Cell_Cycle = c("CDK4", "CDK6", "CCND1", "CCND2", "CCND3", "CCNE1",
                 "RB1", "E2F1", "E2F2", "E2F3"),
  
  # Adhesion/Migration
  Adhesion = c("ITGB1", "ITGB3", "ITGB4", "ITGAV", "CD44", "MCAM",
               "CDH1", "CDH2", "VIM", "FN1"),
  
  # Angiogenesis
  Angiogenesis = c("VEGFA", "VEGFB", "VEGFC", "VEGFR1", "VEGFR2", "VEGFR3",
                   "ANGPT1", "ANGPT2", "TIE2", "HIF1A")
)

# Flatten to get all druggable genes
all_druggable <- unique(unlist(druggable_genes))

cat(paste("  Defined", length(all_druggable), "druggable genes across",
          length(druggable_genes), "families\n"))

# =============================================================================
# 3. IDENTIFY DRUGGABLE TARGETS IN RESISTANCE
# =============================================================================

cat("\n=== 3. IDENTIFYING DRUGGABLE RESISTANCE TARGETS ===\n")

if (!is.null(stage_degs)) {
  
  # Filter for druggable genes
  druggable_resistance_genes <- stage_degs %>%
    filter(gene %in% all_druggable) %>%
    filter(significant) %>%
    arrange(desc(abs(avg_log2FC)))
  
  # Classify by family
  druggable_resistance_genes$family <- sapply(druggable_resistance_genes$gene, function(g) {
    families <- names(druggable_genes)[sapply(druggable_genes, function(x) g %in% x)]
    if (length(families) > 0) families[1] else "Other"
  })
  
  cat(paste("  Found", nrow(druggable_resistance_genes), 
            "druggable genes in resistance signature\n"))
  
  # Top targets
  top_targets <- druggable_resistance_genes %>%
    filter(avg_log2FC > 0.5) %>%
    head(20)
  
  cat("\n  Top 20 Druggable Resistance Targets:\n")
  for (i in 1:min(20, nrow(top_targets))) {
    cat(sprintf("    %2d. %-10s (FC: %5.2f, Family: %s)\n",
                i, top_targets$gene[i], top_targets$avg_log2FC[i], 
                top_targets$family[i]))
  }
  
  # Save results
  write.csv(druggable_resistance_genes,
            file.path("DrugTarget_Results", "Druggable_Resistance_Genes.csv"),
            row.names = FALSE)
  
} else {
  cat("  Stage DEGs not available. Skipping resistance target analysis.\n")
  druggable_resistance_genes <- NULL
}

# =============================================================================
# 4. CELL TYPE-SPECIFIC DRUGGABLE TARGETS
# =============================================================================

cat("\n=== 4. CELL TYPE-SPECIFIC DRUGGABLE TARGETS ===\n")

# Focus on key resistance cell types
resistance_celltypes <- c("Invasive", "Resistant", "Hypoxic", "ImmuneEvasion")
resistance_celltypes <- intersect(resistance_celltypes, unique(GSE164897$celltype))

celltype_druggable_targets <- list()

for (ct in resistance_celltypes) {
  cat(paste("\n  Analyzing:", ct, "\n"))
  
  # Get cells of this type
  ct_cells <- subset(GSE164897, celltype == ct)
  
  # Calculate average expression of druggable genes
  expr_data <- GetAssayData(ct_cells, layer = "data")
  druggable_present <- intersect(all_druggable, rownames(expr_data))
  
  if (length(druggable_present) > 0) {
    avg_expr <- rowMeans(expr_data[druggable_present, ])
    
    # Get top expressed druggable genes
    top_expr <- sort(avg_expr, decreasing = TRUE)[1:20]
    
    # Create data frame
    ct_targets <- data.frame(
      gene = names(top_expr),
      avg_expression = as.numeric(top_expr),
      celltype = ct,
      stringsAsFactors = FALSE
    )
    
    # Add family
    ct_targets$family <- sapply(ct_targets$gene, function(g) {
      families <- names(druggable_genes)[sapply(druggable_genes, function(x) g %in% x)]
      if (length(families) > 0) families[1] else "Other"
    })
    
    celltype_druggable_targets[[ct]] <- ct_targets
    
    cat(paste("    Top targets:", paste(head(ct_targets$gene, 5), collapse = ", "), "\n"))
  }
}

# Combine and save
if (length(celltype_druggable_targets) > 0) {
  all_celltype_targets <- bind_rows(celltype_druggable_targets)
  write.csv(all_celltype_targets,
            file.path("DrugTarget_Results", "CellType_Druggable_Targets.csv"),
            row.names = FALSE)
}

# =============================================================================
# 5. TREATMENT-SPECIFIC TARGETS
# =============================================================================

cat("\n=== 5. TREATMENT-SPECIFIC DRUGGABLE TARGETS ===\n")

if (!is.null(treatment_degs)) {
  
  treatment_targets <- list()
  
  for (trt in names(treatment_degs)) {
    cat(paste("\n  Analyzing:", trt, "\n"))
    
    degs <- treatment_degs[[trt]]
    
    # Filter for druggable genes
    druggable_degs <- degs %>%
      filter(gene %in% all_druggable, significant) %>%
      arrange(desc(abs(avg_log2FC)))
    
    # Add family
    druggable_degs$family <- sapply(druggable_degs$gene, function(g) {
      families <- names(druggable_genes)[sapply(druggable_genes, function(x) g %in% x)]
      if (length(families) > 0) families[1] else "Other"
    })
    
    treatment_targets[[trt]] <- druggable_degs
    
    # Top upregulated targets
    top_up <- druggable_degs %>%
      filter(avg_log2FC > 0.5) %>%
      head(10)
    
    if (nrow(top_up) > 0) {
      cat("    Top upregulated druggable genes:\n")
      for (i in 1:nrow(top_up)) {
        cat(sprintf("      %s (FC: %.2f, %s)\n", 
                    top_up$gene[i], top_up$avg_log2FC[i], top_up$family[i]))
      }
    }
  }
  
  # Save
  all_treatment_targets <- bind_rows(treatment_targets)
  write.csv(all_treatment_targets,
            file.path("DrugTarget_Results", "Treatment_Druggable_Targets.csv"),
            row.names = FALSE)
}

# =============================================================================
# 6. PRIORITIZE TARGETS BY MULTIPLE CRITERIA
# =============================================================================

cat("\n=== 6. PRIORITIZING DRUG TARGETS ===\n")

# Create comprehensive target ranking
if (!is.null(druggable_resistance_genes) && !is.null(treatment_degs)) {
  
  # Get all druggable genes present in data
  all_genes_in_data <- rownames(GSE164897)
  druggable_in_data <- intersect(all_druggable, all_genes_in_data)
  
  target_priority <- data.frame(
    gene = druggable_in_data,
    stringsAsFactors = FALSE
  )
  
  # Add family
  target_priority$family <- sapply(target_priority$gene, function(g) {
    families <- names(druggable_genes)[sapply(druggable_genes, function(x) g %in% x)]
    if (length(families) > 0) families[1] else "Other"
  })
  
  # Criterion 1: Resistance enrichment
  target_priority$resistance_FC <- sapply(target_priority$gene, function(g) {
    if (g %in% druggable_resistance_genes$gene) {
      druggable_resistance_genes$avg_log2FC[druggable_resistance_genes$gene == g]
    } else {
      0
    }
  })
  
  # Criterion 2: Expression level in resistant cells
  resistant_cells <- subset(GSE164897, 
                            celltype %in% c("Resistant", "Invasive", "ImmuneEvasion"))
  expr_data <- GetAssayData(resistant_cells, layer = "data")
  
  target_priority$resistant_expression <- sapply(target_priority$gene, function(g) {
    if (g %in% rownames(expr_data)) {
      mean(expr_data[g, ])
    } else {
      0
    }
  })
  
  # Criterion 3: Specificity to resistant cells
  sensitive_cells <- subset(GSE164897, celltype == "Differentiated")
  expr_sensitive <- GetAssayData(sensitive_cells, layer = "data")
  
  target_priority$specificity_score <- sapply(target_priority$gene, function(g) {
    if (g %in% rownames(expr_data) && g %in% rownames(expr_sensitive)) {
      resistant_expr <- mean(expr_data[g, ])
      sensitive_expr <- mean(expr_sensitive[g, ])
      log2((resistant_expr + 0.1) / (sensitive_expr + 0.1))
    } else {
      0
    }
  })
  
  # Calculate composite priority score
  target_priority <- target_priority %>%
    mutate(
      # Normalize each criterion to 0-1 scale
      resistance_score = (resistance_FC - min(resistance_FC)) / 
                        (max(resistance_FC) - min(resistance_FC) + 0.001),
      expression_score = (resistant_expression - min(resistant_expression)) / 
                        (max(resistant_expression) - min(resistant_expression) + 0.001),
      specificity_norm = (specificity_score - min(specificity_score)) / 
                        (max(specificity_score) - min(specificity_score) + 0.001),
      # Weighted composite score
      priority_score = (resistance_score * 0.4) + 
                      (expression_score * 0.3) + 
                      (specificity_norm * 0.3)
    ) %>%
    arrange(desc(priority_score))
  
  # Top 30 prioritized targets
  top_priorities <- head(target_priority, 30)
  
  cat("\n  Top 30 Prioritized Drug Targets:\n")
  cat("  Rank | Gene       | Family              | Priority | Resistance FC | Expression\n")
  cat("  -----|------------|---------------------|----------|---------------|------------\n")
  for (i in 1:nrow(top_priorities)) {
    cat(sprintf("  %4d | %-10s | %-19s | %8.3f | %13.2f | %10.2f\n",
                i, 
                top_priorities$gene[i],
                substr(top_priorities$family[i], 1, 19),
                top_priorities$priority_score[i],
                top_priorities$resistance_FC[i],
                top_priorities$resistant_expression[i]))
  }
  
  # Save prioritized targets
  write.csv(target_priority,
            file.path("DrugTarget_Results", "Prioritized_Drug_Targets.csv"),
            row.names = FALSE)
}

# =============================================================================
# 7. COMBINATION THERAPY SUGGESTIONS
# =============================================================================

cat("\n=== 7. COMBINATION THERAPY SUGGESTIONS ===\n")

# Suggest combinations based on pathway complementarity
combination_strategies <- data.frame(
  Strategy = character(),
  Target1 = character(),
  Target2 = character(),
  Rationale = character(),
  stringsAsFactors = FALSE
)

# Strategy 1: BRAF + MEK (standard)
combination_strategies <- rbind(combination_strategies, data.frame(
  Strategy = "Standard_of_Care",
  Target1 = "BRAF",
  Target2 = "MEK1/2",
  Rationale = "Dual MAPK pathway inhibition"
))

# Strategy 2: BRAF/MEK + RTK
if (exists("top_priorities")) {
  top_rtks <- top_priorities %>%
    filter(family == "RTKs") %>%
    head(3)
  
  for (i in 1:nrow(top_rtks)) {
    combination_strategies <- rbind(combination_strategies, data.frame(
      Strategy = "MAPK_plus_RTK",
      Target1 = "BRAF/MEK",
      Target2 = top_rtks$gene[i],
      Rationale = paste("Block MAPK reactivation via", top_rtks$gene[i])
    ))
  }
}

# Strategy 3: Targeted + Immunotherapy
combination_strategies <- rbind(combination_strategies, data.frame(
  Strategy = "Targeted_plus_Immuno",
  Target1 = "BRAF/MEK",
  Target2 = "PD1/PDL1",
  Rationale = "Combine targeted therapy with immune checkpoint blockade"
))

# Strategy 4: Epigenetic + Targeted
combination_strategies <- rbind(combination_strategies, data.frame(
  Strategy = "Epigenetic_plus_Targeted",
  Target1 = "BRAF/MEK",
  Target2 = "HDAC/BRD4",
  Rationale = "Reverse epigenetic resistance mechanisms"
))

# Strategy 5: Metabolic + Targeted
combination_strategies <- rbind(combination_strategies, data.frame(
  Strategy = "Metabolic_plus_Targeted",
  Target1 = "BRAF/MEK",
  Target2 = "LDHA/PKM",
  Rationale = "Target metabolic adaptation in resistant cells"
))

cat("\n  Suggested Combination Therapy Strategies:\n\n")
for (i in 1:nrow(combination_strategies)) {
  cat(sprintf("  %d. %s\n", i, combination_strategies$Strategy[i]))
  cat(sprintf("     %s + %s\n", combination_strategies$Target1[i], 
              combination_strategies$Target2[i]))
  cat(sprintf("     Rationale: %s\n\n", combination_strategies$Rationale[i]))
}

write.csv(combination_strategies,
          file.path("DrugTarget_Results", "Combination_Therapy_Strategies.csv"),
          row.names = FALSE)

# =============================================================================
# 8. VISUALIZATIONS
# =============================================================================

cat("\n=== 8. CREATING VISUALIZATIONS ===\n")

# --- 8.1 Target Priority Plot ---
if (exists("top_priorities")) {
  
  p_priority <- ggplot(head(top_priorities, 20), 
                       aes(x = reorder(gene, priority_score), 
                           y = priority_score,
                           fill = family)) +
    geom_col() +
    coord_flip() +
    labs(title = "Top 20 Prioritized Drug Targets",
         x = "Gene",
         y = "Priority Score",
         fill = "Gene Family") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          panel.border = element_rect(fill = NA, color = "black"))
  
  print(p_priority)
  
  ggsave(file.path("DrugTarget_Results", "Top_Prioritized_Targets.pdf"),
         p_priority, width = 10, height = 8)
}

# --- 8.2 Druggable Gene Expression Heatmap ---
if (exists("top_priorities")) {
  
  top_target_genes <- head(top_priorities$gene, 30)
  
  # Get expression across cell types
  expr_matrix <- GetAssayData(GSE164897, layer = "data")[top_target_genes, ]
  
  cell_types_ordered <- sort(unique(GSE164897$celltype))
  avg_expr <- sapply(cell_types_ordered, function(ct) {
    cells <- colnames(GSE164897)[GSE164897$celltype == ct]
    rowMeans(expr_matrix[, cells, drop = FALSE])
  })
  
  # Scale for visualization
  avg_expr_scaled <- t(scale(t(avg_expr)))
  
  # Create heatmap
  pheatmap(avg_expr_scaled,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-2, 2, length.out = 101),
           main = "Top Drug Targets - Expression Across Cell Types",
           fontsize = 10,
           fontsize_row = 8,
           angle_col = 45,
           filename = file.path("DrugTarget_Results", "DrugTarget_Expression_Heatmap.pdf"),
           width = 10,
           height = 12)
}

# --- 8.3 Target Family Distribution ---
if (exists("druggable_resistance_genes")) {
  
  family_counts <- druggable_resistance_genes %>%
    filter(avg_log2FC > 0.5) %>%
    count(family) %>%
    arrange(desc(n))
  
  p_families <- ggplot(family_counts, aes(x = reorder(family, n), y = n, fill = family)) +
    geom_col() +
    coord_flip() +
    labs(title = "Druggable Gene Families in Resistance",
         x = "Gene Family",
         y = "Number of Genes") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          legend.position = "none",
          panel.border = element_rect(fill = NA, color = "black"))
  
  print(p_families)
  
  ggsave(file.path("DrugTarget_Results", "Druggable_Family_Distribution.pdf"),
         p_families, width = 8, height = 6)
}

# =============================================================================
# 9. CLINICAL TRIAL RECOMMENDATIONS
# =============================================================================

cat("\n=== 9. CLINICAL TRIAL RECOMMENDATIONS ===\n\n")

cat("  Based on the analysis, we recommend investigating:\n\n")

if (exists("top_priorities")) {
  # Get top 5 targets
  top_5 <- head(top_priorities, 5)
  
  cat("  HIGH PRIORITY TARGETS:\n")
  for (i in 1:nrow(top_5)) {
    cat(sprintf("    %d. %s (%s family)\n", i, top_5$gene[i], top_5$family[i]))
    cat(sprintf("       - Resistance FC: %.2f\n", top_5$resistance_FC[i]))
    cat(sprintf("       - Expression in resistant cells: %.2f\n", 
                top_5$resistant_expression[i]))
    cat(sprintf("       - Priority score: %.3f\n\n", top_5$priority_score[i]))
  }
}

cat("  RECOMMENDED COMBINATION TRIALS:\n")
for (i in 1:min(3, nrow(combination_strategies))) {
  cat(sprintf("    %d. %s\n", i, combination_strategies$Strategy[i]))
  cat(sprintf("       %s + %s\n\n", combination_strategies$Target1[i],
              combination_strategies$Target2[i]))
}

# =============================================================================
# 10. SAVE SESSION INFO
# =============================================================================

save_session_info(file.path("DrugTarget_Results", "session_info_drugtarget.txt"))

cat("\n=============================================================================\n")
cat("✓ DRUG TARGET IDENTIFICATION ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
cat(paste("\nResults saved in:", file.path(getwd(), "DrugTarget_Results"), "\n"))
cat("\nKey outputs:\n")
cat("  - Druggable resistance genes\n")
cat("  - Cell type-specific targets\n")
cat("  - Treatment-specific targets\n")
cat("  - Prioritized target ranking\n")
cat("  - Combination therapy strategies\n")
cat("  - Clinical trial recommendations\n")
cat("  - Expression heatmaps\n")
cat("=============================================================================\n\n")
