# Fix for Corrupted Packrat GlobalOptions Package
# This script provides solutions for the packrat library corruption issue

message("=== Packrat Library Corruption Fix ===\n")

# Option 1: Reinstall GlobalOptions outside of packrat
message("Option 1: Reinstalling GlobalOptions...")
tryCatch({
  # Remove corrupted package
  remove.packages("GlobalOptions")
  
  # Reinstall from Bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("GlobalOptions", update = FALSE, ask = FALSE)
  
  message("✓ GlobalOptions reinstalled successfully")
}, error = function(e) {
  message("✗ Reinstallation failed:", e$message)
})

# Option 2: Reinstall ComplexHeatmap and dependencies
message("\nOption 2: Reinstalling ComplexHeatmap and all dependencies...")
tryCatch({
  # Remove potentially corrupted packages
  pkgs_to_remove <- c("GlobalOptions", "ComplexHeatmap", "circlize")
  for (pkg in pkgs_to_remove) {
    if (pkg %in% installed.packages()[,"Package"]) {
      remove.packages(pkg)
      message(paste("Removed:", pkg))
    }
  }
  
  # Reinstall fresh
  BiocManager::install(c("GlobalOptions", "circlize", "ComplexHeatmap"), 
                       update = FALSE, ask = FALSE, force = TRUE)
  
  message("✓ ComplexHeatmap and dependencies reinstalled")
}, error = function(e) {
  message("✗ Reinstallation failed:", e$message)
})

# Option 3: Test if ComplexHeatmap works now
message("\nOption 3: Testing ComplexHeatmap...")
test_result <- tryCatch({
  library(ComplexHeatmap)
  library(circlize)
  
  # Create a simple test heatmap
  test_mat <- matrix(rnorm(100), 10, 10)
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  ht <- Heatmap(test_mat, col = col_fun, name = "test")
  
  message("✓ ComplexHeatmap is working correctly!")
  TRUE
}, error = function(e) {
  message("✗ ComplexHeatmap still has issues:", e$message)
  message("\nThe MelanomaStates.R script will automatically use pheatmap or ggplot2 instead.")
  FALSE
})

# Option 4: Disable packrat for this session (if using packrat)
if (exists(".packrat_mode") || file.exists("packrat/packrat.lock")) {
  message("\nOption 4: Packrat detected")
  message("You may want to disable packrat temporarily:")
  message("  packrat::off()")
  message("  # Then reinstall packages")
  message("  # Then run: packrat::on()")
}

# Summary
message("\n=== Summary ===")
if (test_result) {
  message("✓ ComplexHeatmap is working - you can use MelanomaStates.R normally")
} else {
  message("⚠ ComplexHeatmap has issues, but MelanomaStates.R will use alternatives")
  message("  - The script will automatically use pheatmap or ggplot2")
  message("  - All visualizations will still be created")
  message("  - No action needed on your part")
}

message("\n=== Alternative: Use pheatmap directly ===")
message("If you prefer, you can use pheatmap which doesn't have this issue:")
message("  install.packages('pheatmap')")
message("  library(pheatmap)")
