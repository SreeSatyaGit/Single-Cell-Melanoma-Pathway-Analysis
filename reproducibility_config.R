# =============================================================================
# REPRODUCIBILITY CONFIGURATION
# =============================================================================
# This file ensures reproducibility across all analyses in the SCMPA project
# Source this file at the beginning of each analysis script
# =============================================================================

# Set global random seed for reproducibility
GLOBAL_SEED <- 42

# Set seed for R's random number generator
set.seed(GLOBAL_SEED)

# Print reproducibility information
cat("=============================================================================\n")
cat("REPRODUCIBILITY CONFIGURATION LOADED\n")
cat("=============================================================================\n")
cat("Global seed:", GLOBAL_SEED, "\n")
cat("Date:", Sys.Date(), "\n")
cat("Time:", Sys.time(), "\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("=============================================================================\n\n")

# Function to set seed before any stochastic operation
set_analysis_seed <- function(step_name = NULL, seed = GLOBAL_SEED) {
  set.seed(seed)
  if (!is.null(step_name)) {
    cat(sprintf("✓ Seed set to %d for: %s\n", seed, step_name))
  }
}

# Function to save session information
save_session_info <- function(output_file = "session_info.txt") {
  sink(output_file)
  cat("=============================================================================\n")
  cat("SESSION INFORMATION\n")
  cat("=============================================================================\n")
  cat("Date:", as.character(Sys.Date()), "\n")
  cat("Time:", as.character(Sys.time()), "\n\n")
  
  cat("R Version:\n")
  print(R.version)
  cat("\n")
  
  cat("Platform Information:\n")
  print(sessionInfo()$platform)
  cat("\n")
  
  cat("Loaded Packages:\n")
  print(sessionInfo()$otherPkgs)
  cat("\n")
  
  cat("Base Packages:\n")
  print(sessionInfo()$basePkgs)
  cat("\n")
  
  cat("Detailed Session Info:\n")
  print(sessionInfo())
  cat("\n")
  
  cat("Installed Package Versions:\n")
  ip <- as.data.frame(installed.packages()[, c("Package", "Version")])
  print(ip[order(ip$Package), ])
  
  sink()
  cat(sprintf("✓ Session information saved to: %s\n", output_file))
}

# Function to create reproducibility report
create_reproducibility_report <- function(script_name, output_dir = ".") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  report_file <- file.path(output_dir, sprintf("reproducibility_%s_%s.txt", 
                                                gsub("\\.R$", "", script_name), 
                                                timestamp))
  
  sink(report_file)
  cat("=============================================================================\n")
  cat("REPRODUCIBILITY REPORT\n")
  cat("=============================================================================\n")
  cat("Script:", script_name, "\n")
  cat("Date:", as.character(Sys.Date()), "\n")
  cat("Time:", as.character(Sys.time()), "\n")
  cat("Global Seed:", GLOBAL_SEED, "\n")
  cat("Working Directory:", getwd(), "\n\n")
  
  cat("System Information:\n")
  cat("  OS:", Sys.info()["sysname"], "\n")
  cat("  Release:", Sys.info()["release"], "\n")
  cat("  Machine:", Sys.info()["machine"], "\n\n")
  
  cat("R Session Information:\n")
  print(sessionInfo())
  
  sink()
  cat(sprintf("✓ Reproducibility report saved to: %s\n", report_file))
  return(report_file)
}

# Function to verify package versions
verify_package_versions <- function(required_packages = NULL) {
  if (is.null(required_packages)) {
    # Default required packages for SCMPA project
    required_packages <- list(
      "Seurat" = "5.0.0",
      "SCPA" = NULL,  # Any version
      "dyno" = NULL,
      "msigdbr" = NULL,
      "tidyverse" = NULL,
      "ComplexHeatmap" = NULL
    )
  }
  
  cat("\n=============================================================================\n")
  cat("PACKAGE VERSION VERIFICATION\n")
  cat("=============================================================================\n")
  
  all_ok <- TRUE
  for (pkg in names(required_packages)) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      installed_ver <- as.character(packageVersion(pkg))
      required_ver <- required_packages[[pkg]]
      
      if (is.null(required_ver)) {
        cat(sprintf("✓ %s: %s (any version accepted)\n", pkg, installed_ver))
      } else {
        if (package_version(installed_ver) >= package_version(required_ver)) {
          cat(sprintf("✓ %s: %s (required: >= %s)\n", pkg, installed_ver, required_ver))
        } else {
          cat(sprintf("✗ %s: %s (required: >= %s) - UPDATE NEEDED\n", 
                     pkg, installed_ver, required_ver))
          all_ok <- FALSE
        }
      }
    } else {
      cat(sprintf("✗ %s: NOT INSTALLED\n", pkg))
      all_ok <- FALSE
    }
  }
  
  cat("=============================================================================\n\n")
  
  if (!all_ok) {
    warning("Some package version requirements are not met!")
  }
  
  return(all_ok)
}

# Function to set parallel processing seeds
set_parallel_seeds <- function(n_cores = 4, seed = GLOBAL_SEED) {
  # For parallel processing reproducibility
  if (requireNamespace("parallel", quietly = TRUE)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    cat(sprintf("✓ Parallel RNG initialized with seed %d for %d cores\n", seed, n_cores))
  }
}

# Export key variables
REPRODUCIBILITY_CONFIG <- list(
  seed = GLOBAL_SEED,
  date = Sys.Date(),
  r_version = R.version.string,
  platform = R.version$platform
)

cat("✓ Reproducibility configuration loaded successfully\n")
cat("  Use set_analysis_seed('step_name') before stochastic operations\n")
cat("  Use save_session_info() to save package versions\n")
cat("  Use create_reproducibility_report('script_name.R') for detailed reports\n\n")
