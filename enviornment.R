# environment.R
# Install all required packages for RNA-seq analysis
# Run this once before starting the analysis

cat("🔧 Starting package installation...\n")
cat("⏱️  This may take 10-20 minutes\n\n")

# Install CRAN packages
cran_packages <- c(
  "tidyverse",      # Data manipulation & visualization
  "pheatmap",       # Heatmaps
  "RColorBrewer",   # Color palettes
  "ggrepel",        # Label handling in plots
  "plotly",         # Interactive plots
  "DT",             # Interactive tables
  "openxlsx",       # Excel export
  "viridis",        # Color scales
  "gridExtra",      # Arrange multiple plots
  "cowplot",        # Publication-ready plots
  "ggpubr"          # Publication-ready plots
)

cat("📦 Installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("   Installing %s...\n", pkg))
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(sprintf("   ✓ %s already installed\n", pkg))
  }
}

# Install BiocManager
cat("\n📦 Installing BiocManager...\n")
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor packages
bioc_packages <- c(
  "DESeq2",           # Differential expression analysis
  "edgeR",            # Alternative DE tool
  "limma",            # Linear models
  "biomaRt",          # Gene annotation
  "AnnotationDbi",    # Annotation databases
  "org.Hs.eg.db",     # Human gene annotations
  "clusterProfiler",  # Pathway enrichment
  "enrichplot",       # Enrichment visualization
  "DOSE",             # Disease ontology
  "pathview",         # KEGG pathway visualization
  "genefilter",       # Gene filtering
  "vsn",              # Variance stabilization
  "ComplexHeatmap",   # Advanced heatmaps
  "EnhancedVolcano"   # Publication-ready volcano plots
)

cat("\n📦 Installing Bioconductor packages...\n")
cat("   (This is the slower part...)\n\n")

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("   Installing %s...\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("   ✓ %s already installed\n", pkg))
  }
}

# Test loading key libraries
cat("\n🧪 Testing key packages...\n")
test_packages <- c("DESeq2", "tidyverse", "clusterProfiler", 
                   "pheatmap", "EnhancedVolcano")

all_loaded <- TRUE
for (pkg in test_packages) {
  test <- try(library(pkg, character.only = TRUE), silent = TRUE)
  if (inherits(test, "try-error")) {
    cat(sprintf("   ✗ %s failed to load\n", pkg))
    all_loaded <- FALSE
  } else {
    cat(sprintf("   ✓ %s loaded successfully\n", pkg))
  }
}

if (all_loaded) {
  cat("\n✅ All packages installed and tested successfully!\n")
  cat("📊 Ready to start RNA-seq analysis!\n\n")
} else {
  cat("\n⚠️  Some packages failed. Please report errors.\n\n")
}

# Print session info for reproducibility
cat("📋 Session Information:\n")
cat("=" %R% 50, "\n")
sessionInfo()