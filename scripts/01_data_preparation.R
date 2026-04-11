# scripts/01_data_preparation.R
# Create simulated RNA-seq count data for differential expression analysis
# This simulates a realistic colorectal cancer vs normal tissue dataset

library(tidyverse)

cat("📊 Creating simulated RNA-seq dataset...\n\n")

set.seed(123)  # For reproducibility

# Parameters
n_genes <- 20000      # Total number of genes
n_samples <- 12       # Total samples (6 tumor, 6 normal)

# Create sample names
sample_names <- c(
  paste0("Tumor_", 1:6),
  paste0("Normal_", 1:6)
)

# Create gene names (using realistic gene symbols)
# Mix of real-looking gene IDs
gene_names <- paste0("GENE_", sprintf("%05d", 1:n_genes))

cat("🧬 Simulating count matrix...\n")
cat(sprintf("   Genes: %d\n", n_genes))
cat(sprintf("   Samples: %d (%d tumor, %d normal)\n", n_samples, 6, 6))

# Initialize count matrix
counts <- matrix(0, nrow = n_genes, ncol = n_samples)
rownames(counts) <- gene_names
colnames(counts) <- sample_names

# Simulate realistic gene expression levels
# Base expression levels (most genes lowly expressed, few highly expressed)
base_expression <- rnbinom(n_genes, mu = 100, size = 2)

# Gene categories:
# - 16,000 genes (80%): No difference between conditions
# - 2,000 genes (10%): Upregulated in tumor
# - 2,000 genes (10%): Downregulated in tumor

cat("\n🎯 Creating differential expression patterns...\n")
cat("   - 80% genes: No change\n")
cat("   - 10% genes: Upregulated in tumor\n")
cat("   - 10% genes: Downregulated in tumor\n\n")

for (i in 1:n_genes) {
  
  if (i <= 16000) {
    # No significant difference - same distribution for all samples
    counts[i, ] <- rnbinom(n_samples, mu = base_expression[i], size = 2)
    
  } else if (i <= 18000) {
    # Upregulated in tumor (2-8 fold change)
    fold_change <- runif(1, 2, 8)
    counts[i, 1:6] <- rnbinom(6, mu = base_expression[i] * fold_change, size = 2)
    counts[i, 7:12] <- rnbinom(6, mu = base_expression[i], size = 2)
    
  } else {
    # Downregulated in tumor (2-8 fold change)
    fold_change <- runif(1, 2, 8)
    counts[i, 1:6] <- rnbinom(6, mu = base_expression[i] / fold_change, size = 2)
    counts[i, 7:12] <- rnbinom(6, mu = base_expression[i], size = 2)
  }
}

cat("✅ Count matrix created\n\n")

# Create sample metadata
cat("📋 Creating sample metadata...\n")

metadata <- data.frame(
  sample_id = sample_names,
  condition = factor(c(rep("Tumor", 6), rep("Normal", 6)), 
                     levels = c("Normal", "Tumor")),  # Normal as reference
  batch = factor(rep(c("A", "B", "A", "B", "A", "B"), 2)),
  patient_id = factor(c(1:6, 1:6)),
  row.names = sample_names
)

# Add some clinical/demographic variables
metadata$age <- c(65, 58, 72, 55, 68, 61, 64, 59, 70, 56, 67, 62)
metadata$gender <- factor(c("M", "F", "M", "M", "F", "F", 
                            "M", "F", "M", "M", "F", "F"))
metadata$stage <- factor(c("III", "II", "IV", "II", "III", "II",
                           NA, NA, NA, NA, NA, NA))  # Stage only for tumors

cat("✅ Metadata created\n\n")

# Display summary
cat("📊 DATASET SUMMARY\n")
cat("==================\n\n")
cat("Count Matrix:\n")
cat(sprintf("  Dimensions: %d genes × %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("  Total counts: %s\n", format(sum(counts), big.mark = ",")))
cat(sprintf("  Mean counts per gene: %.0f\n", mean(rowSums(counts))))
cat(sprintf("  Median counts per gene: %.0f\n", median(rowSums(counts))))

cat("\nSample Metadata:\n")
print(metadata)

cat("\n\nFirst 5 genes × 6 samples:\n")
print(counts[1:5, 1:6])

# Save data
cat("\n💾 Saving data files...\n")

# Make sure directories exist
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Save as CSV
write.csv(counts, "data/raw/count_matrix.csv", row.names = TRUE)
write.csv(metadata, "data/raw/sample_metadata.csv", row.names = TRUE)

# Also save as RDS for faster loading in R
saveRDS(counts, "data/raw/count_matrix.rds")
saveRDS(metadata, "data/raw/sample_metadata.rds")

cat("\n✅ FILES SAVED:\n")
cat("   📁 data/raw/count_matrix.csv\n")
cat("   📁 data/raw/sample_metadata.csv\n")
cat("   📁 data/raw/count_matrix.rds\n")
cat("   📁 data/raw/sample_metadata.rds\n")

cat("\n🎉 Data preparation complete!\n")
cat("📌 Next step: Quality control analysis\n\n")