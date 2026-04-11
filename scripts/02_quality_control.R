# scripts/02_quality_control.R
# Quality control and exploratory data analysis of RNA-seq data
# This script performs essential QC checks before differential expression analysis

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)

cat("🔍 Starting Quality Control Analysis...\n\n")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("📂 Loading data...\n")

counts <- readRDS("data/raw/count_matrix.rds")
metadata <- readRDS("data/raw/sample_metadata.rds")

cat(sprintf("   ✓ Count matrix: %d genes × %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("   ✓ Metadata: %d samples × %d variables\n\n", nrow(metadata), ncol(metadata)))

# Verify samples match
if (!all(colnames(counts) == rownames(metadata))) {
  stop("❌ Sample names don't match between counts and metadata!")
}

# ============================================================================
# 2. LIBRARY SIZE (Total Counts per Sample)
# ============================================================================

cat("📊 Analyzing library sizes...\n")

# Calculate total counts per sample
lib_sizes <- data.frame(
  sample = colnames(counts),
  total_counts = colSums(counts),
  condition = metadata$condition
)

# Summary statistics
cat("\n   Library Size Summary:\n")
cat(sprintf("   Min:    %s\n", format(min(lib_sizes$total_counts), big.mark = ",")))
cat(sprintf("   Median: %s\n", format(median(lib_sizes$total_counts), big.mark = ",")))
cat(sprintf("   Max:    %s\n", format(max(lib_sizes$total_counts), big.mark = ",")))
cat(sprintf("   Mean:   %s\n", format(mean(lib_sizes$total_counts), big.mark = ",")))

# Plot library sizes
p1 <- ggplot(lib_sizes, aes(x = reorder(sample, total_counts), 
                            y = total_counts/1e6, 
                            fill = condition)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Normal" = "#3498db", "Tumor" = "#e74c3c")) +
  labs(title = "Library Size Distribution",
       subtitle = "Total read counts per sample",
       x = "Sample",
       y = "Total Counts (Millions)",
       fill = "Condition") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "top")

ggsave("results/figures/01_library_sizes.png", p1, width = 10, height = 6, dpi = 300)
cat("   ✓ Saved: results/figures/01_library_sizes.png\n\n")

# ============================================================================
# 3. GENE DETECTION (How many genes detected per sample?)
# ============================================================================

cat("🧬 Analyzing gene detection...\n")

# Count genes with at least 1, 10, and 100 counts per sample
gene_detection <- data.frame(
  sample = colnames(counts),
  genes_detected_1 = colSums(counts >= 1),
  genes_detected_10 = colSums(counts >= 10),
  genes_detected_100 = colSums(counts >= 100),
  condition = metadata$condition
)

cat("\n   Genes Detected (≥10 counts):\n")
cat(sprintf("   Min:    %s\n", format(min(gene_detection$genes_detected_10), big.mark = ",")))
cat(sprintf("   Median: %s\n", format(median(gene_detection$genes_detected_10), big.mark = ",")))
cat(sprintf("   Max:    %s\n", format(max(gene_detection$genes_detected_10), big.mark = ",")))
# Reshape for plotting
gene_det_long <- gene_detection %>%
  select(sample, condition, genes_detected_1, genes_detected_10, genes_detected_100) %>%
  pivot_longer(cols = starts_with("genes"), 
               names_to = "threshold", 
               values_to = "count") %>%
  mutate(threshold = factor(threshold, 
                            levels = c("genes_detected_1", "genes_detected_10", "genes_detected_100"),
                            labels = c("≥1 count", "≥10 counts", "≥100 counts")))

p2 <- ggplot(gene_det_long, aes(x = sample, y = count, fill = condition)) +
  geom_bar(stat = "identity") +
  facet_wrap(~threshold, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Normal" = "#3498db", "Tumor" = "#e74c3c")) +
  labs(title = "Gene Detection Across Samples",
       subtitle = "Number of genes detected at different count thresholds",
       x = "Sample",
       y = "Number of Genes Detected",
       fill = "Condition") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "top",
        strip.background = element_rect(fill = "grey90"))

ggsave("results/figures/02_gene_detection.png", p2, width = 10, height = 10, dpi = 300)
cat("   ✓ Saved: results/figures/02_gene_detection.png\n\n")

# ============================================================================
# 4. COUNT DISTRIBUTION
# ============================================================================

cat("📈 Analyzing count distributions...\n")

# Filter genes with very low counts (keep genes with ≥10 counts in ≥3 samples)
keep <- rowSums(counts >= 10) >= 3
counts_filtered <- counts[keep, ]

cat(sprintf("   Genes before filtering: %d\n", nrow(counts)))
cat(sprintf("   Genes after filtering:  %d\n", nrow(counts_filtered)))
cat(sprintf("   Genes removed: %d (%.1f%%)\n\n", 
            nrow(counts) - nrow(counts_filtered),
            100 * (nrow(counts) - nrow(counts_filtered)) / nrow(counts)))

# Log2 transformation for visualization (add pseudocount to avoid log(0))
log_counts <- log2(counts_filtered + 1)

# Create boxplot of log-transformed counts
log_counts_long <- as.data.frame(log_counts) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "log2_count") %>%
  left_join(metadata %>% rownames_to_column("sample"), by = "sample")

p3 <- ggplot(log_counts_long, aes(x = sample, y = log2_count, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("Normal" = "#3498db", "Tumor" = "#e74c3c")) +
  labs(title = "Distribution of Gene Expression Across Samples",
       subtitle = "Log2(counts + 1) for genes with ≥10 counts in ≥3 samples",
       x = "Sample",
       y = "Log2(Count + 1)",
       fill = "Condition") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "top")

ggsave("results/figures/03_count_distribution.png", p3, width = 12, height = 6, dpi = 300)
cat("   ✓ Saved: results/figures/03_count_distribution.png\n\n")

# ============================================================================
# 5. SAMPLE CORRELATION HEATMAP (IMPROVED VERSION)
# ============================================================================

cat("🔥 Creating sample correlation heatmap...\n")

# Calculate pairwise correlation between samples
sample_cor <- cor(log_counts, method = "pearson")

# Annotation for heatmap
annotation_col <- data.frame(
  Condition = metadata$condition,
  Batch = metadata$batch,
  row.names = rownames(metadata)
)

annotation_colors <- list(
  Condition = c(Normal = "#3498db", Tumor = "#e74c3c"),
  Batch = c(A = "#95a5a6", B = "#34495e")
)

# Create heatmap with better settings
png("results/figures/04_sample_correlation.png", width = 10, height = 9, units = "in", res = 300)
pheatmap(sample_cor,
         annotation_col = annotation_col,
         annotation_row = annotation_col,  # Add row annotation too
         annotation_colors = annotation_colors,
         color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
         main = "Sample-to-Sample Correlation Heatmap",
         fontsize = 10,
         fontsize_row = 9,
         fontsize_col = 9,
         border_color = NA,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",  # Better clustering
         cutree_rows = 2,  # Show 2 main clusters
         cutree_cols = 2,  # Show 2 main clusters
         display_numbers = FALSE,  # Don't show numbers (too crowded)
         legend = TRUE,
         angle_col = 45)
dev.off()

cat("   ✓ Saved: results/figures/04_sample_correlation.png\n\n")
# ============================================================================
# 6. PCA ANALYSIS
# ============================================================================

cat("🎯 Performing PCA analysis...\n")

# Create DESeq2 object for variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ condition
)

# Variance Stabilizing Transformation (better than log for PCA)
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

# Perform PCA
pca_result <- prcomp(t(vsd_mat), scale. = FALSE)

# Calculate variance explained
variance_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

# Create PCA dataframe
pca_data <- data.frame(
  sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  condition = metadata$condition,
  batch = metadata$batch
)

cat(sprintf("   PC1 explains: %.1f%% of variance\n", variance_explained[1]))
cat(sprintf("   PC2 explains: %.1f%% of variance\n\n", variance_explained[2]))

# PCA plot colored by condition
p4 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Normal" = "#3498db", "Tumor" = "#e74c3c")) +
  labs(title = "Principal Component Analysis (PCA)",
       subtitle = "Variance-stabilized transformed data",
       x = paste0("PC1 (", variance_explained[1], "% variance)"),
       y = paste0("PC2 (", variance_explained[2], "% variance)"),
       color = "Condition",
       shape = "Batch") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "right",
        panel.grid.minor = element_blank())

ggsave("results/figures/05_pca_plot.png", p4, width = 10, height = 7, dpi = 300)
cat("   ✓ Saved: results/figures/05_pca_plot.png\n\n")

# ============================================================================
# 7. SAMPLE DISTANCE HEATMAP
# ============================================================================

cat("📏 Creating sample distance heatmap...\n")

# Calculate Euclidean distances between samples
sample_dists <- dist(t(vsd_mat))
sample_dist_matrix <- as.matrix(sample_dists)

# Create heatmap
png("results/figures/06_sample_distances.png", width = 10, height = 9, units = "in", res = 300)
pheatmap(sample_dist_matrix,
         annotation_col = annotation_col,
         annotation_row = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
         main = "Sample-to-Sample Distance Heatmap",
         fontsize = 10,
         fontsize_row = 9,
         fontsize_col = 9,
         border_color = NA)
dev.off()

cat("   ✓ Saved: results/figures/06_sample_distances.png\n\n")

# ============================================================================
# 8. SAVE PROCESSED DATA AND QC SUMMARY
# ============================================================================

cat("💾 Saving processed data...\n")

# Save filtered counts
saveRDS(counts_filtered, "data/processed/counts_filtered.rds")
saveRDS(vsd_mat, "data/processed/vsd_matrix.rds")

# Create QC summary table
qc_summary <- data.frame(
  Sample = rownames(metadata),
  Condition = metadata$condition,
  Batch = metadata$batch,
  Total_Counts = lib_sizes$total_counts,
  Genes_Detected = gene_detection$genes_detected_10,
  Percent_Genes_Detected = round(100 * gene_detection$genes_detected_10 / nrow(counts), 1)
)

write.csv(qc_summary, "results/tables/qc_summary.csv", row.names = FALSE)
cat("   ✓ Saved: data/processed/counts_filtered.rds\n")
cat("   ✓ Saved: data/processed/vsd_matrix.rds\n")
cat("   ✓ Saved: results/tables/qc_summary.csv\n\n")

# Display summary
cat("📊 QC SUMMARY TABLE:\n")
cat("=" %R% 80, "\n")
print(qc_summary)

# ============================================================================
# 9. FINAL QC ASSESSMENT
# ============================================================================

cat("\n\n✅ QUALITY CONTROL COMPLETE!\n")
cat("=" %R% 80, "\n\n")

cat("📋 QC Assessment:\n\n")

# Check for outliers based on library size
lib_size_median <- median(lib_sizes$total_counts)
lib_size_mad <- mad(lib_sizes$total_counts)
outliers <- abs(lib_sizes$total_counts - lib_size_median) > 3 * lib_size_mad

if (any(outliers)) {
  cat("   ⚠️  Potential outliers detected (library size):\n")
  cat(paste("      -", lib_sizes$sample[outliers], collapse = "\n"))
  cat("\n\n")
} else {
  cat("   ✓ No outliers detected based on library size\n")
}

# Check correlation
min_cor <- min(sample_cor[lower.tri(sample_cor)])
if (min_cor < 0.8) {
  cat(sprintf("   ⚠️  Low correlation detected: %.3f (should be >0.8)\n", min_cor))
} else {
  cat(sprintf("   ✓ Good sample correlation (min: %.3f)\n", min_cor))
}

# Check PCA separation
pc1_range_tumor <- range(pca_data$PC1[pca_data$condition == "Tumor"])
pc1_range_normal <- range(pca_data$PC1[pca_data$condition == "Normal"])
pca_separation <- !((pc1_range_tumor[1] > pc1_range_normal[2]) || 
                      (pc1_range_normal[1] > pc1_range_tumor[2]))

if (pca_separation) {
  cat("   ✓ Good separation between conditions in PCA\n")
} else {
  cat("   ⚠️  Conditions overlap in PCA - may affect DE analysis\n")
}

cat("\n📁 Generated Files:\n")
cat("   - 6 publication-quality figures\n")
cat("   - 1 QC summary table\n")
cat("   - 2 processed data files\n")

cat("\n🚀 Next step: Differential expression analysis (Script 03)\n\n")