# scripts/03_differential_expression.R
# Differential expression analysis using DESeq2
# Identifies genes differentially expressed between Tumor and Normal tissue

library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

cat("🧬 DIFFERENTIAL EXPRESSION ANALYSIS WITH DESeq2\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("📂 Loading processed data...\n")

counts_filtered <- readRDS("data/processed/counts_filtered.rds")
metadata <- readRDS("data/raw/sample_metadata.rds")

cat(sprintf("   ✓ Count matrix: %d genes × %d samples\n", nrow(counts_filtered), ncol(counts_filtered)))
cat(sprintf("   ✓ Metadata: %d samples\n\n", nrow(metadata)))

# Verify sample order matches
if (!all(colnames(counts_filtered) == rownames(metadata))) {
  stop("❌ Sample order mismatch!")
}

# ============================================================================
# 2. CREATE DESeq2 OBJECT
# ============================================================================

cat("🔬 Creating DESeq2 dataset...\n")

# Design formula: account for batch effects + condition
# ~ batch + condition means: control for batch, then test condition effect

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ batch + condition
)

cat("   ✓ DESeq2 object created\n")
cat(sprintf("   Design formula: ~ batch + condition\n"))
cat(sprintf("   Reference level: %s (vs %s)\n\n", 
            levels(metadata$condition)[1], 
            levels(metadata$condition)[2]))

# ============================================================================
# 3. RUN DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================

cat("⚡ Running DESeq2 analysis...\n")
cat("   (This may take 1-2 minutes for ~19,779 genes)\n\n")

# Run the DE analysis
dds <- DESeq(dds)

cat("   ✓ DESeq2 analysis complete!\n\n")

# Get results
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))

cat("📊 DESeq2 Results Summary:\n")
summary(res)

# ============================================================================
# 4. EXTRACT AND ANNOTATE RESULTS
# ============================================================================

cat("\n📋 Processing results...\n")

# Convert to data frame
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

# Add significance categories
res_df <- res_df %>%
  mutate(
    significant = case_when(
      is.na(padj) ~ "Not tested",
      padj < 0.05 & abs(log2FoldChange) >= 1 ~ "Significant",
      padj < 0.05 ~ "Padj < 0.05 only",
      TRUE ~ "Not significant"
    ),
    direction = case_when(
      significant == "Significant" & log2FoldChange > 0 ~ "Upregulated",
      significant == "Significant" & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Not DE"
    )
  )

# Summary statistics
n_total <- nrow(res_df)
n_sig <- sum(res_df$significant == "Significant", na.rm = TRUE)
n_up <- sum(res_df$direction == "Upregulated", na.rm = TRUE)
n_down <- sum(res_df$direction == "Downregulated", na.rm = TRUE)

cat(sprintf("\n   Total genes tested: %s\n", format(n_total, big.mark = ",")))
cat(sprintf("   Significant DE genes (padj < 0.05, |LFC| >= 1): %s (%.1f%%)\n", 
            format(n_sig, big.mark = ","),
            100 * n_sig / n_total))
cat(sprintf("      - Upregulated in Tumor: %s\n", format(n_up, big.mark = ",")))
cat(sprintf("      - Downregulated in Tumor: %s\n\n", format(n_down, big.mark = ",")))

# ============================================================================
# 5. VOLCANO PLOT
# ============================================================================

cat("🌋 Creating volcano plot...\n")

# Enhanced volcano plot
p_volcano <- EnhancedVolcano(res_df,
                             lab = res_df$gene,
                             x = 'log2FoldChange',
                             y = 'padj',
                             title = 'Differential Expression: Tumor vs Normal',
                             subtitle = sprintf('%d genes upregulated, %d downregulated', n_up, n_down),
                             pCutoff = 0.05,
                             FCcutoff = 1.0,
                             pointSize = 1.5,
                             labSize = 3.0,
                             labCol = 'black',
                             labFace = 'bold',
                             boxedLabels = TRUE,
                             colAlpha = 0.5,
                             legendPosition = 'right',
                             legendLabSize = 10,
                             legendIconSize = 3.0,
                             drawConnectors = TRUE,
                             widthConnectors = 0.3,
                             colConnectors = 'grey30',
                             max.overlaps = 15,
                             col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                             xlim = c(min(res_df$log2FoldChange, na.rm = TRUE) - 1, 
                                      max(res_df$log2FoldChange, na.rm = TRUE) + 1))

ggsave("results/figures/07_volcano_plot.png", p_volcano, 
       width = 12, height = 10, dpi = 300)

cat("   ✓ Saved: results/figures/07_volcano_plot.png\n\n")

# ============================================================================
# 6. MA PLOT
# ============================================================================

cat("📈 Creating MA plot...\n")

png("results/figures/08_MA_plot.png", width = 10, height = 8, units = "in", res = 300)
plotMA(res, 
       main = "MA Plot: Tumor vs Normal",
       ylim = c(-5, 5),
       colNonSig = "gray60",
       colSig = "red3",
       colLine = "blue")
dev.off()

cat("   ✓ Saved: results/figures/08_MA_plot.png\n\n")

# ============================================================================
# 7. HEATMAP OF TOP DE GENES
# ============================================================================

cat("🔥 Creating heatmap of top DE genes...\n")

# Get top 50 DE genes (25 up, 25 down)
top_genes <- res_df %>%
  filter(significant == "Significant") %>%
  group_by(direction) %>%
  slice_min(order_by = padj, n = 25) %>%
  ungroup() %>%
  pull(gene)

cat(sprintf("   Selected %d top DE genes for heatmap\n", length(top_genes)))

# Get normalized counts for heatmap
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# Subset to top genes
heatmap_data <- vsd_mat[top_genes, ]

# Scale by row (z-score)
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# Annotation
annotation_col <- data.frame(
  Condition = metadata$condition,
  Batch = metadata$batch,
  row.names = rownames(metadata)
)

annotation_colors <- list(
  Condition = c(Normal = "#3498db", Tumor = "#e74c3c"),
  Batch = c(A = "#95a5a6", B = "#34495e")
)

# Create heatmap
png("results/figures/09_heatmap_top_genes.png", 
    width = 10, height = 12, units = "in", res = 300)

pheatmap(heatmap_data_scaled,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         main = "Top 50 Differentially Expressed Genes\n(Z-score normalized)",
         fontsize = 8,
         fontsize_row = 6,
         fontsize_col = 9,
         border_color = NA,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         cutree_cols = 2,
         breaks = seq(-2, 2, length.out = 100))

dev.off()

cat("   ✓ Saved: results/figures/09_heatmap_top_genes.png\n\n")

# ============================================================================
# 8. TOP UPREGULATED AND DOWNREGULATED GENES
# ============================================================================

cat("📊 Identifying top genes...\n\n")

# Top 10 upregulated
top_up <- res_df %>%
  filter(direction == "Upregulated") %>%
  slice_min(order_by = padj, n = 10)

cat("🔺 TOP 10 UPREGULATED GENES (Tumor vs Normal):\n")
cat(paste(rep("-", 80), collapse=""), "\n")
print(top_up %>% 
        select(gene, log2FoldChange, padj) %>%
        mutate(FoldChange = 2^log2FoldChange,
               padj = format(padj, scientific = TRUE, digits = 3)) %>%
        select(gene, FoldChange, log2FoldChange, padj),
      row.names = FALSE)

# Top 10 downregulated
top_down <- res_df %>%
  filter(direction == "Downregulated") %>%
  slice_min(order_by = padj, n = 10)

cat("\n🔻 TOP 10 DOWNREGULATED GENES (Tumor vs Normal):\n")
cat(paste(rep("-", 80), collapse=""), "\n")
print(top_down %>%
        select(gene, log2FoldChange, padj) %>%
        mutate(FoldChange = 2^log2FoldChange,
               padj = format(padj, scientific = TRUE, digits = 3)) %>%
        select(gene, FoldChange, log2FoldChange, padj),
      row.names = FALSE)

# ============================================================================
# 9. FOLD CHANGE DISTRIBUTION
# ============================================================================

cat("\n\n📊 Creating fold change distribution plot...\n")

p_fc_dist <- ggplot(res_df %>% filter(!is.na(padj)), 
                    aes(x = log2FoldChange, fill = significant)) +
  geom_histogram(bins = 60, alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Significant" = "#e74c3c",
                               "Not significant" = "#95a5a6",
                               "Padj < 0.05 only" = "#f39c12")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(title = "Distribution of Log2 Fold Changes",
       subtitle = "Dashed lines: |LFC| = 1 threshold",
       x = "Log2 Fold Change (Tumor vs Normal)",
       y = "Number of Genes",
       fill = "Significance") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "top")

ggsave("results/figures/10_fold_change_distribution.png", p_fc_dist,
       width = 10, height = 7, dpi = 300)

cat("   ✓ Saved: results/figures/10_fold_change_distribution.png\n\n")

# ============================================================================
# 10. P-VALUE DISTRIBUTION
# ============================================================================

cat("📊 Creating p-value distribution plot...\n")

p_pval <- ggplot(res_df %>% filter(!is.na(pvalue)), 
                 aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "#3498db", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "P-value Distribution",
       subtitle = "Well-powered experiment shows enrichment near 0",
       x = "P-value",
       y = "Number of Genes") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("results/figures/11_pvalue_distribution.png", p_pval,
       width = 10, height = 7, dpi = 300)

cat("   ✓ Saved: results/figures/11_pvalue_distribution.png\n\n")

# ============================================================================
# 11. SAVE RESULTS TABLES
# ============================================================================

cat("💾 Saving results tables...\n")

# Full results
write.csv(res_df, "results/tables/DE_results_full.csv", row.names = FALSE)
cat("   ✓ Saved: results/tables/DE_results_full.csv\n")

# Significant genes only
sig_genes <- res_df %>% filter(significant == "Significant")
write.csv(sig_genes, "results/tables/DE_results_significant.csv", row.names = FALSE)
cat(sprintf("   ✓ Saved: results/tables/DE_results_significant.csv (%d genes)\n", 
            nrow(sig_genes)))

# Top 100 genes by adjusted p-value
top100 <- res_df %>% 
  filter(!is.na(padj)) %>%
  slice_min(order_by = padj, n = 100)
write.csv(top100, "results/tables/DE_results_top100.csv", row.names = FALSE)
cat("   ✓ Saved: results/tables/DE_results_top100.csv\n")

# Upregulated genes
up_genes <- res_df %>% filter(direction == "Upregulated")
write.csv(up_genes, "results/tables/DE_genes_upregulated.csv", row.names = FALSE)
cat(sprintf("   ✓ Saved: results/tables/DE_genes_upregulated.csv (%d genes)\n", 
            nrow(up_genes)))

# Downregulated genes
down_genes <- res_df %>% filter(direction == "Downregulated")
write.csv(down_genes, "results/tables/DE_genes_downregulated.csv", row.names = FALSE)
cat(sprintf("   ✓ Saved: results/tables/DE_genes_downregulated.csv (%d genes)\n\n", 
            nrow(down_genes)))

# ============================================================================
# 12. SUMMARY STATISTICS
# ============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("📈 DIFFERENTIAL EXPRESSION SUMMARY\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

cat(sprintf("Total genes analyzed:              %s\n", format(n_total, big.mark = ",")))
cat(sprintf("Significant DE genes:              %s (%.1f%%)\n", 
            format(n_sig, big.mark = ","),
            100 * n_sig / n_total))
cat(sprintf("   - Upregulated (Tumor > Normal):   %s\n", format(n_up, big.mark = ",")))
cat(sprintf("   - Downregulated (Tumor < Normal): %s\n\n", format(n_down, big.mark = ",")))

cat("Significance criteria:\n")
cat("   - Adjusted p-value < 0.05\n")
cat("   - |Log2 Fold Change| >= 1 (≥2-fold change)\n\n")

cat("Fold change range (significant genes):\n")
if (n_sig > 0) {
  fc_range <- sig_genes %>%
    summarise(
      min_lfc = min(log2FoldChange),
      max_lfc = max(log2FoldChange),
      min_fc = 2^min_lfc,
      max_fc = 2^max_lfc
    )
  
  cat(sprintf("   - Upregulation: up to %.1f-fold (log2FC = %.2f)\n", 
              fc_range$max_fc, fc_range$max_lfc))
  cat(sprintf("   - Downregulation: down to %.1f-fold (log2FC = %.2f)\n\n", 
              fc_range$min_fc, fc_range$min_lfc))
}

cat("Generated outputs:\n")
cat("   - 5 publication-quality figures\n")
cat("   - 5 results tables (CSV format)\n")
cat("   - Complete statistical results\n")

cat("\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("✅ DIFFERENTIAL EXPRESSION ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

cat("🚀 Next step: Pathway enrichment analysis (Script 04)\n\n")