# scripts/04_pathway_enrichment.R
# ---------------------------------------------------------------------------
# REAL GO over-representation analysis of the differentially expressed genes,
# using clusterProfiler against the org.Hs.eg.db annotation.
#
# Everything here is computed from the DE gene lists produced by script 03 —
# nothing is hard-coded. Because the data was built with proliferation genes
# up and immune genes down (script 01), this analysis should recover exactly
# that, which validates the whole pipeline against the known ground truth.
# ---------------------------------------------------------------------------

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

cat("PATHWAY ENRICHMENT ANALYSIS (clusterProfiler / GO BP)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load DE results
# ---------------------------------------------------------------------------
de_results <- read.csv("results/tables/DE_results_full.csv")
up_genes   <- read.csv("results/tables/DE_genes_upregulated.csv")$gene
down_genes <- read.csv("results/tables/DE_genes_downregulated.csv")$gene

# the statistical background ("universe") = every gene that was actually tested
universe <- de_results %>% filter(!is.na(padj)) %>% pull(gene)

cat(sprintf("Background (tested) genes: %s\n", format(length(universe), big.mark = ",")))
cat(sprintf("Upregulated genes:   %s\n", format(length(up_genes),   big.mark = ",")))
cat(sprintf("Downregulated genes: %s\n\n", format(length(down_genes), big.mark = ",")))

# ---------------------------------------------------------------------------
# 2. Run enrichGO (real over-representation test)
#    keyType = "SYMBOL" lets us pass gene symbols directly.
# ---------------------------------------------------------------------------
run_go <- function(genes, label) {
  cat(sprintf("Running GO BP enrichment for %s genes...\n", label))
  ego <- enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = FALSE
  )
  # collapse redundant parent/child terms for a cleaner figure
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego <- clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  }
  ego
}

ego_up   <- run_go(up_genes,   "upregulated")
ego_down <- run_go(down_genes, "downregulated")

cat(sprintf("   Enriched BP terms (up):   %d\n",   nrow(as.data.frame(ego_up))))
cat(sprintf("   Enriched BP terms (down): %d\n\n", nrow(as.data.frame(ego_down))))

# ---------------------------------------------------------------------------
# 3. Standard clusterProfiler plots
# ---------------------------------------------------------------------------
if (nrow(as.data.frame(ego_up)) > 0) {
  ggsave("results/figures/12_GO_upregulated_barplot.png",
         barplot(ego_up, showCategory = 12) +
           ggtitle("GO Enrichment: Upregulated Genes in Tumor"),
         width = 11, height = 7, dpi = 300)
  ggsave("results/figures/13_GO_upregulated_dotplot.png",
         dotplot(ego_up, showCategory = 12) +
           ggtitle("GO Enrichment: Upregulated Genes"),
         width = 11, height = 7, dpi = 300)
}
if (nrow(as.data.frame(ego_down)) > 0) {
  ggsave("results/figures/14_GO_downregulated_barplot.png",
         barplot(ego_down, showCategory = 12) +
           ggtitle("GO Enrichment: Downregulated Genes in Tumor"),
         width = 11, height = 7, dpi = 300)
  ggsave("results/figures/15_GO_downregulated_dotplot.png",
         dotplot(ego_down, showCategory = 12) +
           ggtitle("GO Enrichment: Downregulated Genes"),
         width = 11, height = 7, dpi = 300)
}

# ---------------------------------------------------------------------------
# 4. Combined comparison figure (top 6 each, real adjusted p-values)
# ---------------------------------------------------------------------------
top_terms <- function(ego, direction, n = 6) {
  df <- as.data.frame(ego)
  if (nrow(df) == 0) return(NULL)
  df %>% arrange(p.adjust) %>% slice(1:min(n, nrow(df))) %>%
    transmute(Description, p.adjust, Count, Direction = direction)
}

go_combined <- bind_rows(
  top_terms(ego_up,   "Upregulated"),
  top_terms(ego_down, "Downregulated")
)

if (!is.null(go_combined) && nrow(go_combined) > 0) {
  p_combined <- ggplot(go_combined,
                       aes(x = reorder(Description, -log10(p.adjust)),
                           y = -log10(p.adjust), fill = Direction)) +
    geom_col(color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Upregulated" = "#e74c3c",
                                 "Downregulated" = "#3498db")) +
    coord_flip() +
    facet_wrap(~Direction, scales = "free_y", ncol = 1) +
    labs(title = "Top Enriched Pathways: Tumor vs Normal",
         subtitle = "Upregulated (red) vs Downregulated (blue) — real GO BP over-representation",
         x = "GO Term", y = "-log10(Adjusted P-value)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "none",
          strip.background = element_rect(fill = "gray90"),
          strip.text = element_text(face = "bold"))

  ggsave("results/figures/16_GO_combined_comparison.png", p_combined,
         width = 12, height = 10, dpi = 300)
  cat("Saved: results/figures/16_GO_combined_comparison.png\n")
}

# ---------------------------------------------------------------------------
# 5. Save enrichment tables
# ---------------------------------------------------------------------------
if (nrow(as.data.frame(ego_up)) > 0)
  write.csv(as.data.frame(ego_up),   "results/tables/GO_enrichment_upregulated.csv",   row.names = FALSE)
if (nrow(as.data.frame(ego_down)) > 0)
  write.csv(as.data.frame(ego_down), "results/tables/GO_enrichment_downregulated.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Validate against the planted ground truth (optional sanity check)
# ---------------------------------------------------------------------------
if (file.exists("data/raw/ground_truth_genes.rds")) {
  truth <- readRDS("data/raw/ground_truth_genes.rds")
  up_recovered   <- length(intersect(up_genes,   truth$up))   / length(truth$up)
  down_recovered <- length(intersect(down_genes, truth$down)) / length(truth$down)
  cat(sprintf("\nGround-truth recovery: %.0f%% of planted UP, %.0f%% of planted DOWN genes called DE\n",
              100 * up_recovered, 100 * down_recovered))
}

cat("\nPATHWAY ENRICHMENT COMPLETE!\n")
