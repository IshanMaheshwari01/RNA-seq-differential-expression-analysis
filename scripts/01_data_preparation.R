# scripts/01_data_preparation.R
# ---------------------------------------------------------------------------
# Create a SIMULATED but biologically grounded RNA-seq dataset for a
# tumour-vs-normal differential-expression demonstration.
#
# WHY THIS DESIGN:
#   The point of this project is to show a correct, reproducible bulk RNA-seq
#   pipeline (QC -> DESeq2 -> enrichment). To make every step honest, the data
#   is simulated with a KNOWN GROUND TRUTH planted in REAL, annotated human
#   genes:
#       * genes annotated to proliferation / cell-cycle GO terms are made
#         UP-regulated in tumour
#       * genes annotated to immune-response GO terms are made DOWN-regulated
#         in tumour
#   Because the signal lives in real genes, the downstream GO enrichment
#   (script 04) recovers "proliferation up / immune down" FOR REAL, computed
#   from the data — nothing is hard-coded. This is exactly how a pipeline is
#   validated against a known truth.
# ---------------------------------------------------------------------------

library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

set.seed(123)  # reproducibility

cat("Creating biologically grounded simulated RNA-seq dataset...\n\n")

# ---------------------------------------------------------------------------
# 1. Build the gene universe from REAL human gene symbols
# ---------------------------------------------------------------------------
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
all_symbols <- unique(all_symbols)

n_genes <- 20000
universe <- sample(all_symbols, min(n_genes, length(all_symbols)))

# ---------------------------------------------------------------------------
# 2. Pull REAL gene sets to plant the signal in
#    (genes annotated to these GO biological-process terms)
# ---------------------------------------------------------------------------
genes_for_go <- function(go_id) {
  s <- AnnotationDbi::select(org.Hs.eg.db, keys = go_id,
                             keytype = "GOALL", columns = "SYMBOL")$SYMBOL
  unique(na.omit(s))
}

proliferation <- unique(c(
  genes_for_go("GO:0007049"),  # cell cycle
  genes_for_go("GO:0051301"),  # cell division
  genes_for_go("GO:0006260")   # DNA replication
))
immune <- unique(c(
  genes_for_go("GO:0006955"),  # immune response
  genes_for_go("GO:0002376")   # immune system process
))

# keep only genes that are in our universe, and make the two sets disjoint
proliferation <- intersect(proliferation, universe)
immune        <- setdiff(intersect(immune, universe), proliferation)

# plant a manageable, realistic number of DE genes
n_up   <- min(1500, length(proliferation))
n_down <- min(1500, length(immune))
up_genes   <- sample(proliferation, n_up)
down_genes <- sample(immune, n_down)

cat(sprintf("Universe: %d real genes\n", length(universe)))
cat(sprintf("Planted UP (proliferation):   %d genes\n", n_up))
cat(sprintf("Planted DOWN (immune):        %d genes\n", n_down))
cat(sprintf("Null (no change):             %d genes\n\n",
            length(universe) - n_up - n_down))

# ---------------------------------------------------------------------------
# 3. Simulate the count matrix (negative binomial — correct model for counts)
# ---------------------------------------------------------------------------
n_samples <- 12
sample_names <- c(paste0("Tumor_", 1:6), paste0("Normal_", 1:6))

base_expression <- rnbinom(length(universe), mu = 100, size = 2)
base_expression <- pmax(base_expression, 5)
names(base_expression) <- universe

counts <- matrix(0L, nrow = length(universe), ncol = n_samples,
                 dimnames = list(universe, sample_names))

up_set   <- up_genes
down_set <- down_genes
tumor_cols  <- 1:6
normal_cols <- 7:12

for (g in universe) {
  mu <- base_expression[g]
  if (g %in% up_set) {
    fc <- runif(1, 3, 8)
    counts[g, tumor_cols]  <- rnbinom(6, mu = mu * fc, size = 2)
    counts[g, normal_cols] <- rnbinom(6, mu = mu,      size = 2)
  } else if (g %in% down_set) {
    fc <- runif(1, 3, 8)
    counts[g, tumor_cols]  <- rnbinom(6, mu = mu / fc, size = 2)
    counts[g, normal_cols] <- rnbinom(6, mu = mu,      size = 2)
  } else {
    counts[g, ] <- rnbinom(n_samples, mu = mu, size = 2)
  }
}

cat("Count matrix created\n\n")

# ---------------------------------------------------------------------------
# 4. Sample metadata
# ---------------------------------------------------------------------------
metadata <- data.frame(
  sample_id  = sample_names,
  condition  = factor(c(rep("Tumor", 6), rep("Normal", 6)),
                      levels = c("Normal", "Tumor")),   # Normal = reference
  batch      = factor(rep(c("A", "B", "A", "B", "A", "B"), 2)),
  patient_id = factor(c(1:6, 1:6)),
  row.names  = sample_names
)
metadata$age    <- c(65, 58, 72, 55, 68, 61, 64, 59, 70, 56, 67, 62)
metadata$gender <- factor(c("M","F","M","M","F","F","M","F","M","M","F","F"))
metadata$stage  <- factor(c("III","II","IV","II","III","II", NA,NA,NA,NA,NA,NA))

# ---------------------------------------------------------------------------
# 5. Save (also save the ground-truth lists so the pipeline can be validated)
# ---------------------------------------------------------------------------
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write.csv(counts,   "data/raw/count_matrix.csv",   row.names = TRUE)
write.csv(metadata, "data/raw/sample_metadata.csv", row.names = TRUE)
saveRDS(counts,   "data/raw/count_matrix.rds")
saveRDS(metadata, "data/raw/sample_metadata.rds")
saveRDS(list(up = up_genes, down = down_genes),
        "data/raw/ground_truth_genes.rds")

cat("DATASET SUMMARY\n==================\n")
cat(sprintf("  Dimensions: %d genes x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("  Total counts: %s\n", format(sum(counts), big.mark = ",")))
cat("\nFiles saved to data/raw/ (count_matrix, sample_metadata, ground_truth_genes)\n")
cat("\nData preparation complete!\n")

