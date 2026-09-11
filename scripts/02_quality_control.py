"""
02_quality_control.py

Basic QC on the matched TCGA-BRCA tumor/normal count subset produced by
01_data_preparation.py:
  - library size (total counts) per sample, split by condition
  - gene detection rate per sample (fraction of genes with count > 0)
  - low-count gene filtering (kept for DE step)
  - sample-sample correlation (on log-CPM) as a heatmap
  - PCA on log-CPM of filtered genes, colored by condition

Outputs:
  results/tables/qc_library_sizes.csv
  results/tables/filtered_counts.csv          (genes that pass the filter)
  results/figures/01_library_sizes.png
  results/figures/02_pca_tumor_vs_normal.png
  results/figures/03_sample_correlation_heatmap.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

BASE = Path(__file__).resolve().parents[1]
TABLES = BASE / "results" / "tables"
FIGURES = BASE / "results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="talk")

# minimum filtering: gene must have >= 10 total counts across all samples
# AND be detected (count > 0) in at least 20% of samples
MIN_TOTAL_COUNT = 10
MIN_SAMPLE_FRACTION = 0.2


def main():
    counts = pd.read_csv(TABLES / "raw_counts_matched_subset.csv", index_col=0)
    sample_sheet = pd.read_csv(TABLES / "sample_sheet.csv", index_col=0)
    sample_sheet = sample_sheet.loc[counts.columns]

    # ---- library sizes ----
    lib_sizes = counts.sum(axis=0)
    qc = pd.DataFrame(
        {
            "sample_id": counts.columns,
            "library_size": lib_sizes.values,
            "condition": sample_sheet["condition"].values,
            "genes_detected": (counts > 0).sum(axis=0).values,
            "detection_rate": ((counts > 0).sum(axis=0) / counts.shape[0]).values,
        }
    )
    qc.to_csv(TABLES / "qc_library_sizes.csv", index=False)
    print(qc.groupby("condition")[["library_size", "detection_rate"]].describe().T)

    # ---- library size barplot ----
    fig, ax = plt.subplots(figsize=(14, 6))
    order = qc.sort_values(["condition", "sample_id"])["sample_id"]
    colors = qc.set_index("sample_id").loc[order, "condition"].map(
        {"Tumor": "#d62728", "Normal": "#1f77b4"}
    )
    ax.bar(range(len(order)), qc.set_index("sample_id").loc[order, "library_size"] / 1e6, color=colors)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_ylabel("Library size (million reads)")
    ax.set_title("Library sizes per sample (TCGA-BRCA matched tumor/normal subset)")
    handles = [
        plt.Rectangle((0, 0), 1, 1, color="#d62728", label="Tumor"),
        plt.Rectangle((0, 0), 1, 1, color="#1f77b4", label="Normal"),
    ]
    ax.legend(handles=handles)
    fig.tight_layout()
    fig.savefig(FIGURES / "01_library_sizes.png", dpi=150)
    plt.close(fig)

    # ---- gene filtering ----
    total_counts = counts.sum(axis=1)
    frac_detected = (counts > 0).sum(axis=1) / counts.shape[1]
    keep = (total_counts >= MIN_TOTAL_COUNT) & (frac_detected >= MIN_SAMPLE_FRACTION)
    filtered = counts.loc[keep]
    print(f"Genes before filtering: {counts.shape[0]}")
    print(f"Genes after filtering (total count >= {MIN_TOTAL_COUNT}, detected in >= {MIN_SAMPLE_FRACTION*100:.0f}% samples): {filtered.shape[0]}")
    filtered.to_csv(TABLES / "filtered_counts.csv")

    # ---- log-CPM for PCA / correlation (QC visualization only; pydeseq2
    # will work from raw filtered counts directly in script 03) ----
    cpm = filtered.div(filtered.sum(axis=0), axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)

    # ---- PCA ----
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(log_cpm.T.values)
    pca_df = pd.DataFrame(coords, columns=["PC1", "PC2"], index=log_cpm.columns)
    pca_df["condition"] = sample_sheet.loc[pca_df.index, "condition"].values
    pca_df["patient_id"] = sample_sheet.loc[pca_df.index, "patient_id"].values
    pca_df.to_csv(TABLES / "pca_coordinates.csv")

    fig, ax = plt.subplots(figsize=(8, 7))
    sns.scatterplot(
        data=pca_df, x="PC1", y="PC2", hue="condition", style="condition",
        palette={"Tumor": "#d62728", "Normal": "#1f77b4"}, s=90, ax=ax
    )
    var_explained = pca.explained_variance_ratio_ * 100
    ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}% variance)")
    ax.set_title("PCA of log2(CPM+1)\nTCGA-BRCA matched tumor vs normal", fontsize=15)
    fig.tight_layout()
    fig.savefig(FIGURES / "02_pca_tumor_vs_normal.png", dpi=150)
    plt.close(fig)

    # ---- sample correlation heatmap ----
    corr = log_cpm.corr(method="pearson")
    condition_colors = sample_sheet.loc[corr.index, "condition"].map(
        {"Tumor": "#d62728", "Normal": "#1f77b4"}
    )
    g = sns.clustermap(
        corr, cmap="viridis", figsize=(12, 12),
        row_colors=condition_colors, col_colors=condition_colors,
        xticklabels=True, yticklabels=True,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=5)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=5)
    g.figure.suptitle("Sample-sample Pearson correlation (log2 CPM)", y=1.02)
    g.savefig(FIGURES / "03_sample_correlation_heatmap.png", dpi=150)
    plt.close("all")

    print("Saved QC tables and figures.")


if __name__ == "__main__":
    main()
