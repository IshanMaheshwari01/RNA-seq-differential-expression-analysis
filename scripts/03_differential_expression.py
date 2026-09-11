"""
03_differential_expression.py

Runs the differential expression analysis with pydeseq2 (a Python
reimplementation of DESeq2's negative-binomial GLM + Wald test, Muzellec
et al. 2023, Bioinformatics) on the filtered TCGA-BRCA matched tumor/normal
count matrix produced by 02_quality_control.py.

Design: paired by patient. Each patient contributed one tumor and one
normal sample (see 01_data_preparation.py), so the model is

    ~ patient_id + condition

which removes patient-to-patient variation before testing for the tumor vs.
normal effect, similar to a paired t-test / blocked ANOVA design. This is
the design pydeseq2's documentation recommends for matched/paired samples.

Outputs:
  results/tables/DE_full_results.csv        (all genes, DESeq2 stats)
  results/tables/DE_significant.csv         (padj < 0.05 & |log2FC| >= 1)
  results/tables/DE_top20_upregulated.csv
  results/tables/DE_top20_downregulated.csv
  results/figures/04_volcano_plot.png
  results/figures/05_ma_plot.png
  results/figures/06_top50_heatmap.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BASE = Path(__file__).resolve().parents[1]
TABLES = BASE / "results" / "tables"
FIGURES = BASE / "results" / "figures"

PADJ_THRESHOLD = 0.05
LFC_THRESHOLD = 1.0

sns.set_theme(style="whitegrid", context="talk")


def main():
    counts = pd.read_csv(TABLES / "filtered_counts.csv", index_col=0).T  # samples x genes for pydeseq2
    sample_sheet = pd.read_csv(TABLES / "sample_sheet.csv", index_col=0)
    sample_sheet = sample_sheet.loc[counts.index]

    metadata = sample_sheet[["condition", "patient_id"]].copy()
    metadata["condition"] = pd.Categorical(metadata["condition"], categories=["Normal", "Tumor"])
    metadata["patient_id"] = metadata["patient_id"].astype(str)

    print("Sample metadata for design ~ patient_id + condition:")
    print(metadata["condition"].value_counts())
    print(f"{metadata['patient_id'].nunique()} unique patients (paired design)")

    # pydeseq2 requires integer counts
    counts = counts.round().astype(int)

    print("\nFitting pydeseq2 negative-binomial GLM (design: ~patient_id + condition)...")
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design="~patient_id + condition",
        refit_cooks=True,
    )
    dds.deseq2()

    print("Running Wald test for condition Tumor vs Normal...")
    stat_res = DeseqStats(dds, contrast=["condition", "Tumor", "Normal"], alpha=PADJ_THRESHOLD)
    stat_res.summary()

    results = stat_res.results_df.copy()
    results = results.sort_values("padj")
    results.index.name = "gene_id"
    results.to_csv(TABLES / "DE_full_results.csv")
    print(f"Saved DE_full_results.csv ({results.shape[0]} genes)")

    sig = results[(results["padj"] < PADJ_THRESHOLD) & (results["log2FoldChange"].abs() >= LFC_THRESHOLD)].copy()
    sig = sig.sort_values("padj")
    sig.to_csv(TABLES / "DE_significant.csv")
    print(f"Saved DE_significant.csv ({sig.shape[0]} genes with padj<{PADJ_THRESHOLD} & |log2FC|>={LFC_THRESHOLD})")

    up = sig[sig["log2FoldChange"] > 0].sort_values("log2FoldChange", ascending=False)
    down = sig[sig["log2FoldChange"] < 0].sort_values("log2FoldChange")
    up.head(20).to_csv(TABLES / "DE_top20_upregulated.csv")
    down.head(20).to_csv(TABLES / "DE_top20_downregulated.csv")
    print(f"Upregulated in tumor: {up.shape[0]} genes; downregulated in tumor: {down.shape[0]} genes")

    # ---- volcano plot ----
    plot_df = results.dropna(subset=["padj", "log2FoldChange"]).copy()
    plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"].replace(0, 1e-300))
    plot_df["category"] = "Not significant"
    plot_df.loc[(plot_df["padj"] < PADJ_THRESHOLD) & (plot_df["log2FoldChange"] >= LFC_THRESHOLD), "category"] = "Up in tumor"
    plot_df.loc[(plot_df["padj"] < PADJ_THRESHOLD) & (plot_df["log2FoldChange"] <= -LFC_THRESHOLD), "category"] = "Down in tumor"

    fig, ax = plt.subplots(figsize=(9, 8))
    palette = {"Not significant": "#bbbbbb", "Up in tumor": "#d62728", "Down in tumor": "#1f77b4"}
    for cat, color in palette.items():
        sub = plot_df[plot_df["category"] == cat]
        ax.scatter(sub["log2FoldChange"], sub["neg_log10_padj"], s=8, alpha=0.6, color=color, label=f"{cat} (n={len(sub)})")
    ax.axvline(LFC_THRESHOLD, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-LFC_THRESHOLD, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(PADJ_THRESHOLD), color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 fold change (Tumor vs Normal)")
    ax.set_ylabel("-log10(padj)")
    ax.set_title("Volcano plot – TCGA-BRCA matched tumor vs normal")
    ax.legend(loc="upper left", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGURES / "04_volcano_plot.png", dpi=150)
    plt.close(fig)

    # ---- MA plot ----
    ma_df = results.dropna(subset=["baseMean", "log2FoldChange", "padj"]).copy()
    ma_df["log10_baseMean"] = np.log10(ma_df["baseMean"] + 1)
    ma_df["significant"] = ma_df["padj"] < PADJ_THRESHOLD

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(
        ma_df.loc[~ma_df["significant"], "log10_baseMean"],
        ma_df.loc[~ma_df["significant"], "log2FoldChange"],
        s=6, alpha=0.4, color="#bbbbbb", label="Not significant",
    )
    ax.scatter(
        ma_df.loc[ma_df["significant"], "log10_baseMean"],
        ma_df.loc[ma_df["significant"], "log2FoldChange"],
        s=8, alpha=0.7, color="#d62728", label=f"padj < {PADJ_THRESHOLD}",
    )
    ax.axhline(0, color="black", linewidth=1)
    ax.set_xlabel("log10(mean normalized count + 1)")
    ax.set_ylabel("log2 fold change (Tumor vs Normal)")
    ax.set_title("MA plot – TCGA-BRCA matched tumor vs normal")
    ax.legend(loc="upper right", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGURES / "05_ma_plot.png", dpi=150)
    plt.close(fig)

    # ---- heatmap of top 50 DE genes (by padj), z-scored across samples ----
    top50 = sig.head(50).index.tolist()
    if len(top50) > 0:
        norm_counts = pd.DataFrame(
            dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names
        )
        mat = norm_counts[top50].T
        mat_log = np.log2(mat + 1)
        z = mat_log.sub(mat_log.mean(axis=1), axis=0).div(mat_log.std(axis=1), axis=0)

        col_colors = sample_sheet.loc[z.columns, "condition"].map({"Tumor": "#d62728", "Normal": "#1f77b4"})
        g = sns.clustermap(
            z, cmap="vlag", center=0, figsize=(14, 14),
            col_colors=col_colors, xticklabels=True, yticklabels=True,
            col_cluster=True, row_cluster=True,
        )
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=5)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)
        g.figure.suptitle(f"Top {len(top50)} DE genes by padj (z-scored log2 normalized counts)", y=1.02)
        g.savefig(FIGURES / "06_top50_heatmap.png", dpi=150)
        plt.close("all")
        print(f"Saved heatmap for {len(top50)} top DE genes.")
    else:
        print("No significant genes found -- skipping top-50 heatmap.")

    print("Differential expression analysis complete.")


if __name__ == "__main__":
    main()
