"""
04_pathway_enrichment.py

Takes the significant up/down gene sets from 03_differential_expression.py
and runs GO Biological Process enrichment using gseapy against the Enrichr
GO_Biological_Process_2023 library (fetched live from the Enrichr web API,
https://maayanlab.cloud/Enrichr/).

Gene IDs in the count matrix are Ensembl gene IDs (recount3 / Gencode v26
annotation), but Enrichr gene sets are indexed by gene symbol, so IDs are
mapped to symbols first via the MyGene.info REST API
(https://mygene.info/) before submission.

Outputs:
  results/tables/enrichment_up_GO_BP.csv
  results/tables/enrichment_down_GO_BP.csv
  results/figures/07_go_enrichment_up.png
  results/figures/08_go_enrichment_down.png
"""

import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import requests
import seaborn as sns
import gseapy as gp

BASE = Path(__file__).resolve().parents[1]
TABLES = BASE / "results" / "tables"
FIGURES = BASE / "results" / "figures"

GO_LIBRARY = "GO_Biological_Process_2023"
MYGENE_URL = "https://mygene.info/v3/gene"

sns.set_theme(style="whitegrid", context="talk")


def ensembl_to_symbol(ensembl_ids: list[str]) -> dict[str, str]:
    """Batch-map Ensembl gene IDs (versioned, e.g. ENSG00000141510.16) to
    HGNC gene symbols using the MyGene.info REST API."""
    bare_ids = [g.split(".")[0] for g in ensembl_ids]
    mapping = {}
    batch_size = 500
    for i in range(0, len(bare_ids), batch_size):
        batch = bare_ids[i : i + batch_size]
        resp = requests.post(
            MYGENE_URL,
            data={"ids": ",".join(batch), "fields": "symbol", "species": "human"},
            timeout=60,
        )
        resp.raise_for_status()
        for rec in resp.json():
            if "symbol" in rec and "query" in rec:
                mapping[rec["query"]] = rec["symbol"]
        time.sleep(0.5)
    # map back to the original (versioned) ids
    full_mapping = {}
    for orig, bare in zip(ensembl_ids, bare_ids):
        if bare in mapping:
            full_mapping[orig] = mapping[bare]
    return full_mapping


def run_enrichment(gene_symbols: list[str], label: str) -> pd.DataFrame:
    if len(gene_symbols) == 0:
        print(f"No genes to test for {label} -- skipping.")
        return pd.DataFrame()
    enr = gp.enrichr(
        gene_list=gene_symbols,
        gene_sets=[GO_LIBRARY],
        organism="human",
        outdir=None,
        cutoff=1.0,  # keep all returned terms, we filter afterwards
    )
    res = enr.results.sort_values("Adjusted P-value")
    return res


def plot_top_terms(res: pd.DataFrame, title: str, outfile: Path, color: str, top_n: int = 15):
    if res.empty:
        print(f"Skipping plot for '{title}' -- no enrichment results.")
        return
    import numpy as np

    top = res.head(top_n).copy()
    top["neg_log10_padj"] = -np.log10(top["Adjusted P-value"].clip(lower=1e-300))
    top = top.sort_values("neg_log10_padj", ascending=True)

    max_label_len = max(len(t) for t in top["Term"])
    fig_width = max(13, 8 + 0.09 * max_label_len)
    fig, ax = plt.subplots(figsize=(fig_width, max(5, 0.45 * len(top))))
    ax.barh(top["Term"], top["neg_log10_padj"], color=color)
    ax.set_xlabel("-log10(adjusted p-value)")
    ax.set_title(title)
    ax.tick_params(axis="y", labelsize=9)
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    sig = pd.read_csv(TABLES / "DE_significant.csv", index_col=0)
    up_ids = sig[sig["log2FoldChange"] > 0].index.tolist()
    down_ids = sig[sig["log2FoldChange"] < 0].index.tolist()
    print(f"Up-regulated genes (tumor): {len(up_ids)}; down-regulated genes (tumor): {len(down_ids)}")

    all_ids = list(set(up_ids) | set(down_ids))
    print(f"Mapping {len(all_ids)} Ensembl gene IDs to symbols via MyGene.info...")
    id2sym = ensembl_to_symbol(all_ids)
    print(f"Mapped {len(id2sym)}/{len(all_ids)} IDs to symbols.")

    up_symbols = sorted({id2sym[g] for g in up_ids if g in id2sym})
    down_symbols = sorted({id2sym[g] for g in down_ids if g in id2sym})
    print(f"Up symbols: {len(up_symbols)}; Down symbols: {len(down_symbols)}")

    pd.Series(up_symbols, name="symbol").to_csv(TABLES / "up_gene_symbols.csv", index=False)
    pd.Series(down_symbols, name="symbol").to_csv(TABLES / "down_gene_symbols.csv", index=False)

    print(f"Running Enrichr GO Biological Process enrichment ({GO_LIBRARY}) on up-regulated genes...")
    res_up = run_enrichment(up_symbols, "up")
    res_up.to_csv(TABLES / "enrichment_up_GO_BP.csv", index=False)
    print(f"Up-regulated: {res_up.shape[0]} GO BP terms returned; "
          f"{(res_up['Adjusted P-value'] < 0.05).sum() if not res_up.empty else 0} significant at padj<0.05")

    print(f"Running Enrichr GO Biological Process enrichment ({GO_LIBRARY}) on down-regulated genes...")
    res_down = run_enrichment(down_symbols, "down")
    res_down.to_csv(TABLES / "enrichment_down_GO_BP.csv", index=False)
    print(f"Down-regulated: {res_down.shape[0]} GO BP terms returned; "
          f"{(res_down['Adjusted P-value'] < 0.05).sum() if not res_down.empty else 0} significant at padj<0.05")

    plot_top_terms(
        res_up, "Top GO Biological Process terms – genes UP in tumor",
        FIGURES / "07_go_enrichment_up.png", color="#d62728",
    )
    plot_top_terms(
        res_down, "Top GO Biological Process terms – genes DOWN in tumor",
        FIGURES / "08_go_enrichment_down.png", color="#1f77b4",
    )

    print("Pathway enrichment complete.")


if __name__ == "__main__":
    main()
