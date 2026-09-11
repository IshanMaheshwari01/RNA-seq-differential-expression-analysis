"""
01_data_preparation.py

Downloads (or reuses a local cache of) the real TCGA-BRCA gene-level RNA-seq
count matrix from the recount3 project, identifies patient-matched
tumor/normal sample pairs from the TCGA barcode embedded in the recount3
TCGA metadata table, subsamples to a fixed set of pairs for compute
tractability, and writes a clean count matrix + sample sheet to disk for the
downstream QC / DE / enrichment scripts.

Data source
-----------
recount3 (human, TCGA data source, BRCA project), gene-level counts,
Gencode v26 annotation ("G026"), STAR + featureCounts summarized by the
Monorail pipeline used to build recount3. This is NOT simulated data -- it
is the same processed TCGA RNA-seq count data distributed by the recount3
project (https://rna.recount.bio), mirrored on the AWS Open Data registry.

Files pulled (all real, publicly hosted, no login required):
  - Gene counts:
    https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/tcga/gene_sums/CA/BRCA/tcga.gene_sums.BRCA.G026.gz
  - TCGA sample metadata (barcode, sample type, clinical fields):
    https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/tcga/metadata/CA/BRCA/tcga.tcga.BRCA.MD.gz

Why recount3 instead of the Xena GDC hub file named in the original task
notes: the literal Xena S3 URL
(gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.htseq_counts.tsv.gz)
returned HTTP 403 (Access Denied) from this sandbox, and the
gdc.xenahubs.net redirector 302-redirects to the same restricted bucket.
recount3's AWS Open Data mirror was reachable and returned HTTP 200, so it
was used instead, exactly as the task instructions anticipated ("use
whichever source actually works").

Sample type codes (from the GDC barcode spec, 4th barcode field):
  01 = Primary Tumor, 11 = Solid Tissue Normal.
  https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

Units of the count matrix -- important
---------------------------------------
recount3 "gene_sums" files are NOT raw read counts. They are base-pair
coverage sums (the area under the per-base coverage curve, summed over each
gene's exons), inherited from the recount2 processing model. This is
documented in the recount3 Bioconductor package (`transform_counts()`,
`compute_read_counts()`) -- see
https://rna.recount.bio/docs/bioconductor.html and
https://bioconductor.org/packages/release/bioc/manuals/recount3/man/recount3.pdf.

Per the recount3 R source (`compute_read_counts()` in
https://github.com/LieberInstitute/recount3/blob/devel/R/transform_counts.R),
the standard conversion to approximate read counts is:

    read_counts = round( gene_sums / recount_qc.star.average_mapped_length )

i.e. dividing the coverage sum by that sample's average mapped read length
(pulled from the `recount_qc` metadata table). That conversion is applied
below before the counts are handed to pydeseq2, so the "raw_counts_matched_subset.csv"
produced by this script already contains DESeq2-ready approximate integer
read counts, not the unconverted coverage sums.
"""

import gzip
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RANDOM_SEED = 42
N_PAIRS = 25  # subsample size (task suggests 20-30 matched pairs)

BASE = Path(__file__).resolve().parents[1]
DATA_DIR = BASE / "data"
RESULTS_TABLES = BASE / "results" / "tables"
RESULTS_TABLES.mkdir(parents=True, exist_ok=True)

COUNTS_GZ = DATA_DIR / "TCGA_BRCA_gene_sums.gz"
META_TSV = DATA_DIR / "tcga_brca_tcga_meta"  # already gunzipped, tab-sep
QC_TSV = DATA_DIR / "tcga_recount_qc"  # already gunzipped, tab-sep

COUNTS_URL = (
    "https://recount-opendata.s3.amazonaws.com/recount3/release/human/"
    "data_sources/tcga/gene_sums/CA/BRCA/tcga.gene_sums.BRCA.G026.gz"
)
META_URL = (
    "https://recount-opendata.s3.amazonaws.com/recount3/release/human/"
    "data_sources/tcga/metadata/CA/BRCA/tcga.tcga.BRCA.MD.gz"
)
QC_URL = (
    "https://recount-opendata.s3.amazonaws.com/recount3/release/human/"
    "data_sources/tcga/metadata/CA/BRCA/tcga.recount_qc.BRCA.MD.gz"
)


def _download(url: str, dest: Path) -> None:
    import urllib.request

    print(f"Downloading {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)


def ensure_raw_files() -> None:
    """Make sure the raw recount3 files are present locally, downloading
    them from the AWS Open Data mirror if not already cached."""
    if not COUNTS_GZ.exists():
        _download(COUNTS_URL, COUNTS_GZ)
    if not META_TSV.exists():
        meta_gz = DATA_DIR / "tcga_brca_tcga_meta.gz"
        if not meta_gz.exists():
            _download(META_URL, meta_gz)
        with gzip.open(meta_gz, "rt") as fin, open(META_TSV, "w") as fout:
            fout.write(fin.read())
    if not QC_TSV.exists():
        qc_gz = DATA_DIR / "tcga_recount_qc.gz"
        if not qc_gz.exists():
            _download(QC_URL, qc_gz)
        with gzip.open(qc_gz, "rt") as fin, open(QC_TSV, "w") as fout:
            fout.write(fin.read())


def load_metadata() -> pd.DataFrame:
    meta = pd.read_csv(META_TSV, sep="\t", low_memory=False)
    meta["patient_id"] = meta["tcga_barcode"].str[:12]
    # 4th field of the barcode, first 2 chars = sample type code
    meta["sample_type_code"] = meta["tcga_barcode"].str.split("-").str[3].str[:2]
    return meta


def load_qc() -> pd.DataFrame:
    qc = pd.read_csv(QC_TSV, sep="\t", low_memory=False)
    return qc[["external_id", "star.average_mapped_length"]]


def select_matched_pairs(meta: pd.DataFrame, n_pairs: int, seed: int) -> pd.DataFrame:
    tumor = meta[meta["sample_type_code"] == "01"].copy()
    normal = meta[meta["sample_type_code"] == "11"].copy()

    # keep exactly one tumor and one normal sample per patient
    # (a handful of patients have >1 sample of a given type -- take the
    # first occurrence, sorted by external_id for determinism)
    tumor = tumor.sort_values("external_id").drop_duplicates("patient_id", keep="first")
    normal = normal.sort_values("external_id").drop_duplicates("patient_id", keep="first")

    paired_patients = sorted(set(tumor["patient_id"]) & set(normal["patient_id"]))
    print(f"Total patients with both a -01 tumor and -11 normal sample: {len(paired_patients)}")

    rng = random.Random(seed)
    chosen = sorted(rng.sample(paired_patients, min(n_pairs, len(paired_patients))))
    print(f"Randomly selected {len(chosen)} matched pairs (seed={seed}).")

    tumor_sel = tumor[tumor["patient_id"].isin(chosen)].copy()
    normal_sel = normal[normal["patient_id"].isin(chosen)].copy()
    tumor_sel["condition"] = "Tumor"
    normal_sel["condition"] = "Normal"

    sample_sheet = pd.concat([tumor_sel, normal_sel], ignore_index=True)
    sample_sheet = sample_sheet.sort_values(["patient_id", "condition"]).reset_index(drop=True)
    return sample_sheet[
        [
            "external_id",
            "rail_id",
            "tcga_barcode",
            "patient_id",
            "condition",
            "sample_type_code",
            "cgc_sample_sample_type",
        ]
    ]


def load_counts_for_samples(external_ids: list[str]) -> pd.DataFrame:
    """Stream the (large, ~1256-sample) gene count matrix and keep only the
    columns for the selected samples, to avoid loading the whole 315MB
    uncompressed matrix into memory unnecessarily."""
    with gzip.open(COUNTS_GZ, "rt") as f:
        f.readline()  # ##annotation=...
        f.readline()  # ##date.generated=...
        header = f.readline().rstrip("\n").split("\t")

    wanted = set(external_ids)
    keep_idx = [0] + [i for i, c in enumerate(header) if c in wanted]
    keep_cols = [header[i] for i in keep_idx]
    missing = wanted - set(keep_cols)
    if missing:
        raise ValueError(f"{len(missing)} requested sample IDs not found in count matrix header: {list(missing)[:5]}")

    usecols_idx = keep_idx
    df = pd.read_csv(
        COUNTS_GZ,
        sep="\t",
        skiprows=3,
        header=None,
        names=header,
        usecols=usecols_idx,
        index_col=0,
        compression="gzip",
    )
    df.index.name = "gene_id"
    return df


def main():
    ensure_raw_files()
    meta = load_metadata()
    sample_sheet = select_matched_pairs(meta, N_PAIRS, RANDOM_SEED)

    print("Loading counts for selected samples (streaming column subset)...")
    counts = load_counts_for_samples(sample_sheet["external_id"].tolist())
    counts = counts[sample_sheet["external_id"].tolist()]  # enforce order

    # --- convert recount3 base-pair coverage sums to approximate read
    # counts, following the official recount3 R package formula:
    #   read_counts = round(gene_sums / recount_qc.star.average_mapped_length)
    # https://github.com/LieberInstitute/recount3/blob/devel/R/transform_counts.R
    qc = load_qc().set_index("external_id")
    read_lengths = qc.loc[counts.columns, "star.average_mapped_length"]
    print("Converting coverage sums to approximate read counts using average mapped read length per sample...")
    counts = counts.div(read_lengths.values, axis=1).round().astype(int)

    # rename columns to short, readable sample IDs: <patient>_<condition>
    id_map = dict(zip(sample_sheet["external_id"], sample_sheet["patient_id"] + "_" + sample_sheet["condition"]))
    counts.columns = [id_map[c] for c in counts.columns]
    sample_sheet["sample_id"] = sample_sheet["patient_id"] + "_" + sample_sheet["condition"]
    sample_sheet = sample_sheet.set_index("sample_id")

    print("Final count matrix shape (genes x samples):", counts.shape)
    print("Sample sheet:")
    print(sample_sheet["condition"].value_counts())

    counts.to_csv(RESULTS_TABLES / "raw_counts_matched_subset.csv")
    sample_sheet.to_csv(RESULTS_TABLES / "sample_sheet.csv")

    print(f"Saved raw_counts_matched_subset.csv ({counts.shape[0]} genes x {counts.shape[1]} samples)")
    print(f"Saved sample_sheet.csv ({sample_sheet.shape[0]} samples)")


if __name__ == "__main__":
    main()
