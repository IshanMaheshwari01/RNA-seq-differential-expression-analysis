# Data

This folder is where `scripts/01_data_preparation.py` downloads the raw
recount3 files to (about 135 MB total) and writes the matched-pairs subset
it selects from them. The raw downloads and intermediate files are not
committed to this repository to keep it lightweight; run the script below to
fetch them.

```bash
pip install -r requirements.txt
python scripts/01_data_preparation.py
```

This pulls three files directly from the recount3 AWS Open Data mirror (no
authentication required):

- `https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/tcga/gene_sums/CA/BRCA/tcga.gene_sums.BRCA.G026.gz`
- `https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/tcga/metadata/CA/BRCA/tcga.tcga.BRCA.MD.gz`
- `https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/tcga/metadata/CA/BRCA/tcga.recount_qc.BRCA.MD.gz`

and writes `raw_counts_matched_subset.csv` and `sample_sheet.csv` into
`results/tables/`, which are the files the rest of the pipeline (steps
02 to 04) actually reads from. Those two output files are committed in
`results/tables/`, so you can skip this download step entirely and go
straight to `02_quality_control.py` if you just want to reproduce the
downstream analysis without re-pulling the raw data.
