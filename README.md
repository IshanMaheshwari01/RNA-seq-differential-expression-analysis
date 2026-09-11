# TCGA-BRCA Tumor vs. Normal RNA-seq Differential Expression

Differential expression and pathway enrichment analysis on real, publicly
available TCGA-BRCA RNA-seq data, comparing matched tumor and normal breast
tissue from the same patients. Every number in this README, from the sample
counts to the enrichment p-values, comes from actually running the pipeline
in `scripts/` end to end.

**Data**: TCGA-BRCA (breast invasive carcinoma) bulk RNA-seq gene counts,
pulled from the [recount3](https://rna.recount.bio) project's AWS Open Data
mirror, not simulated. Sample selection: 25 patients who each had both a
primary tumor and a matched solid-tissue-normal sample sequenced (50 samples
total), picked at random (seed 42) out of the 113 such patients available in
TCGA-BRCA.

## Why recount3 instead of the Xena hub file

The task brief suggested pulling `TCGA-BRCA.htseq_counts.tsv.gz` straight
from the UCSC Xena GDC hub S3 bucket
(`gdc-hub.s3.us-east-1.amazonaws.com`). That bucket returned `403 Access
Denied` from this environment, and the `gdc.xenahubs.net` redirector just
302s to the same restricted bucket, so that path was a dead end. recount3
was the fallback the brief itself suggested, and it worked on the first
try – its AWS Open Data mirror
(`recount-opendata.s3.amazonaws.com/recount3/release/...`) returned the
files over plain HTTP with no auth. Exact files used:

- Gene counts: `https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/tcga/gene_sums/CA/BRCA/tcga.gene_sums.BRCA.G026.gz`
- Sample metadata (TCGA barcodes, sample type, clinical fields): `https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/tcga/metadata/CA/BRCA/tcga.tcga.BRCA.MD.gz`
- QC metadata (average mapped read length, used for the counts conversion below): `https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/tcga/metadata/CA/BRCA/tcga.recount_qc.BRCA.MD.gz`

recount3 itself sources TCGA sequence data from the NCI Genomic Data
Commons, processed uniformly through the Monorail pipeline (STAR +
featureCounts-style summarization), so this is the same underlying TCGA
sequencing data GDC and Xena also serve, just processed and hosted
differently – see the recount3 paper for the processing details ([Wilks et
al. 2021, *Genome Biology*](https://pmc.ncbi.nlm.nih.gov/articles/PMC8628444/)).

One wrinkle worth flagging explicitly: recount3's `gene_sums` files are
**base-pair coverage sums**, not read counts – this is documented in the
recount3 Bioconductor manual and is easy to miss if you don't read past the
column headers. I converted them to approximate integer read counts using
the same formula the official R package's `compute_read_counts()` uses:
`round(gene_sums / average_mapped_read_length)`, using the per-sample
average mapped read length from the `recount_qc` table. This is done in
`scripts/01_data_preparation.py` before anything touches pydeseq2. If you
skip this step the library sizes come out in the billions instead of tens
of millions, which was in fact how I caught the mistake the first time
through.

## The biological question

Which genes are consistently up- or down-regulated in breast tumor tissue
compared to the adjacent normal tissue from the *same patient*, and what
biological processes do those genes belong to? Using matched pairs instead
of pooling unrelated tumor and normal samples controls for a lot of
inter-patient variation (age, genetic background, batch) that would
otherwise be folded into the "tumor vs normal" signal.

## What I actually got (real numbers, this run)

- 25 matched tumor/normal pairs, 50 samples, from a pool of 113 available
  matched pairs in TCGA-BRCA.
- 63,856 genes in the raw recount3 Gencode v26 annotation; 41,619 genes
  passed filtering (total count ≥ 10 across all samples, and detected in
  ≥ 20% of samples).
- Library sizes after the coverage→reads conversion ranged from about 21.6
  to 101.3 million reads per sample (see `results/tables/qc_library_sizes.csv`).
- pydeseq2, fit with design `~patient_id + condition` (paired), found
  **8,333 significant genes** at padj < 0.05 and |log2FC| ≥ 1 – 3,024 up in
  tumor, 5,309 down in tumor.
- The single largest fold-changes downward in tumor were genes like `LEP`,
  `CIDEA`, `CIDEC`, and `CSN2` – these are adipocyte/mammary secretory genes,
  which makes sense: normal breast tissue is largely fat and secretory
  epithelium, and tumor tissue displaces both. Upward, genes like `MMP13`
  and `COL10A1` (stromal remodeling/collagen) and cell-cycle genes dominate.
- GO Biological Process enrichment (Enrichr `GO_Biological_Process_2023`)
  on the up-regulated set is dominated by mitotic/cell-cycle terms –
  mitotic spindle checkpoint signaling, sister chromatid segregation, that
  kind of thing – which is exactly what you'd expect from a proliferating
  tumor. The down-regulated set enriches for cell-cell adhesion and GPCR
  signaling terms, consistent with loss of normal tissue architecture and
  signaling. 34 GO BP terms were significant (padj<0.05) in the up set,
  12 in the down set. Full term tables are in
  `results/tables/enrichment_up_GO_BP.csv` / `enrichment_down_GO_BP.csv`.

## Pipeline

All four steps are plain Python scripts in `scripts/`, meant to be run in
order:

1. **`01_data_preparation.py`** – downloads the recount3 gene count matrix
   and TCGA metadata tables (or reuses `data/` if already downloaded),
   parses TCGA barcodes to find sample type (field 4: `01` = primary
   tumor, `11` = solid tissue normal) and patient ID (first 12 characters
   of the barcode), finds patients with both, randomly samples 25 pairs
   with a fixed seed, converts coverage sums to approximate read counts,
   and writes `raw_counts_matched_subset.csv` + `sample_sheet.csv`.
2. **`02_quality_control.py`** – library size and gene-detection QC,
   low-count gene filtering, PCA and sample-correlation plots on
   log2(CPM+1).
3. **`03_differential_expression.py`** – fits pydeseq2's negative binomial
   GLM with a paired design (`~patient_id + condition`), runs the Wald test
   for Tumor vs Normal, writes the full and significant DE tables, and
   makes the volcano/MA/top-50-heatmap plots.
4. **`04_pathway_enrichment.py`** – maps significant Ensembl gene IDs to
   gene symbols via the MyGene.info REST API, then runs Enrichr GO
   Biological Process enrichment (via `gseapy.enrichr`) separately on the
   up- and down-regulated gene sets, and plots the top terms.

### Reproducing this

```bash
pip install -r requirements.txt
python scripts/01_data_preparation.py
python scripts/02_quality_control.py
python scripts/03_differential_expression.py
python scripts/04_pathway_enrichment.py
```

Each script reads from and writes to `results/tables/` and
`results/figures/`, so they have to run in order the first time. After
`01_data_preparation.py` has downloaded the recount3 files once into
`data/`, re-running from scratch takes a few minutes total (the DESeq2 fit
on 50 samples × ~42k genes is the slow part, roughly two minutes; the
Enrichr calls also add some wall-clock time since they hit a live API).

## Output files

```
results/tables/
  raw_counts_matched_subset.csv   # genes x samples, converted read counts
  sample_sheet.csv                # sample metadata: patient, condition, barcode
  qc_library_sizes.csv            # per-sample library size + detection rate
  filtered_counts.csv             # counts after low-count gene filtering
  pca_coordinates.csv             # PC1/PC2 per sample
  DE_full_results.csv             # all tested genes, DESeq2 stats
  DE_significant.csv              # padj<0.05 & |log2FC|>=1
  DE_top20_upregulated.csv
  DE_top20_downregulated.csv
  up_gene_symbols.csv / down_gene_symbols.csv
  enrichment_up_GO_BP.csv / enrichment_down_GO_BP.csv

results/figures/
  01_library_sizes.png
  02_pca_tumor_vs_normal.png
  03_sample_correlation_heatmap.png
  04_volcano_plot.png
  05_ma_plot.png
  06_top50_heatmap.png
  07_go_enrichment_up.png
  08_go_enrichment_down.png
```

## Design notes / limitations

- **Subsampling.** TCGA-BRCA has 113 matched tumor/normal pairs; I used 25
  (50 samples) to keep the DESeq2 fit and the Enrichr round-trips fast
  enough to iterate on. The random seed (42) is fixed in
  `01_data_preparation.py` so re-running reproduces the same 25 patients.
  Using all 113 pairs would give more power, particularly for genes with
  modest effect sizes, and would be a reasonable next step if compute time
  weren't a constraint.
- **This is bulk RNA-seq, not single-cell.** Everything reported here is an
  average signal across the entire tissue biopsy – tumor cells, stroma,
  infiltrating immune cells, blood vessels, whatever else was in the
  sample. A gene showing up as "down in tumor" could reflect a real
  drop in expression, or just a change in the tissue's cell-type
  composition (e.g. less adipose tissue in the tumor biopsy than the
  normal biopsy). The strong adipocyte/secretory signature in the
  down-regulated genes is a good illustration of this – it's likely partly
  a composition effect, not purely a per-cell expression change.
- **recount3's counts are coverage sums, not read counts**, as noted above.
  The read-length-based conversion is the standard approach recount3's own
  R package uses, but it's still an approximation – it doesn't correct for
  gene length or GC content the way TPM/RPKM would, and pydeseq2's own
  size-factor normalization inside `DeseqDataSet` handles library-size
  differences on top of that. I didn't additionally apply gene-length
  normalization since DESeq2's model is built around raw counts and its own
  normalization, not length-normalized values.
- **Paired design vs. pure tumor effect.** The `~patient_id + condition`
  design controls for each patient's baseline expression, which is more
  conservative and more defensible than pooling unrelated tumor and normal
  samples, but it does mean the model has a lot of coefficients (25 patient
  terms) for a fairly small sample size. This is standard practice for
  paired designs in DESeq2 but worth knowing if someone asks about the
  model's degrees of freedom.
- **GO enrichment gene ID mapping.** Roughly 18% of Ensembl IDs (1,498 of
  8,333) in the significant gene list didn't resolve to a symbol through
  MyGene.info – mostly non-coding RNAs, readthrough transcripts, and other
  IDs without a clean HGNC symbol – so the enrichment analysis ran on the
  ~6,800 that did resolve. This slightly undercounts the true enrichment
  signal but shouldn't bias which pathways come out on top.
- **Enrichr gene set version.** I used `GO_Biological_Process_2023`
  (Enrichr's most recent GO BP library as of this run). Term names and
  exact p-values would shift slightly with a different GO release.

## References

- Muzellec, B., Teleńczuk, M., Cabeli, V., & Andreux, M. (2023). PyDESeq2: a
  python package for bulk RNA-seq differential expression analysis.
  *Bioinformatics*, 39(11), btad547.
  https://doi.org/10.1093/bioinformatics/btad547 – pydeseq2 GitHub:
  https://github.com/owkin/PyDESeq2
- Fang, Z., Liu, X., & Peltz, G. (2023). GSEApy: a comprehensive package for
  performing gene set enrichment analysis in Python. *Bioinformatics*, 39(1),
  btac757. https://doi.org/10.1093/bioinformatics/btac757 – gseapy docs:
  https://gseapy.readthedocs.io/
- Wilks, C., Zheng, S.C., Chen, F.Y., et al. (2021). recount3: summaries and
  queries for large-scale RNA-seq expression and splicing. *Genome Biology*,
  22, 323. https://doi.org/10.1186/s13059-021-02533-6 – recount3 docs:
  https://rna.recount.bio/
- recount3 data access (AWS Open Data mirror used in this project):
  https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/tcga/
- The Cancer Genome Atlas Program (TCGA), National Cancer Institute /
  Genomic Data Commons: https://www.cancer.gov/ccg/research/genome-sequencing/tcga
  and https://gdc.cancer.gov/ – TCGA barcode structure reference:
  https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
- Chen, Y., et al. MyGene.info gene annotation query service, used here for
  Ensembl-to-symbol mapping: https://mygene.info/
- Enrichr: Kuleshov, M.V., et al. (2016). Enrichr: a comprehensive gene set
  enrichment analysis web server 2016 update. *Nucleic Acids Research*,
  44(W1), W90-97. https://doi.org/10.1093/nar/gkw377 –
  https://maayanlab.cloud/Enrichr/
