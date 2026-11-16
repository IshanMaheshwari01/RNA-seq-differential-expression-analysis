# TCGA-BRCA Tumor vs. Normal RNA-seq Differential Expression

Differential expression and pathway enrichment on public TCGA-BRCA RNA-seq
data, comparing tumor and normal breast tissue taken from the same patients.
All the numbers below come from running the scripts in `scripts/`.

**Data**: TCGA-BRCA bulk RNA-seq gene counts from the [recount3](https://rna.recount.bio)
AWS Open Data mirror. 25 patients who each had both a primary tumor and a
matched solid-tissue-normal sample (50 samples total), picked at random
(seed 42) out of the 113 such patients available in TCGA-BRCA.

## Where the data comes from

I first tried the UCSC Xena GDC hub copy of `TCGA-BRCA.htseq_counts.tsv.gz`,
but the S3 bucket behind it returned 403 Access Denied, so I switched to
recount3's AWS mirror, which serves the same TCGA data (processed by recount3's
Monorail pipeline: STAR plus featureCounts-style counting). Details in the
[recount3 paper (Wilks et al. 2021, Genome Biology)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8628444/).

One thing to watch: recount3's `gene_sums` files are base-pair coverage sums,
not read counts. I converted them to approximate integer read counts with the
same formula the official R package's `compute_read_counts()` uses,
`round(gene_sums / average_mapped_read_length)`, using the per-sample average
mapped read length from the `recount_qc` table. This happens in
`scripts/01_data_preparation.py` before anything touches pydeseq2. If you skip
it the library sizes come out in the billions instead of tens of millions,
which was in fact how I caught the mistake the first time through.

## The biological question

Which genes are consistently up- or down-regulated in breast tumor tissue
compared to the adjacent normal tissue from the same patient, and what
processes do those genes belong to? Matched pairs control for a lot of
inter-patient variation (age, genetic background, batch) that pooling
unrelated samples would fold into the tumor-vs-normal signal.

## Results

- 25 matched pairs, 50 samples; 41,619 of 63,856 genes passed filtering
  (total count >= 10 across samples, detected in >= 20% of samples).
- Library sizes after the coverage-to-reads conversion: about 21.6 to 101.3
  million reads per sample (`results/tables/qc_library_sizes.csv`).
- pydeseq2, design `~patient_id + condition` (paired): **8,333 significant
  genes** at padj < 0.05 and |log2FC| >= 1: 3,024 up in tumor, 5,309 down.
- Biggest fold-changes down in tumor: `LEP`, `CIDEA`, `CIDEC`, `CSN2` - fat
  cell and milk secretion genes, which makes sense since normal breast
  tissue is largely fat and secretory epithelium. Upward, stromal remodeling
  (`MMP13`, `COL10A1`) and cell-cycle genes dominate.
- GO BP enrichment (Enrichr `GO_Biological_Process_2023`) on the up genes is
  mostly cell-cycle terms (mitotic spindle checkpoint, sister chromatid
  segregation); the down set enriches for cell-cell adhesion and GPCR
  signaling. 34 significant terms up, 12 down. Tables in
  `results/tables/enrichment_up_GO_BP.csv` / `enrichment_down_GO_BP.csv`.

## Pipeline

Four plain Python scripts in `scripts/`, run in order:

1. `01_data_preparation.py` - download recount3 files (or reuse `data/`),
   parse TCGA barcodes for sample type and patient ID, sample 25 pairs with
   a fixed seed, convert coverage sums to read counts, write
   `raw_counts_matched_subset.csv` + `sample_sheet.csv`.
2. `02_quality_control.py` - library size and gene-detection QC, gene
   filtering, PCA and correlation plots on log2(CPM+1).
3. `03_differential_expression.py` - pydeseq2 paired design, Wald test for
   Tumor vs Normal, DE tables, volcano/MA/top-50-heatmap plots.
4. `04_pathway_enrichment.py` - map Ensembl IDs to symbols via MyGene.info,
   Enrichr GO BP enrichment on the up and down sets separately.

```bash
pip install -r requirements.txt
python scripts/01_data_preparation.py
python scripts/02_quality_control.py
python scripts/03_differential_expression.py
python scripts/04_pathway_enrichment.py
```

After the first download, a full re-run takes a few minutes (the DESeq2 fit
on 50 samples x ~42k genes is the slow part, roughly two minutes; the
Enrichr calls hit a live API). Tables land in `results/tables/`, figures in
`results/figures/`.

## Limitations and choices I made

- **25 of 113 pairs.** Keeps the DESeq2 fit and Enrichr round-trips fast
  enough to iterate on; the seed is fixed so re-runs reproduce the same 25
  patients. Using all 113 would give more power for modest effect sizes.
- **Bulk RNA-seq averages the whole biopsy** - tumor cells, stroma, immune
  cells, blood vessels. A gene showing up as "down in tumor" could be a real
  per-cell expression change or just a shift in the tissue's cell-type
  composition. The adipocyte/secretory signature in the down genes is likely
  partly a composition effect.
- **The coverage-to-reads conversion is an approximation** - it doesn't
  correct for gene length or GC content, and I didn't add gene-length
  normalization on top since DESeq2's model is built around raw counts and
  its own size-factor normalization.
- **18% of Ensembl IDs** (1,498 of 8,333) didn't resolve to a symbol through
  MyGene.info (mostly non-coding RNAs and readthrough transcripts), so
  enrichment ran on the ~6,800 that did. This undercounts the signal a
  little but shouldn't bias which pathways come out on top.

## References

- Muzellec et al. 2023, PyDESeq2, *Bioinformatics* 39(11), btad547. [doi:10.1093/bioinformatics/btad547](https://doi.org/10.1093/bioinformatics/btad547)
- Fang, Liu & Peltz 2023, GSEApy, *Bioinformatics* 39(1), btac757. [doi:10.1093/bioinformatics/btac757](https://doi.org/10.1093/bioinformatics/btac757)
- Wilks et al. 2021, recount3, *Genome Biology* 22:323. [doi:10.1186/s13059-021-02533-6](https://doi.org/10.1186/s13059-021-02533-6), data at [rna.recount.bio](https://rna.recount.bio/)
- Kuleshov et al. 2016, Enrichr, *Nucleic Acids Research* 44(W1), W90-97. [doi:10.1093/nar/gkw377](https://doi.org/10.1093/nar/gkw377)
- TCGA / NCI Genomic Data Commons: [gdc.cancer.gov](https://gdc.cancer.gov/), barcode structure in the [GDC docs](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/)
