# Tomato scRNA-seq: Virus-Responsive Transporter Candidate Discovery

A single-cell RNA-seq pipeline identifying transporter genes differentially expressed in response to *Tomato chlorosis virus* (ToCV) infection, followed by cross-validated functional prediction combining sequence homology (BLAST/phylogenetics) and 3D structural comparison (AlphaFold/Foldseek).

Source data: Yue et al. 2024, "Single-cell transcriptome landscape elucidates the cellular and developmental responses to tomato chlorosis virus infection in tomato leaf" (*Plant, Cell & Environment*, DOI: 10.1111/pce.14906), GEO accession GSE201931.

## Environment

* Language: Python 3.10 (conda environment `tomato-sc`)
* Libraries: Scanpy, AnnData, Pandas, NumPy, Biopython, scikit-learn, kneed, Matplotlib
* External CLI tools: BLAST+, MAFFT, FastTree (run via WSL/Ubuntu shell, not inside notebooks)
* External services: ColabFold (Google Colab, AlphaFold2) for structure prediction, Foldseek web server for structure-based search

## Getting Started (How to Run)

1. Install dependencies: `pip install -r requirements.txt`
2. Install BLAST+, MAFFT, FastTree (e.g. via `conda install -c bioconda blast mafft fasttree`)
3. Download raw scRNA-seq data from GEO GSE201931 and place under `dataset/raw/`
4. Download reference files (ITAG4.0 annotation/proteome, Arabidopsis TAIR10.1 proteome) and place under `dataset/reference/` and `dataset/blast/`
5. Run notebooks and scripts in order:
`01_preprocessing.ipynb` → `02_pca_clustering_celltype_deg.ipynb` → `03_sequence_extraction.ipynb`
→ `scripts/01_download_arabidopsis_db.sh` → `scripts/02_build_blast_db.sh` → `scripts/03_run_blastp.sh`
→ `04_ortholog_blast_msa.ipynb` (which also calls `scripts/04_run_mafft.sh`, `scripts/05_run_fasttree.sh` mid-pipeline)
→ external: run ColabFold on candidate sequences, then Foldseek web search on resulting structures
→ `05_analysis.ipynb`

## Data

|Dataset|Type|Source|Version|Format|
|-|-|-|-|-|
|Tomato leaf scRNA-seq (H vs ToCV-infected)|Single-cell expression matrix|[GEO GSE201931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201931)|—|10x mtx (barcodes/features/matrix)|
|Tomato gene annotation|Functional description, GO terms|[ITAG4.0](https://solgenomics.net/)|ITAG4.0|TXT, FASTA|
|Arabidopsis reference proteome|Protein sequences|[NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/)|TAIR10.1|FASTA|

## Pipeline

### 1. `01_preprocessing.ipynb` — QC and Normalization

Loads raw CellRanger output, filters low-quality cells (MAD-based outlier detection on total counts/gene counts, percentile cutoff on mitochondrial/chloroplast content), removes doublets (Scrublet), and normalizes expression.

- **Note**: ITAG4.0 GO annotation flags very few mitochondrial (8) and chloroplast (49) genes, likely due to incomplete annotation in this non-model species — the resulting QC metric has limited sensitivity compared to total-count-based filtering.
- **Note**: Scrublet was run on all samples combined rather than per sample; doublet simulation assumes a single homogeneous population, so per-sample execution would be more rigorous.

### 2. `02_pca_clustering_celltype_deg.ipynb` — Clustering, Annotation, and DEG

Clusters the healthy sample (PCA → Leiden), annotates clusters with cell types using marker genes from the source publication's supplementary table, then maps the infected sample onto the same space (`sc.tl.ingest`) and tests differential expression (V vs H)
within each cell type.

- **Note**: `highly_variable_genes` was run with `batch_key="sample_no"` after the data had already been subset to a single sample, so the batch correction had no effect in this run.
- **Note**: `n_neighbors` for the KNN graph reused the PCA elbow component count (13) rather than being set independently via `n_pcs`; the two parameters serve different purposes and were conflated here.
- **Note**: The Vascular cell type is defined almost entirely by a single marker gene (PP2A1), since 3 of 4 reference markers (APL, ACL5, SHR) had no matching tomato gene ID.
- **Note**: Ingest was chosen as a computationally lightweight baseline mapping approach for this pilot study.
- **Note**: DEG results for Guard cell (121 V cells) and Unknown (55 V cells) have reduced statistical power due to small sample size.
- **Note**: Final candidate selection (15 genes) was manual: genes annotated as putative/"-like" were excluded, transporters with clearly resolved function were prioritized by logFC, and ABC transporter family members were retained given their established role in plant defense responses.

### 3. `03_sequence_extraction.ipynb` — Candidate Sequence Retrieval

Maps each candidate gene ID to its ITAG4.0 protein sequence and writes a FASTA file for BLAST search.

### 4. `04_ortholog_blast_msa.ipynb` — Ortholog Inference

Searches candidates against the Arabidopsis proteome (BLAST), builds a multiple sequence alignment and phylogenetic tree per candidate from the top 5 hits, and determines the single closest Arabidopsis ortholog by tree distance.

- **Note**: `evalue=1e-5` and `max_target_seqs=5` (see `scripts/03_run_blastp.sh`) were practical choices — permissive enough to avoid missing true homologs, with the final ortholog resolved more precisely afterward via tree distance rather than raw BLAST score alone. These values were not derived from a parameter sweep.
- **Note**: The Arabidopsis description attached to each candidate in `ortholog_summary.csv` is the single closest ortholog by tree distance, not simply the top BLAST hit. This value is used directly as the sequence-based functional prediction in notebook 05.

### 5. `05_analysis.ipynb` — Structural Cross-Validation

(External) Predicts 3D structures for 13 of 15 candidates via ColabFold (AlphaFold2);
the 2 largest proteins (1205, 1514 aa) were excluded due to GPU memory limits on the free-tier Colab T4 GPU. (External) Searches each structure against the PDB100 database via the Foldseek web server. This notebook parses the Foldseek results, restricts to PDB100 hits (excluding AlphaFold-DB hits, which are largely close homologs from other plant genomes and add little independent information), and compares the top structural match against the sequence-based ortholog prediction from notebook 04.

## Key Findings & Interpretation

Out of 13 candidate transporters evaluated via dual sequence-structure cross-validation:

| Category / Concordance Level | Count | Summary & Interpretation |
| :--- | :---: | :--- |
| **High Confidence (Exact Match)** | **9** | Sequence ortholog & top 3D structural hit align precisely (e.g., Cation efflux ↔ Zinc Transporter YiiP) |
| **Broad Family Match** | **2** | Matches at the broader transporter family level (Nucleotide / Anion transporters) |
| **Structure-Assisted Hypothesis** | **1** | Sequence homology was uncharacterized; 3D structural alignment suggested a plausible transporter fold hypothesis |
| **Low Confidence / Unresolved** | **1** | Weak Foldseek alignment score (treated as unresolved) |

1. **9 candidates** show an exact functional match between the Arabidopsis ortholog and the top Foldseek structural hit (e.g. Cation efflux protein ↔ Zinc Transporter YiiP), indicating high-confidence functional assignments from two independent lines of evidence.
2. **2 candidates** match at the broader transporter-family level (nucleotide and anion transporters respectively), rather than an exact functional match.
3. **1 candidate** (Solyc10g051120.3) showed structural similarity to sugar/metabolite transporter folds (e.g., SWEET-like architecture) via Foldseek. While primary annotations suggest involvement in pyruvate/organic acid transport, this structural match highlights potential alternative substrate specificities or shared structural domain folds that warrant further experimental validation.
4. **1 candidate** (Solyc07g063520.3) had a weak Foldseek evalue and is treated as low-confidence.

## Critical Methodological Evaluation & Future Research Proposal

### 1. Biological / QC Bottlenecks
- **Incomplete organelle gene annotation**: mt/chloroplast QC relied on ITAG4.0 GO terms, which under-annotate these gene sets in tomato.
  *Next step*: use a curated mitochondrial/chloroplast gene ID list if available.
- **Doublet detection across pooled samples**: Scrublet's synthetic-doublet approach assumes a homogeneous population; running it per sample/replicate would better match this assumption.
- **Single-marker cell type definition (Vascular)**: relying on one matched marker gene is a weaker basis for cell-type identity than a multi-marker consensus.
  *Next step*: recover missing marker orthologs (APL, ACL5, SHR) via BLAST rather
  than description-keyword matching.
- **Low-power DEG groups**: Guard cell and Unknown cell-type DEG results should be interpreted cautiously given small sample sizes (121 and 55 V cells respectively).

### 2. Pipeline / Parameter Refinements
- **HVG batch correction**: `batch_key` was applied after the data was already subset to one sample, making it a no-op in this run; batch correction should either be removed here or applied before subsetting if cross-sample HVGs are needed.
- **Neighbor graph parameters**: `n_neighbors` and `n_pcs` should be set independently rather than reusing the PCA elbow component count as the neighbor count.
- **BLAST threshold justification**: `evalue`/`max_target_seqs` were practical defaults rather than optimized; a parameter sensitivity check would strengthen confidence in the ortholog set.
- **AlphaFold scale limits**: 2 of 15 candidates exceeded the free-tier Colab GPU's memory; a higher-memory GPU or domain-split prediction would allow full coverage.

## Repository Structure

```
tomato-sc/
├── notebook/
│   ├── 01_preprocessing.ipynb              # QC, doublet removal, normalization
│   ├── 02_pca_clustering_celltype_deg.ipynb # Clustering, cell-type annotation, DEG
│   ├── 03_sequence_extraction.ipynb         # Candidate sequence retrieval
│   ├── 04_ortholog_blast_msa.ipynb          # BLAST, MSA, tree-based ortholog inference
│   └── 05_analysis.ipynb                    # Sequence vs. structure cross-validation
├── scripts/
│   ├── 01_download_arabidopsis_db.sh
│   ├── 02_build_blast_db.sh
│   ├── 03_run_blastp.sh
│   ├── 04_run_mafft.sh
│   └── 05_run_fasttree.sh
├── dataset/
│   ├── raw/            # Original CellRanger output (not included)
│   ├── reference/       # ITAG4.0 annotation, marker gene table (not included)
│   ├── processed/       # h5ad, candidate tables, ortholog_summary.csv
│   ├── blast/           # BLAST DB, query fasta, blast_results.tsv
│   ├── msa/              # groups/, aligned/
│   └── trees/            # FastTree .nwk output
└── result/
    ├── alphafold/          # AlphaFold output files (not included)
    ├── foldseek/          # Foldseek raw JSON (not included)
    └── analysis/          # ortholog_vs_foldseek.csv
```
