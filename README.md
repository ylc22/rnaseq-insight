# RNAseq Insight 🧬

**A reproducible RNA-seq analysis workflow that turns raw count matrices into interpretable biological signals.**

[![CI](https://github.com/ylc22/rnaseq-insight/actions/workflows/ci.yml/badge.svg)](https://github.com/ylc22/rnaseq-insight/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)

## Why this project

Many RNA-seq examples stop at exploratory plots. **RNAseq Insight** packages the full analysis path into reusable Python modules with QC, normalization, statistics, multiple-testing correction, visualization, tests, and a reproducible demo dataset.

### Pipeline

```text
count matrix
   ↓
sample QC
   ↓
CPM expression filtering
   ↓
log-CPM normalization
   ↓
PCA / sample separation
   ↓
differential expression
   ↓
Benjamini–Hochberg FDR
   ↓
ranked gene report + figures
```

## Highlights

- Library-size and zero-inflation QC
- Low-expression filtering by counts-per-million
- Log-CPM normalization
- PCA with experimental-condition labels
- Welch's t-test differential-expression analysis
- Benjamini–Hochberg false-discovery-rate correction
- Effect-size ranking
- Volcano and PCA visualization
- Synthetic RNA-seq cohort with planted treatment signal
- Unit tests and GitHub Actions CI

## Quickstart

```bash
git clone https://github.com/ylc22/rnaseq-insight.git
cd rnaseq-insight
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_demo.py
```

The demo generates its own input data and writes results to `outputs/`.

## Expected output

The synthetic experiment plants differential expression in a subset of genes. The workflow should recover many of these genes near the top of the ranked report while controlling the false-discovery rate.

Typical artifacts:

```text
outputs/
├── differential_expression.csv
├── pca.png
└── volcano.png
```

## Repository structure

```text
src/rnaseq_insight/   reusable analysis package
scripts/run_demo.py   end-to-end reproducible demo
tests/                unit tests
data/                  generated demo inputs
outputs/               generated reports and figures
```

## Statistical approach

For each gene, the demo compares experimental groups on normalized expression using Welch's t-test, then applies Benjamini–Hochberg correction across genes. The implementation is intentionally transparent so each transformation can be inspected and replaced independently.

For real production studies, the same interfaces could be extended to DESeq2/edgeR outputs, batch-aware models, gene annotation, pathway enrichment, and larger workflow systems such as Nextflow.

## Tech stack

`Python` · `pandas` · `NumPy` · `SciPy` · `scikit-learn` · `matplotlib` · `pytest` · `GitHub Actions`

## Reproducibility

The demo uses a fixed random seed and generates all required data locally; no private or proprietary datasets are included.

## License

MIT
