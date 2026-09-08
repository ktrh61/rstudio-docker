# REBC-THYR driver-conditioned expression analysis

Analysis code for the manuscript on driver-conditioned transcriptomic
differences across radiation-attributability bands in papillary thyroid
carcinoma (REBC-THYR RNA-seq, NCI Genomic Data Commons).

## Layout

Committed code:

- `config.R`, `setup.R` — shared constants (seeds, worker count, thresholds) and paths
- `scripts/` — numbered computation stages (`NNN_<step>.R`); they read `raw/` and
  `processed/` and write `processed/thyr_<step>.rds` (110 and 120 also record
  provenance under `meta/`)
- `lib/` — functions shared by the stages (Brunner–Munzel enumeration with its
  C++ kernel, DEGES–MUREN normalization, gene-set inference, REO panel, plotting)
- `figures/` — one script per figure, `figure_<ID>_<description>.R`, reading
  `processed/` and writing `output/figures/figure_<ID>.png` (300 dpi) and
  `.tif` (600 dpi)
- `tables/` — one script per table or data file, `table_<ID>_<description>.R`
  and `supplementary_data_<n>_<description>.R`, writing `output/tables/<ID>.csv`
  (Table S2 additionally reads `docker/versions.tsv`)
- `tests/testthat` — reference-implementation checks of the in-house statistics
  and the display-item contract
- `Dockerfile`, `docker/` — the qualified execution environment

Working directories (generated at run time, excluded from git except small
committed inputs under `raw/`):

- `raw/` — downloaded GDC counts, clinical table, IREP assigned share, external
  gene lists
- `processed/` — objects passed between stages
- `output/` — figures and tables for the paper
- `meta/` — run provenance (manifests, strand selection, loading metadata)

Table 1 and Table S8 are manuscript text and have no script.

## Environment

- Base image `ubuntu:noble-20260410` (digest-pinned), apt snapshot 2026-04-10
- R 4.5.3 built from source against the reference BLAS/LAPACK 3.12.0
- Bioconductor 3.22; packages from the P3M snapshot of 2026-04-09
  (`docker/versions.tsv`, verified at image build by `docker/verify_environment.R`)
- Gene sets: msigdbr 26.1.0 fetches the pinned MSigDB release (msigdb.2026.1,
  checksum verified by the package) into the R user cache on first use

Reproducibility is scoped to this container; run everything inside it from the
repository root.

## Running

Every script starts with `source("setup.R")`, which defines `paths` and checks
that the working directory is the repository root. Run the stages in
`scripts/` in numeric order, then the scripts in `figures/` and `tables/`:

```sh
docker exec <container> Rscript scripts/010_download_expression.R
# ... through scripts/560_reo_lowmid_confound.R, then figures/*.R and tables/*.R
```

`WORKERS` in `config.R` is part of the reproduction conditions (MUREN gives
each worker its own random stream); the reported runs used 4.

Tests:

```sh
docker exec <container> Rscript -e 'testthat::test_dir("tests/testthat")'
```
