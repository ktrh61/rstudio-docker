# REBC-THYR Analysis Project

RNA-seq analysis of thyroid cancer (REBC-THYR, GDC) for radiation-effect research.

## Directory Structure

Committed (source):

- `scripts/` — analysis pipeline
- `utils/` — shared utility functions

Not committed (runtime working directories):

- `raw/` — downloaded source data
- `processed/` — intermediate objects passed between scripts
- `output/` — figures, tables, and results for the paper
- `meta/` — provenance and run metadata (data source versions, dates, seeds, session info)

The working directories hold only a `.gitkeep`; their contents are generated at
runtime and excluded via `.gitignore`.

## Environment

- **Base image**: `ubuntu:noble-20260410` (digest-pinned)
- **apt snapshot**: 2026-04-10 (matches base image build date; new installs only, no `apt upgrade`)
- **R**: 4.5.3 (built from source)
- **Bioconductor**: 3.22
- **Packages**: P3M snapshot 2026-04-09 (CRAN / Bioconductor)
- **Editor**: VS Code + Claude Code

The container is the intended execution environment; reproducibility is scoped to it.

## Usage

Run from the repository root inside the container. Each script begins with:

```r
source("setup.R")
```

`setup.R` defines `paths` (the four working directories) and verifies the working
directory is the repository root. It does not load packages; each script declares
its own dependencies.

<!-- TODO: build/run instructions, pipeline overview, citation -->