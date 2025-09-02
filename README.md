# REBC-THYR Analysis Project

## Project Structure
- `analysis_v7/`: Current analysis (expanded dataset, 906 samples)  
- `utils/`: Shared utility functions
- `legacy/`: Previous versions (v6)

## Environment
- **Docker**: Bioconductor 3.21 (fixed)
- **R**: 4.5.1

## Quick Start
1. `make build`
2. `make run`  
3. http://localhost:8787
4. In RStudio: `source("analysis_v7/setup.R")`

## Current Analysis (v7)
- Dataset: REBC-THYR from GDC
- Samples: 906 (expanded from v6's 828)
- Branch: v7-data-expansion
