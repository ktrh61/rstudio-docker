# Additional descriptive audit of R_Tumor

This audit reads finalized analysis artifacts. It does not rerun QC, normalization, hypothesis testing, or multiple-testing adjustment. It creates no new p-values and no alternative gene list.

- Existing q<0.10 list: 1765 genes (971 higher and 794 lower in High).
- Among those genes, median absolute signed Brunner–Munzel effect was 0.600 (IQR 0.567–0.644).
- Their median absolute display-only mean log2 fold change was 0.422 (IQR 0.290–0.698).
- PCA used all 15621 tested genes; PC1 and PC2 explained 19.7% and 16.6% of variance.
- In the finalized main cohort, median relative purity was 0.783 in R_Sporadic and 0.822 in R_High; these values are from the main-cohort analysis object and are distinct from the later all-band REO purity diagnostic.
- A candidate adjusted design with group, age at surgery, relative purity, sex, and RET fusion partner used 27 complete cases, had 7 columns, rank 7, and 20 residual degrees of freedom. This establishes numerical estimability only; no adjusted gene model was fitted.
- The existing PC-OD output flagged 0 RET tumor and 0 RET normal samples; the only main-cohort target flag was one B_High tumor. Thus PC-OD did not alter the realized R_Tumor sample set.

See the CSV files for complete summaries and sample-level values.
