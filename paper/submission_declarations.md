# 投稿事務ドラフト(タイトルページ+Additional Information 骨子)

状態: draft(2026-08-21 起草 — 【記入】は研究者記入待ち。BJC「Article formatting」の要求節順に対応。
確定後に投稿版へ転記する。gpt_review/review.md §10 の残作業リストの実体化)。

## Title page

- **Title**: Driver-conditioned transcriptomic differences across radiation-attributability bands in papillary thyroid carcinoma(≤150字 ✓・結論文でない ✓)
- **Authors / Affiliations**:
  1. Kotaro Harakawa (1) *corresponding
  2. Vladimir A. Saenko (1)
  3. Norisato Mitsutake (1)

  Affiliations:
  1. Department of Molecular Oncology and Diagnostic Medicine, Atomic Bomb Disease Institute, Nagasaki University, 1-12-4 Sakamoto, Nagasaki 852-8523, Japan

  (綴りは PubMed 出版名義で照合済み 2026-08-26。所属名は研究者提供の現名称 — 2022 年論文の
  旧表記と異なる。2025-04 改称の新名称 Molecular Oncology and Diagnostic Medicine を採用 —
  研究室ページ注記に基づき研究者確定 2026-08-26。英語ヘッダ未更新のため共著者検収で再確認)
- **Corresponding author + email + ORCID**:
  Kotaro Harakawa; ktrh61@gmail.com; ORCID 0009-0004-1086-8046
  (大学は卒業生・離籍者向けドメインを提供しないため個人メール)
- **Word count**: 投稿時に `paper/make_submission.py` の出力で再計測(2026-08-26 時点: 本文 4,933・Abstract 197)

## Additional Information(BJC 規定順)

### Acknowledgements

【記入】

### Authors' contributions

KH conceived the study, developed and executed the full analysis pipeline, performed all
statistical analyses, and drafted the manuscript. VAS and NM contributed the driver-based
stratification design and the Assigned Share banding approach, including the specification
of the NIH IREP input parameters; VAS additionally proposed the initial concept of a
per-case predictive classification, which motivated the relative-expression-ordering panel
designed and implemented by KH. All authors reviewed and approved the final manuscript.

(貢献実態の研究者指定 2026-08-26: oncogene 設定・AS 分割・IREP パラメータ = VAS/NM、
他全て = KH。文言は共著者ラウンドで本人ら確認)

### Ethics approval and consent to participate

This study is a secondary analysis of publicly available, de-identified data obtained from
the open-access tier of the NCI Genomic Data Commons (project REBC-THYR). No new data
involving human participants were collected by the authors, and no additional ethical
approval was required. Ethical approvals and participant consent for the original data
collection are described in the source study (Morton et al. 2021).
(研究者承認 2026-08-26 — Thyroid 誌投稿時のオールオープンデータ回答の経緯と同型)

### Consent for publication

Not applicable. The manuscript contains no individual person's identifiable data.

### Data availability

案: Gene-level RNA-seq counts and clinical data are openly available from the NCI Genomic
Data Commons (project REBC-THYR; STAR - Counts, open access). The exact 906-file download
manifest (MD5 7defb0c5574453474c67dfac8367a589) is provided with the analysis code
<!-- N-89 -->. Processed analysis objects sufficient to regenerate the reported figures and
tables accompany the paper. 【公開リポジトリの URL/DOI 確定後に記入】

### Code availability

案: Analysis code, versioned inputs, and the date-pinned container build recipe sufficient
to regenerate the reported analyses are available at 【リポジトリ URL/DOI — 公開範囲は
「論文再現に必要なスクリプトのみ」の方針(2026-08-14 研究者決定)に従い確定】.

### Competing interests

The authors declare no competing interests.

### Funding

The authors received no specific funding for this work.

## Supplementary file descriptions(投稿システム入力用 — BJC GTA「各ファイルに ≤50 語の要約」。2026-08-29 起草・共著者レビュー対象)

### Supplementary Material (PDF)

Supplementary Methods and Results, Figures S1–S2, Tables S1–S4, S6 and S7, and Supplementary References (Table S5 and Supplementary Data 1–2 are separate files). Covers data sources, Assigned Share inputs, quality control, normalization, gene-level, omnibus and gene-set inference, REO panel construction, external gene-list comparison, and AI-use verification.

### Table S5 (table_s5_ora_annotation.csv)

Descriptive over-representation annotation of the RET-tumor q<0.10 gene list (18,576 rows). For each family (Hallmark, C2 canonical pathways, C5 GO Biological Process, radiation-curated) and list (higher, lower, combined): set size, observed and expected overlap, one-sided hypergeometric p-value, and family-by-list Benjamini–Hochberg q-value.

### Supplementary Data 1 (supplementary_data_1_gene_level_results.csv)

Complete gene-level results for every tested gene in each of the four contrasts (62,952 rows): contrast, Ensembl identifier, gene symbol, Brunner–Munzel relative effect θ, exact permutation p-value, and Storey q-value. Supports cross-referencing and reanalysis at thresholds other than the reported q<0.10 rule.

### Supplementary Data 2 (supplementary_data_2_set_level_results.csv)

Complete label-permutation gene-set results for all tested sets in each contrast and collection (24,798 rows): set size, enrichment score, normalized enrichment score, sign-conditional permutation p-value, within-collection Benjamini–Hochberg q-value, redundancy annotation, and leading-edge genes.

## AI 開示の転記先メモ

- BJC: Methods(Reproducibility and AI use 節)に短文開示、Supp Methods(Scope and verification of AI use 小節)に機構詳細を記載済み — 追加転記は不要(GTA は Methods 文書化を要求)。
- ERC 転送時: Declarations 欄へ短縮形(numbers_ledger 改訂メモ 2026-08-18 に確定文あり)。
