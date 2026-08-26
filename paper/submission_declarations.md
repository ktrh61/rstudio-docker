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
  1. Department of Molecular Oncology, Atomic Bomb Disease Institute, Nagasaki University, 1-12-4 Sakamoto, Nagasaki 852-8523, Japan

  (綴りは PubMed 出版名義で照合済み 2026-08-26。所属名は研究者提供の現名称 — 2022 年論文の
  旧表記 Radiation Molecular Epidemiology 等と異なるため、現名称であることの最終確認は研究者)
- **Corresponding author + email + ORCID**:
  Kotaro Harakawa; ktrh61@gmail.com; ORCID 0009-0004-1086-8046
  (大学は卒業生・離籍者向けドメインを提供しないため個人メール)
- **Word count**: 投稿時に `paper/make_submission.py` の出力で再計測(2026-08-26 時点: 本文 4,933・Abstract 197)

## Additional Information(BJC 規定順)

### Acknowledgements

【記入】

### Authors' contributions

【記入】(BJC: 全著者をイニシャルで個別に)

### Ethics approval and consent to participate

案(研究者確認要): This study is a secondary analysis of open-access, de-identified data
from the REBC-THYR project obtained through the NCI Genomic Data Commons; no new data were
collected from human participants. 【所属機関の倫理審査の要否・承認番号等を確認して記入】

### Consent for publication

案: Not applicable(個人を特定できる情報・画像は含まれない)。

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

【記入】

### Funding

【記入】(BJC は Acknowledgements 内でも可 — 誌の最新書式で確認)

## AI 開示の転記先メモ

- BJC: Methods(Reproducibility and AI use 節)に短文開示、Supp Methods(Scope and verification of AI use 小節)に機構詳細を記載済み — 追加転記は不要(GTA は Methods 文書化を要求)。
- ERC 転送時: Declarations 欄へ短縮形(numbers_ledger 改訂メモ 2026-08-18 に確定文あり)。
