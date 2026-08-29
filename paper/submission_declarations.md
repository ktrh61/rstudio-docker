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

## 公開リポジトリ・プレプリント運用(2026-08-29 記録 — 共著者の同意待ち。決定後に不採用側を削除)

前提(規定はライブ確認 2026-08-29): BJC・npj Precision Oncology(Nature Portfolio)・ERC(Bioscientifica)の
3 誌ともプレプリント許容(BJC = 投稿時 DOI 開示+本文で参照/npj = DOI+ライセンス開示/Bioscientifica =
投稿時通知+出版後にプレプリントから最終版へリンク)。コードは審査中「求めに応じて利用可能」で足り、
公開義務は出版時。bioRxiv は雑誌出版まで改訂版を同一 DOI 下に投稿可(出版後は不可)、投稿は全著者同意が条件。
判断軸(共著者説明用): 公開データ+公開ソフトの再解析は参入障壁がなく、コホート所有者(Morton ら)が最も
出しやすい当事者/実行コスト低下で希少資源は設計側に移り、防御は非公開でなく日付つき開示/被引用の副次
効果(Fu & Hughey 2019 eLife 等)。制約 = 全著者同意・v1 は永久に残る・優先権は規範であり評価ではない。

### ルート A(推奨 — 共著者同意あり): bioRxiv を雑誌投稿と同日に掲載
1. 公開リポジトリを切り出し(published-repo-scope)→ GitHub 公開 → Zenodo 連携リリース → DOI 取得。
2. bioRxiv に投稿版と同一の全文(本文+Supp+図表)を投稿(Cancer Biology または Genomics、ライセンスは
   全著者で選択)。掲載承認で DOI 確定。
3. 原稿へ反映: Data/Code availability = GitHub URL+Zenodo DOI/プレプリントを本文で参照し文献に収載
   (BJC 規定)/カバーレターの "has not been posted as a preprint" を「posted on bioRxiv (DOI, licence)」へ。
4. 雑誌投稿(同日)。投稿システムでプレプリント DOI・ライセンスを開示。
5. 審査中: 再投稿ごとに改訂版を bioRxiv へ(v2, v3 …)。受理前の最終改訂版まで可。
6. 受理・出版後: Zenodo に最終版タグ/bioRxiv は約 2 週間で出版版へ自動リンク/SharedIt リンクを配布/
   購読ルートなら出版 6 か月後に accepted manuscript を機関リポジトリへ(Springer Nature STM エンバーゴ)。
   転載(BJC 不採択 → ERC/npj)時はプレプリントはそのまま、次誌の通知・開示規定に従う。

### ルート B(フォールバック — 同意なし): Zenodo 予約 DOI+査読専用添付
1. 切り出しは同じ。GitHub は非公開のまま。
2. Zenodo に下書きレコードを作り DOI を予約(Publish しない — 下書きは非公開)。
3. 原稿へ反映: Data/Code availability = 予約 DOI+"to be released on publication"/カバーレターは
   「has not been posted as a preprint」のまま、コード文を「provided for review; released under a reserved
   DOI on publication」へ。
4. 査読用に切り出しアーカイブを投稿システムへ査読専用ファイルとして添付。
5. 受理時: GitHub 公開+Zenodo Publish(予約 DOI がそのまま有効)。出版後は A-6 と同じ。

両ルート共通で今から着手可能: 切り出し・Zenodo アカウント(ORCID 連携)・フェーズ6 再構築確認。

## AI 開示の転記先メモ

- BJC: Methods(Reproducibility and AI use 節)に短文開示、Supp Methods(Scope and verification of AI use 小節)に機構詳細を記載済み — 追加転記は不要(GTA は Methods 文書化を要求)。
- ERC 転送時: Declarations 欄へ短縮形(numbers_ledger 改訂メモ 2026-08-18 に確定文あり)。
