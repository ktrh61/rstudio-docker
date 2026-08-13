# Methods ドラフト(全節)

状態: draft(Claude Code 起草 2026-08-13、研究者が書き直す前提の下書き)。
写像元: scripts/010〜530(各ヘッダ)、diagnostics/(signature_agreement.R、
reo_lowmid_{outliers,purity,confound}.R、external_gene_anchors.R、
deg_ora_annotation.R、age_arm_difference.R)、config.R、lib/。
数値は全て numbers_ledger の N-ID 経由(タグ併記)。設定値・規則の N 行は
セクション P(N-66〜N-76、2026-08-13 追加)として台帳化済み。
D6 較正と ORA 注釈の文面の正本はそれぞれ draft_methods_results_d6.md /
draft_ora_annotation.md — 本稿には接続位置を示すため同一文を再掲し、その旨を注記。

---

## Data sources and expression matrix

<!-- 写像元: 010, 020, 030, 110, 120 -->

Gene-level RNA-seq counts (STAR - Counts, open access) for the REBC-THYR
project were downloaded from the NCI Genomic Data Commons and md5-verified
against the download manifest. Gene lengths (exon-union length per gene) were
derived from the GDC GENCODE v36 reference annotation <!-- N-66 -->. The
clinical table is Data S1 of the REBC-THYR publication (440 cases
<!-- N-11 -->), read without dropping columns or editing values; missing
markers were coerced to NA and columns typed numerically only where every
non-missing value parses.

Expression files were mapped to cases and biospecimen samples through the GDC
API, and a single count assay per sample was assembled into one expression
matrix (58,448 genes × 906 samples after dropping all-zero genes
<!-- N-14 -->). Library strandedness was decided per sample from the STAR
stranded-column totals (ratio rule: the library is called stranded when the
smaller stranded total is at most half the larger, and the larger column is
used; otherwise the unstranded column is used <!-- N-67 -->); TPM and FPKM
columns were discarded at loading.

## Exposure metric: assigned share

<!-- 写像元: 130, 140 -->

Radiation attributability per exposed case was quantified as the NIOSH-IREP
thyroid-model Assigned Share (AS; probability of causation associated with
the expected value of the excess relative risk), computed from each case's
recorded thyroid dose and ages with a fixed input convention (electrons
E > 15 keV; dose in cSv = mGy/10; exposure year 1986; birth year = 1986 −
age at exposure; surgery year = birth year + age at surgery <!-- N-68 -->).
The values are carried as a versioned input file; provenance and the
input-identity audit are recorded in the analysis plan.

Cases were banded on AS with a single pre-fixed rule: unexposed cases
(dose 0) form the Sporadic band; exposed cases fall in Low (0 < AS ≤ 33.3),
Mid (33.3 < AS < 66.6), or High (AS ≥ 66.6); no case sits on a boundary
<!-- N-69 -->. The clinical composition of all 440 cases across driver ×
band × pairing is given in Table 1(仮) <!-- N-11 -->.

## Quality control and analysis cohorts

<!-- 写像元: 210, 220, 230 -->

All eligibility conditions were finalized into explicit per-case flags before
any inference, and every downstream script consumes the same flags. Sample
outliers were screened with principal-component outlier detection on
unnormalized log-CPM within each group × tissue sub-matrix, before purity
estimation so that an anomalous sample cannot contaminate the purity fit.
Relative tumor purity was then estimated per case with ContamDE on
MUREN-normalized paired counts, in one run per driver cohort with both
exposure arms together: ContamDE purities are relative within the set of
pairs estimated jointly, so a pooled run is what gives both arms one common
scale on which a single threshold has a single meaning. The common tumour
reference this assumes rests on the same premise as the driver
stratification itself — tumour expression in both arms is dominated by the
driver's biology — and the tumour-versus-normal contamination axis that
purity measures is a different axis from the exposure contrast under test;
the estimate serves only as a within-cohort relative filter and a
diagnostic covariate.

The main analysis cohort comprises driver-classified (RET or BRAF) cases in
the Sporadic or High band with a paired tumor/normal sample, both tissues
outlier-clean, and pooled relative purity ≥ 0.6 <!-- N-70 -->. The case flow
runs 440 → 248 → 77 → 70 → 69 → 63 <!-- N-08 -->, giving 63 cases: B_High 9,
B_Sporadic 27, R_High 15, R_Sporadic 12 <!-- N-09 -->. The REO training set
is the RET subset of this cohort (27 cases), and the REO evaluation set is
the paired RET tumours of the Low and Mid bands (36 cases: R_Low 17, R_Mid
19 <!-- N-10 -->), deliberately left unfiltered by the outlier and purity
conditions (its own QC is reported as diagnostics, below).

<!-- C-10 -->
Arm age structure is reported descriptively (per-group median and range,
Table 2(仮) <!-- N-12, N-13 -->) together with interval estimates of the
between-arm difference (below); age is not entered as a covariate, because
age at exposure is a component of the AS definition and adjusting for it
would remove part of the exposure metric itself.

## Analysis units

<!-- 写像元: 310 header, lib/units.R -->

All inference is organized in four analysis units — the exposed-vs-sporadic
contrast within each driver stratum and tissue: R_Tumor, R_Normal (RET arm),
B_Tumor, B_Normal (BRAF arm). Each unit is tested within itself; no
cross-unit family-wise inference is performed, and no study-wide FDR is
claimed. In every unit, group x is Sporadic and group y is High, so
effects > 0.5 indicate higher expression in the exposed arm.

## Normalization

<!-- 写像元: 310, lib/norm_deges.R -->

Each unit's count matrix (protein-coding genes, filterByExpr) was normalized
with a DEGES scheme: MUREN normalization alternating with a permutation
Brunner–Munzel screen that removes potential DEGs (BH q < 0.10) from the
scaling-factor estimation, iterated three times (iDEGES) <!-- N-71 -->. The
per-unit tested-gene counts, screen pi0 estimates, iteration convergence
(Jaccard), and resulting scaling-factor ranges are reported in
Supplementary(仮) <!-- N-15 -->.

## Gene-level differential expression

<!-- 写像元: 410 -->

Each gene was tested for a High-vs-Sporadic difference with the exact
permutation Brunner–Munzel test (two-sided; exact enumeration of all
C(n, nx) allocations, so the gene-level p-values need no seed; the saved
label-shuffle index for unit-level nulls uses 9,999 shuffles at seed
19860426 <!-- N-04, N-05 -->). Gene-level inference is the Storey q-value on
the exact p-values with the plug-in pi0 estimate at λ = 0.5, thresholded at
q < 0.10 protocol-wide <!-- N-05 -->; pi0 is reported with
permutation-calibrated uncertainty. A unit-level omnibus permutation test
accompanies the gene table; Higher Criticism (α0 = 0.1) is the pre-specified
primary omnibus statistic, with count- and max-type rows reported
descriptively <!-- N-05 -->. The full rejection curve R(α) is retained so
that results do not depend on displaying a single threshold count.

## Gene-set level inference

<!-- 写像元: 420, lib/gsea_permutation.R, lib/gsea_collections.R -->

Set-level enrichment consumes the whole per-gene ranking (threshold-free; no
DEG-list cut decides what is tested), ranked by tie-averaged normal scores of
the signed Brunner–Munzel statistic. The enrichment score is the weighted
running sum (gseaParam = 1) evaluated at tie-block boundaries only, which on
tie-free input equals the standard GSEA statistic <!-- N-72 -->. The null
reuses the identical label shuffles saved by the gene-level test (9,999 per
unit <!-- N-07 -->), so gene- and set-level results rest on one permutation
null. Inference is the per-set, sign-conditional permutation p-value with
Benjamini–Hochberg adjustment within each collection, q_bh < 0.10
<!-- N-72 -->; no cross-family claim is made.

Four MSigDB collections were tested (msigdbr 26.1.0): Hallmark (50), C2:CP
(3,910; Reactome, WikiPathways, KEGG MEDICUS, BioCarta and PID), C5:GO:BP
(7,538), and a radiation-curated C2:CGP family (28) whose curation rule was
fixed before this inference touched real data <!-- N-29, N-55 -->. Sets
outside the size window 15–500 were excluded <!-- N-72 -->.

<!-- C-06 -->
The complete gene-set inference was calibrated on held-out null replicates —
label permutations pushed through the identical procedure — before being
applied to the real contrasts, and the choice of the set-level FDR procedure
itself was fixed by this calibration prior to any real-data run
(Supplementary Methods) <!-- N-06 -->. A spike-in control (one Hallmark set
× 1.15 in the 9 samples of one arm) served as the sensitivity counterpart
<!-- N-30 -->.
(= draft_methods_results_d6.md 本文 Methods と同一文。正本はあちら。サプリ層も同稿)

<!-- C-14 -->
As a descriptive complement, the discovered R_Tumor list was annotated by
one-sided hypergeometric over-representation against the identical set
universe (up / down / combined lists; universe = the unit's tested genes;
BH within family × list), reported in full as hypothesis-generating
annotation (Supplementary Table SY) <!-- N-59 -->.
(= draft_ora_annotation.md 起草メモの Methods 併置文と同一。正本はあちら)

## Between-arm signature agreement

<!-- 写像元: diagnostics/signature_agreement.R -->

For each tissue, the exposed-vs-sporadic contrast was summarized per gene by
the signed Brunner–Munzel statistic separately in the RET and BRAF arms, and
the two genome-wide profiles compared by Spearman correlation over shared
genes — threshold-free and without any dose assumption. Each coefficient
travels with a permutation-calibrated reference band: labels were shuffled
independently within each unit (9,999 shuffles, per-unit seeds from base
19450809 <!-- N-73 -->) and the shuffled profiles correlated, giving the
null spread of rho when neither arm carries label-aligned structure. The
normal-tissue pair is the pre-specified hypothesis-bearing comparison; the
tumor pair is a descriptive completion. An observed rho outside the band
indicates shared label-aligned structure; by itself it does not identify
that structure as an exposure trace rather than a shared covariate.

## REO panel: construction and out-of-sample evaluation

<!-- 写像元: 510, 520, 530 -->

Relative expression ordering (REO) works on within-sample gene ranks and is
therefore free of between-sample normalization; expression was taken as TPM
(length-normalized) because within-sample comparison across genes requires
it. Candidate pairs were generated in the RET tumour training arms from the
top 500 genes by Brunner–Munzel effect magnitude |effect − 0.5| (a
threshold-free pool; the q < 0.10 DEG set is deliberately not used). A pair
qualifies when its within-sample log2-TPM difference r has a stable sign in
Sporadic (dead zone |r| < log2(1.2); at most one exception among
non-dead-zone samples; 10th percentile of |r| ≥ log2(1.5)) and reverses in
more than 50% but not all of the High samples <!-- N-74 -->, ranked by the
shift in median r. The final panel was selected greedily in rank order,
excluding gene reuse and pairs whose per-sample r profiles correlate at
Spearman ≥ 0.75 with a kept pair, to a target of 10 pairs
<!-- N-75, N-37 -->. From 317 up- and 182 down-genes, 57,694 pairs were
evaluated and 153 passed all criteria <!-- N-35 -->, with median shifts of
1.159–4.700 and reversal rates 0.53–0.87 <!-- N-36 -->; the panel and its
training-derived classification boundary (positive at score > 2, on the R0
scale) are given in Table(仮) <!-- N-37, N-39 -->.

<!-- C-08 の手法側 -->
The finalized panel was then applied, untouched, to the intermediate-band
RET tumours it was never trained on (R_Low, R_Mid <!-- N-10 -->) as a
graded, out-of-sample descriptive check; the single pre-specified test is
one-sided Brunner–Munzel for Mid > Low reversal scores (Monte Carlo, seed
19860426 <!-- N-75 -->). This evaluation does not alter the panel or its
boundary, and the band-wise score profile is reported descriptively without
assuming any dose–response form.

## REO evaluation diagnostics

<!-- 写像元: diagnostics/reo_lowmid_{outliers,purity,confound}.R -->

<!-- C-09 -->
Because the evaluation cohort is deliberately unfiltered, its QC is reported
as diagnostics with no exclusion authority: (i) the training outlier screen
(PC-OD) was mirrored on the R_Low/R_Mid tumours, and the evaluation was run
once on all eligible cases — the screen reports counts only <!-- N-43 -->;
(ii) tumour purity for the evaluation bands was estimated on the same common
scale as training by pooling the whole RET cohort in one ContamDE run
<!-- N-44 -->; <!-- C-12 --> (iii) the association of band with reversal
score controlling for purity was assessed by partial Spearman correlation
with a permutation reference, alongside purity-stratified one-sided
Brunner–Munzel tests of Mid > Low <!-- N-46, N-47 -->.

## Between-arm age difference

<!-- 写像元: diagnostics/age_arm_difference.R -->

<!-- C-15 -->
The between-arm difference in age at surgery (High − Sporadic) was estimated
in each driver stratum by the Hodges–Lehmann median difference and the
rank-based effect P(Sporadic < High) — the same Brunner–Munzel effect
estimator used throughout — each with a 95% percentile bootstrap confidence
interval (within-arm resampling, 9,999 replicates, seed 19450809
<!-- N-63 -->). This is a disclosure of magnitude and uncertainty, not a
confounding test; no p-value is computed. Age at exposure is structurally
not comparable between arms (the Sporadic arm is unexposed) <!-- N-63 -->.

## External anchor cross-reference

<!-- 写像元: diagnostics/external_gene_anchors.R -->

<!-- C-13 -->
Externally validated radiation-associated gene lists (qRT-PCR-validated
cores of Abend 2013, Abend 2012, Dom 2012, and CLIP2 <!-- N-53 -->) were
cross-referenced against each unit's q < 0.10 gene set as a descriptive
membership count (k of n list genes among the unit's tested genes), with no
enrichment statistic: the lists are small, platforms and contrasts differ
across sources, and the reading was fixed symmetrically in advance — any
count is reported as description and no claim moves with the outcome.

## Software, seeds and reproducibility

<!-- 写像元: config.R, tests/, run 保全記録 -->

The publication run executed R 4.5.3 on Ubuntu 24.04 <!-- N-02 --> with the
reference BLAS/LAPACK 3.12.0 <!-- N-03 -->, four workers and 9,999 label
shuffles as the fixed reproduction contract <!-- N-04 -->, from a clean
repository state <!-- N-01 -->. The canonical inference seed is 19860426;
diagnostics draw documented seeds from base 19450809 <!-- N-05, N-06 -->.
The full pipeline was executed independently on two machines: the 1,819 raw
input files agree by md5 <!-- N-52 --> and the pipeline's test suite passes
on both (415 tests, 0 failures <!-- N-51 -->); primary artifacts were
verified identical across machines. Analysis code, versioned inputs, and the
complete decision record (protocol amendments with their timing) accompany
the paper.

---

## 起草メモ(本文に載せない)

- **未タグの数値ゼロを目標に、設定値・規則を N-66〜N-76 として台帳化した**
  (numbers_ledger セクション P。出典は全てコミット済みコードの行参照、計算なし)。
- D6 段落と ORA 段落は正本(draft_methods_results_d6.md / draft_ora_annotation.md)
  からの再掲。編集はあちらで行い、本稿へは機械的に同期する(検査トリガー1で照合)。
- 「Analysis units」節は意図的に A/B 中立・機構言及なしの定義のみ。仮説の文言
  (radiogenic/sporadic の二値観を含む)は研究者領分のため Intro/Discussion 側に置く
  想定で、Methods には書いていない。
- 530 の記述は「graded, descriptive, no dose–response form assumed」で§0.5 第1回の
  水準制約に合わせた。「if the panel captures a radiation-attributable signal…」型の
  予測文(530 ヘッダにある)は Methods に載せず、載せるなら Intro の予告側。
- 純度プーリング(研究者決定 2026-08-13、同日改訂): script ヘッダ・Methods とも
  **理論構造で記述** — 相対尺度の共通化(腕別 run は単一閾値の意味を壊す。腕別で
  先に回して観察された事実)、共通参照仮定は driver 層別設計の前提に乗る、純度軸は
  曝露対比と別軸、役割は相対フィルタ+診断共変量のみ。設計時実測(0.93–0.99 / 0.99)は
  **二重統計のためライセンスに使わない** — run コミット(8eed384)の凍結ヘッダと
  N-76 行(不使用)にのみ保存。
- 手術年齢の共変量非投入の理由文(C-10 段落)は Q-03 の要旨の圧縮。年齢層別診断を
  しない理由の詳細は Discussion/limitation 側(判断点4の残り)に置く想定。
- IREP の入力規約は 130 ヘッダの写し(N-68)。データ提供側の線量そのものの来歴は
  REBC-THYR 原論文への引用で受ける(引用は Intro/Methods 冒頭、書誌は研究者)。
- Table/Supplementary 番号は全て(仮)。図表構成の確定(判断点5)後に振り直す。
- リント自己検査(10項目): (1) unit 横断 FDR 主張なし(「no study-wide FDR is
  claimed」と明文化)/ (2) A/B 中立(機構文言なし)/ (3) 線量依存の仮定なし
  (530 は descriptive・no dose–response form)/ (4) q<0.25 不使用 / (5) 特異性
  確認の文言なし / (6) 探索セル由来なし(ORA は hypothesis-generating と明示)/
  (7) 年齢の p>0.05 型議論なし(推定+CI のみ、p 値なしを明記)/ (8) 全数値
  N タグ済み(新規設定は N-66〜N-76)/ (9) WY-FWER 不使用 / (10) 主張水準は
  手続き記述のみで超えない。
