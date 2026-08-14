# 原稿統合ドラフト(draft_manuscript — 論文へ収束する単一ファイル)

状態: draft(2026-08-14 統合。研究者指示: Methods/Results の対構造の検査可能化と
同期機構の廃止のため、旧4ドラフト — draft_methods / draft_results /
draft_methods_results_d6 / draft_ora_annotation — を本ファイルへ吸収し削除。
内容の変更なし、移設+Discussion 骨組み追加のみ。履歴は git)。
規約: 英語本文が正本。Results の【訳】は研究者確認用の対訳(正本でない)。
数値は N-ID タグ、主張段落は C-ID タグ(検査対象は英語本文)。図表番号は全て(仮)。
各節の由来・決定は末尾の「起草メモ(統合)」。

---

# Introduction

<!-- placeholder: 仮説の文言化(バイナリ観の定式化 — 研究者領分)待ち。
claim_map の「Intro予告」列(現状ほぼ未定)をここで確定する。
予測マップ+解釈規則の論文内提示(Methods 内小表か Intro か)は判断点5と同時に決定。
C-07 の予告(Abend 2013・Ory 2026 の正常組織被曝記憶+「Ory は driver 非層別」の
隙間提示)は claim_map C-07 行に下書きあり -->

# Methods


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
19860426 <!-- N-04, N-05 -->). The Brunner–Munzel contrast is the
protocol-wide two-sample statistic by design: the analysis takes the order
relation as its primitive — the same commitment that selects the rank-based
REO panel and the tie-invariant enrichment statistic below — its effect
P(X<Y) reads directly as an exceedance probability, and at these arm sizes
the permutation distribution is exactly enumerable, so gene-level inference
carries neither a count-model assumption nor Monte Carlo error.
Gene-level inference is the Storey q-value on
the exact p-values with the plug-in pi0 estimate at λ = 0.5, thresholded at
q < 0.10 protocol-wide <!-- N-05 -->. The estimator was fixed a priori from
the working hypothesis and the design rather than tuned: the hypothesis is a
weak signal spread across many genes — small per-gene effects attenuated by
within-arm case mixture — so holding π0 at 1 would build the absence of
exactly that signal into the correction; the conservativeness of the fixed-λ
plug-in requires only marginal uniformity of the null p-values, which the
exact test guarantees regardless of the dependence between genes, whereas
adaptive choices of λ carry guarantees that assume independent or weakly
dependent tests and do not apply to co-expressed genes sharing one label
vector; and under a weak spread signal the alternative p-density is nearly
flat, so raising λ buys little bias reduction at first-order variance cost —
λ = 0.5 is the untuned Storey (2002) default. pi0 is reported with
permutation-calibrated uncertainty. A unit-level omnibus permutation test
accompanies the gene table and carries the unit-level inferential claim —
the question of whether a unit carries any signal is answered here, not by
the size of its gene list; Higher Criticism (α0 = 0.1) is the pre-specified
primary omnibus statistic, chosen a priori for its sensitivity to many weak
effects (Donoho & Jin 2004), with count- and max-type rows reported
descriptively <!-- N-05 -->. The full rejection curve R(α) is retained so
that results do not depend on displaying a single threshold count.

<!-- C-16 -->
No unit-level binary significance label is assigned: unit-level evidence is
reported continuously — the per-gene q-values, the pre-specified omnibus p
and the rejection curve — following the methodological guidance against
dichotomizing evidence near a threshold (Wasserstein & Lazar 2016; Greenland
et al. 2016; Amrhein et al. 2019), with the interpretation of each outcome
pattern pre-assigned before the results were seen rather than left to
post-hoc labeling (interpretation map, Methods).

## Gene-set level inference

<!-- 写像元: 420, 430(ORA 注釈), lib/gsea_permutation.R, lib/gsea_collections.R -->

Set-level enrichment consumes the whole per-gene ranking (threshold-free; no
DEG-list cut decides what is tested), ranked by tie-averaged normal scores of
the signed Brunner–Munzel statistic. The enrichment score is the weighted
running sum (gseaParam = 1) evaluated at tie-block boundaries only, so that
no arbitrary tie-break injects an order the data do not contain; on tie-free
input it equals the standard GSEA statistic <!-- N-72 -->. The null
reuses the identical label shuffles saved by the gene-level test (9,999 per
unit <!-- N-07 -->), so gene- and set-level results rest on one permutation
null. Inference is the per-set, sign-conditional permutation p-value with
Benjamini–Hochberg adjustment within each collection, q_bh < 0.10
<!-- N-72 -->; no cross-family claim is made. Unlike the gene level, π0 is
held at its conservative bound of 1 here: tens to thousands of dependent
p-values riding one collective drift mode leave the plug-in estimate
variance-dominated, with errors pointing anticonservative exactly in the
realizations where the procedure is most fragile.

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
(サプリ層 = Supplementary Methods / Results 節)

<!-- C-14 -->
As a descriptive complement, the discovered R_Tumor list was annotated by
one-sided hypergeometric over-representation against the identical set
universe (up / down / combined lists; universe = the unit's tested genes;
BH within family × list), reported in full as hypothesis-generating
annotation (Supplementary Table SY) <!-- N-59 -->. The two procedures'
p-values refer to different randomness: only the label-permutation null
refers to the experiment actually performed, and experiment-level claims
rest on it; the over-representation q-values describe the discovered list
against a gene-sampling reference (Goeman & Bühlmann 2007).

## Between-arm concordance of the exposure contrast

<!-- 写像元: diagnostics/signature_agreement.R(節題は判断点2決着に従い記述的言い回し。
「signature agreement」は内部語として本文不使用) -->

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
graded, out-of-sample check. Because both training bands entered panel
construction, any separation within them is circular; Low and Mid are the
only bands where the score can be examined out of sample, and the design
fixes the expected direction in advance (the higher assigned-share band
above the lower). The single pre-specified test is therefore the one-sided
Brunner–Munzel comparison of Mid over Low (Monte Carlo, seed 19860426
<!-- N-75 -->). This evaluation does not alter the panel or its
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
<!-- N-44 -->; <!-- C-12 --> (iii) because the reversal score correlates
with purity within the evaluation bands, a graded band profile could in
principle be purity-driven rather than band-driven — the two are separated
by the partial Spearman correlation of band with score conditioning on
purity (permutation reference) and by purity-stratified one-sided
Brunner–Munzel comparisons of Mid over Low <!-- N-46, N-47 -->.

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
verified identical across machines. Analysis code and versioned inputs
sufficient to regenerate the reported analyses accompany the paper.


# Results


## 1. Cohort(C-10, C-15)

<!-- C-10 -->
Of the 440 cases of the REBC-THYR cohort (Morton et al. 2021) <!-- N-11 -->,
the main analysis cohort comprised 63 paired, driver-stratified cases —
9 B_High, 27 B_Sporadic, 15 R_High and 12 R_Sporadic <!-- N-09 --> — reached
through the pre-specified flow of driver classification, band eligibility,
pairing, outlier screening and purity thresholding (440 → 248 → 77 → 70 →
69 → 63 <!-- N-08 -->; Tables 1–2(仮)). The largest reductions are driver
classification (440 → 248) and the restriction to the Sporadic and High
bands (248 → 77) <!-- N-08 -->; driver-specific attrition is shown in the
cohort-flow figure(仮). The REO evaluation set added 36
paired RET tumours of the intermediate bands (17 R_Low, 19 R_Mid
<!-- N-10 -->). In both driver strata the High arm sat somewhat older at
surgery than the Sporadic arm (RET: median 23 [range 14–31] vs 20.5 [14–27]
years; BRAF: 29 [19–39] vs 24 [11–29]) <!-- N-12, N-13 -->.
<!-- C-15 -->
The between-arm age difference is disclosed as an estimate rather than
tested: Hodges–Lehmann +2.5 years (95% bootstrap CI −1.0 to 6.0;
P(Sporadic<High) 0.625 [0.400–0.828]) in the RET stratum <!-- N-64 --> and
+8.0 years (3.0 to 12.0; 0.850 [0.681–0.973]) in the BRAF stratum
<!-- N-65 -->. Age at exposure exists only for exposed cases and admits no
between-arm comparison <!-- N-63 -->.

【訳】REBC-THYR コホート(Morton et al. 2021)の 440 症例 <!-- N-11 --> のうち、
主解析コホートは、driver 分類・帯適格性・ペア有無・外れ値スクリーン・純度閾値という
事前指定のフロー(440 → 248 → 77 → 70 → 69 → 63 <!-- N-08 -->; 表1–2(仮))を経た
63 のペア付き driver 層別症例 — B_High 9・B_Sporadic 27・R_High 15・R_Sporadic 12
<!-- N-09 --> — である。REO 評価セットとして中間帯のペア付き RET 腫瘍 36 例
(R_Low 17・R_Mid 19 <!-- N-10 -->)を加えた。いずれの driver 層でも High 腕は
Sporadic 腕よりやや高い手術時年齢に分布した(RET: 中央値 23 [範囲 14–31] 対
20.5 [14–27] 歳; BRAF: 29 [19–39] 対 24 [11–29])<!-- N-12, N-13 -->。
腕間年齢差は検定でなく推定として開示する: Hodges–Lehmann 中央値差は RET 層で
+2.5 年(95% ブートストラップ CI −1.0〜6.0; P(Sporadic<High) 0.625 [0.400–0.828])
<!-- N-64 -->、BRAF 層で +8.0 年(3.0〜12.0; 0.850 [0.681–0.973])<!-- N-65 -->。
被曝時年齢は被曝症例にのみ定義され、腕間比較は成立しない <!-- N-63 -->。

## 2. Gene-level differential expression(C-01〜C-04, C-16, 410)

<!-- C-01 -->
In R_Tumor — the unit where the pre-specified prediction map expects the
signal — 1,765 of 15,621 tested genes differed between the High and Sporadic
arms at Storey q < 0.10 <!-- N-16, N-15 -->, and the pre-specified unit-level
omnibus supported the presence of signal (Higher Criticism p = 0.0112
<!-- N-20 -->) (Figs. 1–2(仮), Table 3(仮)). <!-- C-02 --> The discovered
genes ran in both directions: 971 higher and 794 lower in the High arm
<!-- N-17 -->. <!-- C-03 --> In R_Normal no gene reached q < 0.10 and the
omnibus lent no support (HC p = 0.3199) <!-- N-16, N-20 -->. <!-- C-04 -->
B_Tumor likewise yielded no discovery (0 genes at q < 0.10; HC p = 0.1815)
<!-- N-16, N-20 -->; under the pre-fixed reading this cell is
direction-agnostic, and its quiet is not read as a specificity control.
<!-- C-16 -->
In B_Normal the evidence is reported as it stands, without a binary label:
one gene crossed the gene-level threshold (BHLHB9, effect 0.967, q = 0.013
<!-- N-22 -->); the pre-specified primary omnibus gave HC p = 0.0773
<!-- N-20 --> while the descriptive max-statistic row reached p = 0.0125
<!-- N-23 -->, and the unit's π0 estimate (0.727) sat below those of the
other two quiet units (0.955, 0.943) <!-- N-19 -->. The reading pre-assigned
to this pattern is taken up in Discussion.

【訳】事前指定の予測マップがシグナルを期待する unit である R_Tumor では、検定対象
15,621 遺伝子のうち 1,765 が Storey q < 0.10 で High 腕と Sporadic 腕の間で発現差を
示し <!-- N-16, N-15 -->、事前指定の unit レベル・オムニバスがシグナルの存在を支持した
(Higher Criticism p = 0.0112 <!-- N-20 -->)(図1–2(仮)・表3(仮))。発見遺伝子は双方向に
分布した: High 腕で高発現 971・低発現 794 <!-- N-17 -->。R_Normal では q < 0.10 の
遺伝子はなく、オムニバスの支持もなかった(HC p = 0.3199)<!-- N-16, N-20 -->。B_Tumor も
同様に発見なし(q < 0.10 は 0 遺伝子; HC p = 0.1815)<!-- N-16, N-20 --> — 事前固定の
読みに従い、このセルは方向不可知であり、その静けさを特異性の対照とは読まない。
B_Normal の証拠は2値ラベルなしにそのまま報告する: 1遺伝子が遺伝子レベル閾値を越え
(BHLHB9、effect 0.967、q = 0.013 <!-- N-22 -->)、事前指定の主オムニバスは HC p = 0.0773
<!-- N-20 -->、記述的な max 統計量行は p = 0.0125 <!-- N-23 -->、この unit の π0 推定値
(0.727)は他の2つの静かな unit(0.955・0.943)より低かった <!-- N-19 -->。このパターンに
事前割当された読みは Discussion で扱う。

## 3. Gene-set level(C-05, C-06, 420 + D6)

<!-- C-05 -->
At the gene-set level, the calibrated test declared no set at q_bh < 0.10 in
any of the 16 unit × collection cells <!-- N-27 --> — Hallmark (50 sets),
C2:CP (3,910; the union of Reactome, WikiPathways, KEGG MEDICUS, BioCarta
and PID), C5:GO:BP (7,538) and the radiation-curated C2:CGP family (28);
6,141–6,242 sets per unit after filtering <!-- N-29, N-55, N-28 -->. The
smallest adjusted value anywhere was q = 0.114, in the B_Tumor × radiation
cell <!-- N-28 --> (complete listing in Supp. Data 1(仮); Table 4(仮)).
<!-- C-06 -->
Under null inputs, the set-level machinery produced at least one discovery
in 102 of 1,600 replicates pooled across the 16 unit × collection cells
(0.064; nominal 0.10 <!-- N-56 -->), with a single disclosed excess
(B_Normal/Hallmark: 0.18, 95% CI 0.110–0.270 <!-- N-25 -->); the spiked set
was recovered at rank 1 of 50 (q 0.0101 <!-- N-31 -->) with no other set at
q < 0.10 <!-- N-32 --> (per-cell calibration in Supp. Tab. 1(仮)).

【訳】遺伝子セットレベルでは、較正済み検定は 16 の unit × collection セルのいずれでも
q_bh < 0.10 のセットを宣言しなかった <!-- N-27 --> — Hallmark(50 セット)・C2:CP
(3,910; Reactome・WikiPathways・KEGG MEDICUS・BioCarta・PID の統合)・C5:GO:BP(7,538)・
放射線キュレーションの C2:CGP ファミリー(28)、フィルタ後は unit あたり 6,141–6,242
セット <!-- N-29, N-55, N-28 -->。全セルで最小の調整値は B_Tumor × radiation セルの
q = 0.114 だった <!-- N-28 -->(全結果は Supp. Data 1(仮)・表4(仮))。帰無入力の下では、
セットレベル機構は 16 セル合算 1,600 レプリケート中 102 で1つ以上の発見を生じ
(0.064; 名目 0.10 <!-- N-56 -->)、開示済みの超過は1セルのみ(B_Normal/Hallmark: 0.18、
95% CI 0.110–0.270 <!-- N-25 -->)。spike-in セットは 50 セット中 rank 1(q 0.0101
<!-- N-31 -->)で回収され、それ以外に q < 0.10 のセットはなかった <!-- N-32 -->
(セル別較正は Supp. Tab. 1(仮))。

## 4. Composition of the discovered list(C-14, 430)

<!-- C-14 -->
As a descriptive annotation of the discovered list (hypothesis-generating;
Methods), the 794 genes lower in the High arm were strongly concentrated in
proliferation, cell-cycle and DNA-repair programs (E2F targets 46/199, G2M
checkpoint 41/198, Reactome DNA repair 68/322 <!-- N-61, N-62 -->), a single
theme that extends into the radiation-curated family, whose leading flagged
sets are themselves cell-cycle genes responding to irradiation (46 of 126 in
the down list, expected 6.4 <!-- N-60 -->). The 971 genes higher in the High
arm showed no such concentration in any curated family <!-- N-59 --> (full
results in Supp. Tab. 2(仮)).

【訳】発見済みリストの記述的注釈(仮説生成; Methods)として、High 腕で低発現の 794
遺伝子は増殖・細胞周期・DNA 修復プログラムに強く集中し(E2F targets 46/199・
G2M checkpoint 41/198・Reactome DNA repair 68/322 <!-- N-61, N-62 -->)、この単一テーマは
放射線キュレーション・ファミリーにも及ぶ — そこでフラグが立った上位セット自体が
照射に応答する細胞周期遺伝子である(down リストで 126 中 46、期待 6.4 <!-- N-60 -->)。
High 腕で高発現の 971 遺伝子には、どのファミリーにもそのような集中はなかった
<!-- N-59 -->(全結果は Supp. Tab. 2(仮))。

## 5. Between-arm concordance of the exposure contrast(C-07, C-17)

<!-- C-07 -->
The pre-specified between-arm comparison in normal tissue — Spearman
correlation, across arms, of the per-gene signed statistics of the exposure
contrast — gave rho = +0.376 over 15,459 shared genes, inside its
label-shuffle reference band ([−0.46, +0.46]; two-sided p = 0.1199)
<!-- N-34 -->. With no within-unit signal in R_Normal <!-- N-16 -->, the
pre-fixed reading applies: whether the two normal-tissue contrasts share an
exposure trace is not identifiable here. <!-- C-17 --> The symmetric
tumour-pair comparison, computed as a design completion, is reported in
Supp. Tab. 3(仮) <!-- N-33 --> and taken up as hypothesis-generating in
Discussion.

【訳】正常組織における事前指定の腕間比較 — 曝露対比の遺伝子別符号付き統計量を
腕間で Spearman 相関 — は、共有 15,459 遺伝子で rho = +0.376 となり、ラベル
シャッフルの参照帯([−0.46, +0.46]; 両側 p = 0.1199)の内側だった <!-- N-34 -->。
R_Normal に unit 内シグナルがない <!-- N-16 --> ため、事前固定の読みが適用される:
2つの正常組織対比が曝露痕跡を共有するか否かは、ここでは識別できない。対称補完として
計算した腫瘍ペアの比較は Supp. Tab. 3(仮)に報告し <!-- N-33 -->、Discussion で
仮説生成として扱う。

## 6. REO grading(C-08, C-09, C-12, 510–530)

<!-- C-08 -->
The 10-pair REO panel separated its training arms as designed (all 12
R_Sporadic negative, 13 of 15 R_High positive; boundary at score > 2
<!-- N-38, N-37 -->; panel composition in Table(仮) <!-- N-39 -->). Applied
unchanged to the intermediate bands, the per-case reversal score rose in
band order — medians 0 / 1 / 4 / 6 for Sporadic / Low / Mid / High
<!-- N-41 --> (Fig. 3(仮)) — with R_Low classified 9 negative / 8 positive
and R_Mid 8 / 11 <!-- N-42 -->. The one comparison available out of sample
(Mid vs Low; direction pre-specified, Methods) gave one-sided
Brunner–Munzel p = 0.1127, effect P(Low<Mid) = 0.616 <!-- N-40 -->. The
graded profile is reported as a descriptive observation; no dose–response
form is assumed or claimed.
<!-- C-09 -->
The outlier screen mirrored from training flagged no case in either
evaluation band (0 of 17 and 0 of 19 <!-- N-43 -->); the screen carries no
exclusion authority, and the evaluation was run once on all eligible cases.
<!-- C-12 -->
Within the evaluation bands the reversal score correlated with tumour purity
(pooled Spearman +0.538 <!-- N-45 -->), while the band–score correlation was
+0.142 and, conditioned on purity, +0.146 (partial Spearman; one-sided
permutation p = 0.2162) <!-- N-48, N-46 -->; band and purity themselves were
nearly uncorrelated (+0.036 <!-- N-48 -->). Purity-stratified comparisons
are given in Supplementary <!-- N-47 -->.

【訳】10 ペアの REO パネルは訓練腕を設計どおり分離した(R_Sporadic は 12 例全て
negative、R_High は 15 例中 13 が positive; 境界は score > 2 <!-- N-38, N-37 -->;
パネル構成は表(仮)<!-- N-39 -->)。パネルを変更せず中間帯に適用すると、症例別の
逆転スコアは帯順に上昇した — 中央値は Sporadic / Low / Mid / High で 0 / 1 / 4 / 6
<!-- N-41 -->(図3(仮))— 分類は R_Low が negative 9 / positive 8、R_Mid が 8 / 11
<!-- N-42 -->。out-of-sample で検査可能な唯一の比較(Mid 対 Low; 方向は事前指定、
Methods 参照)は片側 Brunner–Munzel p = 0.1127、効果 P(Low<Mid) = 0.616 だった
<!-- N-40 -->。この段階的プロファイルは記述的観察として報告し、線量反応の形は
仮定も主張もしない。
訓練側から鏡映した外れ値スクリーンはどちらの評価帯でも該当例を出さなかった
(17 例中 0・19 例中 0 <!-- N-43 -->)— スクリーンは除外権限を持たず、評価は適格
全例で一度だけ実施した。評価帯内では逆転スコアは腫瘍純度と相関し(pooled Spearman
+0.538 <!-- N-45 -->)、帯–スコア相関は +0.142、純度を条件付けると +0.146(偏 Spearman;
片側置換 p = 0.2162)<!-- N-48, N-46 -->、帯と純度自体はほぼ無相関だった(+0.036
<!-- N-48 -->)。純度層別の比較は Supplementary に示す <!-- N-47 -->。


# Discussion

<!-- 執筆状態: 実体1段落(C-14)+骨組み。本文化は Results の研究者書き直しと
仮説文言化の後。見出しと順序は仮 — 各見出しは claim_map の「Disc受け」列に対応 -->

## 主結果の受け(C-01〜C-05)【骨組みのみ】

<!-- C-01〜C-05 の Disc 受け。未起草 -->

## B_Normal — 事前固定マップの条件適用(C-16)【骨組みのみ】

<!-- C-16: 規則5の条件適用(共有腺痕跡で説明困難 → 交絡第一・年齢筆頭 —
N-65 + Coclet の正常組織側。B_Tumor null との対比で「年齢は正常組織に効き
腫瘍に効かない」に束ね Q-13 (iv) と相互強化)。D6 超過 0.18(N-25)の受けもここ。
「q<0.10 陽性論」への受けは Q-16 -->

## 年齢(C-15、Q-13 (i)〜(iv))【骨組みのみ】

<!-- Q-13 の4層((i) 組織区画の訂正+Vriens / (ii) 帯順序の反転 / (iii) 純度の逆行 /
(iv) 内部対照)+ C-15 の開示推定(N-64, N-65)。Q-03(共変量非投入)の詳細もここか
limitation -->

## リスト構成注釈の読み(C-14)


<!-- C-14 -->
The over-representation annotation of the R_Tumor list is read under a
gene-sampling null that ignores co-expression and is anticonservative
relative to the calibrated set-level inference; its flags are therefore
candidates supported only under that weaker null. Three observations bound
what it can mean. First, the flags across all four families express one
theme — relative attenuation of proliferation, cell-cycle and DNA-repair
programs in the High arm — rather than several independent lines of
evidence; the radiation-curated flags in particular ride the same cell-cycle
content. Second, the calibrated set-level test ranked the same Hallmark sets
first without crossing its threshold <!-- N-28 -->, so the annotation and
the calibrated inference disagree about evidence, not about ordering. Third,
whether this attenuation reflects exposure-associated biology or differences
in age structure and latency between the arms cannot be decided here; the
arms' age difference is, however, disclosed with interval estimates
(Hodges–Lehmann +2.5 years [−1.0, 6.0] in this stratum <!-- N-64 -->), and a
published age contrast an order of magnitude larger produced no genome-wide
differential expression and no set-level findings in papillary thyroid
carcinoma (Vriens et al. 2011). No claim of the primary analyses rests on
this annotation.


## 腕間 concordance — tumor ペア(C-17)【骨組みのみ】

<!-- C-17: driver 横断の曝露関連成分と整合/共有共変量(純度筆頭・年齢)と識別不能/
核仮説の検定には不使用(B.7 整合)。隣接断片 = B_Tumor×radiation min q 0.114(N-28)、
いずれも較正済み閾値未達を併記。N-33 -->

## 正常組織の対照 — Abend/Ory の受け(C-07)【骨組みのみ】

<!-- C-07 Disc 受け(claim_map 行に下書き): Ory 2026 対照段落 — 共有記憶の主張を
driver 層別下で検証可能にする装置がこの比較であり、結果は識別不能に終わった -->

## Limitations(C-11 ほか)【骨組みのみ】

<!-- C-11: A/B 非識別性の明文化(§0.5 第6回)。ほか: ウェット追試不能・REO の
中間データ依存(§0.5 第1回の控えめ制約)、被曝時年齢の腕間比較不能(N-63) -->

# Supplementary Methods / Results

## Supplementary Methods — Calibration of the gene-set level inference

<!-- C-06 -->
Because gene-set q-values inherit the dependence structure of the expression
matrix — genes are co-expressed, and sets overlap — nominal false discovery
rates at the set level cannot be taken on faith. We therefore measured the
operating characteristics of the entire set-level procedure on held-out null
replicates before it was applied to the real contrasts, and we report that
measurement alongside the results it calibrates.

For each analysis unit, exposure labels were shuffled B + R times with a seed
independent of the one used for inference (base seed 19450809; B = 9999,
R = 100 <!-- N-06 -->), and the block enrichment statistic was computed once
for all shuffles. The first B shuffles form a shared null pool; each of the
R remaining shuffles was then treated as a pseudo-observation and pushed
through exactly the inference applied to the real data — normalized enrichment
against the pool, per-set permutation p-values, Benjamini–Hochberg adjustment
within each collection, and the pre-specified threshold q < 0.10
<!-- N-06 -->. Under label exchange the pseudo-observations are exchangeable
with the pool, so for every unit × collection cell the share of replicates
yielding at least one discovery estimates P(≥1 false discovery) under the
complete null, where FDR and family-wise error coincide. We report this share
with an exact binomial confidence interval and the mean discovery count per
replicate (Supplementary Table SX). Because the R replicates share one null
pool, their discovery indicators are weakly positively correlated and the
binomial interval is accordingly somewhat narrow; we record this rather than
correct it.

The calibration also fixed the choice of the set-level inference itself,
before any real-data run. The tail-ratio FDR estimated from pooled normalized
enrichment scores — the inference originally specified — was measured
miscalibrated on these data (pooled P(≥1) 0.140, and 0.221 for a
restandardized variant, against the nominal 0.10; worst cell 0.44
<!-- N-58 -->), and was replaced by the per-set permutation p with
within-collection Benjamini–Hochberg adjustment, which measured 0.045 in the
same pre-real-data setting <!-- N-57 -->.

## Supplementary Results — 較正表の読み(Supp Table SX に併記する注記)

<!-- C-06 -->
Across the 16 unit × collection cells, the share of null replicates producing
at least one discovery at q < 0.10 ranged from 0.01 to 0.18 <!-- N-24 -->,
and 13 of 16 cells were at or below the nominal 0.10. Two further cells in
the radiation-curated family straddled it (B_Tumor 0.10, CI 0.049–0.176;
B_Normal 0.12, CI 0.064–0.200 <!-- N-26 -->). Together with the spike-in
recovery (NES 2.28, p 0.0002 <!-- N-31 -->), the two controls bound the
machinery from both sides: null inputs do not generate discoveries beyond the
nominal level, and a planted coherent signal of modest size is detected. The
one calibration excess (B_Normal/Hallmark) is disclosed and is taken into
account wherever set-level results for that cell are read.


---

# 引用文献(仮書式 — 投稿書式の整合は後工程)

本文が現に引くもの:

- Morton LM, et al. Radiation-related genomic profile of papillary thyroid
  carcinoma after the Chernobyl accident. Science 2021;372:eabg2538.
  doi:10.1126/science.abg2538(コホート出典 — Results §1)
- Storey JD. A direct approach to false discovery rates. J R Stat Soc B
  2002;64:479–498(λ=0.5 の無チューニング既定 — Methods)
- Donoho D, Jin J. Higher criticism for detecting sparse heterogeneous
  mixtures. Ann Stat 2004;32:962–994. doi:10.1214/009053604000000265(HC — Methods)
- Goeman JJ, Bühlmann P. Analyzing gene expression data in terms of gene
  sets: methodological issues. Bioinformatics 2007;23:980–987.
  doi:10.1093/bioinformatics/btm051(sampling-model 軸 — Methods)
- Wasserstein RL, Lazar NA. The ASA statement on p-values. Am Stat
  2016;70:129–133. doi:10.1080/00031305.2016.1154108(非2値化 — Methods)
- Greenland S, Senn SJ, Rothman KJ, Carlin JB, Poole C, Goodman SN,
  Altman DG. Eur J Epidemiol 2016;31:337–350. doi:10.1007/s10654-016-0149-3(同上)
- Amrhein V, Greenland S, McShane B. Scientists rise up against statistical
  significance. Nature 2019;567:305–307. doi:10.1038/d41586-019-00857-9(同上)
- Vriens MR, et al. Cancer 2011;117:259–267. doi:10.1002/cncr.25369
  (年齢対比 — Discussion C-14 段落)

Methods が言及し正式引用が必要になる見込み(執筆時に確定):

- Subramanian A, et al. PNAS 2005;102:15545–15550. doi:10.1073/pnas.0506580102(GSEA)
- Brunner E, Munzel U. The nonparametric Behrens-Fisher problem. Biom J
  2000;42:17–25(BM 検定)
- MSigDB/msigdbr・ContamDE・TCC(DEGES)・MUREN・edgeR ほか実装の原典(執筆時に確定)

Intro/Discussion で必要になる見込み:

- Abend M, et al. 2012 (PLoS One)・2013 (Br J Cancer); Dom G, et al. 2012
  (Br J Cancer); Hess J, et al. 2011 (PNAS; CLIP2)— 外部アンカー(C-13)
- Ory C, et al. Sci Rep 2026. doi:10.1038/s41598-026-53030-4(正常組織対照 — C-07)
- Coclet J, et al. Clin Endocrinol 1989. doi:10.1111/j.1365-2265.1989.tb01290.x
  (正常甲状腺の細胞回転 — Q-13 (i))

---

# 起草メモ(統合 — 本文に載せない)

- **統合の記録(2026-08-14、研究者指示)**: 旧4ドラフトを本ファイルへ吸収し削除。
  各段落の実体は本ファイルの1箇所のみとなり、正本/再掲の同期機構は廃止。
  以下のサブ節は旧ファイルの起草メモをそのまま移設したもの(歴史込み — 「正本はあちら」
  等の記述は統合前の状態を指す)。

## Methods 由来(旧 draft_methods.md)


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
  しない理由の詳細は Discussion/limitation 側(Q-03/Q-13 — 判断点4は決着済み)に置く想定。
- **選定事由の記載方針(2026-08-14 確定)**: 分野の既定から外れる選択で査読者が
  「なぜ?」と聞くと予想される箇所に1文の事由を書く。慣行どおりの選択には書かない。
- **深度配分の原則(2026-08-14 確定、Mid>Low の往復から一般化)**: Results に載る全ての
  検定・read について「それが何を検定するのか」を Methods が1文で持ち、Results は
  結果+Methods 参照に留める。この原則で4点を追補: HC の役割文(unit 水準の主張を担う —
  Q-16 の水準論を本文に可視化)/ ORA の濃淡開示一文(Q-15 決着の実装)/ 節題を
  「Between-arm concordance of the exposure contrast」へ(判断点2の用語決定)/
  REO 診断 (iii) の動機文(純度駆動の可能性の分離)。
  この基準で4点を追記: BM 検定の採用(順序関係を原始とするコミットメント+完全枚挙)/
  HC の選定(多数の弱い効果への感度)/ tie-block ES の理由(順序の注入禁止)/
  セットレベル π0=1 の非対称(依存下で plug-in が分散支配・反保守側に誤る)。
  n_perm 9,999 は慣行的値のため本文で事由を書かない(D4 の導出は手元)。
- λ=0.5 の選定事由を 410 節に追記(2026-08-14、研究者指示): D1 の批准済み導出の圧縮
  (π0=1 は仮説の検出チャネルを塞ぐ/固定 λ の保守性は帰無 p の周辺一様性のみで成立/
  適応的 λ は独立性仮定で適用外/平坦な対立密度下で低 λ が MSE 優位)。弱拡散前提は
  **選定事由(設計時の作業仮説)として書く**のであって形の主張ではない — Q-15 (2)
  (形に賭けない)と整合。「small per-gene effects attenuated by within-arm case
  mixture」の文言はバイナリ観の反映だが、仮説の正式な文言化(Intro、研究者領分)と
  執筆時に整合させること。
- 非2値化の一文(C-16)を 410 節末尾に追加(2026-08-14、判断点1決着)。引用は
  ASA 声明系3本(DOI は claim_map C-16 の根拠列)。unit の陽性/非陽性ラベルは全 unit で
  不使用 — C-01 は両水準結合(DEG かつ HC)、C-03/C-04 は「検出されなかった」の記述形。
  文言検査時に positive/negative 型のラベル語が紛れていないか確認する。
- 決定履歴語の掃引(2026-08-14、研究者指示): 「protocol amendment」「decision record」への
  言及を本文から削除(査読者は決定の行き来を考慮しない。再現性文は「code and versioned
  inputs sufficient to regenerate」に限定 — 公開対象の選別は研究者裁量、§0.5b 明確化と整合)。
- **予測マップ+解釈規則は論文内に提示が必要**(正本=論文 — 外部計画への参照で済ませない。
  非2値化の一文が「interpretation map, Methods」を指すため)。置き場所(Methods 内の
  小表 or Intro)は図表構成(判断点5)と同時に確定。
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

## Results 由来(旧 draft_results.md)


- 英語版が正本、【訳】は研究者確認用の対訳(タグは両方に付けたが検査対象は英語版)。
- 再掲2段落の正本: C-06 = draft_methods_results_d6.md、C-14 = draft_ora_annotation.md。
  編集は正本側で行い本稿へ同期(トリガー1で照合)。
- **C-14 正本の数値取り違えを修正**(2026-08-14): radiation 例示が down の k(46/126)に
  combined の期待値(14.2)を付けていた → down の期待 6.4 に修正(N-60)。正本側も修正済み。
- 丸めの規約: rho は3桁表示(+0.376 = N-34 の +0.3756、帯 [−0.46, +0.46] = 同
  [−0.4615, +0.4580])。他の数値は台帳の桁のまま。凍結時に表示桁を台帳へ追記予定。
- B_Normal(C-16)は評語なしの連続量記述で、読み(規則5条件適用)は Discussion 送り。
  「他の2つの静かな unit より低い π0」は N-19 の数値の並置で、比較検定ではない。
- tumor 側 concordance(C-17)は本文ポインタ+Supp. Tab. 3(仮)+Discussion 送り
  (判断点2決着どおり)。
- ORA の可視性は現状「短報1段落」— 共著者交渉で拡張可(C-14 の可視性無制約)。
  拡張時は正本(draft_ora_annotation)側を先に改稿。
- **フロー提示の方針(2026-08-14 研究者承諾)**: 駆動遺伝子別の途中経過(N-08 の
  RET 73→…→27 / BRAF 175→…→36)は本文でなく**コホートフロー図(両層並記)**が担う。
  本文は合算チェーン+最大削減2段の半文+フロー図参照まで。図の最終形は作図時
  (判断点5)に確定。
- **Mid>Low 検定の提示(2026-08-14 研究者指摘への対応)**: 検定の問い(訓練帯は循環・
  訓練外は Low/Mid のみ・方向は設計が事前指定 → 過適合でないことの out-of-sample 検査)は
  **Methods の深度**として Methods の REO 評価段落に記載。Results は近隣と同密度の
  1句+Methods 参照のみ(深度アンマッチの回避)。検定自体は維持(事前指定 read の
  非掲載は選択的報告になるため。役割は「勾配の不確かさの開示」で、感度解析ではない)。
- **REO ペアの生物学(2026-08-14 研究者決定)**: 本文・Results ではペア名を挙げず表参照
  のみ(器械としての提示)。共著者は生物学的解釈を好む傾向 — 求めがあれば Discussion/
  共著者ラウンドで扱う余地は残すが、その際は P9/P10 が再正規化で入れ替わった経緯
  (選定順位の端の揺れ)を必ず添えること。
- 図表番号は全て(仮)。判断点5(図表構成)確定後に振り直し+図表台帳と照合。
- リント自己検査(10項目): (1) unit 横断 FDR 主張なし / (2) A/B 中立(機構文言なし)/
  (3) 線量依存の仮定なし(C-08 に明示の否定文)/ (4) q<0.25 不使用 / (5) B_Tumor の
  「特異性確認」文言なし(明示の否定文)/ (6) 探索セル産なし(ORA は hypothesis-
  generating 明示、C-17 は Disc 送り)/ (7) 年齢は推定+CI のみ(検定なしを明示)/
  (8) 全数値 N タグ済み / (9) WY 不使用 / (10) 主張水準は C 行の水準列どおり
  (主張 = C-01 のみ、他は記述的観察・仮説生成)。

## D6 二層(旧 draft_methods_results_d6.md)


- 二層の分割方針: 本文には「較正済み・pooled 0.064・超過1セル開示・spike-in 回収」の
  4事実のみ(査読者がサプリを開く前に見えるべき情報)。機構・CI 注意・選定経緯はサプリ。
- 「0.045」は**採用時測定**(2026-08-08、開発 B=999、実データ前)で、本番凍結値は
  16セル表+pooled 0.064(N-56)。前者は選定の時系列防御、後者が論文の較正値。
- WY-FWER 0.112 は N-57 に参考として台帳化したが本文・サプリとも不使用(B.2)。
- Supplementary Results 最終文の受け先: **確定(2026-08-14、判断点1決着)** — Discussion の
  B_Normal 段落(C-16 の規則5条件適用段落)の一文で受ける。B_Normal の扱いは2値ラベルなしの
  連続量記述+規則5条件適用(claim_map C-16、Q-16)。
- リント自己検査: unit 横断 FDR 主張なし / A/B 文言なし / 線量依存なし / q<0.25 不使用 /
  全数値 N タグ済み / C-06 タグ済み / WY 本文不使用。

## ORA 注釈(旧 draft_ora_annotation.md)


- **数値修正(2026-08-14、Results 起草時の照合で検出)**: Results 文の radiation 例示が
  down リストの k(46/126)に combined の期待値(14.2)を付けていた取り違えを修正 —
  正しくは down リストの k=46/126・期待 6.4(N-60。combined は k=50/126・期待 14.2)。
- **2026-08-13 来歴訂正・根拠差し替え(研究者承諾 — claim_map C-14 改訂メモ参照)**:
  ORA は Thyroid 原プロトコルの設計解析の復元であり、GSEA 主・ORA 副の濃淡の根拠は
  「事前固定」ではなく G&B 2007 の sampling-model 軸+本文での一文開示。開示文の例
  (Methods 併置文の直後を想定): "The two procedures' p-values refer to different
  randomness: only the label-permutation null refers to the experiment actually
  performed, and experiment-level claims rest on it; the over-representation
  q-values describe the discovered list against a gene-sampling reference."
  可視性は無制約 — 本稿の「短報2文」は**下限**であり、Results 段落規模・図・
  Abstract 言及への拡張は共著者交渉で可(水準ラベル携行のみ非交渉)。
  転落基準(ORA q を実験レベルの単独主張に使わない)は Q-15 (3b)。
  濃淡開示の一文は **draft_methods の C-14 段落に配置済み**(2026-08-14。
  文言は上の例文と同一、G&B 引用付き)。

- Methods 併置用の1文(420 Methods の末尾を想定):
  "As a descriptive complement, the discovered R_Tumor list was annotated by
  one-sided hypergeometric over-representation against the identical set
  universe (up / down / combined lists; universe = the unit's tested genes;
  BH within family × list), reported in full as hypothesis-generating
  annotation (Supplementary Table SY)." <!-- C-14, N-59 -->
- ガードレール3点(C-14 に記録済み): 単一テーマの明示 / 年齢・潜伏期の候補説明の
  名指し / 420 の同順位・閾値未達の併記。
- 「420 の同順位」の具体: 較正済み検定でも H ファミリー最上位は同じセット群
  (SPERMATOGENESIS p 0.0036 q 0.179、KRAS_SIGNALING_UP・E2F_TARGETS・
  G2M_CHECKPOINT が NES 負の上位 — thyr_enrichment_test.rds; min q は N-28)。
  本文で個別 NES を引く場合は N 行の追加が必要(現状は N-28 の min q のみ台帳化)。
- 年齢・潜伏期の受け: **判断点4決着により確定(2026-08-14 反映)** — 年齢層別診断は
  実施せず(Q-03・撤回済み3案は復活させない)、開示された推定(C-15、N-64)+文献側
  (Vriens 2011、Q-13 (i))で受ける。Discussion の [ ] は埋め済み(N-64 タグ)。
  潜伏期は年齢と不可分のため候補説明としての名指しのまま。
- Supp Table SY = Supp.Tab.2(仮)(図表台帳)。D6 の Supp Table SX と番号整合は
  図表構成確定時に振り直し。
- リント自己検査: 仮説生成の水準明示あり / A/B 中立(機構主張なし)/ 線量依存なし /
  q<0.25 不使用 / unit 横断 FDR なし / 全数値 N タグ(NES 個別引用は保留)/
  C-14 タグ済み。
