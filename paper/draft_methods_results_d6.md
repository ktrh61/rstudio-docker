# D6(セットレベル推論の帰無較正)Methods / Results ドラフト — 二層構成

状態: draft(Claude Code 起草 2026-08-12、二層化 同日。研究者が書き直す前提の下書き)。
写像元: diagnostics/gsea_null_calibration.R(ヘッダ含む)、scripts/420_test_gene_sets.R、
diagnostics/gsea_spikein_control.R。数値は全て numbers_ledger の N-ID 経由(タグ併記)。
構成: 本文は「較正の存在と4事実」のみ(約120語)、機構の説明はサプリへ。
Supp Table SX = D6 16 セル表(図表台帳 Supp.Tab.1(仮)、N-24)。

---

## 本文 Methods(圧縮版)

<!-- C-06 -->
The complete gene-set inference was calibrated on held-out null replicates —
label permutations pushed through the identical procedure — before being
applied to the real contrasts, and the choice of the set-level FDR procedure
itself was fixed by this calibration prior to any real-data run
(Supplementary Methods) <!-- N-06 -->. A spike-in control (one Hallmark set
× 1.15 in the 9 samples of one arm) served as the sensitivity counterpart
<!-- N-30 -->.

## 本文 Results(圧縮版)

<!-- C-06 -->
Under null inputs, the set-level machinery produced at least one discovery in
102 of 1,600 replicates pooled across the 16 unit × collection cells (0.064;
nominal 0.10 <!-- N-56 -->), with a single disclosed excess (B_Normal /
Hallmark: 0.18, 95% CI 0.110–0.270 <!-- N-25 -->); the spiked set was
recovered at rank 1 of 50 (q 0.0101 <!-- N-31 -->) with no other set at
q < 0.10 <!-- N-32 --> (per-cell calibration in Supplementary Table SX).

---

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
same pre-real-data setting <!-- N-57 -->. The replacement is recorded as a
protocol amendment with the state of knowledge at the time of the change; the
superseded measurement logs are retained with the analysis code.

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

## 起草メモ(本文に載せない)

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
