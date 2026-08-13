# ORA 注釈(C-14)Results / Discussion ドラフト

状態: draft(Claude Code 起草 2026-08-12、研究者が書き直す前提)。
写像元: diagnostics/deg_ora_annotation.R(+同 output/ の log・rds)。
数値は numbers_ledger N-59〜N-62 経由。水準 = 仮説生成(C-14 で実行前に固定)。
Methods 側の記述(手法1文)は 420/D6 の Methods に併置する想定で本稿には含めない
(必要文言は起草メモ参照)。

---

## 本文 Results(短報、2文)

<!-- C-14 -->
As a descriptive annotation of the discovered list (hypothesis-generating;
Methods), the 794 genes lower in the High arm were strongly concentrated in
proliferation, cell-cycle and DNA-repair programs (E2F targets 46/199, G2M
checkpoint 41/198, Reactome DNA repair 68/322 <!-- N-61, N-62 -->), a single
theme that extends into the radiation-curated family, whose leading flagged
sets are themselves cell-cycle genes responding to irradiation (46/126,
expected 14.2 in the combined list <!-- N-60 -->). The 971 genes higher in
the High arm showed no such concentration in any curated family
<!-- N-59 --> (full results in Supplementary Table SY).

## Discussion(仮説生成段落、4要素)

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
in age structure and latency between the arms cannot be decided here
[年齢の受け — 判断点4の決定後に確定]. No claim of the primary analyses
rests on this annotation.

---

## 起草メモ(本文に載せない)

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
- 年齢・潜伏期の受け(第2回・Q-03・判断点4): 年齢層別診断を追補計算で実施するか、
  記述+limitation で受けるかの決定待ち。決定後に Discussion の [ ] を埋める。
- Supp Table SY = Supp.Tab.2(仮)(図表台帳)。D6 の Supp Table SX と番号整合は
  図表構成確定時に振り直し。
- リント自己検査: 仮説生成の水準明示あり / A/B 中立(機構主張なし)/ 線量依存なし /
  q<0.25 不使用 / unit 横断 FDR なし / 全数値 N タグ(NES 個別引用は保留)/
  C-14 タグ済み。
