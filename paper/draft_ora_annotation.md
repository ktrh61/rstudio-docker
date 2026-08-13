# ORA 注釈(C-14)Results / Discussion ドラフト

状態: draft(Claude Code 起草 2026-08-12、研究者が書き直す前提)。
写像元: scripts/430_annotate_deg_ora.R(2026-08-13 に 430 編入。正準出力 =
processed/thyr_deg_ora_annotation.rds + output/430_annotate_deg_ora.log。
旧 diagnostics 版と同値確認済み)。
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
