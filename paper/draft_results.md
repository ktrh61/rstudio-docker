# Results ドラフト(全節、英日対訳)

状態: draft(Claude Code 起草 2026-08-14。研究者が両版を確認し、最終的に自分の文章として
書き直す前提の下書き)。
写像元: claim_map C-01〜C-10, C-12, C-14〜C-17(全て批准済み決着に基づく draft 行)。
数値は全て numbers_ledger の N-ID 経由(タグ併記)。英語版が正本、日本語訳は確認用の対訳。
D6 較正文(C-06)と ORA 注釈文(C-14)はそれぞれ draft_methods_results_d6.md /
draft_ora_annotation.md が正本 — 本稿には接続位置を示すため同一文を再掲(その旨注記)。
図表番号は全て(仮)。引用文献は著者年表記(書式整合は後工程)。

---

## 1. Cohort(C-10, C-15)

<!-- C-10 -->
Of the 440 cases of the REBC-THYR cohort (Morton et al. 2021) <!-- N-11 -->,
the main analysis cohort comprised 63 paired, driver-stratified cases —
9 B_High, 27 B_Sporadic, 15 R_High and 12 R_Sporadic <!-- N-09 --> — reached
through the pre-specified flow of driver classification, band eligibility,
pairing, outlier screening and purity thresholding (440 → 248 → 77 → 70 →
69 → 63 <!-- N-08 -->; Tables 1–2(仮)). The REO evaluation set added 36
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
(この段落 = draft_methods_results_d6.md 本文 Results と同一文。正本はあちら)

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
(この段落 = draft_ora_annotation.md 本文 Results と同一文(数値修正 2026-08-14 反映)。
正本はあちら)

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
and R_Mid 8 / 11 <!-- N-42 -->. The single pre-specified out-of-sample test
gave one-sided Brunner–Munzel p = 0.1127 for Mid > Low, with effect
P(Low<Mid) = 0.616 <!-- N-40 -->. The graded profile is reported as a
descriptive observation; no dose–response form is assumed or claimed.
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
<!-- N-42 -->。唯一の事前指定 out-of-sample 検定である Mid > Low の片側
Brunner–Munzel は p = 0.1127、効果 P(Low<Mid) = 0.616 だった <!-- N-40 -->。
この段階的プロファイルは記述的観察として報告し、線量反応の形は仮定も主張もしない。
訓練側から鏡映した外れ値スクリーンはどちらの評価帯でも該当例を出さなかった
(17 例中 0・19 例中 0 <!-- N-43 -->)— スクリーンは除外権限を持たず、評価は適格
全例で一度だけ実施した。評価帯内では逆転スコアは腫瘍純度と相関し(pooled Spearman
+0.538 <!-- N-45 -->)、帯–スコア相関は +0.142、純度を条件付けると +0.146(偏 Spearman;
片側置換 p = 0.2162)<!-- N-48, N-46 -->、帯と純度自体はほぼ無相関だった(+0.036
<!-- N-48 -->)。純度層別の比較は Supplementary に示す <!-- N-47 -->。

---

## 引用文献(Results 本文で必要になるもの — 書式は仮)

- Morton LM, et al. Radiation-related genomic profile of papillary thyroid
  carcinoma after the Chernobyl accident. Science 2021;372:eabg2538.
  doi:10.1126/science.abg2538(コホート出典 — §1 で引用)

(Results 本文で必要な文献は現状この1本。Storey・Brunner–Munzel・GSEA・HC・ASA 系・
G&B・Donoho-Jin は Methods 側、Vriens・Abend・Ory・Coclet は Discussion/Intro 側で引用)

---

## 起草メモ(本文に載せない)

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
- 図表番号は全て(仮)。判断点5(図表構成)確定後に振り直し+図表台帳と照合。
- リント自己検査(10項目): (1) unit 横断 FDR 主張なし / (2) A/B 中立(機構文言なし)/
  (3) 線量依存の仮定なし(C-08 に明示の否定文)/ (4) q<0.25 不使用 / (5) B_Tumor の
  「特異性確認」文言なし(明示の否定文)/ (6) 探索セル産なし(ORA は hypothesis-
  generating 明示、C-17 は Disc 送り)/ (7) 年齢は推定+CI のみ(検定なしを明示)/
  (8) 全数値 N タグ済み / (9) WY 不使用 / (10) 主張水準は C 行の水準列どおり
  (主張 = C-01 のみ、他は記述的観察・仮説生成)。
