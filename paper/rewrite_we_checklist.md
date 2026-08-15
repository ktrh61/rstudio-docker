# 書き直し用チェックリスト: We で立てるべき判断文(足場・使い捨て)

位置づけ: 執筆足場。書き直し(受動→能動の使い分け)の道標であり、規則や決定の記録先ではない。
使い終わったら削除してよい(履歴は git)。2026-08-14 作成(Methods 監査直後のコンテキストで採録)。
運用(2026-08-14 確定): 書き直し一通り完了時に本リストと照合(消し込み+★項目の生存確認+
必要ならアンカー張り替え)— **トリガー2(投稿前一括検査)の前段**として実施し、通過後に削除する。

方針(確認済み): 最終形は混合 — **判断・選択・主張スコープの文だけ We の能動**、手続きの連鎖は受動のまま。
npj PO(Nature 系)は能動推奨、BJC/ERC は能動可。

凡例: ★ = 批准済みリント/claim_map が文言まで保護している文(言い回しを変える場合は検査を通すこと)。
アンカーは現行文の冒頭数語(行番号は書き直しで動くため使わない)。

## Methods

### Exposure metric

- [ ] "Radiation attributability per exposed case was quantified..." — 曝露指標として IREP AS を採用した選択(線量そのものでなく)。例: We quantified radiation attributability with the NIH IREP thyroid model as the AS...
- [ ] "Cases were banded on AS with a single pre-fixed rule..." — 帯規則の事前固定。例: We banded cases with a single rule fixed in advance...

### Quality control and analysis cohorts

- [ ] "Sample outliers were screened before purity estimation, so that..." — スクリーニングを純度推定より前に置く順序判断。例: We screened outliers before estimating purity so that...
- [ ] "...in one run per driver cohort with both groups (Sporadic and High) together:" — contamDE を群プールで1回実行する判断+その理由(共通尺度)。例: We estimated purity in one pooled run per driver cohort because...
- [ ] "The main analysis cohort comprises... (n = 63)" — コホート定義(採用条件の束)。例: We defined the main analysis cohort as...
- [ ] "...deliberately left unfiltered by the outlier and purity conditions" — 評価セットを意図的に未濾過とする判断。例: We deliberately left the evaluation set unfiltered...

### Analysis contrasts

- [ ] ★ "All inference is organized in four analysis contrasts... no cross-contrast family-wise inference is performed, and no study-wide FDR is claimed" — 推論スコープの宣言(リント(1)が「no study-wide FDR is claimed」の明文化を要求 — We 化しても否定宣言そのものは残すこと)。例: We organize all inference in four contrasts... and make no study-wide FDR claim.
- [ ] "Throughout, two-sample comparison uses the Brunner–Munzel statistic..." — プロトコル全体の原理コミットメント(最重要の判断文; 2026-08-14 に Gene-level から Analysis contrasts へ格上げ移動済み)。例: We take the order relation as the protocol's primitive...

### Normalization

- [ ] "Each contrast's count matrix... was normalized with a DEGES scheme... reimplemented in-house" — 自作再実装+BM スクリーン置換という手法選択。例: We normalized each matrix with a DEGES scheme we reimplemented, substituting...

### Gene-level differential expression

- [ ] "The estimator was fixed a priori from the working hypothesis and the design rather than tuned" — λ=0.5 の事前固定(選定事由の核)。例: We fixed the estimator a priori...

### Contrast-level omnibus inference(2026-08-15 に Gene-level から分割)

- [ ] "Higher Criticism (α0 = 0.1) is the pre-specified primary omnibus statistic, chosen..." — 主オムニバスの事前指定。例: We pre-specified Higher Criticism as...
- [ ] ※Intro 起草後の再訪: 主張担持文("The contrast-level inferential claim is carried by...")の文言を Intro 予告列とハーモナイズ(2026-08-15 に "label-aligned signal" 形へ精密化済み — 曝露帰属はオムニバスでなく解釈マップの領分、の線を保つこと)
- [ ] "The full rejection curve R(α) is retained so that..." — 単一閾値表示に依存させない判断。例: We retain the full rejection curve so that...
- [ ] ★ "No contrast-level binary significance label is assigned... pre-assigned before the results were seen" — 非2値化+解釈の事前割当(C-16)。例: We assign no binary label... We pre-assigned the interpretation of each pattern...

### Gene-set level inference

- [ ] "Set-level enrichment consumes the whole per-gene ranking (threshold-free...)" — 閾値フリー設計の選択。例: We test set-level enrichment on the whole ranking...
- [ ] "...evaluated at tie-block boundaries only, so that no arbitrary tie-break injects an order..." — tie-block 統計量の選定事由。例: We evaluate the running sum at tie-block boundaries only, so that...
- [ ] "Unlike the gene level, π0 is held at its conservative bound of 1 here:" — 遺伝子/セット水準で推定量を非対称にする判断。例: Unlike the gene level, we hold π0 at 1 here because...
- [ ] "The complete gene-set inference was calibrated on held-out null replicates... fixed by this calibration prior to any real-data run" — 較正先行+FDR 手続きの較正による固定(D6)。例: We calibrated the complete inference before touching the real contrasts, and fixed the FDR procedure by that calibration.
- [ ] "...whose curation rule was fixed before this inference touched real data" — 放射線ファミリーのキュレーション規則の事前固定。例: We fixed the curation rule before...
- [ ] ★ "The two procedures' p-values refer to different randomness: ... experiment-level claims rest on it" — 実験レベル主張の置き場所の宣言(Q-15 決着の一文開示)。例: We rest experiment-level claims only on the label-permutation null...

### Between-stratum concordance

- [ ] "The normal-tissue pair is the pre-specified hypothesis-bearing comparison; the tumor pair is a descriptive completion." — 仮説担持/補完の役割割当。例: We pre-specified the normal-tissue pair as the hypothesis-bearing comparison...
- [ ] "An observed rho outside the interval indicates... it does not identify that structure as..." — 読みの上限の事前固定。例: We do not read a rho outside the interval as identifying...

### REO panel

- [ ] "...the q < 0.10 DEG set is deliberately not used" — DEG 集合を使わない判断。例: We deliberately did not use the q < 0.10 set...
- [ ] "The single pre-specified test is therefore the one-sided BM comparison of Mid over Low" — 唯一の事前指定検定の宣言。例: We pre-specified a single test: ...
- [ ] ★ "This evaluation does not alter the panel or its boundary, and the band-wise score profile is reported descriptively without assuming any dose–response form." — 無改変適用+線量反応非仮定(リント(3))。例: We applied the panel untouched and report the profile descriptively, assuming no dose–response form.

### REO evaluation diagnostics

- [ ] "...its QC is reported as diagnostics with no exclusion authority" — 診断に除外権限を与えない判断。例: We report the evaluation QC as diagnostics with no exclusion authority...

### Candidate confounders(旧 Between-group age difference — 2026-08-15 に3候補構成へ再構成し Analysis contrasts 直後へ移設)

- [ ] "Three candidate confounders of the exposure contrasts follow from the design and are disclosed rather than adjusted for..." — 交絡フェーズの宣言(設計からの導出+調整でなく開示。「pre-specified」は 2026-08-15 撤回 — 暦の先行を主張しない)。例: We name three candidate confounders that follow from the design, and disclose rather than adjust for them...
- [ ] ★ "this is a disclosure of magnitude and uncertainty, not a confounding test, and no p-value is computed" — 開示であって検定でないという宣言(リント(7))。例: We disclose magnitude and uncertainty and compute no p-value...
- [ ] ★ "Age is not entered as a covariate, because..." — 共変量非投入の判断(Q-03 の受け; QC 節から本節へ移設済み)。例: We did not enter age as a covariate because...

### External anchor cross-reference

- [ ] ★ "...the reading was fixed symmetrically in advance — any count is reported as description and no claim moves with the outcome" — 対称読みの事前固定(C-13)。例: We fixed the reading symmetrically in advance...

### Software, seeds and reproducibility

- [ ] "...four workers and 9,999 label shuffles as the fixed reproduction contract" — 再現契約の宣言(任意 — 手続き文のままでも可)。例: We fix four workers and 9,999 shuffles as the reproduction contract.

## Discussion(書き直し時に We を検討する宣言文)

- [ ] ★ "Throughout this section, readings drawn from the interpretation map (Methods) were assigned to outcome patterns before the results were seen; interpretations that go beyond the map are labelled as hypothesis-generating." (Disc §1) — 事前固定の一括宣言(2026-08-15 改稿で各節の反復をこの1文に集約。スコープは「マップ由来の読み」に限定 — C-17 の仮説生成枠は対象外という線を保つこと)。例: We assigned the readings ... before seeing the results...
- [ ] ★ "no claim ... moves with these counts / rests on it / rests on this annotation" 系(C-13・C-17・C-14)— 主張非依存の宣言。言い回しの散らしは 2026-08-15 改稿で実施済み(moves with / rests on it / rests on this annotation)。全て必要な文 — 削らない。
- [ ] ★ "instead the pattern receives the reading the map assigned to it in advance — confounding first" (B_Normal, C-16) — 事前固定マップの条件適用(2026-08-15 改稿で短縮形へ — 実体は不変)。
- [ ] ★ "Whether the two strata share a normal-tissue exposure trace is not identifiable in this cohort" (C-07) — 識別不能の宣言(2026-08-15 改稿で前置き句を削除 — 実体は不変)。

## 受動のままで良い代表例(変換不要の確認用)

- データ取得・照合の連鎖("Counts ... were downloaded and md5-verified"、"Expression files were mapped...")
- 帯割付・フロー・正規化・検定の実行文("Cases were banded..."、"Each gene was tested...")
- 実測開示("all 906 libraries were stranded (reverse)"、"no case sits on a boundary")
