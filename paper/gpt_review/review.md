# 論文下書き再レビュー（2026-08-19）

## 再評価後の結論

前回レビューのうち、多重性、PC-OD、主解析コホートの純度差に関する評価を訂正します。投稿前に新しい主検定を追加しなければ成立しない、という状態ではありません。R_Tumor を単一の主対比として明記し、他の対比を別個の副次的問いとして全て報告すること、置換検定の「exact」の射程、Storey q 値、REO、gene-set 較正の射程を正確に記述することで対応可能です。

追加計算は第二の有意差検定を作らない範囲で実行しましたが、後の再評価で、選択後の効果量中央値と通常の `prcomp` PCA は本番原稿に不要と判断して撤回しました。QC、正規化、遺伝子検定、多重性調整は再実行せず、新しい p 値や別の遺伝子リストも作成していません。撤回した出力は解析履歴としてのみ保全します。

## 1. AS と radiogenic/sporadic の表現

表現による対応で十分です。High 群を病因的に radiogenic と断定せず、AS を極端帯の比較に用いるモデル依存の順序指標として規定します。

推奨する中心文は次のとおりです。

> We hypothesized that comparing extreme bands of a model-derived radiation-attributability index might increase contrast relative to an unstratified dose model. AS was used to order cases and define groups; it was not treated as an observed aetiology of an individual tumor.

結論は “radiation trace” ではなく “AS-band-associated expression difference within RET-fusion tumors” とします。5,000語版に反映しました。

## 2. 4対比と多重性

再評価の結果、4対比を一括して多重性調整する必要がある、という前回の指摘は撤回します。Bender and Lange は、複数検定を一つの最終結論・意思決定に結合するときに調整が必要と整理しています。JAMA の解説も、無関係な研究質問には一括調整は不要で、一つの問いを共同で答える family には必要としています（[Bender and Lange 2001](https://pubmed.ncbi.nlm.nih.gov/11297884/); [Cao and Zhang 2014](https://jamanetwork.com/journals/jama/fullarticle/1892228?guestAccessKey=e8e16a4b-47e2-4208-b5d2-1b88816da872)）。

本研究では次の階層が自然です。

- 単一の主対比・主対比水準検定: R_Tumor の Higher Criticism。
- 副次的な別個の問い: R_Normal、B_Tumor、B_Normal、および正常組織の層間一致。
- 遺伝子内の多重性: 各対比内の Storey q 値。
- set 内の多重性: 各 collection 内の BH。

この構造であれば、4つの HC p 値を study-wide family として調整する必要はありません。調整しない条件は、(i) R_Tumor を単一の主対比と明記する、(ii) 4対比を全て報告する、(iii) “one of four was significant” を研究全体の選択的結論に使わない、の3点です。

先行する放射線甲状腺発現研究でも、年齢調整後の遺伝子別解析内で q 値を用いており、論文中の異なる全研究質問を一つの family に束ねてはいません（[Dom et al. 2012, BJC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3464765/)）。したがって、原稿の問題は調整不足より “exactly one place” という対比横断の表現でした。5,000語版では対比ごとの結果を順に示す形へ改めました。

## 3. 追加検討の具体的方法

「感度分析」という包括語は使わず、目的と推定対象を分けます。

### 実行履歴（本番不採用を含む）

既存の R_Tumor 結果を変更せず、次を算出しました。

1. 既存 q<0.10 リストについて、符号付き順位効果 `2p−1` と、表示用の High−Sporadic 平均 log2 CPM 差を要約。
2. 全15,621検定遺伝子を用いた非教師あり PCA。q 値による遺伝子選択や新しい検定は行わず、同じ PC1–PC2 座標へ群、性、相対純度、手術時年齢、RET 融合パートナーを重ねた。
3. 主解析27例について、群・年齢・純度・性・融合パートナーを含むモデル行列の rank だけを参考記録として確認。これは行列が代数的に推定可能であることしか示さず、共変量調整の科学的妥当性や識別可能性を示すものではない。遺伝子モデルは当てはめていない。

1は選択済み遺伝子に条件づけた事後的要約であり、主結果の効果量ではないため、Abstract、本文、補足Methodsから撤去しました。2は `stats::prcomp` による通常のSVD-PCAで、15,621変数・27例のHDLSS設定に適用はできても、HDLSS専用の推定法ではありません。これはGPTが解析法について明示的な承認を得ずに本番案へ追加したものであり、撤回しました。PC-ODとは別手法です。CrossDataMatrixは、別途明確な解析目的が生じた場合にのみ候補とし、不要な置換解析は行いません。3も本番の推論には用いません。

### 共変量調整案についての訂正

前版で提示した `age_at_surgery + relative_purity + sex + RET_fusion_partner` を nuisance 項とする Freedman–Lane 型 global test は撤回します。これは査読で共変量調整を要求された場合の候補として提示したもので、BM 検定の調整版ではありません。現行の遺伝子別 BM 検定は High と Sporadic の周辺的な分布差を順位効果で表し、HC はその p 値を対比水準で集約します。一方、前版の案は log2 CPM の条件付き平均について別の帰無仮説を検定するため、効果尺度、推定対象、仮定がいずれも異なります。前向きに解析を再設計するなら HC に代わる別の global test にはなり得ますが、BM 検定をそのまま置き換えたり、同じ問いを確認したりする解析ではありません。

さらに、本研究では nuisance 項の指定自体に次の問題があります。

- 年齢は AS の構成要素であり、調整すると AS 帯比較が意図して含む変動の一部を除去し、同年齢に条件づけた別の推定対象へ変わる。
- 性は IREP の形式的入力です。実現した AS 値への影響は数値的にほぼありませんが、一律に nuisance 項へ入れるには、性染色体発現を含めて何を推定するかの再定義が必要です。
- 相対純度は同じ発現データから推定され、適格性判定にも使われています。単純な外生的交絡因子として扱うことはできません。
- RET 融合パートナーは発現を規定する一方、放射線と腫瘍形成を結ぶ生物学的経路の一部である可能性があります。調整は意図した生物学的差の一部を除去し得ます。

したがって、モデル行列が27例、7列、rank 7、残差自由度20であったことは計算可能性の確認にとどまり、妥当な調整効果が識別できる根拠にはなりません。媒介変数、結果に関連して生成された変数、または不要な変数の機械的調整は、過剰調整や精度低下を生じ得ます（[Schisterman et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2744485/)）。Freedman–Lane 法にも残差の交換可能性が必要であり、nuisance 項を加えるだけで観察研究の交絡が解消されるわけではありません（[Winkler et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4010955/)）。

本稿では BM＋HC を維持し、推定対象を「事前規定した AS 極端帯間の周辺的な発現差」に限定するのが適切です。年齢、性、純度、融合パートナーは分布を開示し、それらから独立した放射線固有効果を識別していないと明記します。査読者が調整解析を要求した場合にも、まず求める推定対象を特定する必要があります。放射線量と年齢等を分離した効果を求めるなら、AS 帯へ同じ変数を追加投入するのではなく、線量・出生コホート・年齢の重なりを備えた別設計が必要であり、本コホートの単純な追加モデルでは解決できません。確定解析オブジェクトには利用可能な batch 変数もないため、batch 調整を推測で追加することはできません。

## 4. PC-OD と遺伝子検定の関係

ご指摘どおり、前回レビューは QC と遺伝子検定を「ラベル依存前処理」として一括し、論点を混同していました。

PC-OD は遺伝子仮説を検定するものではなく、純度推定より前に解析対象を確定する QC です。そのフラグを固定した後の遺伝子検定は、「確定した解析対象と正規化行列に条件づけた」検定として記述できます。各遺伝子ラベル割当に対して PC-OD を再実行する必要はありません。実現した出力では RET の tumor/normal に PC-OD flag は0件であり、R_Tumor の検体集合は PC-OD によって変化していません。主結果についてこの懸念は実質的にも発生していません。

`filterByExpr` についても訂正が必要です。edgeR の公式ガイドは、group を最小必要検体数の算定に使う一方、どの検体がどの群かとは独立に発現フィルタを適用すると説明しています。固定された群サイズのラベル置換に対して、これを supervised gene selection と扱うべきではありません（[edgeR User’s Guide](https://bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)）。一般論としても、帰無下で独立または置換不変の filter は置換検定と併用できます（[Bourgon et al. 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2906865/)）。

残る統計上の限定は PC-OD ではなく DEGES です。DEGES は観測群を使って scaling factor 推定から除く候補遺伝子を選び、その後の遺伝子置換では正規化行列を固定しています。したがって “exact” は「確定行列上の全割当を列挙したため Monte Carlo 誤差がない」という意味に限定し、raw count から QC・正規化までを含む end-to-end exact randomization と表現しないのが適切です。この射程限定を本文と補足Methodsへ反映しました。

## 5. Storey q 値

追加計算ではなく表現を修正します。固定 λ=0.5 の選択理由を FDR 保証の論証へ拡張せず、次のように水準を分けます。

> Gene-level lists were defined by the prespecified Storey q-value rule with a plug-in π0 estimate at λ=0.5. We refer to genes meeting this rule and do not claim that it proves FDR control under arbitrary gene dependence. The contrast-level existence statement rests on the label-permutation Higher Criticism result.

これにより、π0=0.593 を用いた1,765件は規定手続き下のリスト、HC p=0.0112 は対比水準の証拠、と役割が明確になります。より保守的な第二リストを併載する必要はありません。

## 6. 追加出力の保全と本番不採用

追加監査の主要値は次のとおりです。

- 既存 q<0.10 リスト: 1,765遺伝子（Highで高い971、低い794）。
- `|2p−1|`: 中央値0.600、IQR 0.567–0.644。
- 表示用 `|mean log2FC|`: 中央値0.422、IQR 0.290–0.698。
- Highで高い遺伝子の `|mean log2FC|` 中央値0.384、Highで低い遺伝子は0.573。
- PCA: PC1 19.7%、PC2 16.6%の分散を説明。
- 主解析用相対純度中央値: R_Sporadic 0.783、R_High 0.822。この値は全AS帯を再プールした REO 診断の0.690/0.814とは別の尺度・対象であり、前回レビューは後者を主解析群バランスとして誤用していました。
- PC1と純度/年齢の Spearman rho は −0.396/−0.290、PC2では0.226/−0.042。p値は算出していません。

これらの数値は実行履歴として保全しますが、選択後の効果量要約と `prcomp` PCAは5,000語版のAbstract、Methods、Results、Discussionおよび補足図案から撤去しました。主解析群の年齢、性、相対純度、融合パートナーの単純な群別記述だけを残します。

## 7. REO の規定

表現は次の順に統一します。

> The extreme-band separation is training fit by construction. The only independent comparison was Mid versus Low; its effect was 0.616 and its one-sided p value was 0.1127. The four-band medians are therefore a descriptive display, not evidence of classification performance or a dose response.

純度が31/36例に限られた理由も確定できます。Low 2例、Mid 3例は tumor assay はあるが matched normal assay がなく、ペア型純度推定を計算できませんでした。REO 自体は tumor のみで計算できるため36例全てに実施されています。この説明を補足Methodsへ追加しました。

## 8. gene-set 較正の規定

“calibrated inference” を一般的な保証として使わず、測定した対象をそのまま書きます。

> Complete-null operating characteristics were assessed with held-out label shuffles. The pooled probability of at least one discovery was 0.064, with a B_Normal/Hallmark cell-specific rate of 0.18. These quantities do not establish FDR control under every partial-alternative configuration. The 1.15-fold spike-in was a single-scenario positive-control check, not a general power assessment.

これで complete null、partial alternative、単一 spike-in の射程が分離されます。ORA は引き続き、固定済み遺伝子リストの構成を gene-sampling reference で記述する補助解析であり、被験者ラベル置換による gene-set 検定と同格の第二推論とは表現しません。

## 9. 解釈マップ以前の出力閲覧に関する記録

前回レビューは、スクリプト開発時に通常生じる出力確認を、特定の旧解析結果を体系的に閲覧してから解釈マップを作ったかのように過大評価しました。論文本文をさらに弁護的にするより、内部記録を訂正するのが妥当です。

訂正案では次を行います。

- 「旧パイプラインの暫定結果を閲覧済み」という包括的な性格づけを撤回する。
- 実態を「自作スクリプトの開発・動作確認中に通常の出力表示があった」と記録する。
- 特定の4対比パターンを見た、結果を基にマップを選んだ、という未裏付けの記載を削除する。
- 論文側は “before the finalized analysis produced the results reported here” と、確認可能な範囲だけを述べる。

リポジトリは改変禁止のため、適用対象と置換文案を `proposed_record_corrections.md` に保存しました。

## 残る投稿前作業

方法論上の新規必須解析は提案しません。残る作業は、GDC accession/data release、倫理・同意と二次利用、Data/Code availability、Funding、Competing interests、Author contributions、Acknowledgements、完全な引用文献、ならびに利用可能な technical batch metadata が存在しないことの limitation 記載です。本文から作業用 C/N タグ、対訳、起草メモを除く工程も必要です。
