# 論文下書き再レビュー（現状反映版、2026-08-19）

## この文書の位置づけ

本レビューは `paper/gpt_review/draft_manuscript_5000_words.md` の編集判断を記録するものであり、`paper/` 直下の正本や解析計画を置き換えません。以下では、解析出力から確認できる事実、原稿上の表現判断、将来の査読対応案を区別します。過去版に含まれていた未承認解析の結果、過度に強い事前指定の表現、REOを単一の「独立比較」へ還元する説明は撤回します。

## 現在の総合評価

本稿の中心結果は、RET fusion-positive PTC層内でHigh-AS腫瘍とdose-zero腫瘍の発現分布が異なったことです。この結果はAS帯との関連を示しますが、放射線固有の発現痕跡を識別しません。ASが被ばく時年齢をリスク修飾因子として取り込むこと自体は欠点ではありませんが、その役割と年齢・出生コホートが発現に直接関連する経路を本デザイン内では分離できません。相対純度、融合パートナー、性、利用できないtechnical batch情報も、放射線固有の解釈を限定する未解決の構造として残ります。

投稿前に別の主検定や共変量モデルを追加しなければ論文が成立しない、とは判断しません。必要なのは、各推論の対象と限界を揃えること、事務情報と完全な参考文献を補うこと、図表・補足資料を本文と一致させることです。

## 1. ASと病因の表現

ASは症例を順序づけて帯に分けるモデル由来指標として扱い、個々の腫瘍がradiogenicかsporadicかを観測した変数として扱いません。したがって、中心結論は “AS-band-associated expression difference within RET fusion-positive PTCs” の範囲に置きます。radiation-specific trace、radiation marker、probability of radiogenic originといった病因的表現には進めません。

acuteのexposure-rate設定は主解析で維持しますが、post-Chernobylの内部ヨウ素131被曝の時間経過を物理的に表すとは主張しません。本文では事前に固定した計算規約として定義し、33.3と66.6をその数値尺度上の操作的境界として扱います。Zurnadzhy et al. 2022とBogdanova et al. 2022はacuteの妥当性を独立に検証する根拠ではなく、同じ研究プログラム内での計算規約の来歴としてSupplementary Methodsでのみ位置づけます。

既存の簡易的な入力規約監査では、chronic算定は被曝症例の順序を保存し、ペアごとの順位逆転はありませんでした。これはASの順序構造が安定していることを示しますが、acuteの物理的妥当性や、同一の数値境界で帯所属が不変であることまでは示しません。したがって、本稿はacute尺度で固定した帯を報告し、chronicで同じ発現結果が得られるという未確認の主張は行いません。

NCOA4-RETとCCDC6-RETについては、過去の名称との対応をIntroductionで一度示す現在の方針が適切です。RET/PTC3はNCOA4-RET、RET/PTC1はCCDC6-RETに対応し、その後はgene-partner表記に統一します。

## 2. 4対比と多重性

4つの対比を必ず一つのstudy-wide familyとして調整すべき、という以前の指摘は撤回します。ただし「調整が不要であることが先行研究によって証明された」とも書けません。familyの範囲は、本稿がどの問いを一つの結論にまとめるかによって決まります（[Bender and Lange 2001](https://pubmed.ncbi.nlm.nih.gov/11297884/); [Cao and Zhang 2014](https://jamanetwork.com/journals/jama/fullarticle/1892228)）。

現行稿では、RET-tumorを単一のprimary comparisonとし、RET-normal、BRAF-tumor、BRAF-normalを別のsecondary questionsとして全て報告します。この階層を維持し、“one of four was significant” のような選択的な研究全体の結論を作らない限り、4つのHigher Criticism p値を一括調整しない説明は可能です。遺伝子のStorey q値は各対比内、gene-setのBH調整は各collection内に限定され、study-wide FDRは主張しません。これは編集上採用したfamily定義であり、査読者が別のfamilyを要求する可能性までは否定しません。

本文では解析コード名 `R_Tumor` などをIntroductionへ出さず、組織とdriverを通常の語で記述します。コード名はSupplementary Methodsで実装との対応を示す場合に限ります。

## 3. 共変量調整と追加解析

年齢、相対純度、性、RET融合パートナーを一括してnuisance項に入れる追加モデルは採用しません。被ばく時年齢は放射線リスクの修飾因子としてASに意図的に取り込まれており、これを遺漏したnuisance変数とは扱いません。AS構成に使った年齢入力に条件づければ、現在の周辺的なAS帯比較から年齢条件付きの異なる推定対象へ移ります。さらに、純度は発現から推定され適格性判定にも使われ、融合パートナーは生物学的経路の一部である可能性があります。これらは共通のnuisance解釈を持たず、一括調整しても「放射線固有効果」が識別されるわけではありません。そこで本稿は、AS構成の妥当性と、AS帯に整列する年齢・出生コホートの発現への直接関連を分離できないという識別限界を区別します。

過去にGPTが追加した次の計算は、本稿の根拠から除外します。

- 選択済みq<0.10遺伝子における効果量中央値
- `stats::prcomp`による通常のPCAとPC相関
- 候補共変量モデル行列のrank確認

これらは主解析を検証せず、選択後要約、HDLSS専用でない未承認PCA、代数的rankの確認にすぎません。Abstract、本文、Supplementary Methodsには使用しません。CrossDataMatrixも実行していません。PCAに答えるべき明確な科学的問いが現在ないため、単なる手法の置換としてCrossDataMatrixを追加する必要もありません。将来必要になった場合は、目的、対象行列、評価基準を先に定めた別解析として扱うべきです。

## 4. PC-OD、正規化、置換検定

PC-ODは遺伝子仮説検定ではなく、純度推定前に行った症例レベルQCです。確定したQCフラグに条件づけて遺伝子検定を行う現在の解析では、ラベル割当ごとにPC-ODを再実行するとは規定しません。実現した主RET解析ではtumor・normalともPC-OD flagは0件で、QCによるRET症例の除外はありませんでした。

`filterByExpr`は固定された群サイズから必要検体数を定めますが、どの症例がどちらの群かを用いた効果方向の遺伝子選択ではありません。一方、DEGESは観測群を用いて正規化因子を推定し、その後の遺伝子別置換では正規化行列を固定しています。したがって “exact” は、確定済み正規化行列上で全ラベル割当を列挙しMonte Carlo誤差がない、という範囲に限定します。raw countからQC・正規化を含めて毎回再計算するend-to-end exact randomizationとは表現しません。

## 5. Storey q値とHigher Criticism

1,765遺伝子は、固定λ=0.5のplug-in π0推定を用いる所定のStorey q<0.10規則で得たリストです。この事実と、任意の遺伝子依存構造の下でFDR制御が証明されたという主張は同じではありません。本文では前者だけを述べます。

RET-tumorのHigher Criticism p=0.0112は、対比全体にラベルと整列した発現構造があるかを問う結果です。遺伝子リストの大きさとomnibus evidenceは役割が異なり、相互の代用品として扱いません。

## 6. Gene-set解析とORA

被験者ラベル置換に基づくgene-set解析では、q<0.10を満たすsetはありませんでした。一方、固定済みRET-tumor遺伝子リストに対するORAでは、High-ASで低発現側にcell-cycle、proliferation、DNA-repair関連の注釈が集まりました。これらを一つの相関した生物学的テーマとして報告することには情報価値があります。

両者は同じ解析の陽性・陰性判定ではありません。前者は症例ラベルの交換を参照し、後者は検定対象遺伝子を母集団とするgene-sampling referenceを使います。ORAを独立した確認または第二の主検定とは呼ばず、gene-set側のnull resultだけでORAの記述的知見を消去もしません。

complete-null評価のpooled rate 0.064とcell-specific maximum 0.18は、実行したcomplete-null入力における挙動です。あらゆるpartial-alternative設定でのFDR保証または一般的なcalibrationとは表現しません。1.15-fold spike-inも単一条件のpositive-control checkに限ります。

## 7. REOの正しい位置づけ

REOにはbinary readoutとgraded readoutがあります。dose-zeroとHigh-ASはパネルの構築帯であり、境界はdose-zero例の最大スコアから設定されました。この2帯で全12 dose-zero例がnegative、15 High-AS例中13例がpositiveとなったことは、固定パネルを他帯へ適用する前提となるconstruction fitを示します。構築例がpair selectionに使われ、dose-zero例が境界設定にも使われているため、不偏な分類性能ではありませんが、結果として不要という意味でもありません。

固定した境界をLow-ASとMid-ASへ変更なく適用すると、Lowはnegative 9／positive 8、Midはnegative 8／positive 11でした。graded scoreの帯別中央値はdose-zero／Low／Mid／Highで0／1／4／6でした。これは固定スコアで観察された4帯の記述的プロファイルです。

事前に定めた片側Brunner–Munzel比較はMidとLowだけを用い、P(Low<Mid)=0.616、p=0.1127でした。この検定は2帯間の方向的な確率順序を評価するもので、4帯の線形性、単調性、dose-response formを検定しません。したがって、p=0.1127だけでconstruction fitや4帯表示を無効とすることも、0／1／4／6の中央値だけで線形性を主張することも不適切です。また、パネルとスコア方向が構築帯から定義されているため、Mid–Low比較や4帯表示を“independent validation”とは呼びません。

純度推定が得られたのは中間帯36例中31例です。残るLow 2例とMid 3例にはtumor assayがある一方、ペア型純度推定に必要なmatched normal assayがありません。REOはtumorのみで計算できるため36例全てを含みます。純度との関連解析は除外権限を持たない診断であり、純度から独立したband effectを証明しません。

## 8. 解釈フレームワークの時制

確認できるのは、現在報告する確定解析が結果を生成する前に解釈フレームワークを定めたことです。開発中には実データを用いた通常の出力確認がありました。したがって、“prospectively registered”、“before any analysis”、“before any results were seen”、“blinded”とは書きません。

本文の時制は次の範囲に限定します。

> before the finalized analysis produced the results reported here

開発時の通常の出力閲覧を「旧4対比パターンを体系的に確認してからマップを選んだ」と記録した過去の説明は過大でした。詳細な訂正案は `proposed_record_corrections.md` に分離されています。

## 9. 残る投稿前作業

新規の必須統計解析は提案しません。残る主要作業は、GDC accession/data release、倫理・同意と二次利用、Data/Code availability、Funding、Competing interests、Author contributions、Acknowledgements、完全な引用文献、本文と図表・Supplementary Dataの数値および名称の照合です。利用可能なtechnical batch metadataがないことはlimitationとして明記します。
