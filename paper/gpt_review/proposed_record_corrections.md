# 内部記録の訂正案（リポジトリ未適用）

作成日: 2026-08-19

## 訂正の趣旨

現行の一部台帳は、自作スクリプトの開発・動作確認中に通常の出力表示があったことを、「旧パイプラインの暫定結果を体系的に閲覧してから解釈マップを固定した」と過大に性格づけています。特定の4対比パターンを結果として取得・検討したことや、それを根拠に解釈マップを選んだことは、今回確認された事実ではありません。

適切な記録は次の範囲です。

- 自作スクリプトの開発中、実データによる通常の動作確認と出力表示はあった。
- これは正式解析または体系的な暫定解析結果のレビューとは位置づけない。
- 解釈フレームワークは、現在報告する確定解析を実行してその結果を得る前に定めた。
- “prospectively registered”, “before any data analysis”, “blinded” など、事実以上に強い主張はしない。

以下は適用候補であり、リポジトリファイルには反映していません。

## 1. `paper/objections_ledger.md` の Q-18

現行の「前世代実装の出力（R_Tumor 有意・B_Tumor 静か・B_Normal 異常高）を閲覧」という列挙、およびそれを弁護する4段構成を削除し、次の趣旨へ置換します。

> **回答方針** — 本研究のスクリプト開発中には、実データを用いた通常の動作確認と出力表示があった。ただし、特定の4対比パターンを正式または体系的な暫定結果として取得・検討したうえで解釈マップを選定した、という記録は過大であり撤回する。解釈フレームワークは、現在の確定解析が本稿で報告する結果を生成する前に定めた。したがって紙面では “before the finalized analysis produced the results reported here” と記述する。これは prospective registration、全データ解析前、blind review を意味しない。問われた場合は、開発中に実データによる動作確認があったことと、確定解析結果を基にフレームワークを選んだものではないことを、この範囲で回答する。

根拠欄から、未確認の旧4対比パターンと N-16/N-22/N-23 の接続を外します。状態は研究者確認後に ratified とするのが適切です。

## 2. `paper/numbers_ledger.md` の 2026-08-18 閲覧歴エントリ

373–390行相当の長い弁護記録を、次の訂正記録へ置換します。

> 2026-08-19 訂正: 2026-08-18 の「暫定結果閲覧歴」記録は、スクリプト開発時の通常の出力確認を、旧解析結果の体系的閲覧として過大に性格づけていたため撤回する。確認できる事実は、実データを用いた開発・動作確認中に出力表示があったこと、および現在報告する確定解析の結果生成前に解釈フレームワークを定めたことである。特定の旧4対比パターンを取得・検討したとの列挙は削除する。紙面の時制は “before the finalized analysis produced the results reported here” とし、prospectively specified / before any analysis / blinded へ強化しない。

この訂正は数値の変更ではなく、開発来歴の性格づけの修正です。

## 3. `paper/claim_map.md` の改訂メモ

260–262行相当の「旧パイプライン暫定結果の閲覧歴」という記述を、次へ変更します。

> 時制は、確定解析が本稿で報告する結果を生成する前にフレームワークを定めたことを表す形へ統一する。開発中の通常の出力表示を、体系的な暫定結果閲覧としては記録しない。

268行相当の「閲覧歴に関する保留事項の決着」は、「開発来歴の過大な性格づけを numbers_ledger 2026-08-19 訂正および Q-18 で修正」に変更します。

## 4. `paper/rewrite_we_checklist.md`

4箇所の凍結文言を次へ変更します。

> before the finalized analysis produced the results reported here

チェック項目は次のとおりです。

- “prospectively registered”, “before any data analysis”, “blinded”, “before any results were seen” へ強化しない。
- “before the reported results existed” という存在論的な表現より、どの解析・どの結果を指すかが明確な上記文を使う。
- 開発中の通常出力確認を追加開示する必要はないが、尋ねられた場合に否定もしない。

## 5. `paper/draft_manuscript.md` の4箇所に対する表現案

内部記録訂正と整合させる場合の置換案です。

1. Introduction:

   > We evaluated four contrasts under an interpretation framework set before the finalized analysis produced the results reported here.

2. Methods:

   > The interpretation framework for the four contrasts is given in Table 2; readings for the possible outcome patterns were assigned before the finalized analysis produced the results reported here.

3. Discussion:

   > Readings drawn from the interpretation framework were assigned before the finalized analysis produced the results reported here; interpretations beyond that framework are labelled as hypothesis-generating.

4. Table 2 caption:

   > Interpretation framework for the four exposure contrasts, set before the finalized analysis produced the results reported here.

“pre-specified” を残す場合は、確定解析に対する事前設定という限定で用い、臨床試験登録またはデータ非閲覧を含意させないようにします。
