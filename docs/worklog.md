# worklog

常時更新される唯一の生きた文書。確定した設計選択は [decisions](decisions.md)、詳細論拠・実施記録は [records](records.md)、恒常前提は [profile](profile.md) に置く。ここには「まだ動いている事実」だけを残す。

---

## A. 未確定課題（オープンな問い）

- `status()$data_release` の返り値構造の確認（01 改修時、実装時に `status()` を1回実行して確認）。
- リリース番号取得失敗時の挙動（fail-fast 停止か、番号のみ取得不可で続行か。後者は広義フォールバックで方針と衝突）。
- 環境情報記録スクリプトの番号・命名・位置（[D-009](decisions.md#d-009-環境情報記録は専用スクリプトで自動出力git-追跡) の昇格に関与）。
- MUREN 参照選択（`include_self=FALSE` / `refs_cap`）の理論的再検討と是正方針（設計未決。詳細分析は records 5-1）。
- contamDE の MUREN 置換の妥当性（5-1 に従属）。
- gls.cpp 採用の結果影響の最終評価。
- RcppArmadillo/OpenBLAS 廃止の可否（[D-016](decisions.md#d-016-blas-スレッド1固定) に関与、廃止しない可能性も残す）。
- 開発リポジトリの品質是正の範囲（[D-008](decisions.md#d-008-開発公開分業) に従属）。
- 実在スクリプト集合と連番の実態確認（00→13 は暫定。欠番・重複・archive 済みの切り分け、utils・setup の現役集合との対応を 4-2 冒頭で棚卸しする）。
- **本番（ThinkPad）イメージ構築時に再判断する点**（開発の便法を本番へ黙って漏らさない。[D-020](decisions.md#d-020-開発環境の-blas-構成openblas-pthreadupdate-alternatives-方式) の趣旨。4-1b は OMEN 開発環境の成果であり、以下は本番で明示的に決め直す）:
    - BLAS: 開発の OpenBLAS pthread（alternatives 方式）→ 本番で参照 BLAS に揃えるか（[D-020](decisions.md#d-020-開発環境の-blas-構成openblas-pthreadupdate-alternatives-方式)）。
    - ライブラリ配置: `/opt/r-extra-lib` のホストマウント → イメージ内に焼き込むか（records「ubuntu ベース転換と 4-1b 成果」5）。
    - 作業ユーザー: noble 標準 `ubuntu`（UID 1000）流用 → 本番名に調整するか（同 7）。
    - ベース固定: 日付タグ `noble-20260410` → digest 明示（`sha256:c4a8…`）に上げるか（[D-019](decisions.md#d-019-ベースイメージubuntunoble-20260410-r-453-自前ソースビルドrstudio-廃止)）。
- **決定間の相互参照リンクの整合性検査（今後の課題・現時点は運用注意）**: 採番/supersede 時に手書きの `decisions.md#d-0XX-…` リンクが静かにドリフトする。自動検査（チェッカ/フック/CI）を検討したが、参照数が二桁前半で費用対効果が薄く、チェッカ自身がレンダラのアンカー生成規則を写し取る必要があるため見送り。当面は supersede/採番時に grep 目視で突き合わせる。

---

## B. Pending-verification の昇格条件（decisions と連動）

- [D-017](decisions.md#d-017-contamde-の-p-値調整を-qvaluebh-是正)：4-2 で該当スクリプト是正・結果影響評価。
- [D-016](decisions.md#d-016-blas-スレッド1固定)：4-2 で結果影響の最終評価（実装済み・評価待ち）。
- [D-003](decisions.md#d-003-パッケージ版固定p3m-日付スナップショット-2026-04-09)：install 実行時に `cran/2026-04-09` の全依存有効性を最終検証。4-1b で明示49（依存含め総216）の導入が P3M 4/09 で完走したことを確認（残: 4-2 実行時に判明し得る追加依存）。

---

## C. 手順の現在地

- **[完了]** 作業同期プロトコル確定。`rebuild/r453` 作成（ThinkPad 2026-07-05／OMEN 2026-07-09、いずれも 87a75a0 土台）。4-1 環境の器の初期検証（`rocker/rstudio:4.5.3` で P3M 2026-04-09 動作・主要パッケージ版確認。器は [D-019](decisions.md#d-019-ベースイメージubuntunoble-20260410-r-453-自前ソースビルドrstudio-廃止) で ubuntu ベース＋ R 4.5.3 自前ソースビルドへ転換〔RStudio 廃止〕）。
- **[運用]** 開発環境の二拠点化：4-2 以降の開発検証は OMEN（速度優先・必要パッケージ数が未確定なため）で行い、最終出力段階で ThinkPad を使う。ThinkPad が絶対的 production である位置づけは不変（開発の場を OMEN に置くのみ。本番イメージ構築・最終出力・正典は ThinkPad）。OMEN は外部アクセス可能な状態。
- **[進行]** docs/ 文書運用移行。
- **[進行／一部完了]** 4-1b 環境実体化（OMEN、2026-07-09）。`ubuntu:noble-20260410` ＋ R 4.5.3 自前ソースビルドの開発用コンテナ `rebc-r453-dev` を構築・検証（明示49導入完走・主要版一致。材料は records「ubuntu ベース転換と 4-1b 成果」、決定は [D-019](decisions.md#d-019-ベースイメージubuntunoble-20260410-r-453-自前ソースビルドrstudio-廃止)／[D-018](decisions.md#d-018-apt-層の日付固定snapshotubuntucom-2026-04-10level-2)／[D-020](decisions.md#d-020-開発環境の-blas-構成openblas-pthreadupdate-alternatives-方式)）。残: 暫定 Dockerfile 化（この材料を元に別途）、`contamde_gls.cpp` のコンパイル・実行確認（次段）、4-2 単独実行で判明し得る追加 system-lib。確定後に ThinkPad で本番イメージ構築（三層構成 (B) 発動）。
- **[残点]** `snapshot.ubuntu.com` の取得可否。4/09 で 503 を観測したが、apt 固定日を **4/10 へ移行**（[D-018](decisions.md#d-018-apt-層の日付固定snapshotubuntucom-2026-04-10level-2) が [D-004](decisions.md#d-004-apt-層の日付固定snapshotubuntucom-2026-04-09level-2) を supersede）、4-1b で 4/10 HTTPS 到達を実測（ca-certificates は Verify-Peer 一時無効化で初回導入）。503 の恒常性は依然未判定。システムライブラリは因果確立分＋`libssl-dev`（今回初判明）まで確定、最終確定は 4-2 依存分析。
- **[未着手]**
    - 4-2 スクリプト依存順の作り直し（実在集合と連番を冒頭で棚卸し。00→13 は暫定）。
    - 4-3 環境記録スクリプト新設。
    - 4-4 通し検証（整えたスクリプト群が依存順に完走し、妥当な中間・最終出力を産むことの確認。旧結果値の再現を目標としない＝現結果保持の圧力は外している）。
    - 4-5 凍結・記録。

前提：4-2 は 4-1b の OMEN 開発用コンテナ上で行う。スクリプト単独実行で判明した依存（パッケージ・システムライブラリ・構築時インストール）は暫定 Dockerfile に反映して育てる。

**4-2 各スクリプトの手順テンプレート**（1本ずつ順守）:

1. パス参照を新構造（フラット）へ書き直す（setup も utils も統一規則で）。
2. 遺物との切り分け（現役依存のみ残す）。
3. 版表記・手書き日付・スクリプト間相互参照を削除（Git に委譲）。
4. fail-fast 確認（暗黙フォールバック・握りつぶすエラー・存在を仮定した外部依存がないか）。
5. R 4.5.3 で単独実行し、期待される中間出力を確認。
6. 追跡可能な単位（原則 1 スクリプト）でコミット。構造・パス変更とロジック修正はコミットを分ける。

---

## D. 現役の棚卸し事実（方向性のみ）

- データ入出力パスは `paths$*` 経由で規則的。setup.R の paths 定義変更で一括追随し、フラット化は安全側。
- utils は一本化し source を統一する。現役 source は4ファイル（`utils_improved.R`, `norm_improved.R`, `with_openblas_threads.R`, `contamde_purity_functions.R`）。配置は新構造決定に従属（utils だけ先行して寄せない）。
- 依存チェーン（00→13、07=deges 正規化／08=pca 可視化、11=reo_pair_selection、14 なし）。各実入出力は 4-2 の単独実行で確定する。

---

## E. 未処理事項

- `v7-data-expansion` のローカル残存削除（OMEN、reconciliation 残処理）。
- main の整理（v6 の `461429c` で止まった幹を正系列へ合流させるか）。
- macOS の 00 以外の untracked 掃除（macOS 再開時）。
- 9月WIP（`archive/wip-v7-data-expansion-2025-09-17` タグ）の中身の最終判定。
- `data/raw/gdc/DR40/downloads` と `manifest_rebc_thyr_20250902_031018.tsv` の正体確認。
