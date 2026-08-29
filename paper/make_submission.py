#!/usr/bin/env python3
# 投稿版生成(派生物 — 正本は不変)
#
# 正本 = paper/draft_manuscript.md / paper/supplementary_material.md
#   (author-year 引用+N/C タグ = 照合可能状態を維持する。本スクリプトは読むだけ)
# 出力 = output/submission/
#   manuscript_submission.md : タグ・メタ記述を除去し、本文引用を初出順の
#     Vancouver 数値 [n] へ変換、References を引用台帳の書誌で番号順に生成
#   supplementary_submission.md : タグ・メタ記述の除去のみ(SI は author-year のまま —
#     BJC は SI を submitted のまま掲載するため番号体系を本文と分離しない)
#   warnings.txt : 未変換候補・台帳不一致・除去タグ数の全数レポート
#
# claim_map「紐付け機構」(投稿版の生成時に機械的に除去する)の実装。
# 書式の残作業(組版時): [n] の上付き化・6著者超の et al. 短縮確認・Word 化。

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "output" / "submission"
OUT.mkdir(parents=True, exist_ok=True)

ORG_AUTHORS = ["R Core Team"]  # 多語の団体著者

CJK = re.compile(r"[　-ヿ㐀-鿿＀-￯]")


def load_ledger():
    """引用台帳の書誌列 → (著者キー, 年) 索引。著者キーは 'Surname et al.' / 'Surname and Surname' / 'Surname' / 団体名。"""
    text = (ROOT / "paper" / "citation_ledger.md").read_text(encoding="utf-8")
    index = {}
    for line in text.split("\n"):
        if not line.startswith("| ") or line.startswith("| 書誌") or line.startswith("| ---"):
            continue
        biblio = line.split("|")[1].strip()
        m = re.match(r"([A-ZÀ-Þ][\w'’-]*)", biblio)
        # 出版年の抽出は構造を優先する(表題中の年 — 例「1985 NIH Tables」— を拾わない):
        # (1) 巻号の直前「2008;95」 (2) ソフトウェア型「(2025)」 (3) 「; 2003」型 (4) 最初の年
        ym = (re.search(r"\b((?:19|20)\d\d[a-z]?)(?=;\d)", biblio)
              or re.search(r"\(((?:19|20)\d\d[a-z]?)\)", biblio)
              or re.search(r";\s((?:19|20)\d\d[a-z]?)\b", biblio)
              or re.search(r"\b((?:19|20)\d\d[a-z]?)\b", biblio))
        if not m or not ym:
            continue
        year = ym.group(1)
        for org in ORG_AUTHORS:
            if biblio.startswith(org):
                index.setdefault((org, year), biblio)
                break
        else:
            surname = m.group(1)
            # 第2著者の姓(and 形式の照合用)
            # ジャーナル型「Davis S,」「Davis S.」に加えソフトウェア型「Davis S (2025)」も許容
            second = re.search(r"^[A-ZÀ-Þ][\w'’-]+ [A-Z]{1,3}, ([A-ZÀ-Þ][\w'’-]+) [A-Z]{1,3}(?=[,.]| \()", biblio)
            # 同姓同年の衝突は先勝ち(setdefault)— and 形式キーが後続行を区別する
            # (例: Morgan et al. 2025 = SummarizedExperiment / Morgan and Davis 2025 = GDC)
            index.setdefault((surname, year), biblio)
            if second:
                index.setdefault((f"{surname}+{second.group(1)}", year), biblio)
    return index


CITE_PART = re.compile(
    r"(?:^|;\s*|\(\s*)"
    r"(?P<pre>[^();]*?)"
    r"(?P<auth>(?:R Core Team)|(?:[A-ZÀ-Þ][\w'’-]+(?: et al\.| and [A-ZÀ-Þ][\w'’-]+)?))"
    r"[,]? (?P<years>(?:19|20)\d\d[a-z]?(?:,\s*(?:19|20)\d\d[a-z]?)*)"
    r"\s*(?=$|;|\))"
)


def auth_key(auth):
    if auth in ORG_AUTHORS:
        return auth
    m = re.match(r"([A-ZÀ-Þ][\w'’-]+) and ([A-ZÀ-Þ][\w'’-]+)$", auth)
    if m:
        return f"{m.group(1)}+{m.group(2)}"
    return auth.split(" ")[0]


def convert(text, index, order, warnings, do_convert=True):
    """括弧内・地の文の author-year を [n] へ。order は初出順の共有リスト。"""

    def assign(auth, year):
        key = (auth_key(auth), year)
        if key not in index:
            # and 形式が第2著者索引に無い場合は筆頭姓で再試行
            key = (auth.split(" ")[0], year)
        if key not in index:
            warnings.append(f"台帳不一致: {auth} {year}")
            return None
        biblio = index[key]
        if biblio not in order:
            order.append(biblio)
        return order.index(biblio) + 1

    def paren_repl(m):
        inner = m.group(1)
        if not re.search(r"\b(19|20)\d\d[a-z]?\b", inner):
            return m.group(0)
        parts = [p.strip() for p in inner.split(";")]
        nums, keep = [], []
        for p in parts:
            pm = re.match(
                r"^(?:(?P<pre>[^,]*?[;,]\s*)?)"
                r"(?P<auth>(?:R Core Team)|(?:[A-ZÀ-Þ][\w'’-]+(?: et al\.| and [A-ZÀ-Þ][\w'’-]+)?))"
                r",? (?P<years>(?:19|20)\d\d[a-z]?(?:,\s*(?:19|20)\d\d[a-z]?)*)$", p)
            if pm:
                for y in re.split(r",\s*", pm.group("years")):
                    n = assign(pm.group("auth"), y)
                    if n:
                        nums.append(n)
                    else:
                        keep.append(p)
                pre = (pm.group("pre") or "").strip().rstrip(";,")
                if pre:
                    keep.insert(0, pre)
            else:
                keep.append(p)
        if not nums:
            return m.group(0)
        numtxt = "[" + ",".join(str(n) for n in sorted(set(nums))) + "]"
        if keep:
            return "(" + "; ".join(keep) + ") " + numtxt
        return numtxt

    def narrative_repl(m):
        auth, years = m.group(1), m.group(2)
        nums = []
        for y in re.split(r",\s*", years):
            n = assign(auth, y)
            if n:
                nums.append(n)
        if not nums:
            return m.group(0)
        return auth + " [" + ",".join(str(n) for n in sorted(set(nums))) + "]"

    if do_convert:
        # 地の文: "Morton et al. (2021)" / "Abend et al. (2012, 2013)"
        text = re.sub(
            r"\b((?:R Core Team)|(?:[A-ZÀ-Þ][\w'’-]+(?: et al\.| and [A-ZÀ-Þ][\w'’-]+)?)) "
            r"\(((?:19|20)\d\d[a-z]?(?:,\s*(?:19|20)\d\d[a-z]?)*)\)",
            narrative_repl, text)
        # 括弧内
        text = re.sub(r"\(([^()]*)\)", paren_repl, text)
    return text


def clean_biblio(b):
    """書誌から台帳の管理注記を除く(「、DOI なし」・末尾の「 (PMID n)」)。"""
    b = b.replace("、DOI なし", "")
    return re.sub(r"\s*\(PMID \d+\)", "", b).strip()


def strip_meta(text, drop_sections):
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    lines = text.split("\n")
    out, skip = [], False
    for line in lines:
        if line.startswith("## "):
            skip = any(line[3:].startswith(d) for d in drop_sections)
        if skip:
            continue
        out.append(line)
    text = "\n".join(out)
    # 状態行(日本語メタ)と孤立空行の整理
    text = "\n".join(l for l in text.split("\n") if not (CJK.search(l) and not l.startswith("|")))
    text = re.sub(r"\n{3,}", "\n\n", text)
    text = re.sub(r"  +", " ", text)          # タグ除去痕の連続スペース
    text = re.sub(r" ([.,;)])", r"\1", text)  # 句読点前の孤立スペース
    text = re.sub(r"^[ \t]+", "", text, flags=re.M)  # タグ除去痕の行頭スペース
    text = re.sub(r"[ \t]+$", "", text, flags=re.M)  # タグ除去痕の行末スペース
    return text


def main():
    index = load_ledger()
    warnings, order = [], []

    main_src = (ROOT / "paper" / "draft_manuscript.md").read_text(encoding="utf-8")
    body = strip_meta(main_src, drop_sections=["引用文献", "起草メモ"])
    body = convert(body, index, order, warnings, do_convert=True)

    refs = ["## References", ""]
    # 台帳の管理注記(「、DOI なし」・「(PMID n)」)は投稿用表示から除去する
    # (BJC の Vancouver 書式は DOI までで PMID を含まない — フェーズ2 整形 2026-08-29)
    refs += [f"{i+1}. {clean_biblio(b)}" for i, b in enumerate(order)]
    submission = body.rstrip("\n") + "\n\n" + "\n".join(refs) + "\n"
    (OUT / "manuscript_submission.md").write_text(submission, encoding="utf-8")

    supp_src = (ROOT / "paper" / "supplementary_material.md").read_text(encoding="utf-8")
    supp = strip_meta(supp_src, drop_sections=[])
    # SI は独立文書として自己完結させる: 本文と同じ照合機構で Supp 内の author-year
    # 引用を収集し(変換結果は捨てる)、書誌をアルファベット順で末尾に付す。
    # 台帳に書誌のない引用は既存の「台帳不一致」警告に乗る。
    supp_order = []
    convert(supp, index, supp_order, warnings, do_convert=True)
    if supp_order:
        supp = (supp.rstrip("\n") + "\n\n## Supplementary References\n\n"
                + "\n\n".join(clean_biblio(b)
                              for b in sorted(supp_order, key=str.lower)) + "\n")
    (OUT / "supplementary_submission.md").write_text(supp, encoding="utf-8")

    # 残存 author-year 候補(未変換の検出 — References 節より前の本文のみ)
    body_only = submission[: submission.index("## References")]
    residual = [x.group(0) for x in re.finditer(
        r"[A-ZÀ-Þ][\w'’-]+(?: et al\.| and [A-ZÀ-Þ][\w'’-]+)? \(?(?:19|20)\d\d[a-z]?\)?", body_only)]
    residual = [r for r in residual if not re.match(r"^(Data|Table|Figure|Fig)", r)]
    report = [
        f"変換済み引用(一意): {len(order)}",
        f"除去タグ数: main {len(re.findall(r'<!--', main_src))} / supp {len(re.findall(r'<!--', supp_src))}",
        f"警告: {len(warnings)}",
        *[f"  - {w}" for w in warnings],
        f"未変換の author-year 残存候補: {len(residual)}",
        *[f"  - {r}" for r in residual],
    ]
    (OUT / "warnings.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report[:6]))
    print(f"出力: {OUT}")


if __name__ == "__main__":
    sys.exit(main())
