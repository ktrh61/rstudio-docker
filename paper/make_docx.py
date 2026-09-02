#!/usr/bin/env python3
# 投稿用 Word 生成(派生物 — 正本は不変)
#
# 入力 = output/submission/manuscript_submission.md(make_submission.py の出力)
# 出力 = output/submission/manuscript_submission.docx
#
# BJC 書式(GTA ライブ照合 2026-08-25)を機械注入する:
#   (1) BJC 節順に組み替え: タイトルページ → Abstract → Background → Methods → Results →
#       Discussion → Additional Information → References → Figure legends → Tables → Figures
#       (改ページつき。表は投稿形の列名・脚注、添字は subscript 記法で表示層整形)
#   (2) 引用番号 [n] を上付きへ(pandoc の ^…^ 記法へ前処理)
#   (3) 1.5 行間(styles.xml の Normal に w:line=360)
#   (4) 全行番号(sectPr に lnNumType)
#   (5) 全ページ番号(PAGE フィールドのフッタを新設)
# 変換後、docx 本文のトークン集合を入力 md と照合し内容同一性を検査する。
#
# カバーレター(paper/cover_letter.md と参考訳 cover_letter_ja.md)も同じタグで書簡体裁の
# docx にする(行番号・頁番号なし、管理メモは除去)。
# 日本語参考訳(paper/manuscript_ja.md / paper/supplementary_ja.md)も同時に
# Word 化し、英語版と同一時刻タグで対応づける(冒頭に対応する英語版ファイル名を
# 自動記載)。和文書体 = 本文 BIZ UD明朝 Medium 11pt+Times New Roman、見出し
# BIZ UDゴシック太字(研究者裁定 2026-08-27 — 共有は PDF 前提でフォント埋め込み
# されるため他環境の置換懸念なし、画面可読性を優先)。行番号は付けない。
#
# pandoc は環境変数 PANDOC で指定(既定 pandoc 3.6.4 を想定。バイナリは
# リポジトリ外に置く — コミットされるのは本スクリプトのみ)。

import os
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SUB = ROOT / "output" / "submission"
PANDOC = os.environ.get("PANDOC", "pandoc")

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"

# ヒト遺伝子記号のイタリック対象(研究者裁定 2026-08-25。PTC 単体は疾患略語で対象外)
GENES = ("RET|BRAF|CCDC6|NCOA4|CLIP2|BHLHB9|S100A10|TESC|EHD4|"
         "ATP5MF|MRPL52|NTHL1|URM1|USE1|PXDN|P2RY1|PLK2|PTC1|PTC3")

# 閲覧用の図埋め込み(研究者裁定 2026-08-27 — 表と同じ一体ファイル方針の延長。
# 投稿時は BJC 規定どおり図を個別 TIFF/PDF で分離するため閲覧モード専用の一方向変換)
FIGS = [
    ("Figure 1", "図 1", "fig_cohort_flow.png"),
    ("Figure 2", "図 2", "fig_gene_bm_evidence.png"),
    ("Figure 3", "図 3", "fig_reo_grading.png"),
    ("Figure S1", "図 S1", "fig_ma_gene_bm.png"),
    ("Figure S2", "図 S2", "fig_d6_calibration.png"),
]
FIG_PATH = {en: ROOT / "output" / "figures" / f for en, _, f in FIGS}

# 投稿形の改ページ(pandoc raw openxml — 閲覧 docx の組版のみ)
PAGE_BREAK = "\n\n```{=openxml}\n<w:p><w:r><w:br w:type=\"page\"/></w:r></w:p>\n```\n\n"

PG_MAR = ('<w:pgMar w:top="1417" w:right="1417" w:bottom="1417" w:left="1417" '
          'w:header="709" w:footer="709" w:gutter="0"/>')
FOOTER_REF = ('<w:footerReference xmlns:r="http://schemas.openxmlformats.org/officeDocument/'
              '2006/relationships" w:type="default" r:id="rIdFooterPg"/>')


def section_break(landscape):
    """直前までの区間を閉じるセクション区切り(A4、向きは引数)。フッタ(頁番号)と
    行番号は各区間に明示(Word はセクションごとに保持するため)。"""
    pg = ('<w:pgSz w:w="16838" w:h="11906" w:orient="landscape"/>' if landscape
          else '<w:pgSz w:w="11906" w:h="16838"/>')
    return ("\n\n```{=openxml}\n<w:p><w:pPr><w:sectPr>" + FOOTER_REF
            + '<w:type w:val="nextPage"/>' + pg + PG_MAR
            + '<w:lnNumType w:countBy="1" w:restart="continuous"/>'
            + "</w:sectPr></w:pPr></w:p>\n```\n\n")


def embed_figures(text):
    """(日本語版) 図凡例行の直下へ対応 PNG を幅 160mm(A4 余白内)で挿入する。"""
    for en, ja, fname in FIGS:
        img = f"![]({ROOT / 'output' / 'figures' / fname}){{width=160mm}}"
        for label in (en, ja):
            text = re.sub(rf"(^\*\*{label} \|[^\n]*\n)",
                          lambda m, i=img: m.group(1) + "\n" + i + "\n\n",
                          text, count=1, flags=re.M)
    return text


def typeset_symbols(text):
    """投稿形の添字(表示層のみ — 正本は素の ASCII 表記を維持): π0 → π₀、α0、log2/log10、
    C(n,nx) → C(n, n_X)。pandoc の subscript 記法へ。"""
    text = text.replace("π0", "π~0~").replace("α0", "α~0~")
    text = re.sub(r"\blog(2|10)(?!\d)", r"log~\1~", text)
    text = text.replace("C(n,nx)", "C(n, n~X~)")
    text = re.sub(r"\bnx\b", "n~X~", text)
    return text


def _num(v, k):
    return f"{float(v):.{k}f}".replace("-", "−")


def _sci(v):
    m, e = f"{float(v):.2e}".split("e")
    e = int(e)
    return f"{m} × 10^{'−' if e < 0 else ''}{abs(e)}^"


COLLECTION = {"H": "Hallmark", "C2:CP": "C2 canonical pathways",
              "C5:GO:BP": "C5 GO Biological Process",
              "C2:CGP:radiation": "Radiation-curated (C2:CGP)"}

# 表の投稿形整形(表示層 — 凍結 CSV の列名・値は不変。研究者指示 2026-08-29: 共著者
# レビュー前に「未完成に見える」列名を最終形へ)。丸めは本文の桁に合わせる。
TABLE_SPECS = {
    "tab_case_characteristics.csv": dict(
        headers={"group": "Group", "n": "n", "pool": "Stratum × band pool, n (paired)",
                 "female": "Female", "male": "Male",
                 "age_surgery": "Age at surgery, years, median [range]",
                 "age_exposure": "Age at exposure, years, median [range]",
                 "driver_detail": "Designated driver, n"},
        rows={"R_Sporadic": "R_Sporadic (dose-zero)", "R_Low": "R_Low (Low-AS)",
              "R_Mid": "R_Mid (Mid-AS)", "R_High": "R_High (High-AS)",
              "B_Sporadic": "B_Sporadic (dose-zero)", "B_High": "B_High (High-AS)"}),
    "tab_gene_level_summary.csv": dict(
        headers={"unit": "Contrast", "n_tested": "Genes tested", "pi0": "π0",
                 "deg_q10": "Genes at q<0.10", "up": "Higher in High-AS",
                 "down": "Lower in High-AS", "min_p_exact": "Minimum exact p",
                 "hc_p": "Higher Criticism p"},
        fmt={"min_p_exact": _sci}),
    "table_s1_normalization_diagnostics.csv": dict(
        headers={"contrast": "Contrast", "reference_group": "Reference group",
                 "high_group": "High-AS group", "n_reference": "n (reference)",
                 "n_high": "n (High-AS)", "protein_coding_genes": "Protein-coding genes",
                 "genes_after_filterByExpr": "Genes after filterByExpr",
                 "deges_iterations": "DEGES iterations",
                 "final_screen_pi0": "Screen π0 (iteration 3)",
                 "final_jaccard": "Jaccard index (iteration 3)",
                 "norm_factor_min": "Normalization factor, minimum",
                 "norm_factor_max": "Normalization factor, maximum"}),
    "table_s2_software_versions.csv": dict(headers={"package": "Package", "version": "Version"}),
    "table_s3_gene_set_summary.csv": dict(
        headers={"unit": "Contrast", "collection": "Collection",
                 "n_sets": "Sets tested (15–500 genes)", "min_q_bh": "Minimum BH q"},
        values={"collection": COLLECTION}),
    "table_s4_complete_null_calibration.csv": dict(
        headers={"unit": "Contrast", "collection": "Collection", "m_sets": "Sets tested",
                 "replicates": "Pseudo-observations",
                 "n_any_discovery": "With ≥1 discovery, n", "p_any": "Proportion",
                 "mean_discoveries": "Mean discoveries",
                 "max_discoveries": "Maximum discoveries"},
        values={"collection": COLLECTION},
        merge=[("95% Clopper–Pearson interval", ["ci_lo", "ci_hi"],
                lambda lo, hi: f"{_num(lo, 3)}–{_num(hi, 3)}")]),
    "table_s6_between_stratum_concordance.csv": dict(
        headers={"pair": "Tissue", "units": "Contrasts", "n_shared_genes": "Shared genes",
                 "rho": "Spearman ρ", "p_two_sided": "Two-sided shuffle p",
                 "n_perm": "Shuffles"},
        fmt={"rho": lambda v: _num(v, 3)},
        values={"pair": {"normal": "Normal tissue", "tumor": "Tumor tissue"}},
        merge=[("Central 95% shuffle interval", ["interval_lo", "interval_hi"],
                lambda lo, hi: f"{_num(lo, 2)} to {_num(hi, 2)}")]),
    "table_s7_reo_panel.csv": dict(
        headers={"pair_id": "Pair", "up": "Higher-expression gene (Ensembl)",
                 "up_name": "Symbol", "down": "Lower-expression gene (Ensembl)",
                 "down_name": "Symbol",
                 "median_diff": "Absolute shift in median within-sample log2-TPM difference",
                 "reversal_rate": "High-AS reversal rate",
                 "r0_q10": "10th percentile of absolute difference, dose-zero"},
        fmt={"median_diff": lambda v: _num(v, 3), "reversal_rate": lambda v: _num(v, 2),
             "r0_q10": lambda v: _num(v, 3)},
        italic=["up_name", "down_name"]),
}


def fmt_table(path, spec):
    """凍結 CSV を投稿形の Markdown 表へ(列名の付け替え・値の表示整形のみ)。"""
    import csv
    rows = [r for r in csv.reader(open(path, encoding="utf-8"))]
    head, body = rows[0], rows[1:]
    idx = {c: i for i, c in enumerate(head)}
    merges = spec.get("merge", [])
    merged = {c for _, cols, _ in merges for c in cols}
    cols = []
    for c in head:
        if c in merged:
            for disp, mcols, fn in merges:
                if c == mcols[0]:
                    cols.append((disp, lambda r, mcols=mcols, fn=fn: fn(*[r[idx[x]] for x in mcols])))
            continue
        f = spec.get("fmt", {}).get(c, lambda v: v)
        vmap = spec.get("values", {}).get(c, {})
        ital = c in spec.get("italic", [])

        def cell(r, c=c, f=f, vmap=vmap, ital=ital):
            v = r[idx[c]]
            if v == "NA":
                return "—"
            v = vmap.get(v, f(v))
            if c == head[0]:
                v = spec.get("rows", {}).get(v, v)
            if ital:
                v = f"*{v}*"
            return v.replace("|", "\\|")
        cols.append((spec.get("headers", {}).get(c, c), cell))
    out = ["| " + " | ".join(d.replace("|", "\\|") for d, _ in cols) + " |",
           "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(fn(r) for _, fn in cols) + " |" for r in body]
    return "\n".join(out)


def csv_md_table(path):
    """凍結 CSV を Markdown 表へ(整形仕様があれば投稿形、無ければ素のまま)。"""
    spec = TABLE_SPECS.get(Path(path).name)
    if spec:
        return fmt_table(path, spec)
    import csv
    rows = [[c.replace("|", "\\|") for c in r]
            for r in csv.reader(open(path, encoding="utf-8"))]
    out = ["| " + " | ".join(rows[0]) + " |",
           "|" + "---|" * len(rows[0])]
    out += ["| " + " | ".join(r) + " |" for r in rows[1:]]
    return "\n".join(out)


def _declarations():
    return (ROOT / "paper" / "submission_declarations.md").read_text(encoding="utf-8")


def author_block():
    """submission_declarations.md の記入値から著者ブロックを描画する(一方向)。"""
    d = _declarations()
    sec = d.split("- **Authors / Affiliations**:", 1)[1]
    authors = re.findall(
        r"^\s+\d+\. (.+?) \((\d+(?:,\s*\d+)*)\)( \*corresponding)?\s*$", sec, re.M)
    aff_zone = sec.split("Affiliations:", 1)[1]
    affs = re.findall(r"^\s+(\d+)\. ([A-Z][^\n]+)$", aff_zone, re.M)
    corr = re.search(r"\*\*Corresponding author[^\n]*\n\s+([^\n]+)", d).group(1)
    names = ", ".join(
        n + "^" + a.replace(" ", "") + "^" + ("\\*" if c else "")
        for n, a, c in authors)
    lines = [names, ""]
    lines += [f"^{n}^ {a}" for n, a in affs]
    lines += ["", "\\* Correspondence: " + corr]
    return "\n".join(lines)


def title_page():
    """BJC のタイトルページ(題名・全著者と所属・対応著者の e-mail・ORCID)。"""
    d = _declarations()
    title = re.search(r"\*\*Title\*\*: ([^(\n]+)", d).group(1).strip()
    return "# " + title + "\n\n" + author_block()


CJK = re.compile(r"[　-ヿ㐀-鿿＀-￯]")


def additional_information():
    """declarations の Additional Information を BJC 規定順で描画する。
    日本語管理注記(CJK 行)は除去、未記入節(【記入】)は明示プレースホルダ、
    URL/DOI の未確定括弧も明示プレースホルダへ置換。"""
    d = _declarations()
    zone = d.split("## Additional Information", 1)[1].split("\n## ", 1)[0]
    zone = zone.split("\n", 1)[1]
    out = ["## Additional Information"]
    for m in re.finditer(r"^### (.+?)$\n(.*?)(?=^### |\Z)", zone, re.S | re.M):
        head, body = m.group(1).strip(), m.group(2)
        body = re.sub(r"<!--.*?-->", "", body, flags=re.S)
        # 管理接頭辞「案: 」は CJK 行除去より前に外す(先頭行の英文が行ごと脱落する不具合の修正 2026-09-02)
        body = re.sub(r"^案: ", "", body, flags=re.M)
        body = re.sub(r"【リポジトリ URL/DOI[^】]*】", "[URL/DOI to be added at publication]", body)
        body = re.sub(r"【公開リポジトリの URL/DOI[^】]*】", "", body)
        lines = [l for l in body.split("\n") if not CJK.search(l)]
        body = re.sub(r"\n{3,}", "\n\n", "\n".join(lines)).strip()
        if not body or "【記入】" in body:
            body = "[To be completed before submission]"
        out += [f"### {head}", "", body, ""]
    return "\n".join(out).rstrip()


def _finish(text):
    """英語版共通の仕上げ: 変異表記の上付き・遺伝子イタリック・添字。"""
    text = text.replace("BRAF V600E", "*BRAF*^V600E^")
    text = re.sub(rf"(?<![\w*])({GENES})(?![\w*])", r"*\1*", text)
    return typeset_symbols(text)


def _legend_paras(section):
    return [p for p in section.split("\n\n") if p.strip() and not p.startswith("## ")]


SUPP_TABLE_FILES = {
    "S1": "table_s1_normalization_diagnostics.csv",
    "S2": "table_s2_software_versions.csv",
    "S3": "table_s3_gene_set_summary.csv",
    "S4": "table_s4_complete_null_calibration.csv",
    "S5": None,  # ORA(18,576 行)は別ファイル供給
    "S6": "table_s6_between_stratum_concordance.csv",
    "S7": "table_s7_reo_panel.csv",
}
PORTRAIT_TEXT_MM = 160  # A4 縦置きの本文幅(余白 25 mm)
LANDSCAPE_OVER_MM = 140  # 自然幅がこれを超えたら横置き(見出しの折返しで窮屈になる手前で切替)


def table_width_mm(md, char_mm=1.8, pad_mm=2.5, wrap_cap=16):
    """Markdown 表の自然幅(mm)を見積もる: 各列 = max(折返し不能な最長セル, 見出し)。
    数値・ID(空白なし)は折り返せないので実長、文章セルは wrap_cap 字で折り返す前提、
    見出しは 2 行までの折返しを許容(最長語と全長の半分の大きい方)。10 pt 前提。"""
    rows = [[c.strip() for c in r.strip("|").split("|")]
            for r in md.split("\n") if r.startswith("|") and not r.startswith("|---")]
    head, body = rows[0], rows[1:]
    total = 0.0
    for j, h in enumerate(head):
        cells = [r[j] for r in body if j < len(r)]
        cell_w = max([len(c) if not re.search(r"\s", c) else min(len(c), wrap_cap) for c in cells] or [0])
        words = h.split() or [""]
        head_w = max(max(len(w) for w in words), (len(h) + 1) // 2)
        total += max(cell_w, head_w) * char_mm + pad_mm
    return total


def is_landscape(block):
    """表ブロック(キャプション+Markdown 表)の向き。表が無ければ縦置き。"""
    m = re.search(r"(?:^\|.*\n)+", block, re.M)
    return bool(m) and table_width_mm(m.group(0)) > LANDSCAPE_OVER_MM


def orient_blocks(blocks, heading):
    """1 頁 1 表で並べ、向きが変わる境目にセクション区切りを置く。blocks は
    (label, block) のリスト。返り値は連結済み文字列(末尾で縦置きへ戻す)。"""
    parts, landscape = [], False
    for k, (label, block) in enumerate(blocks):
        want = is_landscape(block)
        parts.append(section_break(landscape) if want != landscape else PAGE_BREAK)
        landscape = want
        parts.append((heading + "\n\n" if k == 0 else "") + block)
        print(f"  {label}: ~{table_width_mm(re.search(r'(?:^\|.*\n)+', block, re.M).group(0)) if re.search(r'(?:^\|.*\n)+', block, re.M) else 0:.0f} mm -> {'landscape' if want else 'portrait'}")
    parts.append(section_break(True) if landscape else PAGE_BREAK)
    return "".join(parts)
SUPP_DATA_FILES = {"1": "supplementary_data_1_gene_level_results.csv",
                   "2": "supplementary_data_2_set_level_results.csv"}


def supp_preprocess(text):
    """Supplementary Material の投稿形整形: 表紙 → 本文 → Supplementary References →
    図(凡例+図、1 ページ 1 図)→ 表(キャプション+表、1 ページ 1 表)→ Data のキャプション。"""
    d = _declarations()
    title = re.search(r"\*\*Title\*\*: ([^(\n]+)", d).group(1).strip()
    cover = "# Supplementary Material\n\n**" + title + "**\n\n" + author_block() + PAGE_BREAK
    text = re.sub(r"^-{3,}\s*$", "", text, flags=re.M)
    leg = re.search(r"^## Supplementary figure and table legends$.*?(?=^## |\Z)", text, re.S | re.M).group(0)
    text = text.replace(leg, "")
    refs = re.search(r"^## Supplementary References$.*?(?=^## |\Z)", text, re.S | re.M).group(0)
    text = text.replace(refs, "")
    paras = _legend_paras(leg)
    supp_files = ROOT / "paper" / "gpt_review" / "supplementary_files"
    figs, tabs, data = [], [], []
    for p in paras:
        if p.startswith("**Figure S"):
            n = re.match(r"\*\*(Figure S\d)", p).group(1)
            figs.append(p + "\n\n" + f"![]({FIG_PATH[n]}){{width=160mm}}")
        elif p.startswith("**Table S"):
            tag = re.match(r"\*\*Table (S\d)", p).group(1)
            fname = SUPP_TABLE_FILES[tag]
            if fname:
                tabs.append((tag, p + "\n\n" + csv_md_table(supp_files / fname)))
            else:
                tabs.append((tag, p + "\n\n*Provided as a separate file (table_s5_ora_annotation.csv).*"))
        elif p.startswith("**Supplementary Data"):
            n = re.match(r"\*\*Supplementary Data (\d)", p).group(1)
            data.append(p + f"\n\n*Provided as a separate file ({SUPP_DATA_FILES[n]}).*")
    # 補足表は 1 頁 1 表。向きは表の自然幅(セル内容+見出し)で決め、境目にセクション区切り。
    assembled = (cover + text.rstrip("\n") + "\n\n" + refs.rstrip("\n")
                 + PAGE_BREAK + "## Supplementary figures\n\n" + PAGE_BREAK.join(figs)
                 + orient_blocks([(f"Table {tag}", b) for tag, b in tabs], "## Supplementary tables")
                 + "## Supplementary data\n\n" + "\n\n".join(data)
                 + PAGE_BREAK + supp_file_descriptions() + "\n")
    return _finish(assembled)


def supp_file_descriptions():
    """投稿システムに入力する各補足ファイルの ≤50 語要約(declarations 記載)を、
    共著者レビュー用に Supp 末尾へ描画する。"""
    d = _declarations()
    zone = d.split("## Supplementary file descriptions", 1)[1].split("\n## ", 1)[0]
    out = ["## Supplementary file descriptions (submission-system summaries, ≤50 words each)"]
    for m in re.finditer(r"^### (.+?)$\n\n(.+?)$", zone, re.M):
        out += ["", f"**{m.group(1)}.** {m.group(2).strip()}"]
    return "\n".join(out)


def preprocess(text):
    """本文の投稿形整形(BJC GTA 節順): タイトルページ → Abstract → Background → Methods →
    Results → Discussion → Additional Information → References → Figure legends →
    Tables(1 ページ 1 表)→ Figures(1 ページ 1 図)。引用 [n] は上付き。"""
    text = re.sub(r"^-{3,}\s*$", "", text, flags=re.M)
    text = re.sub(r"^\*\*Title:\*\*[^\n]*\n", lambda m: title_page() + PAGE_BREAK, text,
                  count=1, flags=re.M)
    leg = re.search(r"^## Figure legends and table captions$.*?(?=^## )", text,
                    flags=re.S | re.M).group(0)
    text = text.replace(leg, "")
    refs = re.search(r"^## References$.*?(?=^## |\Z)", text, flags=re.S | re.M).group(0)
    text = text.replace(refs, "")
    # References の番号はリスト化させず静的テキストとして残し(再採番事故の予防)、
    # 各文献を 1 件 1 段落にする
    fixed = re.sub(r"^(\d+)\. ", r"\1\\. ", refs, flags=re.M)
    fixed = re.sub(r"\n(?=\d+\\\. )", "\n\n", fixed)
    paras = _legend_paras(leg)
    fig_legends = [p for p in paras if p.startswith("**Figure ")]
    blocks, cur = [], None
    for p in paras:
        if p.startswith("**Table "):
            cur = [p]
            blocks.append(cur)
        elif p.startswith("**Figure "):
            cur = None
        elif cur is not None:
            cur.append(p)
    # Table 2・3 の実体(凍結 CSV)をキャプション直下へ(脚注は CSV 表の後ろに残る)
    for b in blocks:
        if b[0].startswith("**Table 2 |"):
            b.insert(1, csv_md_table(ROOT / "output" / "tables" / "tab_case_characteristics.csv"))
        if b[0].startswith("**Table 3 |"):
            b.insert(1, csv_md_table(ROOT / "output" / "tables" / "tab_gene_level_summary.csv"))
    tables_md = orient_blocks(
        [(re.match(r"\*\*(Table \d)", b[0]).group(1), "\n\n".join(b)) for b in blocks], "## Tables")
    figures_md = PAGE_BREAK.join(
        f"**{n}**\n\n![]({FIG_PATH[n]}){{width=160mm}}" for n in ("Figure 1", "Figure 2", "Figure 3"))
    text = (text.rstrip("\n") + "\n\n" + additional_information() + "\n\n" + fixed.rstrip("\n")
            + PAGE_BREAK + "## Figure legends\n\n" + "\n\n".join(fig_legends)
            + tables_md
            + "## Figures\n\n" + figures_md + "\n")
    text = re.sub(r"\[(\d+(?:,\d+)*)\]", r"^\1^", text)  # 本文中の [1] / [1,2] を上付き化
    return _finish(text)


def letter_preprocess(text, ja=False, pair=None):
    """カバーレター(書簡体裁): md 冒頭の管理メモ(見出し・状態・確認事項の箇条書き)を
    除き、区切り線以降の本文だけを渡す。和文版は対応版の一行を先頭に付す。"""
    body = text.split("\n---\n", 1)[1] if "\n---\n" in text else text
    body = body.strip("\n") + "\n"
    if ja and pair:
        body = pair + "\n\n" + body
    return body


def ja_preprocess(text, en_stem, tag, commit, src_label):
    """日本語参考訳の閲覧用整形: 対応する英語版 docx の名前を冒頭に記載し
    (同一時刻タグで同時生成される対)、遺伝子イタリック・V600E 上付き・添字を
    英語版と同じ規則で適用する。引用 [n] は素のまま(英語版 References 対応)。"""
    text = re.sub(r"^-{3,}\s*$", "", text, flags=re.M)
    pair = (f"**対応版**: 英語版 {en_stem}_{tag}.docx(同時生成の対)。"
            f"ソース: {src_label} @{commit}")
    text = re.sub(r"^(# [^\n]+\n)", lambda m: m.group(1) + "\n" + pair + "\n",
                  text, count=1)
    text = text.replace("BRAF V600E", "*BRAF*^V600E^")
    text = re.sub(rf"(?<![\w*])({GENES})(?![\w*])", r"*\1*", text)
    return typeset_symbols(embed_figures(text))


def patch_docx(path, ja=False, letter=False):
    """styles.xml(1.5 行間)・document.xml(行番号・フッタ参照)・
    footer1.xml(ページ番号)を注入する。ja=True では和文書体
    (本文 游明朝+Times New Roman・見出し 游ゴシック・両端揃え)とし、
    行番号は付けない。"""
    zin = zipfile.ZipFile(path)
    parts = {n: zin.read(n) for n in zin.namelist()}
    zin.close()

    styles = parts["word/styles.xml"].decode("utf-8")
    normal = re.search(
        r'<w:style [^>]*w:styleId="Normal"[^>]*>.*?</w:style>', styles, re.S
    ).group(0)
    # 和文版は固定行送り 19pt(atLeast)にする — 和文書体は縦メトリクスが大きく、
    # 倍率指定(auto)では見かけの行間が過大になるため。揃えは左のまま
    # (両端揃えは長い英字トークンで文字間が間延びする)。
    spacing = ('<w:spacing w:after="120" w:line="380" w:lineRule="atLeast"/>'
               if ja else
               '<w:spacing w:after="160" w:line="360" w:lineRule="auto"/>')
    if "<w:pPr>" in normal:
        new_normal = re.sub(r"<w:spacing[^/]*/>", "", normal)
        new_normal = new_normal.replace("<w:pPr>", "<w:pPr>" + spacing, 1)
    else:
        # スキーマ上 pPr は style 末尾側 — 閉じタグ直前に挿入
        new_normal = normal.replace(
            "</w:style>", "<w:pPr>" + spacing + "</w:pPr></w:style>"
        )
    styles = styles.replace(normal, new_normal)
    # 原稿慣例の書体へ: 英語版は本文 Times New Roman 12pt、和文版は
    # BIZ UD明朝 Medium 11pt(欧文 Times New Roman — 英語版と同一の欧文見えを保つ)。
    # テーマ経由の解決(日本語環境では Aptos に落ちる)を全経路で遮断するため、
    # eastAsia も含む明示指定を docDefaults・テーマにも適用する。
    ea = "BIZ UD明朝 Medium" if ja else "Times New Roman"
    tnr = ('<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" '
           f'w:eastAsia="{ea}" w:cs="Times New Roman"/>')
    styles = re.sub(r"<w:rFonts [^/]*Theme[^/]*/>", tnr, styles)
    # 派生スタイル(BodyText 等)の spacing が Normal の行間を要素ごと上書きする —
    # w:line を持たない全 spacing に同じ行送りを明示付与
    rule, pitch = ("atLeast", "380") if ja else ("auto", "360")
    styles = re.sub(r'<w:spacing (?![^/>]*w:line=)([^/>]*)/>',
                    rf'<w:spacing \1 w:line="{pitch}" w:lineRule="{rule}"/>', styles)
    normal2 = re.search(
        r'<w:style [^>]*w:styleId="Normal"[^>]*>.*?</w:style>', styles, re.S
    ).group(0)
    bsz = "22" if ja else "24"  # 和文 11pt / 英文 12pt
    styles = styles.replace(
        normal2,
        normal2.replace(
            "</w:style>",
            "<w:rPr>" + tnr + f'<w:sz w:val="{bsz}"/><w:szCs w:val="{bsz}"/></w:rPr>'
            "</w:style>",
        ),
    )
    # 見出し: 英語版は本文と同サイズの太字(節)/太字イタリック(小節)。
    # 和文版は 明朝本文+ゴシック見出し の標準対(14/12.5/11pt 太字、イタリック不使用)
    if ja:
        hfont = ('<w:rFonts w:ascii="BIZ UDゴシック" w:hAnsi="BIZ UDゴシック" '
                 'w:eastAsia="BIZ UDゴシック" w:cs="BIZ UDゴシック"/>')
        hsizes = (("Heading1", "28"), ("Heading2", "25"), ("Heading3", "22"))
    else:
        hfont = tnr
        hsizes = (("Heading1", "24"), ("Heading2", "24"), ("Heading3", "24"))
    for hid, size in hsizes:
        m2 = re.search(
            rf'<w:style [^>]*w:styleId="{hid}"[^>]*>.*?</w:style>', styles, re.S
        )
        if not m2:
            continue
        h = m2.group(0)
        h2 = re.sub(r"<w:rFonts[^/]*/>", "", h)
        h2 = re.sub(r"<w:color[^/]*/>", "", h2)
        h2 = re.sub(r"<w:sz[^/]*/>", "", h2)
        h2 = re.sub(r"<w:szCs[^/]*/>", "", h2)
        ins = (hfont + '<w:b/>'
               + ('<w:i/>' if hid == "Heading3" and not ja else '')
               + '<w:color w:val="000000"/>'
               f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/>')
        if "<w:rPr>" in h2:
            h2 = h2.replace("<w:rPr>", "<w:rPr>" + ins, 1)
        else:
            h2 = h2.replace("</w:style>", "<w:rPr>" + ins + "</w:rPr></w:style>")
        styles = styles.replace(h, h2)
    if ja:
        # 冒頭メタ情報(箇条書き = 両訳ともヘッダブロックのみが該当。pandoc は
        # タイトなリストへ Compact を割り当てる)は 9pt 灰色・行送り 14pt に
        # 圧縮し、本文より控えめに見せる
        m3 = re.search(
            r'<w:style [^>]*w:styleId="Compact"[^>]*>.*?</w:style>', styles, re.S)
        if m3:
            lp = m3.group(0)
            lp2 = lp.replace('w:line="380" w:lineRule="atLeast"',
                             'w:line="280" w:lineRule="atLeast"')
            lp2 = re.sub(r"<w:sz[^/]*/>", "", lp2)
            lp2 = re.sub(r"<w:szCs[^/]*/>", "", lp2)
            ins = (tnr + '<w:color w:val="595959"/>'
                   '<w:sz w:val="18"/><w:szCs w:val="18"/>')
            if "<w:rPr>" in lp2:
                lp2 = lp2.replace("<w:rPr>", "<w:rPr>" + ins, 1)
            else:
                lp2 = lp2.replace("</w:style>",
                                  "<w:rPr>" + ins + "</w:rPr></w:style>")
            styles = styles.replace(lp, lp2)
    parts["word/styles.xml"] = styles.encode("utf-8")

    doc = parts["word/document.xml"].decode("utf-8")
    # 表のセル本文は 10pt・黒を run に明示(段落スタイル Compact の継承・和文版の
    # メタ箇条書き用 9pt 灰色の影響を受けない。表の見出し行の太字は pandoc 既定)
    def _tbl_runs(m):
        tbl = m.group(0)
        tbl = re.sub(r"<w:rPr>(?![^<]*<w:sz )", '<w:rPr><w:sz w:val="20"/><w:szCs w:val="20"/>'
                     '<w:color w:val="000000"/>', tbl)
        tbl = re.sub(r"<w:r>(?!<w:rPr>)", '<w:r><w:rPr><w:sz w:val="20"/><w:szCs w:val="20"/>'
                     '<w:color w:val="000000"/></w:rPr>', tbl)
        # 表内は 1 倍行送り(本文の 1.5 行間は表には不要 — 縦長の表を 1 頁に収める)
        tbl = re.sub(r'<w:pStyle w:val="Compact" ?/>',
                     '<w:pStyle w:val="Compact"/><w:spacing w:after="0" w:line="240" w:lineRule="auto"/>',
                     tbl)
        return tbl
    doc = re.sub(r"<w:tbl>.*?</w:tbl>", _tbl_runs, doc, flags=re.S)
    # 列幅は内容に合わせて自動調整(pandoc は均等固定幅を出し、狭い列で語が分断される)。
    # 罫線は学術表の三線式に揃える: 表の上下に実線(見出し下線は pandoc 既定)
    doc = doc.replace('<w:tblLayout w:type="fixed" />', '<w:tblLayout w:type="autofit" />')
    doc = re.sub(r'(<w:tblW [^>]*/>)',
                 r'\1<w:tblBorders><w:top w:val="single" w:sz="8" w:space="0" w:color="000000"/>'
                 r'<w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000"/></w:tblBorders>',
                 doc)
    # 文書末尾の sectPr(本文中に挿入したセクション区切りではなく最終区間)を対象にする
    sect = re.findall(r"<w:sectPr[^>]*>(?:(?!<w:sectPr).)*?</w:sectPr>", doc, re.S)[-1]
    new_sect = sect
    if "pgSz" not in new_sect:
        # pandoc の sectPr は最小構成 — A4・余白 2.5cm を明示してロケール依存を排し、
        # 英語版のみ続けて lnNumType(スキーマ順: pgSz → pgMar → lnNumType)
        ln = "" if (ja or letter) else '<w:lnNumType w:countBy="1" w:restart="continuous"/>'
        new_sect = new_sect.replace(
            "</w:sectPr>",
            '<w:pgSz w:w="11906" w:h="16838"/>'
            '<w:pgMar w:top="1417" w:right="1417" w:bottom="1417" w:left="1417" '
            'w:header="709" w:footer="709" w:gutter="0"/>'
            + ln + '</w:sectPr>',
        )
    footer_ref = '<w:footerReference xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" w:type="default" r:id="rIdFooterPg"/>'
    if "footerReference" not in new_sect and not letter:
        new_sect = re.sub(r"(<w:sectPr[^>]*>)", r"\1" + footer_ref, new_sect)
    parts["word/document.xml"] = doc.replace(sect, new_sect).encode("utf-8")

    parts["word/footer1.xml"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        f'<w:ftr xmlns:w="{W}"><w:p><w:pPr><w:jc w:val="center"/></w:pPr>'
        '<w:r><w:fldChar w:fldCharType="begin"/></w:r>'
        '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
        '<w:r><w:fldChar w:fldCharType="end"/></w:r></w:p></w:ftr>'
    ).encode("utf-8")

    # 画像の自動圧縮を文書設定で無効化(Word の「ファイル内のイメージを圧縮しない」に相当。
    # 文書単位の設定なので生成時に埋め込み、開くたびの再設定を不要にする。CT_Settings の
    # 順序: clrSchemeMapping の後・decimalSymbol の前)
    st = parts["word/settings.xml"].decode("utf-8")
    if "doNotAutoCompressPictures" not in st:
        flag = "<w:doNotAutoCompressPictures/>"
        st = (st.replace("<w:decimalSymbol", flag + "<w:decimalSymbol", 1)
              if "<w:decimalSymbol" in st else st.replace("</w:settings>", flag + "</w:settings>"))
    parts["word/settings.xml"] = st.encode("utf-8")

    ct = parts["[Content_Types].xml"].decode("utf-8")
    if "footer1.xml" not in ct:
        ct = ct.replace(
            "</Types>",
            '<Override PartName="/word/footer1.xml" ContentType="application/'
            'vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/>'
            "</Types>",
        )
    parts["[Content_Types].xml"] = ct.encode("utf-8")

    rels = parts["word/_rels/document.xml.rels"].decode("utf-8")
    if "rIdFooterPg" not in rels:
        rels = rels.replace(
            "</Relationships>",
            '<Relationship Id="rIdFooterPg" Type="http://schemas.openxml'
            'formats.org/officeDocument/2006/relationships/footer" '
            'Target="footer1.xml"/></Relationships>',
        )
    parts["word/_rels/document.xml.rels"] = rels.encode("utf-8")

    # テーマの既定書体も明示書体へ(テーマ参照の取り残し対策)
    if "word/theme/theme1.xml" in parts:
        theme = parts["word/theme/theme1.xml"].decode("utf-8")
        theme = re.sub(r'(<a:latin typeface=")[^"]*(")',
                       r"\1Times New Roman\2", theme)
        if ja:
            theme = re.sub(r'(<a:ea typeface=")[^"]*(")',
                           r"\1BIZ UD明朝 Medium\2", theme)
        parts["word/theme/theme1.xml"] = theme.encode("utf-8")

    zout = zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED)
    for name, data in parts.items():
        zout.writestr(name, data)
    zout.close()


def tokens_md(text):
    text = re.sub(r"```\{=openxml\}.*?```", " ", text, flags=re.S)  # 改ページ(非テキスト)
    text = re.sub(r"!\[[^\]]*\]\([^)]*\)(\{[^}]*\})?", " ", text)  # 画像構文は非テキスト
    text = re.sub(r"\^([^\^\s]+)\^", r" \1 ", text)  # 上付き(引用番号・V600E・10^−6)
    text = re.sub(r"~([^~\s]+)~", r" \1 ", text)  # 添字(π0・log2・n_X)
    text = re.sub(r"[#*|`^~]", " ", text)
    return re.findall(r"[A-Za-z0-9]+", text)


def tokens_docx(path):
    doc = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8")
    body = " ".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", doc))
    body = (body.replace("&lt;", "<").replace("&gt;", ">")
            .replace("&amp;", "&").replace("&quot;", '"'))
    return re.findall(r"[A-Za-z0-9]+", body)


def build(src_name, out_name, prep, ja=False, letter=False):
    src_path = src_name if isinstance(src_name, Path) else SUB / src_name
    src = src_path.read_text(encoding="utf-8")
    pre = prep(src)
    md = SUB / out_name.replace(".docx", "_input.md")
    md.write_text(pre, encoding="utf-8")
    out = SUB / out_name
    subprocess.run(
        [PANDOC, str(md), "-f", "markdown+superscript+subscript+raw_attribute",
         "-o", str(out)],
        check=True,
    )
    patch_docx(out, ja=ja, letter=letter)
    a, b = sorted(tokens_md(pre)), sorted(tokens_docx(out))
    from collections import Counter
    diff = Counter(a) - Counter(b) | Counter(b) - Counter(a)
    doc = zipfile.ZipFile(out).read("word/document.xml").decode("utf-8")
    styles = zipfile.ZipFile(out).read("word/styles.xml").decode("utf-8")
    media = len([n for n in zipfile.ZipFile(out).namelist()
                 if n.startswith("word/media/")])
    checks = {"頁番号": ("rIdFooterPg" not in doc) if letter else ("rIdFooterPg" in doc),
              "図埋め込み": media == pre.count("![](")}
    if letter:
        checks["行番号なし"] = "lnNumType" not in doc
    elif ja:
        checks["行送り19pt"] = 'w:line="380" w:lineRule="atLeast"' in styles
        checks["行番号なし"] = "lnNumType" not in doc
        checks["和文書体"] = ("BIZ UD明朝 Medium" in styles
                          and "BIZ UDゴシック" in styles)
    else:
        checks["1.5行間"] = 'w:line="360"' in styles
        checks["行番号"] = "lnNumType" in doc
        if "## Additional Information" in pre and "## References" in pre:
            checks["節順"] = (pre.index("## Additional Information") < pre.index("## References")
                            < pre.index("## Figure legends"))
    print(f"{out_name}: token差 {sum(diff.values())} 件"
          + (f" {dict(list(diff.items())[:5])} " if diff else " ")
          + " ".join(f"{k}:{'OK' if v else 'NG'}" for k, v in checks.items()))
    if diff or not all(checks.values()):
        sys.exit(1)
    return out


def main():
    import time
    print(subprocess.run([PANDOC, "--version"], capture_output=True,
                         text=True).stdout.splitlines()[0])
    # Word/OneDrive のキャッシュ・書き戻しを避けるため毎回新しいファイル名で渡す。
    # 同一タグ = 同一実行の 4 ファイル(英語版と日本語訳の対応はタグ一致で判別)
    tag = time.strftime("%Y%m%d_%H%M")
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                            capture_output=True, text=True, cwd=ROOT).stdout.strip()
    outs = [
        build("manuscript_submission.md", "manuscript_submission.docx", preprocess),
        build("supplementary_submission.md", "supplementary_submission.docx",
              supp_preprocess),
        build(ROOT / "paper" / "manuscript_ja.md", "manuscript_ja.docx",
              lambda t: ja_preprocess(t, "manuscript_submission", tag, commit,
                                      "paper/manuscript_ja.md"), ja=True),
        build(ROOT / "paper" / "supplementary_ja.md", "supplementary_ja.docx",
              lambda t: ja_preprocess(t, "supplementary_submission", tag, commit,
                                      "paper/supplementary_ja.md"), ja=True),
        # カバーレター(編集部宛の独立文書 — 原稿には混ぜない。書簡体裁: 行番号・頁番号なし)
        build(ROOT / "paper" / "cover_letter.md", "cover_letter.docx",
              lambda t: letter_preprocess(t), letter=True),
        build(ROOT / "paper" / "cover_letter_ja.md", "cover_letter_ja.docx",
              lambda t: letter_preprocess(t, ja=True, pair=(
                  f"**対応版**: 英語版 cover_letter_{tag}.docx(同時生成の対)。"
                  f"ソース: paper/cover_letter_ja.md @{commit}")), ja=True, letter=True),
    ]
    # 投稿システムの Cover Letter 欄への貼り付け用プレーンテキスト(書式記号なし)
    letter_txt = SUB / "cover_letter.txt"
    subprocess.run([PANDOC, str(SUB / "cover_letter_input.md"), "-f", "markdown+superscript",
                    "-t", "plain", "--wrap=none", "-o", str(letter_txt)], check=True)
    dest = Path("/mnt/c/Users/kotaro/OneDrive/論文関連（説明用資料含）/word_check")
    if dest.parent.exists():
        dest.mkdir(exist_ok=True)
        for out in outs:
            view = dest / f"{out.stem}_{tag}.docx"
            shutil.copy2(out, view)
            print(f"閲覧用コピー: {view}")
        shutil.copy2(letter_txt, dest / f"cover_letter_{tag}.txt")
        print(f"貼り付け用テキスト: {dest / f'cover_letter_{tag}.txt'}")


if __name__ == "__main__":
    main()
