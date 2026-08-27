#!/usr/bin/env python3
# 投稿用 Word 生成(派生物 — 正本は不変)
#
# 入力 = output/submission/manuscript_submission.md(make_submission.py の出力)
# 出力 = output/submission/manuscript_submission.docx
#
# BJC 書式(GTA ライブ照合 2026-08-25)を機械注入する:
#   (1) 図表凡例節を References の後ろへ移動(BJC: legends は References 後の別ページ)
#   (2) 引用番号 [n] を上付きへ(pandoc の ^…^ 記法へ前処理)
#   (3) 1.5 行間(styles.xml の Normal に w:line=360)
#   (4) 全行番号(sectPr に lnNumType)
#   (5) 全ページ番号(PAGE フィールドのフッタを新設)
# 変換後、docx 本文のトークン集合を入力 md と照合し内容同一性を検査する。
#
# 日本語参考訳(paper/manuscript_ja.md / paper/supplementary_ja.md)も同時に
# Word 化し、英語版と同一時刻タグで対応づける(冒頭に対応する英語版ファイル名を
# 自動記載)。和文書体 = 本文 游明朝 10.5pt+Times New Roman、見出し 游ゴシック
# 太字、両端揃え。行番号は付けない(BJC 要件は英語版のみ)。
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


def csv_md_table(path):
    """凍結 CSV を Markdown 表へ(閲覧用の一方向レンダリング — 列名・値は素のまま)。"""
    import csv
    rows = [[c.replace("|", "\\|") for c in r]
            for r in csv.reader(open(path, encoding="utf-8"))]
    out = ["| " + " | ".join(rows[0]) + " |",
           "|" + "---|" * len(rows[0])]
    out += ["| " + " | ".join(r) + " |" for r in rows[1:]]
    return "\n".join(out)


def author_block():
    """submission_declarations.md の記入値から閲覧用の著者ブロックを描画する
    (一方向レンダリング — タイトルページは投稿時に別ファイル化)。"""
    d = (ROOT / "paper" / "submission_declarations.md").read_text(encoding="utf-8")
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


CJK = re.compile(r"[　-ヿ㐀-鿿＀-￯]")


def additional_information():
    """declarations の Additional Information を BJC 規定順で描画する。
    日本語管理注記(CJK 行)は除去、未記入節(【記入】)は見出しごと省略、
    URL/DOI の未確定括弧は明示プレースホルダへ置換。"""
    d = (ROOT / "paper" / "submission_declarations.md").read_text(encoding="utf-8")
    zone = d.split("## Additional Information", 1)[1].split("\n## ", 1)[0]
    zone = zone.split("\n", 1)[1]
    out = ["## Additional Information"]
    for m in re.finditer(r"^### (.+?)$\n(.*?)(?=^### |\Z)", zone, re.S | re.M):
        head, body = m.group(1).strip(), m.group(2)
        body = re.sub(r"<!--.*?-->", "", body, flags=re.S)
        body = re.sub(r"【リポジトリ URL/DOI[^】]*】", "[URL/DOI to be added at publication]", body)
        body = re.sub(r"【公開リポジトリの URL/DOI[^】]*】", "", body)
        lines = [l for l in body.split("\n") if not CJK.search(l)]
        body = re.sub(r"\n{3,}", "\n\n", "\n".join(lines)).strip()
        body = re.sub(r"^案: ", "", body)
        if not body or "【記入】" in body:
            continue  # 未確定節(Acknowledgements 等)は記入後の再生成で自動的に現れる
        out += [f"### {head}", "", body, ""]
    return "\n".join(out).rstrip()


def supp_preprocess(text):
    """Supplementary Material の閲覧用整形: 表紙情報を付与し、遺伝子イタリック・
    V600E 上付き・区切り線除去を本文と同じ規則で適用する。追加情報が必要になったら
    この表紙ブロックに足す。"""
    d = (ROOT / "paper" / "submission_declarations.md").read_text(encoding="utf-8")
    title = re.search(r"\*\*Title\*\*: ([^(\n]+)", d).group(1).strip()
    cover = ("**Supplementary Material for:** " + title + "\n\n" + author_block()
             + "\n\n")
    text = re.sub(r"^-{3,}\s*$", "", text, flags=re.M)
    # 小表(S1・S2・S4〜S8)は出荷用コピーの実体をキャプション直下へ埋め込む。
    # S3(18,577 行)と Data 1/2 は物理的に埋め込み不能 — キャプションのみ(別ファイル参照)
    supp_files = ROOT / "paper" / "gpt_review" / "supplementary_files"
    for tag, fname in [
        ("S1", "table_s1_cohort_composition.csv"),
        ("S2", "table_s2_normalization_diagnostics.csv"),
        ("S4", "table_s4_reo_panel.csv"),
        ("S5", "table_s5_software_versions.csv"),
        ("S6", "table_s6_gene_set_summary.csv"),
        ("S7", "table_s7_complete_null_calibration.csv"),
        ("S8", "table_s8_between_stratum_concordance.csv"),
    ]:
        tbl = csv_md_table(supp_files / fname)
        text = re.sub(rf"(^\*\*Table {tag} \|[^\n]*\n)",
                      lambda m, t=tbl: m.group(1) + "\n" + t + "\n\n",
                      text, count=1, flags=re.M)
    text = text.replace("BRAF V600E", "*BRAF*^V600E^")
    text = re.sub(rf"(?<![\w*])({GENES})(?![\w*])", r"*\1*", text)
    return cover + text


def preprocess(text):
    """凡例節を References 後へ移動し、引用 [n] を上付き記法へ。"""
    # Markdown の区切り線(---)は原稿の可視要素ではない — 罫線化させない
    text = re.sub(r"^-{3,}\s*$", "", text, flags=re.M)
    # タイトル直下に著者ブロック(declarations 記入値の一方向描画)
    blk = author_block()
    text = re.sub(r"(\*\*Title:\*\*[^\n]*\n)",
                  lambda m: m.group(1) + "\n" + blk + "\n\n", text, count=1)
    # 閲覧用: Table 2・3 の実体(凍結 CSV)をキャプション直下へ埋め込む
    # (研究者決定 2026-08-26: 共著者版は一体ファイル — 版ズレ回避・コメント集約。
    #  投稿時は 3 表とも個別ファイルへ分離するため、この埋め込みは閲覧モード専用)
    t2 = csv_md_table(ROOT / "output" / "tables" / "tab_case_characteristics.csv")
    text = re.sub(r"(\*\*Table 2 \| Case characteristics\.\*\*[^\n]*\n)",
                  r"\1\n" + t2.replace("\\", "\\\\") + "\n\n", text, count=1)
    t3 = csv_md_table(ROOT / "output" / "tables" / "tab_gene_level_summary.csv")
    text = re.sub(r"(\*\*Table 3 \| Gene-level results\.\*\*[^\n]*\n)",
                  r"\1\n" + t3.replace("\\", "\\\\") + "\n\n", text, count=1)
    m = re.search(r"^## Figure legends and table captions$.*?(?=^## )", text,
                  flags=re.S | re.M)
    legends = m.group(0)
    text = text.replace(legends, "")
    # BJC 節順: 本文 → References → Additional Information → 凡例・表
    text = (text.rstrip("\n") + "\n\n" + additional_information()
            + "\n\n" + legends.rstrip("\n") + "\n")
    # References の番号はリスト化させず静的テキストとして残し(再採番事故の予防)、
    # 各文献を 1 件 1 段落にする(空行区切り — 連続行は 1 段落に融合してしまうため)
    refs = re.search(r"^## References$.*?(?=^## |\Z)", text, flags=re.S | re.M)
    fixed = re.sub(r"^(\d+)\. ", r"\1\\. ", refs.group(0), flags=re.M)
    fixed = re.sub(r"\n(?=\d+\\\. )", "\n\n", fixed)
    text = text.replace(refs.group(0), fixed)
    # 本文中の [1] / [1,2] を上付き化
    text = re.sub(r"\[(\d+(?:,\d+)*)\]", r"^\1^", text)
    # 変異表記は BJC 現行慣行(2024–25 掲載例)の上付き形へ(正本は不変 — 派生層の組版)
    text = text.replace("BRAF V600E", "*BRAF*^V600E^")
    # ヒト遺伝子記号のイタリック(解析ラベル中の RET/BRAF も含む — 研究者裁定 2026-08-25。
    # 既にイタリック化済みトークンは再包止め)
    text = re.sub(rf"(?<![\w*])({GENES})(?![\w*])", r"*\1*", text)
    return text


def ja_preprocess(text, en_stem, tag, commit, src_label):
    """日本語参考訳の閲覧用整形: 対応する英語版 docx の名前を冒頭に記載し
    (同一時刻タグで同時生成される対)、遺伝子イタリック・V600E 上付きを
    英語版と同じ規則で適用する。引用 [n] は素のまま(英語版 References 対応)。"""
    text = re.sub(r"^-{3,}\s*$", "", text, flags=re.M)
    pair = (f"**対応版**: 英語版 {en_stem}_{tag}.docx(同時生成の対)/"
            f"ソース {src_label} @{commit}")
    text = re.sub(r"^(# [^\n]+\n)", lambda m: m.group(1) + "\n" + pair + "\n",
                  text, count=1)
    text = text.replace("BRAF V600E", "*BRAF*^V600E^")
    text = re.sub(rf"(?<![\w*])({GENES})(?![\w*])", r"*\1*", text)
    return text


def patch_docx(path, ja=False):
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
    # 和文版は固定行送り 18pt(atLeast)にする — 游書体は縦メトリクスが大きく、
    # 倍率指定(auto)では見かけの行間が過大になるため。揃えは左のまま
    # (両端揃えは長い英字トークンで文字間が間延びする)。
    spacing = ('<w:spacing w:after="120" w:line="360" w:lineRule="atLeast"/>'
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
    # 游明朝 10.5pt(欧文 Times New Roman — 英語版と同一の欧文見えを保つ)。
    # テーマ経由の解決(日本語環境では Aptos に落ちる)を全経路で遮断するため、
    # eastAsia も含む明示指定を docDefaults・テーマにも適用する。
    ea = "游明朝" if ja else "Times New Roman"
    tnr = ('<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" '
           f'w:eastAsia="{ea}" w:cs="Times New Roman"/>')
    styles = re.sub(r"<w:rFonts [^/]*Theme[^/]*/>", tnr, styles)
    # 派生スタイル(BodyText 等)の spacing が Normal の行間を要素ごと上書きする —
    # w:line を持たない全 spacing に同じ行送りを明示付与
    rule = "atLeast" if ja else "auto"
    styles = re.sub(r'<w:spacing (?![^/>]*w:line=)([^/>]*)/>',
                    rf'<w:spacing \1 w:line="360" w:lineRule="{rule}"/>', styles)
    normal2 = re.search(
        r'<w:style [^>]*w:styleId="Normal"[^>]*>.*?</w:style>', styles, re.S
    ).group(0)
    bsz = "21" if ja else "24"  # 和文 10.5pt / 英文 12pt
    styles = styles.replace(
        normal2,
        normal2.replace(
            "</w:style>",
            "<w:rPr>" + tnr + f'<w:sz w:val="{bsz}"/><w:szCs w:val="{bsz}"/></w:rPr>'
            "</w:style>",
        ),
    )
    # 見出し: 英語版は本文と同サイズの太字(節)/太字イタリック(小節)。
    # 和文版は 明朝本文+ゴシック見出し の標準対(14/12/10.5pt 太字、イタリック不使用)
    if ja:
        hfont = ('<w:rFonts w:ascii="游ゴシック" w:hAnsi="游ゴシック" '
                 'w:eastAsia="游ゴシック" w:cs="游ゴシック"/>')
        hsizes = (("Heading1", "28"), ("Heading2", "24"), ("Heading3", "21"))
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
            lp2 = lp.replace('w:line="360" w:lineRule="atLeast"',
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
    sect = re.search(r"<w:sectPr[^>]*>.*?</w:sectPr>", doc, re.S).group(0)
    new_sect = sect
    if "pgSz" not in new_sect:
        # pandoc の sectPr は最小構成 — A4・余白 2.5cm を明示してロケール依存を排し、
        # 英語版のみ続けて lnNumType(スキーマ順: pgSz → pgMar → lnNumType)
        ln = "" if ja else '<w:lnNumType w:countBy="1" w:restart="continuous"/>'
        new_sect = new_sect.replace(
            "</w:sectPr>",
            '<w:pgSz w:w="11906" w:h="16838"/>'
            '<w:pgMar w:top="1417" w:right="1417" w:bottom="1417" w:left="1417" '
            'w:header="709" w:footer="709" w:gutter="0"/>'
            + ln + '</w:sectPr>',
        )
    footer_ref = '<w:footerReference xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" w:type="default" r:id="rIdFooterPg"/>'
    if "footerReference" not in new_sect:
        new_sect = re.sub(r"(<w:sectPr[^>]*>)", r"\1" + footer_ref, new_sect)
    parts["word/document.xml"] = doc.replace(sect, new_sect).encode("utf-8")

    parts["word/footer1.xml"] = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        f'<w:ftr xmlns:w="{W}"><w:p><w:pPr><w:jc w:val="center"/></w:pPr>'
        '<w:r><w:fldChar w:fldCharType="begin"/></w:r>'
        '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
        '<w:r><w:fldChar w:fldCharType="end"/></w:r></w:p></w:ftr>'
    ).encode("utf-8")

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
            theme = re.sub(r'(<a:ea typeface=")[^"]*(")', r"\1游明朝\2", theme)
        parts["word/theme/theme1.xml"] = theme.encode("utf-8")

    zout = zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED)
    for name, data in parts.items():
        zout.writestr(name, data)
    zout.close()


def tokens_md(text):
    text = re.sub(r"\^(\d[\d,]*)\^", r" \1 ", text)
    text = re.sub(r"[#*|`^]", " ", text)
    return re.findall(r"[A-Za-z0-9]+", text)


def tokens_docx(path):
    doc = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8")
    body = " ".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", doc))
    body = (body.replace("&lt;", "<").replace("&gt;", ">")
            .replace("&amp;", "&").replace("&quot;", '"'))
    return re.findall(r"[A-Za-z0-9]+", body)


def build(src_name, out_name, prep, ja=False):
    src_path = src_name if isinstance(src_name, Path) else SUB / src_name
    src = src_path.read_text(encoding="utf-8")
    pre = prep(src)
    md = SUB / out_name.replace(".docx", "_input.md")
    md.write_text(pre, encoding="utf-8")
    out = SUB / out_name
    subprocess.run(
        [PANDOC, str(md), "-f", "markdown+superscript", "-o", str(out)],
        check=True,
    )
    patch_docx(out, ja=ja)
    a, b = sorted(tokens_md(pre)), sorted(tokens_docx(out))
    from collections import Counter
    diff = Counter(a) - Counter(b) | Counter(b) - Counter(a)
    doc = zipfile.ZipFile(out).read("word/document.xml").decode("utf-8")
    styles = zipfile.ZipFile(out).read("word/styles.xml").decode("utf-8")
    checks = {"頁番号": "rIdFooterPg" in doc}
    if ja:
        checks["行送り18pt"] = 'w:line="360" w:lineRule="atLeast"' in styles
        checks["行番号なし"] = "lnNumType" not in doc
        checks["和文書体"] = "游明朝" in styles and "游ゴシック" in styles
    else:
        checks["1.5行間"] = 'w:line="360"' in styles
        checks["行番号"] = "lnNumType" in doc
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
    ]
    dest = Path("/mnt/c/Users/kotaro/OneDrive/論文関連（説明用資料含）/word_check")
    if dest.parent.exists():
        dest.mkdir(exist_ok=True)
        for out in outs:
            view = dest / f"{out.stem}_{tag}.docx"
            shutil.copy2(out, view)
            print(f"閲覧用コピー: {view}")


if __name__ == "__main__":
    main()
