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


def preprocess(text):
    """凡例節を References 後へ移動し、引用 [n] を上付き記法へ。"""
    # Markdown の区切り線(---)は原稿の可視要素ではない — 罫線化させない
    text = re.sub(r"^-{3,}\s*$", "", text, flags=re.M)
    m = re.search(r"^## Figure legends and table captions$.*?(?=^## )", text,
                  flags=re.S | re.M)
    legends = m.group(0)
    text = text.replace(legends, "")
    text = text.rstrip("\n") + "\n\n" + legends.rstrip("\n") + "\n"
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
    # PTC 単体は疾患略語のため対象外、既にイタリック化済みトークンは再包止め)
    genes = ("RET|BRAF|CCDC6|NCOA4|CLIP2|BHLHB9|S100A10|TESC|EHD4|"
             "ATP5MF|MRPL52|NTHL1|URM1|USE1|PTC1|PTC3")
    text = re.sub(rf"(?<![\w*])({genes})(?![\w*])", r"*\1*", text)
    return text


def patch_docx(path):
    """styles.xml(1.5 行間)・document.xml(行番号・フッタ参照)・
    footer1.xml(ページ番号)を注入する。"""
    zin = zipfile.ZipFile(path)
    parts = {n: zin.read(n) for n in zin.namelist()}
    zin.close()

    styles = parts["word/styles.xml"].decode("utf-8")
    normal = re.search(
        r'<w:style [^>]*w:styleId="Normal"[^>]*>.*?</w:style>', styles, re.S
    ).group(0)
    spacing = '<w:spacing w:after="160" w:line="360" w:lineRule="auto"/>'
    if "<w:pPr>" in normal:
        new_normal = re.sub(r"<w:spacing[^/]*/>", "", normal)
        new_normal = new_normal.replace("<w:pPr>", "<w:pPr>" + spacing, 1)
    else:
        # スキーマ上 pPr は style 末尾側 — 閉じタグ直前に挿入
        new_normal = normal.replace(
            "</w:style>", "<w:pPr>" + spacing + "</w:pPr></w:style>"
        )
    styles = styles.replace(normal, new_normal)
    # 原稿慣例の書体へ: 本文 Times New Roman 12pt、見出しは同書体・黒・太字。
    # テーマ経由の解決(日本語環境では Aptos に落ちる)を全経路で遮断するため、
    # eastAsia も含む明示指定を docDefaults・テーマにも適用する。
    tnr = ('<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" '
           'w:eastAsia="Times New Roman" w:cs="Times New Roman"/>')
    styles = re.sub(r"<w:rFonts [^/]*Theme[^/]*/>", tnr, styles)
    # 派生スタイル(BodyText 等)の spacing が Normal の行間を要素ごと上書きする —
    # w:line を持たない全 spacing に 1.5 行間を明示付与
    styles = re.sub(r'<w:spacing (?![^/>]*w:line=)([^/>]*)/>',
                    r'<w:spacing \1 w:line="360" w:lineRule="auto"/>', styles)
    normal2 = re.search(
        r'<w:style [^>]*w:styleId="Normal"[^>]*>.*?</w:style>', styles, re.S
    ).group(0)
    styles = styles.replace(
        normal2,
        normal2.replace(
            "</w:style>",
            "<w:rPr>" + tnr + '<w:sz w:val="24"/><w:szCs w:val="24"/></w:rPr>'
            "</w:style>",
        ),
    )
    # 見出しは本文と同サイズの太字(節)/太字イタリック(小節)— 投稿原稿の慣例形
    for hid, size in (("Heading1", "24"), ("Heading2", "24"), ("Heading3", "24")):
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
        ins = (tnr + '<w:b/>' + ('<w:i/>' if hid == "Heading3" else '')
               + '<w:color w:val="000000"/>'
               f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/>')
        if "<w:rPr>" in h2:
            h2 = h2.replace("<w:rPr>", "<w:rPr>" + ins, 1)
        else:
            h2 = h2.replace("</w:style>", "<w:rPr>" + ins + "</w:rPr></w:style>")
        styles = styles.replace(h, h2)
    parts["word/styles.xml"] = styles.encode("utf-8")

    doc = parts["word/document.xml"].decode("utf-8")
    sect = re.search(r"<w:sectPr[^>]*>.*?</w:sectPr>", doc, re.S).group(0)
    new_sect = sect
    if "lnNumType" not in new_sect:
        # pandoc の sectPr は最小構成 — A4・余白 2.5cm を明示してロケール依存を排し、
        # 続けて lnNumType(スキーマ順: pgSz → pgMar → lnNumType)
        new_sect = new_sect.replace(
            "</w:sectPr>",
            '<w:pgSz w:w="11906" w:h="16838"/>'
            '<w:pgMar w:top="1417" w:right="1417" w:bottom="1417" w:left="1417" '
            'w:header="709" w:footer="709" w:gutter="0"/>'
            '<w:lnNumType w:countBy="1" w:restart="continuous"/></w:sectPr>',
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

    # テーマの既定書体も Times New Roman へ(テーマ参照の取り残し対策)
    if "word/theme/theme1.xml" in parts:
        theme = parts["word/theme/theme1.xml"].decode("utf-8")
        theme = re.sub(r'(<a:latin typeface=")[^"]*(")',
                       r"\1Times New Roman\2", theme)
        parts["word/theme/theme1.xml"] = theme.encode("utf-8")

    zout = zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED)
    for name, data in parts.items():
        zout.writestr(name, data)
    zout.close()


def tokens_md(text):
    text = re.sub(r"\^(\d[\d,]*)\^", r"\1", text)
    text = re.sub(r"[#*|`^]", " ", text)
    return re.findall(r"[A-Za-z0-9]+", text)


def tokens_docx(path):
    doc = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8")
    body = " ".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", doc))
    body = (body.replace("&lt;", "<").replace("&gt;", ">")
            .replace("&amp;", "&").replace("&quot;", '"'))
    return re.findall(r"[A-Za-z0-9]+", body)


def main():
    src = (SUB / "manuscript_submission.md").read_text(encoding="utf-8")
    pre = preprocess(src)
    md = SUB / "manuscript_docx_input.md"
    md.write_text(pre, encoding="utf-8")

    out = SUB / "manuscript_submission.docx"
    ver = subprocess.run([PANDOC, "--version"], capture_output=True,
                         text=True).stdout.splitlines()[0]
    subprocess.run(
        [PANDOC, str(md), "-f", "markdown+superscript", "-o", str(out)],
        check=True,
    )
    patch_docx(out)

    a, b = sorted(tokens_md(pre)), sorted(tokens_docx(out))
    from collections import Counter
    diff = Counter(a) - Counter(b) | Counter(b) - Counter(a)
    doc = zipfile.ZipFile(out).read("word/document.xml").decode("utf-8")
    checks = {
        "lnNumType(行番号)": "lnNumType" in doc,
        "footerReference(頁番号)": "rIdFooterPg" in doc,
        "spacing360(1.5行間)": b'w:line="360"'
        in zipfile.ZipFile(out).read("word/styles.xml"),
        "上付き引用": "<w:vertAlign" in doc,
    }
    print(ver)
    print(f"出力: {out}")
    print(f"内容照合: token差 {sum(diff.values())} 件"
          + (f" {dict(list(diff.items())[:5])}" if diff else ""))
    for k, v in checks.items():
        print(f"  {k}: {'OK' if v else 'NG'}")
    if not all(checks.values()):
        sys.exit(1)

    dest = Path("/mnt/c/Users/kotaro/OneDrive/論文関連（説明用資料含）/word_check")
    if dest.parent.exists():
        dest.mkdir(exist_ok=True)
        # Word/OneDrive のキャッシュ・書き戻しを避けるため毎回新しいファイル名で渡す
        import time
        tag = time.strftime("%H%M")
        view = dest / f"manuscript_submission_{tag}.docx"
        shutil.copy2(out, view)
        print(f"閲覧用コピー: {view}")


if __name__ == "__main__":
    main()
