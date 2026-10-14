#!/usr/bin/env python3
# curate_external_gene_lists.py
#
# Transcription of published multivariate gene signatures from the publishers'
# supplementary files into raw/external/*.csv (the 9-column schema read by
# scripts/460_external_gene_anchors.R). This is a curation record, not a
# pipeline step: the published pipeline consumes the committed CSVs only.
#
# Input : raw/external/source/  (git-ignored; publisher supplementary files
#         placed there by hand, file names as distributed)
#           41598_2026_53030_MOESM1_ESM.pdf   Ory 2026, Sci Rep 16:53030,
#             Supplementary Information (Tables S2, S4, S6); CC BY.
# Output: raw/external/ory2026_gene_signatures.csv  (--write), otherwise the
#         transcription is compared with the committed CSV and must be
#         identical (row order included); any difference is printed and the
#         exit status is 1.
#
# Symbol handling: symbol_published is copied verbatim from the table.
# symbol_current / ensembl_id / note apply the explicit ANNOTATIONS table
# below; everything else is carried unchanged, and resolution against the
# study's own GENCODE v36 annotation happens in scripts/460 (unresolved
# symbols are printed there, never dropped silently).
#
# Environment: Python 3 with pdfplumber (transcribed with pdfplumber 0.11.10).
import argparse, csv, re, sys, io
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "raw" / "external" / "source"
OUT = ROOT / "raw" / "external"
COLS = ["signature", "tissue", "mapping", "doi", "source_table",
        "symbol_published", "symbol_current", "ensembl_id", "note"]
# Curation decisions carried in the committed CSV (published symbol ->
# symbol_current, ensembl_id, note). Everything not listed is copied verbatim.
ANNOTATIONS = {
    "SPEN-AS1": ("SPEN-AS1", "", "not present in GENCODE v36 annotation"),
    "MCTS2":    ("MCTS2", "", "not present in GENCODE v36 annotation"),
    "SNORD47":  ("SNORD47", "", "not present in GENCODE v36 annotation"),
    "SMIM11":   ("SMIM11A", "", "GENCODE v36 symbol; source Entrez Gene 54065"),
    "HERC2P7":  ("HERC2P7", "ENSG00000281909", "explicit Ensembl mapping resolves duplicated GENCODE symbol"),
}

ORY2026 = dict(
    file="41598_2026_53030_MOESM1_ESM.pdf", doi="10.1038/s41598-026-53030-4",
    # (start marker, stop marker, signature id, tissue, mapping): mapping is
    # "shared" for the shared-tissue list and "direct" for tissue-specific lists.
    tables=[("Table S2:", "Table S3:", "ory2026_shared", "shared", "shared"),
            ("Table S4:", "Table S5:", "ory2026_normal", "normal", "direct"),
            ("Table S6:", "Table S7:", "ory2026_tumor", "tumor", "direct")])

ROW = re.compile(r"^(\S+)\s+(\d+)\s+[Cc]omponent\s+(\d)\s+(-?[\d.,]+)\s+(High/Low|Low/High)\s+([\d.,]+)\s+(HD|LD) vs NE$")
CONT = re.compile(r"^(-?[\d.,]+)\s+(High/Low|Low/High)\s+([\d.,]+)\s+(HD|LD) vs NE$")
SPLIT_A = re.compile(r"^(\S+-)\s+(\d+)\s+(\S+)\s+([Cc]omponent.*)$")            # "MORF4L2-" / "340544" / "AS1 component ..."
SPLIT_B = re.compile(r"^(\S+-)\s+([A-Za-z0-9]+)\s+(\d+)\s+([Cc]omponent.*)$")   # "DNAJC27- AS1 729723 Component ..."
SKIP = re.compile(r"^(=====|\d{1,2}$|Table S|Adjusted|Name ID|p-value|131$)")


def pdf_text(path):
    import pdfplumber
    with pdfplumber.open(str(path)) as pdf:
        return "\n".join(f"\n===== PAGE {i} =====\n" + (p.extract_text() or "")
                         for i, p in enumerate(pdf.pages, 1))


def parse_ory2026(text):
    out = []
    for start, stop, sig, tissue, mapping in ORY2026["tables"]:
        a = text.rindex(start)
        b = text.index(stop, a)
        lines = [l.strip() for l in text[a:b].split("\n") if l.strip() and not SKIP.match(l.strip())]
        buf, cur, n = "", None, 0
        for l in lines:
            cand = (buf + " " + l).strip() if buf else l
            m = SPLIT_A.match(cand)
            if m:
                cand = f"{m.group(1)}{m.group(3)} {m.group(2)} {m.group(4)}"
            m = SPLIT_B.match(cand)
            if m:
                cand = f"{m.group(1)}{m.group(2)} {m.group(3)} {m.group(4)}"
            m = ROW.match(cand)
            if m:
                sym = m.group(1)
                cur, n, buf = sym, n + 1, ""
                cur_sym, ensembl, note = ANNOTATIONS.get(sym, (sym, "", ""))
                out.append(dict(signature=sig, tissue=tissue, mapping=mapping, doi=ORY2026["doi"],
                                source_table=start.rstrip(":"), symbol_published=sym,
                                symbol_current=cur_sym, ensembl_id=ensembl, note=note))
                continue
            if CONT.match(cand) and cur is not None:
                buf = ""
                continue
            buf = cand
            if len(buf.split()) > 12:
                sys.exit(f"unparsed line in {sig}: {buf}")
        if buf:
            sys.exit(f"unparsed tail in {sig}: {buf}")
        print(f"{sig}: {n} genes transcribed from {start.rstrip(':')}")
    return out


def to_csv_text(rows):
    s = io.StringIO()
    w = csv.DictWriter(s, fieldnames=COLS, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return s.getvalue()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true", help="overwrite the committed CSV instead of comparing")
    args = ap.parse_args()
    rows = parse_ory2026(pdf_text(SRC / ORY2026["file"]))
    new = to_csv_text(rows)
    target = OUT / "ory2026_gene_signatures.csv"
    if args.write:
        target.write_text(new)
        print("wrote", target)
        return
    old = target.read_text()
    if old == new:
        print(f"identical to committed {target.name} ({len(rows)} rows)")
    else:
        import difflib
        sys.stdout.writelines(difflib.unified_diff(old.splitlines(True), new.splitlines(True), "committed", "transcribed"))
        sys.exit(1)


if __name__ == "__main__":
    main()
