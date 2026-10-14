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
#           41598_2026_53030_MOESM1_ESM.pdf      Ory 2026, Sci Rep 16:53030,
#             Supplementary Information, Tables S2/S4/S6 (CC BY)
#           100205_erc100205_sm13_table_3.pdf    Ory 2011, Endocr Relat Cancer
#             18:193, Supplementary Table 3 (322-probe discriminating
#             signature; expression = log mean sporadic - log mean
#             radiation-induced, so a negative value is higher expression in
#             radiation-induced tumors)
#           pone.0023581.s006.doc                Ugolin 2011, PLoS One 6:e23581,
#             Table S3 on the journal page (the file is titled "Supplementary
#             table 2" inside): 106-gene / 109-probe signature discriminating
#             sporadic from post-Chernobyl PTCs (GEO GSE3950, Detours 2007);
#             expression = log mean post-Chernobyl - log mean sporadic, so a
#             positive value is higher expression in post-Chornobyl tumors.
#             Word 97 binary; text pieces are read as 8-bit and UTF-16 runs.
#         raw/reference/gencode.v36.annotation.gtf  (gene_name universe used
#             to choose symbol_current; final resolution happens in 460)
# Output: raw/external/ory2026_gene_signatures.csv
#         raw/external/ory2011_gene_signature.csv
#         raw/external/ugolin2011_gene_signature.csv
#         Without --write the transcription is compared with the committed CSV
#         and must be identical (row order included); differences are printed
#         and the exit status is 1.
#
# Symbol handling: symbol_published is copied verbatim. symbol_current is the
# GENCODE v36 gene_name under which the gene is carried: the published symbol
# when GENCODE v36 has it, otherwise the current HGNC/NCBI symbol listed in
# RENAMES when GENCODE v36 has that; otherwise the published symbol is kept
# and the note says so (460 prints such entries as unresolved, never drops
# them silently). Probes that hybridise two paralogues (two-line cells in the
# Ory 2011 table) give one row per gene sharing the probe's direction; probes
# without a gene symbol (ESTs, immunoglobulin loci, clones) are kept with the
# accession as symbol_published so that the listed size stays 322 probes.
#
# Environment: Python 3 with pdfplumber (transcribed with pdfplumber 0.11.10).
import argparse, csv, re, sys, io
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "raw" / "external" / "source"
OUT = ROOT / "raw" / "external"
GTF = ROOT / "raw" / "reference" / "gencode.v36.annotation.gtf"
COLS = ["signature", "tissue", "mapping", "doi", "source_table",
        "symbol_published", "symbol_current", "ensembl_id", "note"]

# ---------------------------------------------------------------- Ory 2026 --
ANNOTATIONS_2026 = {  # published symbol -> (symbol_current, ensembl_id, note)
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
ROW26 = re.compile(r"^(\S+)\s+(\d+)\s+[Cc]omponent\s+(\d)\s+(-?[\d.,]+)\s+(High/Low|Low/High)\s+([\d.,]+)\s+(HD|LD) vs NE$")
CONT26 = re.compile(r"^(-?[\d.,]+)\s+(High/Low|Low/High)\s+([\d.,]+)\s+(HD|LD) vs NE$")
SPLIT_A = re.compile(r"^(\S+-)\s+(\d+)\s+(\S+)\s+([Cc]omponent.*)$")            # "MORF4L2-" / "340544" / "AS1 component ..."
SPLIT_B = re.compile(r"^(\S+-)\s+([A-Za-z0-9]+)\s+(\d+)\s+([Cc]omponent.*)$")   # "DNAJC27- AS1 729723 Component ..."
SKIP26 = re.compile(r"^(=====|\d{1,2}$|Table S|Adjusted|Name ID|p-value|131$)")

# ---------------------------------------------------------------- Ory 2011 --
ORY2011 = dict(file="100205_erc100205_sm13_table_3.pdf", doi="10.1677/ERC-10-0205",
               signature="ory2011_tumor", tissue="tumor", mapping="direct",
               source_table="Supplementary Table 3", n_probes=322, n_higher=137, n_lower=185)
# published symbol -> (current symbol, provenance of the rename)
RENAMES_2011 = {
    "ATP5G1": ("ATP5MC1", "HGNC:841"), "ATPGD1": ("CARNS1", "HGNC:29268"), "C16orf75": ("RMI2", "HGNC:28349"),
    "C19orf20": ("TPGS1", "HGNC:25058"), "C9orf11": ("EQTN", "HGNC:1359"), "C9orf3": ("AOPEP", "HGNC:1361"),
    "CD97": ("ADGRE5", "HGNC:1711"), "DCI": ("ECI1", "HGNC:2703"), "FAIM3": ("FCMR", "HGNC:14315"),
    "FAM105B": ("OTULIN", "HGNC:25118"), "FAM108C1": ("ABHD17C", "HGNC:26925"), "GRAMD3": ("GRAMD2B", "HGNC:24911"),
    "GRASP": ("TAMALIN", "HGNC:18707"), "H2AFV": ("H2AZ2", "HGNC:20664"), "HIST1H2AB": ("H2AC4", "HGNC:4734"),
    "IMPAD1": ("BPNT2", "HGNC:26019"), "KIAA0355": ("GARRE1", "HGNC:29016"), "LRMP": ("IRAG2", "HGNC:6690"),
    "MINA": ("RIOX2", "HGNC:19441"), "MOBKL2A": ("MOB3A", "HGNC:29802"), "MRVI1": ("IRAG1", "HGNC:7237"),
    "NcRAN": ("SNHG16", "HGNC:44352 alias"), "ODZ2": ("TENM2", "HGNC:29943"), "RABL4": ("IFT27", "HGNC:18626"),
    "SC4MOL": ("MSMO1", "HGNC:10545"), "SFRS15": ("SCAF4", "HGNC:19304"),
    "LOC285758": ("MROCKI", "NCBI Gene 285758"), "LOC342346": ("C16orf96", "NCBI Gene 342346"),
    "LOC644246": ("KANSL1-AS1", "NCBI Gene 644246"), "LOC387790": ("GVQW3", "NCBI Gene 387790, now 100506127"),
}
# One probe lists its symbol inside a multi-word name cell ("NcRAN RNA for
# non-coding RNA expressed in aggressive neuroblastoma, long form"); the
# accession identifies it and the symbol is set explicitly.
SPECIAL_SYMBOLS_2011 = {"BC042949": "NcRAN"}
VAL11 = re.compile(r"^(-?\d+,\d+)\s+(\d+,\d+)(?:\s+\((.*)\))?$")
FULL11 = re.compile(r"^(\S+)\s+(.+?)\s+(-?\d+,\d+)\s+(\d+,\d+)(?:\s+\((.*)\))?$")
ACCVAL11 = re.compile(r"^(\S+)\s+(-?\d+,\d+)\s+(\d+,\d+)(?:\s+\((.*)\))?$")
ACCSYM11 = re.compile(r"^(\S+)\s+(\S+)$")
FOOT11 = re.compile(r"^(Downloaded|from|journals\..*|by|guest|on|\d{1,2}|October|2026|Accession|Symbol or name|number value|Expression)$")
SYMBOL11 = re.compile(r"^[A-Za-z][A-Za-z0-9\-\._]*$")


# --------------------------------------------------------------- Ugolin 2011 --
UGOLIN2011 = dict(file="pone.0023581.s006.doc", doi="10.1371/journal.pone.0023581",
                  signature="ugolin2011_tumor", tissue="tumor", mapping="direct",
                  source_table="Table S3 (file pone.0023581.s006.doc)", n_probes=109, n_genes=106)
RENAMES_UGOLIN = {
    "ATP5E": ("ATP5F1E", "HGNC:838"), "ATP5H": ("ATP5PD", "HGNC:845"), "CD97": ("ADGRE5", "HGNC:1711"),
    "PKM2": ("PKM", "HGNC:9021"), "RAD51L1": ("RAD51B", "HGNC:9822"), "SHFM1": ("SEM1", "HGNC:10845"),
    "TCEB1": ("ELOC", "HGNC:11617"), "TMEM178": ("TMEM178A", "HGNC:28517"), "RNU1A": ("RNU1-1", "HGNC:10120"),
    "LOC203274": ("LINC00537", "NCBI Gene 203274"),
}
# One row carries the NCBI locus identifier in the accession cell with the
# name "Hypothetical protein LOC203274" and no symbol; the identifier is used
# as the published symbol so that the NCBI rename above can apply.
SPECIAL_SYMBOLS_UGOLIN = {"LOC203274": "LOC203274"}
NUMU = re.compile(r"^-?\d+,\d+$")
ACCU = re.compile(r"^(N[MR]_\d+|X[MR]_\d+|[A-Z]{1,3}\d{5,6}|[A-Z]{2}\d{6})$")


def doc_text_runs(path):
    """Printable text runs of a Word 97 binary file (8-bit and UTF-16LE pieces), in file order."""
    b = open(path, "rb").read()
    runs = [(m.start(), m.group().decode("ascii")) for m in re.finditer(rb"[\x20-\x7e]{2,}", b)]
    runs += [(m.start(), m.group().decode("utf-16le")) for m in re.finditer(rb"(?:[\x20-\x7e\xa0-\xff]\x00){2,}", b)]
    return [t.strip() for _, t in sorted(runs) if t.strip()]


def parse_ugolin2011(names):
    toks = doc_text_runs(SRC / UGOLIN2011["file"])
    toks = toks[toks.index("Acc") + 1:]
    if "PAGE" in toks:                      # header/footer field after the last row
        toks = toks[:toks.index("PAGE")]
    probes, cur = [], None
    for t in toks:
        if ACCU.match(t):
            if cur:
                probes.append(cur)
            cur = dict(acc=t, text=[], nums=[])
        elif cur is not None:
            (cur["nums"] if NUMU.match(t) else cur["text"]).append(t)
    if cur:
        probes.append(cur)
    probes = [p for p in probes if len(p["nums"]) >= 2]
    if len(probes) != UGOLIN2011["n_probes"]:
        sys.exit(f"Ugolin 2011: {len(probes)} probes parsed, expected {UGOLIN2011['n_probes']}")
    # merged gene cells: a probe without text is a second probe of the preceding gene
    genes = []
    for p in probes:
        if not p["text"]:
            genes[-1]["probes"].append(p)
            continue
        genes.append(dict(text=p["text"], probes=[p]))
    if len(genes) != UGOLIN2011["n_genes"]:
        sys.exit(f"Ugolin 2011: {len(genes)} genes parsed, expected {UGOLIN2011['n_genes']}")
    out = []
    for g in genes:
        exprs = [float(p["nums"][0].replace(",", ".")) for p in g["probes"]]
        if len({e > 0 for e in exprs}) != 1:
            sys.exit(f"Ugolin 2011: probes of one gene disagree in sign: {g}")
        direction = "higher_in_post_chornobyl" if exprs[0] > 0 else "lower_in_post_chornobyl"
        accs = "/".join(p["acc"] for p in g["probes"])
        notes = [f"published_direction={direction}", f"accession {accs}"]
        if len(g["probes"]) > 1:
            notes.append(f"{len(g['probes'])} probes for this gene (merged cells in source)")
        sym = g["text"][0] if len(g["text"]) >= 2 and SYMBOL11.match(g["text"][0]) else ""
        sym = SPECIAL_SYMBOLS_UGOLIN.get(g["probes"][0]["acc"], sym)
        if not sym:
            key = g["probes"][0]["acc"]
            out.append(dict(signature=UGOLIN2011["signature"], tissue=UGOLIN2011["tissue"], mapping=UGOLIN2011["mapping"],
                            doi=UGOLIN2011["doi"], source_table=UGOLIN2011["source_table"], symbol_published=key,
                            symbol_current=key, ensembl_id="",
                            note="; ".join(notes[:1]) + f"; no gene symbol in source ({' '.join(g['text'])})"))
            continue
        cur = sym
        if sym not in names:
            if sym in RENAMES_UGOLIN and RENAMES_UGOLIN[sym][0] in names:
                cur = RENAMES_UGOLIN[sym][0]; notes.append(f"renamed ({RENAMES_UGOLIN[sym][1]})")
            elif sym in RENAMES_UGOLIN:
                notes.append(f"current symbol {RENAMES_UGOLIN[sym][0]} ({RENAMES_UGOLIN[sym][1]}); not present in GENCODE v36 annotation")
            else:
                notes.append("not present in GENCODE v36 annotation")
        out.append(dict(signature=UGOLIN2011["signature"], tissue=UGOLIN2011["tissue"], mapping=UGOLIN2011["mapping"],
                        doi=UGOLIN2011["doi"], source_table=UGOLIN2011["source_table"], symbol_published=sym,
                        symbol_current=cur, ensembl_id="", note="; ".join(notes)))
    keys = [o["symbol_published"] for o in out]
    dup = sorted({k for k in keys if keys.count(k) > 1})
    if dup:
        sys.exit(f"duplicate symbol_published in ugolin2011_tumor: {dup}")
    higher = sum(1 for o in out if "higher_in" in o["note"])
    resolved = sum(1 for o in out if o["symbol_current"] in names)
    print(f"ugolin2011_tumor: {len(probes)} probes, {len(out)} genes ({higher} higher / {len(out) - higher} lower in post-Chornobyl PTCs; {resolved} with a GENCODE v36 symbol)")
    return out


def pdf_pages(path):
    import pdfplumber
    with pdfplumber.open(str(path)) as pdf:
        return [(p.extract_text() or "") for p in pdf.pages]


def gencode_names():
    names = set()
    with open(GTF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 9)
            if len(parts) > 8 and parts[2] == "gene":
                m = re.search(r'gene_name "([^"]+)"', parts[8])
                if m:
                    names.add(m.group(1))
    return names


def parse_ory2026(pages):
    text = "\n".join(f"\n===== PAGE {i} =====\n{t}" for i, t in enumerate(pages, 1))
    out = []
    for start, stop, sig, tissue, mapping in ORY2026["tables"]:
        a = text.rindex(start)
        b = text.index(stop, a)
        lines = [l.strip() for l in text[a:b].split("\n") if l.strip() and not SKIP26.match(l.strip())]
        buf, cur, n = "", None, 0
        for l in lines:
            cand = (buf + " " + l).strip() if buf else l
            m = SPLIT_A.match(cand)
            if m:
                cand = f"{m.group(1)}{m.group(3)} {m.group(2)} {m.group(4)}"
            m = SPLIT_B.match(cand)
            if m:
                cand = f"{m.group(1)}{m.group(2)} {m.group(3)} {m.group(4)}"
            m = ROW26.match(cand)
            if m:
                sym = m.group(1)
                cur, n, buf = sym, n + 1, ""
                cur_sym, ensembl, note = ANNOTATIONS_2026.get(sym, (sym, "", ""))
                out.append(dict(signature=sig, tissue=tissue, mapping=mapping, doi=ORY2026["doi"],
                                source_table=start.rstrip(":"), symbol_published=sym,
                                symbol_current=cur_sym, ensembl_id=ensembl, note=note))
                continue
            if CONT26.match(cand) and cur is not None:
                buf = ""
                continue
            buf = cand
            if len(buf.split()) > 12:
                sys.exit(f"unparsed line in {sig}: {buf}")
        if buf:
            sys.exit(f"unparsed tail in {sig}: {buf}")
        print(f"{sig}: {n} genes transcribed from {start.rstrip(':')}")
    return out


def parse_ory2011(pages, names):
    """One record per probe: dict(expr, p, pieces=[(token, rest), ...])."""
    probes = []
    for pno, text in enumerate(pages[:7], 1):
        lines = [l.strip() for l in text.split("\n") if l.strip() and not FOOT11.match(l.strip())]
        if pno == 1:
            lines = lines[next(i for i, l in enumerate(lines) if FULL11.match(l)):]
        if "References" in lines:
            lines = lines[:lines.index("References")]
        pend, i = [], 0
        while i < len(lines):
            l = lines[i]
            nxt = lines[i + 1] if i + 1 < len(lines) else ""
            nxt2 = lines[i + 2] if i + 2 < len(lines) else ""
            m = FULL11.match(l)
            if m:
                acc, name, e, p, refs = m.groups()
                pieces = pend + [(acc, name)]
                pend = []
                # a lower-case fragment on the next line continues this name (3-line cell)
                if nxt and nxt[0].islower() and not FULL11.match(nxt) and not VAL11.match(nxt):
                    pieces.append(("", nxt)); i += 1
                probes.append(dict(page=pno, expr=e, p=p, refs=refs or "", pieces=pieces)); i += 1; continue
            m = ACCVAL11.match(l)
            if m:
                acc, e, p, refs = m.groups()
                probes.append(dict(page=pno, expr=e, p=p, refs=refs or "", pieces=pend + [(acc, "")])); pend = []; i += 1; continue
            m = VAL11.match(l)
            if m:
                e, p, refs = m.groups()
                pieces = list(pend); pend = []
                # bottom half of a two-line paralogue cell: "ACC SYM" not followed by a values line
                if ACCSYM11.match(nxt) and not FULL11.match(nxt) and not VAL11.match(nxt2):
                    pieces.append(ACCSYM11.match(nxt).groups()); i += 1
                probes.append(dict(page=pno, expr=e, p=p, refs=refs or "", pieces=pieces)); i += 1; continue
            pend.append(tuple(l.split(" ", 1)) if " " in l else (l, ""))
            i += 1
        if pend:
            sys.exit(f"unparsed tail on page {pno}: {pend}")
    n = len(probes)
    higher = sum(1 for r in probes if float(r["expr"].replace(",", ".")) < 0)
    print(f"ory2011_tumor: {n} probes transcribed ({higher} higher / {n - higher} lower in radiation-induced tumours)")
    if (n, higher, n - higher) != (ORY2011["n_probes"], ORY2011["n_higher"], ORY2011["n_lower"]):
        sys.exit("Ory 2011 probe count or direction split does not match the paper (322; 137/185)")
    out = []
    for r in probes:
        direction = "higher_in_radiation_induced" if float(r["expr"].replace(",", ".")) < 0 else "lower_in_radiation_induced"
        accs = [t for t, rest in r["pieces"] if rest == "" and t] + [t for t, rest in r["pieces"] if rest and SYMBOL11.match(rest) and t]
        genes = [rest for t, rest in r["pieces"] if rest and SYMBOL11.match(rest)]
        for t, rest in r["pieces"]:
            if t in SPECIAL_SYMBOLS_2011 and SPECIAL_SYMBOLS_2011[t] not in genes:
                genes.append(SPECIAL_SYMBOLS_2011[t]); accs.append(t)
        desc = " ".join(x for t, rest in r["pieces"] for x in ([t, rest] if (rest and not SYMBOL11.match(rest)) else []) if x).strip()
        acc_note = "accession " + "/".join(dict.fromkeys(accs)) if accs else ""
        if not genes:
            key = accs[0] if accs else desc.split()[0]
            if key == "Unknown":  # three probes mapped only to a chromosome in the source
                key = "Unknown_" + desc.split()[-1]
            out.append(dict(signature=ORY2011["signature"], tissue=ORY2011["tissue"], mapping=ORY2011["mapping"],
                            doi=ORY2011["doi"], source_table=ORY2011["source_table"], symbol_published=key,
                            symbol_current=key, ensembl_id="",
                            note=f"published_direction={direction}; no gene symbol in source ({desc or 'unnamed probe'})"))
            continue
        for g in genes:
            notes = [f"published_direction={direction}", acc_note] if acc_note else [f"published_direction={direction}"]
            if len(genes) > 1:
                notes.append("probe hybridises " + " and ".join(genes) + " (two-line cell in source)")
            cur = g
            if g not in names:
                if g in RENAMES_2011 and RENAMES_2011[g][0] in names:
                    cur = RENAMES_2011[g][0]; notes.append(f"renamed ({RENAMES_2011[g][1]})")
                elif g in RENAMES_2011:
                    notes.append(f"current symbol {RENAMES_2011[g][0]} ({RENAMES_2011[g][1]}); not present in GENCODE v36 annotation")
                else:
                    notes.append("not present in GENCODE v36 annotation")
            out.append(dict(signature=ORY2011["signature"], tissue=ORY2011["tissue"], mapping=ORY2011["mapping"],
                            doi=ORY2011["doi"], source_table=ORY2011["source_table"], symbol_published=g,
                            symbol_current=cur, ensembl_id="", note="; ".join(notes)))
    keys = [o["symbol_published"] for o in out]
    dup = sorted({k for k in keys if keys.count(k) > 1})
    if dup:
        sys.exit(f"duplicate symbol_published in ory2011_tumor: {dup}")
    resolved = sum(1 for o in out if o["symbol_current"] in names)
    print(f"ory2011_tumor: {len(out)} rows ({resolved} with a GENCODE v36 symbol)")
    return out


def to_csv_text(rows):
    s = io.StringIO()
    w = csv.DictWriter(s, fieldnames=COLS, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return s.getvalue()


def emit(rows, target, write):
    new = to_csv_text(rows)
    if write:
        target.write_text(new)
        print("wrote", target)
        return True
    if not target.exists():
        print(f"{target.name} does not exist yet; run with --write {target.stem}")
        return False
    old = target.read_text()
    if old == new:
        print(f"identical to committed {target.name} ({len(rows)} rows)")
        return True
    import difflib
    sys.stdout.writelines(difflib.unified_diff(old.splitlines(True), new.splitlines(True), "committed", "transcribed"))
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", nargs="*", default=None, metavar="LIST",
                    help="write these CSVs (ory2026_gene_signatures, ory2011_gene_signature, ugolin2011_gene_signature) instead of comparing")
    args = ap.parse_args()
    write = set(args.write) if args.write is not None else set()
    names = gencode_names()
    ok = True
    ok &= emit(parse_ory2026(pdf_pages(SRC / ORY2026["file"])), OUT / "ory2026_gene_signatures.csv",
               "ory2026_gene_signatures" in write)
    ok &= emit(parse_ory2011(pdf_pages(SRC / ORY2011["file"]), names), OUT / "ory2011_gene_signature.csv",
               "ory2011_gene_signature" in write)
    ok &= emit(parse_ugolin2011(names), OUT / "ugolin2011_gene_signature.csv",
               "ugolin2011_gene_signature" in write)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
