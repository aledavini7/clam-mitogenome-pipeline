#!/usr/bin/env python3
"""Build a concise, prioritization-oriented HTML summary for nuclear WES."""

import argparse
import csv
import html
import math
from pathlib import Path
from urllib.parse import unquote


# This watchlist is a display-prioritization aid only, not a clinical or
# variant-pathogenicity classification.
REVIEW_GENE_WATCHLIST = {
    "TP53", "KRAS", "NRAS", "HRAS", "BRAF", "EGFR", "ERBB2", "PIK3CA",
    "PTEN", "APC", "RB1", "BRCA1", "BRCA2", "NF1", "NF2", "CDKN2A",
    "STK11", "MET", "ALK", "RET", "ROS1", "KIT", "PDGFRA", "IDH1",
    "IDH2", "FGFR1", "FGFR2", "FGFR3", "POLE", "MSH2", "MSH6", "MLH1",
    "PMS2", "NPM1", "FLT3", "JAK2", "NOTCH1", "TSC1", "TSC2", "SMAD4",
    "CTNNB1", "MYC", "CDH1", "VHL", "KEAP1", "ARID1A", "ARID2", "KMT2D",
    "KMT2C", "SETD2", "GATA3", "MAP2K1", "MAP2K2", "RHOA", "RUNX1",
    "SF3B1", "SRSF2", "U2AF1", "DNMT3A", "TET2", "ASXL1", "CUX1",
    "ABL1", "ABRAXAS1", "ACVR1", "AKT1", "AKT2", "AKT3", "AMER1", "AR",
    "ARAF", "ARHGAP35", "ARID2", "ASXL2", "ATR", "ATRX", "B2M", "BAP1",
    "BCL10", "BCL11A", "BCL2", "BCL6", "BCL7A", "BCL9", "BIRC3", "BLM",
    "BMPR1A", "BRAF", "BRIP1", "CARD11", "CASP8", "CBL", "CCND1", "CCND2",
    "CCND3", "CCNE1", "CD274", "CD28", "CD58", "CDC73", "CDH1", "CDK4",
    "CDK6", "CDK12", "CDKN1B", "CDKN2B", "CEBPA", "CHD2", "CHD4", "CHEK1",
    "CHEK2", "CIC", "CREBBP", "CSF1R", "CTNNA1", "CUL3", "DAXX", "DDX41",
    "DDX3X", "DICER1", "DIS3", "DNMT3B", "EED", "EP300", "ERBB3", "ERBB4",
    "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ESR1", "ETV6", "FAM46C", "FANCA",
    "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FAS", "FBXW7", "FGFR4",
    "FH", "FLCN", "FOXL2", "GATA1", "GATA2", "GNA11", "GNAQ", "GNAS",
    "H3-3A", "H3-3B", "HIST1H3B", "HNF1A", "HOXB13", "IKZF1", "IL7R", "JAK1",
    "JAK3", "KDM5C", "KDM6A", "KDR", "KLF4", "KMT2A", "KMT2B", "KMT2D",
    "LZTR1", "MALT1", "MAP2K1", "MAP2K2", "MAPK1", "MDM2", "MDM4", "MED12",
    "MEN1", "MGA", "MITF", "MLH1", "MTOR", "MYCN", "MYCL", "NBN", "NFE2L2",
    "NOTCH2", "NOTCH3", "NSD1", "NSD2", "PALB2", "PAX5", "PBRM1", "PDGFRB",
    "PIK3CB", "PIK3R1", "POLD1", "POT1", "PPP2R1A", "PRDM1", "PRKAR1A",
    "PTCH1", "PTPN11", "RAD21", "RAD51B", "RAD51C", "RAD51D", "RARA", "RBM10",
    "REL", "RHOA", "RIT1", "RNF43", "RUNX2", "RUNX3", "SETBP1", "SH2B3",
    "SMAD2", "SMAD3", "SMARCA4", "SMARCB1", "SMARCE1", "SMC1A", "SMC3",
    "SOCS1", "SOCS3", "SOX2", "SPOP", "SRC", "STAG2", "STAT3", "STAT5B",
    "SUFU", "TERT", "TGFBR2", "TNFAIP3", "TNFRSF14", "TP63", "TSC1", "TSC2",
    "VAV1", "WHSC1", "WT1", "ZRSR2", "ZBTB7A", "ZMYM3", "ZNF217",
}
# Broaden review priority beyond oncology to include clinically relevant immune,
# hematologic, and constitutional-disease genes that may matter in WES samples.
REVIEW_GENE_WATCHLIST.update({
    "ADA", "AIRE", "ATM", "BCL11B", "BLNK", "BTK", "CARMIL2", "CASP10",
    "CD40", "CD40LG", "CD79A", "CD79B", "CTLA4", "CYBB", "DOCK8", "ELANE",
    "FASLG", "FOXP3", "G6PC3", "GATA2", "HAX1", "ICOS", "IFNAR1", "IFNAR2",
    "IFNGR1", "IFNGR2", "IL10RA", "IL10RB", "IL12B", "IL12RB1", "IL17RA",
    "IL17RC", "IL2RA", "IL2RG", "IL21R", "IL23R", "IL6ST", "IRAK4", "IRF3",
    "IRF7", "IRF9", "JAGN1", "LAT", "LIG4", "LIPA", "LYST", "MALT1", "MCM4",
    "MPO", "MYD88", "NBN", "NFKB1", "NFKB2", "NFKBIA", "ORAI1", "PGM3",
    "PIK3CD", "PRF1", "RAC2", "RAG1", "RAG2", "RELA", "RFX5", "RFXANK",
    "RFXAP", "SAMD9", "SAMD9L", "SH2D1A", "SLC29A3", "SPPL2A", "STAT1",
    "STIM1", "STK4", "STX11", "STXBP2", "TBX1", "TNFRSF13B", "TNFRSF13C",
    "TNFRSF14", "TRAC", "TRAF3IP2", "TYK2", "UNC13D", "VPS13B", "WAS", "ZAP70",
    "ABCA4", "APOE", "ATP7B", "COL1A1", "COL1A2", "DMD", "FBN1", "GBA1",
    "HFE", "LDLR", "MYH7", "PAH", "PKD1", "PKD2", "RYR1", "SCN1A", "TTR",
})
CORE_DRIVER_WATCHLIST = {
    "TP53", "KRAS", "NRAS", "BRAF", "EGFR", "ERBB2", "PIK3CA", "PTEN",
    "APC", "RB1", "BRCA1", "BRCA2", "NF1", "CDKN2A", "STK11", "ALK",
    "RET", "MET", "IDH1", "IDH2", "FLT3", "NPM1", "JAK2",
}

PROTEIN_ALTERING = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
    "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site",
    "Nonstop_Mutation", "Translation_Start_Site",
}
COMMON_AF_THRESHOLD = 0.01
SEVERITY = {
    "Frame_Shift_Del": 0,
    "Frame_Shift_Ins": 0,
    "Nonsense_Mutation": 1,
    "Nonstop_Mutation": 1,
    "Translation_Start_Site": 1,
    "Splice_Site": 1,
    "In_Frame_Del": 2,
    "In_Frame_Ins": 2,
    "Missense_Mutation": 3,
}


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def integer(value):
    value = number(value)
    return int(value) if value is not None else None


def escape(value):
    return html.escape(str(value if value not in (None, "") else "—"))


def fmt_depth(value):
    return f"{value:,.1f}×" if value is not None else "Not available"


def variant_key(row):
    return tuple(row.get(key, "") for key in (
        "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"
    ))


def allele_fraction(row):
    depth = number(row.get("t_depth"))
    alt = number(row.get("t_alt_count"))
    if depth and alt is not None and depth > 0:
        return alt / depth
    return None


def compact_allele(value, max_bases=5):
    value = str(value)
    return value if len(value) <= max_bases else value[:max_bases] + "..."


def unique_variants(rows):
    seen = set()
    result = []
    for row in rows:
        key = variant_key(row)
        if key not in seen:
            seen.add(key)
            result.append(row)
    return result


def rank_key(row, common_af_threshold=COMMON_AF_THRESHOLD):
    gene = row.get("Hugo_Symbol", "").upper()
    consequence = row.get("Variant_Classification", "")
    af = number(row.get("gnomAD_AF"))
    vaf = allele_fraction(row)
    # Watchlist gene first, then predicted consequence severity and database
    # evidence. Missing AF is deliberately not treated as rare.
    return (
        0 if gene == "TP53" else 1 if gene in CORE_DRIVER_WATCHLIST else 2 if gene in REVIEW_GENE_WATCHLIST else 3,
        SEVERITY.get(consequence, 9),
        1 if af is not None and af >= common_af_threshold else 0,
        0 if row.get("COSMIC_ID", "").strip() not in ("", ".") else 1,
        af if af is not None else float("inf"),
        -(vaf if vaf is not None else -1),
        gene,
        integer(row.get("Start_Position")) or 0,
    )


def population_flag(value, threshold):
    af = number(value)
    if af is None:
        return "Not annotated"
    if af >= threshold:
        return f"Common (≥{threshold:g})"
    return f"Below {threshold:g}"


def choose_primary_segments(cnv_dir):
    paths = sorted(cnv_dir.glob("*.seg"))
    primary = [p for p in paths if ".bintest." not in p.name and ".call.call." not in p.name]
    return primary or paths


def read_segments(paths):
    segments = []
    for path in paths:
        try:
            for row in read_tsv(path):
                chrom = row.get("chrom", row.get("chromosome", ""))
                start = integer(row.get("loc.start", row.get("start")))
                end = integer(row.get("loc.end", row.get("end")))
                mean = number(row.get("seg.mean", row.get("log2")))
                marks = integer(row.get("num.mark", row.get("probes")))
                if start is None or end is None or mean is None:
                    continue
                segments.append({
                    "source": path.name,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "length": max(0, end - start + 1),
                    "mean": mean,
                    "marks": marks,
                })
        except (OSError, csv.Error):
            continue
    return segments


def target_coverage(cnv_dir):
    paths = sorted(cnv_dir.glob("*.targetcoverage.cnn"))
    if not paths:
        return None
    depths = []
    for path in paths:
        try:
            for row in read_tsv(path):
                depth = number(row.get("depth"))
                start = integer(row.get("start"))
                end = integer(row.get("end"))
                if depth is not None and start is not None and end is not None and end > start:
                    depths.append((depth, end - start))
        except (OSError, csv.Error):
            continue
    if not depths:
        return None
    total_bases = sum(length for _, length in depths)
    return {
        "mean": sum(depth * length for depth, length in depths) / total_bases,
        "min": min(depth for depth, _ in depths),
        "max": max(depth for depth, _ in depths),
        "bins": len(depths),
    }


def table(headers, rows, empty_message, classes=None):
    head = "".join(f"<th>{escape(value)}</th>" for value in headers)
    body = []
    for row in rows:
        cells = "".join(f"<td>{value}</td>" for value in row)
        body.append(f'<tr class="{classes(row) if classes else ""}">{cells}</tr>')
    if not body:
        body = [f'<tr><td colspan="{len(headers)}" class="empty">{escape(empty_message)}</td></tr>']
    return f"<div class=\"table-wrap\"><table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--calling-mode", choices=("tumor_only", "paired", "germline"), required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--maf", required=True)
    parser.add_argument("--cnv-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--common-af-threshold", type=float, default=COMMON_AF_THRESHOLD,
                        help="Flag population AF at/above this value as common (default: 0.01).")
    parser.add_argument("--top-variants", type=int, default=100,
                        help="Legacy compatibility option; all protein-altering variants are shown.")
    args = parser.parse_args()
    if not 0 <= args.common_af_threshold <= 1:
        parser.error("--common-af-threshold must be between 0 and 1")
    if args.top_variants < 1:
        parser.error("--top-variants must be at least 1")

    sample_id = html.escape(args.sample_id)
    maf_rows = read_tsv(args.maf)
    variants = unique_variants(maf_rows)
    coding = [row for row in variants if row.get("Variant_Classification", "") in PROTEIN_ALTERING]
    coding.sort(key=lambda row: rank_key(row, args.common_af_threshold))
    high_impact = sum(row.get("Variant_Classification") in {
        "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation",
        "Translation_Start_Site", "Splice_Site"
    } for row in variants)
    missense = sum(row.get("Variant_Classification") == "Missense_Mutation" for row in variants)
    gnomad_count = sum(number(row.get("gnomAD_AF")) is not None for row in variants)
    common_count = sum(
        (number(row.get("gnomAD_AF")) or 0) >= args.common_af_threshold
        for row in variants
    )

    cnv_dir = Path(args.cnv_dir)
    seg_paths = choose_primary_segments(cnv_dir)
    segments = read_segments(seg_paths)
    extreme_segments = [segment for segment in segments if abs(segment["mean"]) > 5]
    candidate_segments = [
        segment for segment in segments
        if abs(segment["mean"]) >= 0.3
        and abs(segment["mean"]) <= 5
        and segment["length"] >= 5_000_000
        and (segment["marks"] is None or segment["marks"] >= 10)
    ]
    candidate_segments.sort(key=lambda segment: (-segment["length"], -abs(segment["mean"])))
    reference_mode = "paired/reference CNVkit baseline" if list(cnv_dir.glob("paired_reference.cnn")) else "flat CNVkit baseline (single-sample exploratory)"
    coverage = target_coverage(cnv_dir)

    variant_table_rows = []
    for row_index, row in enumerate(coding):
        gene = row.get("Hugo_Symbol", "")
        gene_upper = gene.upper()
        priority = (
            "TP53 focus" if gene_upper == "TP53"
            else "Core driver" if gene_upper in CORE_DRIVER_WATCHLIST
            else "Review-priority gene" if gene_upper in REVIEW_GENE_WATCHLIST
            else "—"
        )
        vaf = allele_fraction(row)
        vaf_text = f"{vaf * 100:.2f}%" if vaf is not None else "—"
        counts = "—"
        if row.get("t_alt_count") or row.get("t_depth"):
            counts = f"{row.get('t_alt_count', '—')}/{row.get('t_depth', '—')}"
        full_alt = row.get("Tumor_Seq_Allele2", "")
        loc = f"{row.get('Chromosome', '')}:{row.get('Start_Position', '')} {row.get('Reference_Allele', '')}>{compact_allele(full_alt)}"
        full_loc = f"{row.get('Chromosome', '')}:{row.get('Start_Position', '')} {row.get('Reference_Allele', '')}>{full_alt}"
        variant_table_rows.append([
            f"<strong>{escape(gene)}</strong>", escape(priority),
            f'<span title="{escape(full_loc)}">{escape(loc)}</span>',
            escape(row.get("Variant_Classification", "")), escape(unquote(row.get("HGVSc", ""))),
            escape(unquote(row.get("HGVSp", ""))),
            escape(vaf_text), escape(counts), escape(row.get("gnomAD_AF", "")),
            escape(population_flag(row.get("gnomAD_AF", ""), args.common_af_threshold)),
            escape(row.get("COSMIC_ID", "")),
        ])

    cnv_table_rows = []
    for segment in candidate_segments[:15]:
        direction = "Gain-like" if segment["mean"] > 0 else "Loss-like"
        cnv_table_rows.append([
            escape(segment["chrom"]), f"{segment['start']:,}–{segment['end']:,}",
            f"{segment['length'] / 1_000_000:.1f} Mb", f"{segment['mean']:+.2f}", direction,
            str(segment["marks"] if segment["marks"] is not None else "—"),
        ])

    summary_cards = [
        ("PASS variant records", f"{len(variants):,}", "Unique chromosome-position-allele records in the MAF"),
        ("Protein-altering candidates", f"{len(coding):,}", f"{missense:,} missense; {high_impact:,} truncating/splice/start/stop"),
        ("gnomAD frequency", f"{gnomad_count:,}/{len(variants):,} annotated", f"{common_count:,} at or above AF {args.common_af_threshold:g}; missing AF is unknown"),
        ("Non-synonymous calls", f"{len(coding):,}", "Protein-altering variants in the MAF"),
        ("Mean on-target depth", fmt_depth(coverage["mean"]) if coverage else "Not available", f"{coverage['bins']:,} target bins" if coverage else "Target coverage file absent from supplied outputs"),
        ("Target-bin depth range", f"{fmt_depth(coverage['min'])} – {fmt_depth(coverage['max'])}" if coverage else "Not available", "Minimum to maximum target-bin mean depth"),
    ]
    cards = "".join(
        f'<div class="card"><div class="label">{escape(label)}</div><div class="value">{escape(value)}</div><div class="hint">{escape(hint)}</div></div>'
        for label, value, hint in summary_cards
    )
    variant_table_html = table(
        ["Gene", "Priority", "Variant", "Consequence", "HGVS cDNA", "HGVS protein", "VAF", "Alt/Depth", "gnomAD AF", "Population frequency", "COSMIC"],
        variant_table_rows,
        "No protein-altering variants found in the MAF.",
        classes=lambda row: "watchlist" if any(
            marker in row for marker in ("TP53 focus", "Core driver", "Review-priority gene")
        ) else "",
    )
    # Keep the full protein-altering callset visible and scrollable; filters only
    # narrow rows when the user enters a search or selects a consequence.
    variant_rows_html = []
    for idx, cells in enumerate(variant_table_rows):
        rendered = "".join(f"<td>{cell}</td>" for cell in cells)
        row_class = "watchlist" if cells[1] in ("TP53 focus", "Core driver", "Review-priority gene") else ""
        source_row = coding[idx]
        full_alt = source_row.get("Tumor_Seq_Allele2", "")
        searchable = " ".join(str(cell).replace("<strong>", "").replace("</strong>", "") for cell in cells) + " " + full_alt
        variant_rows_html.append(f'<tr class="{row_class.strip()}" data-search="{escape(searchable).lower()}">{rendered}</tr>')
    variant_table_html = f'''<div class="table-controls">
<label>Search all protein-altering variants <input id="variant-search" type="search" placeholder="Gene, HGVS, consequence, coordinate"></label>
<label>Consequence <select id="consequence-filter"><option value="">All consequences</option>{"".join(f"<option>{escape(consequence)}</option>" for consequence in sorted({row.get("Variant_Classification", "") for row in coding}))}</select></label>
<span id="variant-count">Showing all {len(coding):,} protein-altering variants</span>
</div><div class="table-wrap"><table id="nuclear-variants"><thead><tr>{"".join(f"<th>{escape(header)}</th>" for header in ["Gene", "Priority", "Variant", "Consequence", "HGVS cDNA", "HGVS protein", "VAF", "Alt/Depth", "gnomAD AF", "Population frequency", "COSMIC"])}</tr></thead><tbody>{"".join(variant_rows_html) if variant_rows_html else f'<tr><td colspan="11" class="empty">No protein-altering variants found in the MAF.</td></tr>'}</tbody></table></div>'''
    cnv_table_html = table(
        ["Chromosome", "Interval", "Length", "Segment log2", "Direction", "Bins"],
        cnv_table_rows,
        "No broad candidate segments passed the report thresholds. This does not exclude focal CNVs.",
    )
    cnv_status = (
        f"{len(segments):,} primary segments from {len(seg_paths)} segment file(s); "
        f"{len(candidate_segments):,} broad candidate(s) shown; {len(extreme_segments):,} extreme outlier(s) flagged."
        if segments else "No CNVkit segment files were found."
    )
    document = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>CLAM nDNA summary: {sample_id}</title>
<style>
:root{{color-scheme:light;--ink:#16202a;--muted:#667085;--line:#d9e2ec;--blue:#9a3412;--panel:#fff;--bg:#f5f7fb;--amber:#fff7df;--red:#fff0ee}}
*{{box-sizing:border-box}} body{{margin:0;font:15px/1.45 Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:var(--ink);background:var(--bg);letter-spacing:0}}
header{{padding:30px 36px 22px;background:#7c2d12;color:#fff}} header h1{{margin:0;font-size:30px;font-weight:760;letter-spacing:0}} header p{{margin:8px 0 0;color:#ffedd5}}
main{{max-width:1220px;margin:0 auto;padding:24px 24px 44px}} h2{{margin:30px 0 8px;font-size:20px}} p{{margin:5px 0;color:var(--muted)}} .meta{{margin-bottom:20px}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:20px 0}}
.card{{padding:15px;border:1px solid var(--line);border-radius:10px;background:white;min-height:112px}}
.label{{font-size:12px;text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}} .value{{font-size:24px;font-weight:700;margin:5px 0}} .hint{{font-size:12px;color:var(--muted)}}
.notice{{border-left:4px solid #d97706;background:var(--amber);padding:12px 16px;margin:18px 0;border-radius:4px}} .notice ul{{margin:4px 0 0;padding-left:20px}} .notice li{{margin:4px 0}}
.table-wrap{{overflow:auto;max-height:70vh;border:1px solid var(--line);border-radius:8px;margin:10px 0 20px}} table{{border-collapse:collapse;width:100%;font-size:13px;white-space:nowrap}} th,td{{border-bottom:1px solid var(--line);padding:8px 10px;text-align:left}} #nuclear-variants td:first-child{{white-space:nowrap}} th{{position:sticky;top:0;background:#edf2f7;z-index:1}} tr:last-child td{{border-bottom:0}} #nuclear-variants tbody tr:not(.watchlist){{background:#f8fafc}} .watchlist{{background:#eff6ff}} .empty{{padding:16px;color:var(--muted)}}
.table-controls{{display:flex;gap:14px;align-items:end;flex-wrap:wrap;margin:12px 0}} .table-controls label{{display:flex;flex-direction:column;gap:4px;font-size:12px;color:var(--muted)}} .table-controls input,.table-controls select{{padding:8px;border:1px solid var(--line);border-radius:6px;min-width:220px}} .table-controls span{{font-size:12px;color:var(--muted)}}
.small{{font-size:12px;color:var(--muted)}} .badge{{display:inline-block;border-radius:99px;background:#dbeafe;color:#1e40af;padding:2px 8px;font-size:11px;font-weight:650}}
@media print{{body{{max-width:none;padding:0}} .table-wrap{{overflow:visible}} table{{white-space:normal}} th{{position:static}}}}
</style></head><body>
<header><h1>CLAM Nuclear DNA Report</h1><p>Whole-exome sequencing (WES) &nbsp;|&nbsp; Sample {sample_id}</p></header>
<main>
<div class="meta"><p><b>Sample:</b> {sample_id} &nbsp; <b>Calling mode:</b> {escape(args.calling_mode)} &nbsp; <b>CNV baseline:</b> {escape(reference_mode)}</p>
 </div>
<div class="cards">{cards}</div>
<h2>Protein-altering candidates to review</h2>
{variant_table_html}
<script>
(() => {{
 const table = document.getElementById('nuclear-variants'); if (!table) return;
 const rows = [...table.querySelectorAll('tbody tr')];
 const search = document.getElementById('variant-search');
 const consequence = document.getElementById('consequence-filter');
 const count = document.getElementById('variant-count');
 function update() {{
   const query = search.value.trim().toLowerCase(), effect = consequence.value;
   let shown = 0;
 for (const row of rows) {{
     const matches = (!query || row.dataset.search.includes(query)) && (!effect || row.cells[3]?.textContent === effect);
     row.style.display = matches ? '' : 'none';
     if (row.style.display !== 'none') shown++;
   }}
   count.textContent = (query || effect ? 'Showing matching ' : 'Showing all ') + shown.toLocaleString() + ' of {len(coding):,} protein-altering variants';
 }}
 search.addEventListener('input', update); consequence.addEventListener('change', update);
}})();
</script>
<h2>CNV overview</h2><p>{escape(cnv_status)}</p>
<p>Major-event screen: absolute segment log2 ratio ≥ 0.30, interval ≥ 5 Mb, at least 10 bins when the field is available; segments with absolute log2 ratio &gt; 5 are flagged rather than summarized as events.</p>
{cnv_table_html}
<p class="small">CNVkit baseline: {escape(reference_mode)}. Segment deviations are not absolute copy-number calls. Review with matched/pooled normals and inspect the underlying CNVkit plots before interpreting.</p>
</main></body></html>"""
    Path(args.output).write_text(document, encoding="utf-8")


if __name__ == "__main__":
    main()
