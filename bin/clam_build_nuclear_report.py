#!/usr/bin/env python3
"""Build a compact per-sample HTML summary for the nuclear WES workflow."""

import argparse
import csv
import html
from pathlib import Path


def read_tsv(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--calling-mode", choices=("tumor_only", "paired", "germline"), required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--maf", required=True)
    parser.add_argument("--cnv-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    variants = read_tsv(args.maf)
    cnv_dir = Path(args.cnv_dir)
    segments = []
    for segment_file in sorted(cnv_dir.glob("*.seg")):
        try:
            segments.extend(read_tsv(segment_file))
        except (OSError, csv.Error):
            continue
    variant_rows = "".join(
        "<tr>" + "".join(f"<td>{html.escape(row.get(col, ''))}</td>" for col in
                            ("Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "HGVSp", "gnomAD_AF", "COSMIC_ID")) + "</tr>"
        for row in variants
    ) or '<tr><td colspan="9">No PASS variants in the annotated MAF.</td></tr>'
    cnv_rows = "".join(
        "<tr>" + "".join(f"<td>{html.escape(row.get(key, ''))}</td>" for key in ("ID", "chrom", "loc.start", "loc.end", "seg.mean")) + "</tr>"
        for row in segments
    ) or '<tr><td colspan="5">No CNVkit segment table found.</td></tr>'
    mode_note = (
        "Germline-style calls identify variants present in this sample; they do not establish constitutional inheritance. "
        "In T-cell DNA, some calls may be acquired somatic mutations."
        if args.calling_mode == "germline"
        else "Somatic call mode: " + args.calling_mode.replace("_", "-") + "."
    )
    document = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>CLAM nDNA report: {html.escape(args.sample_id)}</title>
<style>body{{font:15px system-ui,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#222}}table{{border-collapse:collapse;width:100%;margin:1rem 0 2rem}}th,td{{border:1px solid #ccc;padding:.4rem;text-align:left}}th{{background:#f2f4f7}}code{{word-break:break-all}}</style></head>
<body><h1>Nuclear WES report</h1><p><b>Sample:</b> {html.escape(args.sample_id)}</p><p><b>Calling mode:</b> {html.escape(args.calling_mode)}</p>
<p><b>Interpretation note:</b> {html.escape(mode_note)}</p>
<p><b>Annotated VCF:</b> <code>{html.escape(Path(args.vcf).name)}</code></p><p><b>PASS variants in MAF:</b> {len(variants)}</p>
<h2>Variants</h2><table><thead><tr><th>Gene</th><th>Chromosome</th><th>Position</th><th>Ref</th><th>Alt</th><th>Class</th><th>Protein change</th><th>gnomAD AF</th><th>COSMIC</th></tr></thead><tbody>{variant_rows}</tbody></table>
<h2>CNVkit segments</h2><p>Segments reported: {len(segments)}. Review the accompanying CNVkit files and interpret against suitable sample controls/reference.</p>
<table><thead><tr><th>ID</th><th>Chromosome</th><th>Start</th><th>End</th><th>Segment mean</th></tr></thead><tbody>{cnv_rows}</tbody></table>
<p>This report is a workflow summary and is not a clinical interpretation.</p></body></html>"""
    Path(args.output).write_text(document, encoding="utf-8")


if __name__ == "__main__":
    main()
