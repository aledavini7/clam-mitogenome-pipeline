#!/usr/bin/env python3
import argparse
import csv
import html
import statistics
from collections import Counter
from pathlib import Path


CONFIDENCE_ORDER = {"HIGH": 0, "MEDIUM": 1, "LOW": 2, "": 3}
MT_FEATURES = [
    (1, 576, "D-loop"), (577, 647, "MT-TF"), (648, 1601, "MT-RNR1"),
    (1602, 1670, "MT-TV"), (1671, 3229, "MT-RNR2"), (3230, 3304, "MT-TL1"),
    (3307, 4262, "MT-ND1"), (4263, 4331, "MT-TI"), (4329, 4400, "MT-TQ"),
    (4402, 4469, "MT-TM"), (4470, 5511, "MT-ND2"), (5512, 5579, "MT-TW"),
    (5587, 5655, "MT-TA"), (5657, 5729, "MT-TN"), (5761, 5826, "MT-TC"),
    (5826, 5891, "MT-TY"), (5904, 7445, "MT-CO1"), (7446, 7514, "MT-TS1"),
    (7518, 7585, "MT-TD"), (7586, 8269, "MT-CO2"), (8295, 8364, "MT-TK"),
    (8366, 8572, "MT-ATP8"), (8527, 9207, "MT-ATP6"), (9207, 9990, "MT-CO3"),
    (9991, 10058, "MT-TG"), (10059, 10404, "MT-ND3"), (10405, 10469, "MT-TR"),
    (10470, 10766, "MT-ND4L"), (10760, 12137, "MT-ND4"), (12138, 12206, "MT-TH"),
    (12207, 12265, "MT-TS2"), (12266, 12336, "MT-TL2"), (12337, 14148, "MT-ND5"),
    (14149, 14673, "MT-ND6"), (14674, 14742, "MT-TE"), (14747, 15887, "MT-CYB"),
    (15888, 15953, "MT-TT"), (15956, 16023, "MT-TP"), (16024, 16569, "D-loop"),
]


def mitochondrial_gene(position):
    pos = as_int(position)
    if pos is None:
        return "Unknown"
    hits = [name for start, end, name in MT_FEATURES if start <= pos <= end]
    return " / ".join(dict.fromkeys(hits)) if hits else "Intergenic"


def as_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def as_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def pct(value, digits=1):
    if value is None:
        return "NA"
    return f"{value:.{digits}f}%"


def num(value, digits=1):
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}f}"


def esc(value):
    return html.escape("" if value is None else str(value))


def short_text(value, limit=100):
    value = "" if value is None else str(value)
    if len(value) <= limit:
        return value
    return value[: limit - 1] + "..."


def pretty_analysis_mode(value):
    labels = {
        "wgs_numt_correction": "with_NUMT_correction",
        "rcrs_only": "without_NUMT_correction",
    }
    return labels.get(value, value)


def pretty_filter_value(value):
    if value in ("", "."):
        return "not evaluated"
    return value


def read_tsv(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def parse_vcf(path):
    counts = Counter()
    records = 0
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            records += 1
            counts[fields[6] or "."] += 1
    return {"records": records, "filters": counts}


def coverage_stats(path, mt_length):
    depths = [0] * mt_length
    observed_contigs = Counter()
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            pos = as_int(fields[1])
            depth = as_int(fields[2])
            if pos is None or depth is None:
                continue
            observed_contigs[fields[0]] += 1
            if 1 <= pos <= mt_length:
                depths[pos - 1] = depth

    covered = sum(1 for value in depths if value > 0)
    total_depth = sum(depths)
    if mt_length:
        mean_depth = total_depth / mt_length
    else:
        mean_depth = None

    return {
        "contigs": observed_contigs,
        "covered_bases": covered,
        "mean_depth": mean_depth,
        "median_depth": statistics.median(depths) if depths else None,
        "min_depth": min(depths) if depths else None,
        "max_depth": max(depths) if depths else None,
        "breadth_1x": covered / mt_length * 100 if mt_length else None,
        "breadth_30x": sum(1 for value in depths if value >= 30) / mt_length * 100 if mt_length else None,
        "breadth_100x": sum(1 for value in depths if value >= 100) / mt_length * 100 if mt_length else None,
    }


def parse_haplogroup(path):
    rows = read_tsv(path)
    if not rows:
        return {"haplogroup": "NA", "quality": "NA", "sample": "NA"}
    row = rows[0]
    return {
        "sample": row.get("SampleID", "NA"),
        "haplogroup": row.get("Haplogroup", "NA"),
        "quality": row.get("Quality", "NA"),
    }


def summarize_variants(rows):
    confidence = Counter(row.get("confidence_tier", "") for row in rows)
    callers = Counter(row.get("callers", "") for row in rows)
    heteroplasmy = [
        as_float(row.get("consensus_heteroplasmy_pct_0_100"))
        for row in rows
        if as_float(row.get("consensus_heteroplasmy_pct_0_100")) is not None
    ]
    return {
        "count": len(rows),
        "confidence": confidence,
        "callers": callers,
        "heteroplasmy": heteroplasmy,
    }


def summarize_mitomap(rows, enabled):
    if not enabled:
        return {"enabled": False, "match": Counter(), "type": Counter()}
    return {
        "enabled": True,
        "match": Counter(row.get("mitomap_match", "") for row in rows),
        "type": Counter(row.get("mitomap_match_type", "") for row in rows),
    }


def summarize_numt_exploration(summary_rows, contig_rows, window_rows, enabled):
    if not enabled:
        return {
            "enabled": False,
            "summary": {},
            "contigs": [],
            "windows": [],
        }

    summary = summary_rows[0] if summary_rows else {}

    contigs = sorted(
        contig_rows,
        key=lambda row: as_int(row.get("alignments")) or 0,
        reverse=True,
    )
    windows = sorted(
        window_rows,
        key=lambda row: (
            as_int(row.get("unique_read_ids")) or 0,
            as_int(row.get("alignments")) or 0,
        ),
        reverse=True,
    )

    return {
        "enabled": True,
        "summary": summary,
        "contigs": contigs,
        "windows": windows,
    }


def summarize_bam_qc(rows):
    by_stage = {}
    for row in rows:
        stage = row.get("stage", "unknown")
        by_stage[stage] = {
            "total": as_int(row.get("total_alignments")),
            "mapped": as_int(row.get("mapped_alignments")),
            "primary": as_int(row.get("primary_alignments")),
            "properly_paired": as_int(row.get("properly_paired_alignments")),
            "duplicates": as_int(row.get("duplicate_alignments")),
            "mt_contig": as_int(row.get("mt_contig_alignments")),
        }
    return by_stage


def vcf_filter_summary(vcf_summary):
    if not vcf_summary["filters"]:
        return "No records"
    return ", ".join(
        f"{pretty_filter_value(key)}: {value}"
        for key, value in sorted(vcf_summary["filters"].items())
    )


def bar(label, value, maximum, color="#2563eb"):
    width = 0 if not maximum else max(0, min(100, value / maximum * 100))
    return f"""
    <div class="bar-row">
      <div class="bar-label">{esc(label)}</div>
      <div class="bar-track"><span style="width:{width:.1f}%; background:{color};"></span></div>
      <div class="bar-value">{esc(value)}</div>
    </div>
    """


def metric_card(title, value, subtitle=""):
    return f"""
    <article class="metric-card">
      <div class="metric-title">{esc(title)}</div>
      <div class="metric-value">{esc(value)}</div>
      <div class="metric-subtitle">{esc(subtitle)}</div>
    </article>
    """


def simple_table(headers, rows, css_class="compact-table"):
    header_html = "".join(f"<th>{esc(header)}</th>" for header in headers)
    row_html = []
    for row in rows:
        row_html.append("<tr>" + "".join(f"<td>{esc(value)}</td>" for value in row) + "</tr>")
    if not row_html:
        row_html.append(f"<tr><td colspan=\"{len(headers)}\">No rows available</td></tr>")
    return f"""
    <table class="{css_class}">
      <thead><tr>{header_html}</tr></thead>
      <tbody>{''.join(row_html)}</tbody>
    </table>
    """


def top_variant_rows(rows, limit=None):
    selected = sorted(rows, key=lambda row: (
        CONFIDENCE_ORDER.get(row.get("confidence_tier", ""), 3),
        as_int(row.get("position")) if as_int(row.get("position")) is not None else float("inf"),
        row.get("ref", ""), row.get("alt", ""),
        0 if row.get("callers", "") == "both" else 1,
    ))
    if limit is not None:
        selected = selected[:limit]
    table_rows = []
    for row in selected:
        table_rows.append(
            [
                row.get("position", ""),
                row.get("ref", ""),
                row.get("alt", ""),
                row.get("callers", ""),
                row.get("confidence_tier", ""),
                row.get("coverage_depth", ""),
                row.get("mutect2_depth", ""),
                row.get("mutect2_alt_depth", ""),
                row.get("mutect2_heteroplasmy_pct_0_100", ""),
                row.get("mutserve_depth", ""),
                row.get("mutserve_alt_depth", ""),
                row.get("mutserve_heteroplasmy_pct_0_100", ""),
                row.get("consensus_heteroplasmy_pct_0_100", ""),
                row.get("mitomap_match", ""),
                row.get("mitomap_match_type", ""),
                row.get("mitomap_gbfreq", ""),
                short_text(row.get("confidence_reason", ""), 60),
            ]
        )
    return table_rows


def mitomap_frequency_display(raw_counts, match, match_type, denominator):
    if match != "YES":
        return "NA"
    if not raw_counts:
        return "NA"

    values = []
    for item in str(raw_counts).split("||"):
        token = item.strip()
        if not token:
            continue
        try:
            count = float(token)
        except ValueError:
            continue
        if count < 0:
            continue
        pct_value = count / denominator * 100
        label = f"{pct_value:.3f}% ({count:g}/{denominator:,})"
        if label not in values:
            values.append(label)
    if not values:
        return "NA"
    if match_type != "exact":
        return "Position only"
    return "; ".join(values)


def compact_allele(value, max_bases=5):
    value = str(value)
    return value if len(value) <= max_bases else value[:max_bases] + "..."


def interactive_variant_table(rows, gb_denominator):
    headers = ["Position", "Ref", "Alt", "Gene / region", "Callers", "Confidence",
               "BAM coverage", "Mutect2 depth", "Mutect2 ALT", "Mutect2 het %",
               "mutserve depth", "mutserve ALT", "mutserve het %", "Mean caller het %",
               "MITOMAP", "MITOMAP match", f"MITOMAP frequency (% FL; n/{gb_denominator:,})",
               "Reason"]
    variants = top_variant_rows(rows)
    out_rows = []
    for idx, values in enumerate(variants):
        confidence, callers, match = values[4], values[3], values[13]
        cls = {"HIGH": "confidence-high", "MEDIUM": "confidence-medium", "LOW": "confidence-low"}.get(confidence, "confidence-unscored")
        if callers != "mutect2,mutserve":
            cls += " single-caller"
        frequency_display = mitomap_frequency_display(values[15], values[13], values[14], gb_denominator)
        searchable = " ".join(str(value) for value in values).lower()
        gene = mitochondrial_gene(values[0])
        # Keep the region beside the allele so a reader sees the affected gene
        # without scrolling across the annotation columns.
        display_values = values[:2] + [compact_allele(values[2])] + [gene] + values[3:15] + [frequency_display, values[16]]
        searchable = (searchable + " " + gene + " " + frequency_display).lower()
        cells = "".join(
            f'<td title="Full ALT: {esc(values[2])}">{esc(value)}</td>' if cell_index == 2
            else f'<td>{esc(value)}</td>'
            for cell_index, value in enumerate(display_values)
        )
        out_rows.append(f'<tr class="{cls}" data-search="{esc(searchable)}" data-confidence="{esc(confidence)}" data-callers="{esc(callers)}" data-mitomap="{esc(match)}" data-gene="{esc(gene)}">{cells}</tr>')
    genes = sorted({mitochondrial_gene(row[0]) for row in variants})
    return f'''<div class="table-controls">
<label>Search variants <input id="mt-search" type="search" placeholder="Position, allele, caller, MITOMAP ID"></label>
<label>Confidence <select id="mt-confidence"><option value="">All tiers</option><option>HIGH</option><option>MEDIUM</option><option>LOW</option><option value="other">Unscored</option></select></label>
<label>Caller support <select id="mt-callers"><option value="">All caller support</option><option value="mutect2,mutserve">Both callers</option><option value="mutect2">Mutect2 only</option><option value="mutserve">mutserve only</option></select></label>
<label>MITOMAP <select id="mt-mitomap"><option value="">All MITOMAP states</option><option>YES</option><option>NO</option></select></label>
<label>Gene / region <select id="mt-gene"><option value="">All mtDNA genes / regions</option>{''.join(f'<option>{esc(gene)}</option>' for gene in genes)}</select></label>
<span id="mt-count">Showing all {len(variants):,} variants</span></div>
<div class="table-wrap"><table id="mt-variants"><thead><tr>{''.join(f'<th>{esc(h)}</th>' for h in headers)}</tr></thead><tbody>{''.join(out_rows) if out_rows else f'<tr><td colspan="{len(headers)}">No mitochondrial variants found.</td></tr>'}</tbody></table></div>'''


def top_numt_contig_rows(rows, limit=10):
    table_rows = []
    for row in rows[:limit]:
        table_rows.append(
            [
                row.get("nuclear_contig", ""),
                row.get("alignments", ""),
                row.get("high_mapq_alignments", ""),
                row.get("min_pos", ""),
                row.get("max_pos", ""),
                row.get("mean_mapq", ""),
            ]
        )
    return table_rows


def top_numt_window_rows(rows, limit=10):
    table_rows = []
    for row in rows[:limit]:
        table_rows.append(
            [
                row.get("nuclear_contig", ""),
                row.get("window_start", ""),
                row.get("window_end", ""),
                row.get("alignments", ""),
                row.get("unique_read_ids", ""),
                row.get("high_mapq_alignments", ""),
                row.get("mean_mapq", ""),
            ]
        )
    return table_rows


def heteroplasmy_bins(values):
    bins = [
        ("0-1%", 0, 1),
        ("1-5%", 1, 5),
        ("5-20%", 5, 20),
        ("20-80%", 20, 80),
        ("80-100%", 80, 101),
    ]
    counts = []
    for label, low, high in bins:
        counts.append((label, sum(1 for value in values if low <= value < high)))
    return counts


def write_qc_summary(path, metrics):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            writer.writerow([key, value])


def build_html(args, data):
    final_cov = data["coverage_final"]
    before_cov = data["coverage_before"]
    filtered = data["filtered"]
    all_calls = data["all"]
    mutect2_vcf = data["mutect2_vcf"]
    mutserve_vcf = data["mutserve_vcf"]
    mutect2_haplo = data["mutect2_haplo"]
    mutserve_haplo = data["mutserve_haplo"]
    mitomap = data["mitomap"]
    bam_qc = data["bam_qc"]
    numt_exploration = data["numt_exploration"]
    analysis_mode_label = pretty_analysis_mode(args.analysis_mode)

    high = filtered["confidence"].get("HIGH", 0)
    medium = filtered["confidence"].get("MEDIUM", 0)
    final_stage = "after" if as_bool(args.run_numt_correction) else "final"
    final_bam = bam_qc.get(final_stage, {})
    before_bam = bam_qc.get("before", {})

    cards = [
        metric_card("NUMT correction", "YES" if as_bool(args.run_numt_correction) else "NO", ""),
        metric_card("Final Mean Depth", num(final_cov["mean_depth"], 1), f"Median {num(final_cov['median_depth'], 1)}"),
        metric_card("Breadth >=30x", pct(final_cov["breadth_30x"], 1), f"Covered bases: {final_cov['covered_bases']}"),
        metric_card("Filtered Calls", filtered["count"], f"HIGH {high} | MEDIUM {medium}"),
        metric_card("Mutect2 VCF Records", mutect2_vcf["records"], vcf_filter_summary(mutect2_vcf)),
        metric_card("Mutserve VCF Records", mutserve_vcf["records"], vcf_filter_summary(mutserve_vcf)),
    ]

    if mitomap["enabled"]:
        cards.append(
            metric_card(
                "MITOMAP Matches",
                mitomap["match"].get("YES", 0),
                f"Exact {mitomap['type'].get('exact', 0)} | Position {mitomap['type'].get('position', 0)}",
            )
        )
    else:
        cards.append(metric_card("MITOMAP", "Not requested", "Provide --mitomap_variant_table to enable"))

    if final_bam:
        cards.append(metric_card("Final BAM Alignments", final_bam.get("total", "NA"), f"Mapped {final_bam.get('mapped', 'NA')}"))

    if numt_exploration["enabled"]:
        numt_summary = numt_exploration["summary"]
        cards.append(
            metric_card(
                "NUMT-like Reads",
                numt_summary.get("mitoscape_rejected_read_ids", "NA"),
                f"{numt_summary.get('mitoscape_rejected_fraction_pct', 'NA')}% of pre-MitoScape read IDs",
            )
        )

    confidence_max = max([1] + list(all_calls["confidence"].values()))
    confidence_bars = "".join(
        bar(label, all_calls["confidence"].get(label, 0), confidence_max, color)
        for label, color in [("HIGH", "#15803d"), ("MEDIUM", "#ca8a04"), ("LOW", "#dc2626")]
    )

    caller_max = max([1] + list(all_calls["callers"].values()))
    caller_bars = "".join(
        bar(label or "unlabelled", value, caller_max, "#4f46e5")
        for label, value in sorted(all_calls["callers"].items())
    )

    heteroplasmy_counts = heteroplasmy_bins(filtered["heteroplasmy"])
    heteroplasmy_max = max([1] + [count for _, count in heteroplasmy_counts])
    heteroplasmy_bars = "".join(bar(label, count, heteroplasmy_max, "#0891b2") for label, count in heteroplasmy_counts)

    coverage_rows = [
        ["Before NUMT correction" if as_bool(args.run_numt_correction) else "Final", num(before_cov["mean_depth"], 1), num(before_cov["median_depth"], 1), before_cov["min_depth"], before_cov["max_depth"], pct(before_cov["breadth_30x"], 1), pct(before_cov["breadth_100x"], 1)],
    ]
    if as_bool(args.run_numt_correction):
        coverage_rows.append(["After NUMT correction", num(final_cov["mean_depth"], 1), num(final_cov["median_depth"], 1), final_cov["min_depth"], final_cov["max_depth"], pct(final_cov["breadth_30x"], 1), pct(final_cov["breadth_100x"], 1)])

    bam_rows = []
    for stage, values in bam_qc.items():
        bam_rows.append(
            [
                stage,
                values.get("total", "NA"),
                values.get("mapped", "NA"),
                values.get("primary", "NA"),
                values.get("properly_paired", "NA"),
                values.get("duplicates", "NA"),
                values.get("mt_contig", "NA"),
            ]
        )

    haplo_rows = [
        ["Mutect2", mutect2_haplo["sample"], mutect2_haplo["haplogroup"], mutect2_haplo["quality"]],
        ["mutserve", mutserve_haplo["sample"], mutserve_haplo["haplogroup"], mutserve_haplo["quality"]],
    ]

    vcf_rows = [
        ["Mutect2 raw", mutect2_vcf["records"], vcf_filter_summary(mutect2_vcf)],
        ["Mutect2 filtered", data["filtered_mutect2_vcf"]["records"], vcf_filter_summary(data["filtered_mutect2_vcf"])],
        ["mutserve", mutserve_vcf["records"], vcf_filter_summary(mutserve_vcf)],
    ]

    output_links = [
        ["BAMs", "../bams/"],
        ["Coverage", "../coverage/"],
        ["Variants", "../variants/"],
        ["Annotation TSVs", "../annotation/"],
        ["MAFs", "../annotation/mafs/"],
        ["MITOMAP TSVs", "../annotation/mitomap/" if mitomap["enabled"] else "not requested"],
        ["NUMT-like exploration", "../numt_exploration/" if numt_exploration["enabled"] else "not requested"],
        ["Haplogroups", "../haplogroups/"],
        ["QC TSVs", "../qc/"],
    ]

    numt_section = ""
    if numt_exploration["enabled"]:
        numt_summary = numt_exploration["summary"]
        numt_summary_rows = [
            ["Pre-MitoScape unique read IDs", numt_summary.get("pre_mitoscape_unique_read_ids", "NA")],
            ["Post-MitoScape unique read IDs", numt_summary.get("post_mitoscape_unique_read_ids", "NA")],
            ["MitoScape-rejected read IDs", numt_summary.get("mitoscape_rejected_read_ids", "NA")],
            ["Rejected fraction (%)", numt_summary.get("mitoscape_rejected_fraction_pct", "NA")],
            ["Rejected nuclear alignments", numt_summary.get("rejected_nuclear_alignments", "NA")],
            ["High-MAPQ nuclear alignments", numt_summary.get("rejected_nuclear_high_mapq_alignments", "NA")],
            ["Interpretation label", numt_summary.get("interpretation", "NA")],
        ]
        numt_section = f"""
    <section class="panel">
      <h2>NUMT-like Read Exploration</h2>
      <div class="notice">This section describes MitoScape-rejected reads and their nuclear remapping pattern. It is an exploratory signal, not a confirmed NUMT insertion call.</div>
      {simple_table(["Metric", "Value"], numt_summary_rows)}
      <h2>Top Nuclear Contigs</h2>
      {simple_table(["Contig", "Alignments", "High-MAPQ", "Min pos", "Max pos", "Mean MAPQ"], top_numt_contig_rows(numt_exploration["contigs"]))}
      <h2>Top Nuclear Windows</h2>
      {simple_table(["Contig", "Window start", "Window end", "Alignments", "Unique read IDs", "High-MAPQ", "Mean MAPQ"], top_numt_window_rows(numt_exploration["windows"]))}
    </section>
        """

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CLAM Report - {esc(args.sample_id)}</title>
  <style>
    :root {{
      --ink: #16202a;
      --muted: #667085;
      --line: #d9e2ec;
      --panel: #ffffff;
      --bg: #f5f7fb;
      --accent: #1d4ed8;
    }}
    body {{
      margin: 0;
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--ink);
      background: var(--bg);
      letter-spacing: 0;
    }}
    header {{
      padding: 30px 36px 22px;
      background: #0f172a;
      color: #fff;
    }}
    header h1 {{
      margin: 0;
      font-size: 30px;
      font-weight: 760;
      letter-spacing: 0;
    }}
    header p {{
      margin: 8px 0 0;
      color: #cbd5e1;
    }}
    main {{
      max-width: 1220px;
      margin: 0 auto;
      padding: 24px 24px 44px;
    }}
    section {{
      margin-top: 20px;
    }}
    h2 {{
      margin: 0 0 12px;
      font-size: 20px;
    }}
    .metrics {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
      gap: 12px;
    }}
    .metric-card, .panel {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
    }}
    .metric-card {{
      padding: 15px 16px;
      min-height: 88px;
    }}
    .metric-title {{
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
      font-weight: 700;
    }}
    .metric-value {{
      margin-top: 8px;
      font-size: 28px;
      line-height: 1.05;
      font-weight: 760;
    }}
    .metric-subtitle {{
      margin-top: 8px;
      font-size: 13px;
      color: var(--muted);
    }}
    .grid-2 {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 16px;
    }}
    .panel {{
      padding: 16px;
      overflow-x: auto;
    }}
    .bar-row {{
      display: grid;
      grid-template-columns: 120px minmax(120px, 1fr) 52px;
      align-items: center;
      gap: 10px;
      margin: 10px 0;
      font-size: 13px;
    }}
    .bar-track {{
      height: 10px;
      border-radius: 999px;
      background: #e5e7eb;
      overflow: hidden;
    }}
    .bar-track span {{
      display: block;
      height: 100%;
      border-radius: 999px;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 13px;
    }}
    th, td {{
      padding: 8px 9px;
      border-bottom: 1px solid #e5e7eb;
      text-align: left;
      vertical-align: top;
    }}
    th {{
      color: #475569;
      background: #f8fafc;
      font-size: 12px;
      text-transform: uppercase;
    }}
    .table-wrap {{ overflow-x: auto; max-height: 70vh; }}
    .table-controls {{ display: flex; gap: 12px; align-items: end; flex-wrap: wrap; margin: 12px 0; }}
    .table-controls label {{ display: flex; flex-direction: column; gap: 4px; color: var(--muted); font-size: 12px; }}
    .table-controls input, .table-controls select {{ padding: 8px; border: 1px solid var(--line); border-radius: 6px; min-width: 165px; background: white; }}
    .table-controls span {{ color: var(--muted); font-size: 12px; }}
    .confidence-high {{ background: #ecfdf5; }} .confidence-medium {{ background: #fffbeb; }} .confidence-low {{ background: #fef2f2; }}
    .confidence-unscored {{ background: #f1f5f9; }} .single-caller td:nth-child(5) {{ font-weight: 700; color: #9a3412; }}
    #mt-variants td:nth-child(4) {{ white-space: nowrap; }}
    .notice {{
      background: #fff7ed;
      border: 1px solid #fed7aa;
      color: #7c2d12;
      border-radius: 8px;
      padding: 12px 14px;
      margin: 8px 0;
    }}
    .links a {{
      color: var(--accent);
      text-decoration: none;
      font-weight: 650;
    }}
  </style>
</head>
<body>
  <header>
    <h1>CLAM Mitochondrial Report</h1>
    <p>Sample {esc(args.sample_id)} | Reference: {esc(args.reference_set)}</p>
  </header>
  <main>
    <section class="metrics">
      {''.join(cards)}
    </section>

    <section class="grid-2">
      <div class="panel">
        <h2>Confidence Tiers</h2>
        {confidence_bars}
      </div>
      <div class="panel">
        <h2>Caller Support</h2>
        {caller_bars}
      </div>
      <div class="panel">
        <h2>Filtered Heteroplasmy</h2>
        {heteroplasmy_bars}
      </div>
      <div class="panel">
        <h2>Haplogroups</h2>
        {simple_table(["Caller", "HaploGrep sample", "Haplogroup", "Quality"], haplo_rows)}
      </div>
    </section>

    <section class="panel">
      <h2>Coverage</h2>
      {simple_table(["Stage", "Mean depth", "Median depth", "Min", "Max", "Breadth >=30x", "Breadth >=100x"], coverage_rows)}
    </section>

    <section class="panel">
      <h2>BAM QC</h2>
      {simple_table(["Stage", "Total", "Mapped", "Primary", "Proper pairs", "Duplicates", args.mt_contig], bam_rows)}
    </section>

    <section class="panel">
      <h2>Variant Calling</h2>
      {simple_table(["File", "Records", "Filters"], vcf_rows)}
    </section>

    <section class="panel">
      <h2>Mitochondrial variants</h2>
      {interactive_variant_table(data["top_rows"], args.mitomap_gb_denominator)}
    </section>

    {numt_section}

    <section class="panel links">
      <h2>Output Index</h2>
      {simple_table(["Category", "Relative path"], output_links)}
    </section>
  </main>
  <script>
  (() => {{
    const table = document.getElementById('mt-variants'); if (!table) return;
    const rows = [...table.querySelectorAll('tbody tr')];
    const search = document.getElementById('mt-search'), confidence = document.getElementById('mt-confidence');
    const callers = document.getElementById('mt-callers'), mitomap = document.getElementById('mt-mitomap'), gene = document.getElementById('mt-gene'), count = document.getElementById('mt-count');
    function update() {{
      const q = search.value.trim().toLowerCase(), tier = confidence.value, caller = callers.value, match = mitomap.value, geneFilter = gene.value;
      let shown = 0;
      rows.forEach((row, i) => {{
        const scored = row.dataset.confidence || '';
        const ok = (!q || row.dataset.search.includes(q)) && (!tier || (tier === 'other' ? !['HIGH','MEDIUM','LOW'].includes(scored) : scored === tier)) && (!caller || row.dataset.callers === caller) && (!match || row.dataset.mitomap === match) && (!geneFilter || row.dataset.gene === geneFilter);
        row.style.display = ok ? '' : 'none';
        if (row.style.display !== 'none') shown++;
      }});
      count.textContent = `Showing ${{shown.toLocaleString()}} of {len(data["top_rows"]):,} variants`;
    }}
    [search, confidence, callers, mitomap, gene].forEach(el => el.addEventListener(el.tagName === 'INPUT' ? 'input' : 'change', update));
  }})();
  </script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description="Build a self-contained CLAM HTML report.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--analysis-mode", required=True)
    parser.add_argument("--input-type", required=True)
    parser.add_argument("--reference-set", required=True)
    parser.add_argument("--run-numt-correction", required=True)
    parser.add_argument("--mt-length", type=int, default=16569)
    parser.add_argument("--mt-contig", default="MT")
    parser.add_argument("--bam-qc", required=True)
    parser.add_argument("--coverage-before", required=True)
    parser.add_argument("--coverage-final", required=True)
    parser.add_argument("--mutect2-vcf", required=True)
    parser.add_argument("--filtered-mutect2-vcf", required=True)
    parser.add_argument("--mutserve-vcf", required=True)
    parser.add_argument("--mutserve-summary", required=True)
    parser.add_argument("--mutect2-haplogroups", required=True)
    parser.add_argument("--mutserve-haplogroups", required=True)
    parser.add_argument("--mutect2-tsv", required=True)
    parser.add_argument("--mutserve-tsv", required=True)
    parser.add_argument("--all-tsv", required=True)
    parser.add_argument("--confidence-filtered-tsv", required=True)
    parser.add_argument("--has-mitomap", required=True)
    parser.add_argument("--mitomap-confidence-tsv", required=True)
    parser.add_argument("--mitomap-all-tsv", required=True)
    parser.add_argument("--mitomap-gb-denominator", type=int, default=65286,
                        help="Full-length GenBank sequence count for this MITOMAP snapshot (default: 65,286 for the 2026-05-22 table).")
    parser.add_argument("--has-numt-exploration", required=True)
    parser.add_argument("--numt-summary-tsv", required=True)
    parser.add_argument("--numt-contig-tsv", required=True)
    parser.add_argument("--numt-window-tsv", required=True)
    parser.add_argument("--output-html", required=True)
    parser.add_argument("--output-qc", required=True)
    args = parser.parse_args()
    if args.mitomap_gb_denominator < 1:
        parser.error("--mitomap-gb-denominator must be a positive integer")

    all_rows = read_tsv(args.all_tsv)
    filtered_rows = read_tsv(args.confidence_filtered_tsv)
    mitomap_enabled = as_bool(args.has_mitomap)
    mitomap_rows = read_tsv(args.mitomap_all_tsv) if mitomap_enabled else []
    top_rows = mitomap_rows if mitomap_enabled else all_rows
    numt_exploration_enabled = as_bool(args.has_numt_exploration)
    numt_summary_rows = read_tsv(args.numt_summary_tsv) if numt_exploration_enabled else []
    numt_contig_rows = read_tsv(args.numt_contig_tsv) if numt_exploration_enabled else []
    numt_window_rows = read_tsv(args.numt_window_tsv) if numt_exploration_enabled else []

    data = {
        "bam_qc": summarize_bam_qc(read_tsv(args.bam_qc)),
        "coverage_before": coverage_stats(args.coverage_before, args.mt_length),
        "coverage_final": coverage_stats(args.coverage_final, args.mt_length),
        "mutect2_vcf": parse_vcf(args.mutect2_vcf),
        "filtered_mutect2_vcf": parse_vcf(args.filtered_mutect2_vcf),
        "mutserve_vcf": parse_vcf(args.mutserve_vcf),
        "mutect2": summarize_variants(read_tsv(args.mutect2_tsv)),
        "mutserve": summarize_variants(read_tsv(args.mutserve_tsv)),
        "all": summarize_variants(all_rows),
        "filtered": summarize_variants(filtered_rows),
        "mitomap": summarize_mitomap(mitomap_rows, mitomap_enabled),
        "numt_exploration": summarize_numt_exploration(
            numt_summary_rows,
            numt_contig_rows,
            numt_window_rows,
            numt_exploration_enabled,
        ),
        "mutect2_haplo": parse_haplogroup(args.mutect2_haplogroups),
        "mutserve_haplo": parse_haplogroup(args.mutserve_haplogroups),
        "top_rows": top_rows,
    }

    final_cov = data["coverage_final"]
    qc_metrics = {
        "sample_id": args.sample_id,
        "analysis_mode": args.analysis_mode,
        "analysis_mode_label": pretty_analysis_mode(args.analysis_mode),
        "input_type": args.input_type,
        "reference_set": args.reference_set,
        "numt_correction": args.run_numt_correction,
        "final_mean_depth": num(final_cov["mean_depth"], 3),
        "final_median_depth": num(final_cov["median_depth"], 3),
        "final_breadth_30x_pct": pct(final_cov["breadth_30x"], 3),
        "final_breadth_100x_pct": pct(final_cov["breadth_100x"], 3),
        "all_variant_count": data["all"]["count"],
        "confidence_filtered_variant_count": data["filtered"]["count"],
        "high_confidence_count": data["filtered"]["confidence"].get("HIGH", 0),
        "medium_confidence_count": data["filtered"]["confidence"].get("MEDIUM", 0),
        "mitomap_enabled": mitomap_enabled,
        "mitomap_exact_matches": data["mitomap"]["type"].get("exact", 0),
        "mitomap_position_matches": data["mitomap"]["type"].get("position", 0),
        "numt_exploration_enabled": numt_exploration_enabled,
        "numt_like_rejected_read_ids": data["numt_exploration"]["summary"].get("mitoscape_rejected_read_ids", "NA"),
        "numt_like_rejected_fraction_pct": data["numt_exploration"]["summary"].get("mitoscape_rejected_fraction_pct", "NA"),
        "numt_like_high_mapq_nuclear_alignments": data["numt_exploration"]["summary"].get("rejected_nuclear_high_mapq_alignments", "NA"),
        "mutect2_haplogroup": data["mutect2_haplo"]["haplogroup"],
        "mutserve_haplogroup": data["mutserve_haplo"]["haplogroup"],
    }

    write_qc_summary(args.output_qc, qc_metrics)
    html_text = build_html(args, data)
    with open(args.output_html, "w", encoding="utf-8") as handle:
        handle.write(html_text)


if __name__ == "__main__":
    main()
