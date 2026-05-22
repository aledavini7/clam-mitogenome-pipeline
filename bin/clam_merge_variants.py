#!/usr/bin/env python3
import argparse
import csv
import math
from collections import defaultdict


def split_number_list(value):
    if value in (None, "", "."):
        return []

    values = []
    for item in value.replace("[", "").replace("]", "").split(","):
        item = item.strip()
        if item in ("", "."):
            values.append(None)
            continue
        try:
            values.append(float(item))
        except ValueError:
            values.append(None)
    return values


def as_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def normalize_af(value):
    if value is None:
        return None
    if 1.0 < value <= 100.0:
        return value / 100.0
    return value


def parse_info(info_text):
    info = {}
    if info_text in ("", "."):
        return info

    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def parse_sample(format_text, sample_text):
    if format_text in ("", ".") or sample_text in ("", "."):
        return {}

    keys = format_text.split(":")
    values = sample_text.split(":")
    return {key: values[idx] if idx < len(values) else None for idx, key in enumerate(keys)}


def pick(values, idx):
    return values[idx] if idx < len(values) else None


def parse_vcf(path, caller):
    calls = {}

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom, pos_text, _id, ref, alt_text, qual, filt, info_text = fields[:8]
            fmt = fields[8] if len(fields) > 8 else "."
            sample = fields[9] if len(fields) > 9 else "."
            info = parse_info(info_text)
            sample_data = parse_sample(fmt, sample)

            alts = alt_text.split(",")
            ad_values = split_number_list(sample_data.get("AD"))
            af_values = split_number_list(sample_data.get("AF"))
            dp = as_int(sample_data.get("DP")) or as_int(info.get("DP"))
            dp_source = "DP" if dp is not None else ""

            for alt_idx, alt in enumerate(alts):
                if alt in ("", ".", "*") or alt.startswith("<"):
                    continue

                ref_depth = as_int(pick(ad_values, 0))
                alt_depth = as_int(pick(ad_values, alt_idx + 1))
                ref_depth_source = "AD" if ref_depth is not None else ""
                alt_depth_source = "AD" if alt_depth is not None else ""
                af = normalize_af(as_float(pick(af_values, alt_idx)))

                if af is None and dp and alt_depth is not None and dp > 0:
                    af = alt_depth / dp

                if dp is None and ad_values:
                    numeric_ad = [as_int(value) for value in ad_values]
                    numeric_ad = [value for value in numeric_ad if value is not None]
                    dp = sum(numeric_ad) if numeric_ad else None
                    dp_source = "AD_sum" if dp is not None else ""

                if alt_depth is None and af is not None and dp is not None:
                    alt_depth = round(af * dp)
                    alt_depth_source = "AFxDP"

                if ref_depth is None and dp is not None and alt_depth is not None:
                    ref_depth = max(dp - alt_depth, 0)
                    ref_depth_source = "DP-alt"

                key = (chrom, int(pos_text), ref, alt)
                calls[key] = {
                    "caller": caller,
                    "filter": filt,
                    "qual": qual,
                    "dp": dp,
                    "ref_depth": ref_depth,
                    "alt_depth": alt_depth,
                    "dp_source": dp_source,
                    "ref_depth_source": ref_depth_source,
                    "alt_depth_source": alt_depth_source,
                    "af": af,
                }

    return calls


def parse_coverage(path):
    coverage = {}

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom, pos_text, depth_text = fields[:3]
            depth = as_int(depth_text)
            if depth is not None:
                coverage[(chrom, int(pos_text))] = depth

    return coverage


def wilson_interval(alt_depth, total_depth, z=1.96):
    if alt_depth is None or total_depth is None or total_depth <= 0:
        return None, None

    phat = max(0.0, min(1.0, alt_depth / total_depth))
    denominator = 1.0 + z * z / total_depth
    center = (phat + z * z / (2.0 * total_depth)) / denominator
    margin = (
        z
        * math.sqrt((phat * (1.0 - phat) + z * z / (4.0 * total_depth)) / total_depth)
        / denominator
    )
    return max(0.0, center - margin), min(1.0, center + margin)


def fmt_float(value, digits=6):
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def fmt_percent(value):
    if value is None:
        return ""
    return f"{value * 100.0:.3f}"


def is_pass(filter_value):
    return filter_value in (None, "", ".", "PASS")


def confidence_tier(callers, calls, consensus_depth, alt_depth, af_values, args):
    reasons = []

    if consensus_depth is None:
        reasons.append("missing_depth")
    elif consensus_depth < args.medium_depth:
        reasons.append(f"depth<{args.medium_depth}")

    if alt_depth is None:
        reasons.append("missing_alt_depth")
    elif alt_depth < args.medium_alt_reads:
        reasons.append(f"alt_reads<{args.medium_alt_reads}")

    non_pass_filters = [
        f"{name}_filter={calls[name]['filter']}" for name in callers if not is_pass(calls[name]["filter"])
    ]
    reasons.extend(non_pass_filters)

    if len(callers) == 2 and len(af_values) == 2:
        discordance = abs(af_values[0] - af_values[1])
        if discordance > args.max_af_discordance:
            reasons.append(f"caller_af_discordance>{args.max_af_discordance}")
    else:
        discordance = None

    if (
        not non_pass_filters
        and consensus_depth is not None
        and consensus_depth >= args.high_depth
        and alt_depth is not None
        and alt_depth >= args.high_alt_reads
        and (discordance is None or discordance <= args.max_af_discordance)
    ):
        if len(callers) == 2:
            return "HIGH", "both_callers_pass_high_depth_alt_support"
        return "HIGH", f"{callers[0]}_pass_high_depth_alt_support"

    if (
        consensus_depth is not None
        and consensus_depth >= args.medium_depth
        and alt_depth is not None
        and alt_depth >= args.medium_alt_reads
        and not non_pass_filters
        and (discordance is None or discordance <= args.max_af_discordance)
    ):
        return "MEDIUM", "single_or_partial_support_with_acceptable_depth"

    return "LOW", ",".join(reasons) if reasons else "limited_support"


def write_table(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--mutect2-vcf", required=True)
    parser.add_argument("--mutserve-vcf", required=True)
    parser.add_argument("--coverage", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--mutect2-output")
    parser.add_argument("--mutserve-output")
    parser.add_argument("--confidence-filtered-output")
    parser.add_argument("--high-depth", type=int, default=100)
    parser.add_argument("--medium-depth", type=int, default=30)
    parser.add_argument("--high-alt-reads", type=int, default=10)
    parser.add_argument("--medium-alt-reads", type=int, default=3)
    parser.add_argument("--max-af-discordance", type=float, default=0.10)
    args = parser.parse_args()

    caller_calls = defaultdict(dict)
    for key, call in parse_vcf(args.mutect2_vcf, "mutect2").items():
        caller_calls[key]["mutect2"] = call
    for key, call in parse_vcf(args.mutserve_vcf, "mutserve").items():
        caller_calls[key]["mutserve"] = call

    coverage = parse_coverage(args.coverage)

    all_fieldnames = [
        "sample_id",
        "chrom",
        "position",
        "ref",
        "alt",
        "callers",
        "coverage_depth",
        "consensus_depth",
        "consensus_alt_depth",
        "consensus_af_fraction",
        "consensus_heteroplasmy_pct_0_100",
        "heteroplasmy_ci95_low_pct_0_100",
        "heteroplasmy_ci95_high_pct_0_100",
        "confidence_tier",
        "confidence_reason",
        "mutect2_filter",
        "mutect2_af_fraction",
        "mutect2_heteroplasmy_pct_0_100",
        "mutect2_depth",
        "mutect2_ref_depth",
        "mutect2_alt_depth",
        "mutserve_filter",
        "mutserve_af_fraction",
        "mutserve_heteroplasmy_pct_0_100",
        "mutserve_depth",
        "mutserve_ref_depth",
        "mutserve_alt_depth",
    ]

    caller_fieldnames = [
        "sample_id",
        "caller",
        "chrom",
        "position",
        "ref",
        "alt",
        "filter",
        "af_fraction",
        "heteroplasmy_pct_0_100",
        "depth",
        "ref_depth",
        "alt_depth",
        "coverage_depth",
        "confidence_tier",
        "confidence_reason",
        "depth_source",
        "ref_depth_source",
        "alt_depth_source",
    ]

    all_rows = []
    mutect2_rows = []
    mutserve_rows = []

    for key in sorted(caller_calls, key=lambda item: (item[0], item[1], item[2], item[3])):
        chrom, pos, ref, alt = key
        calls = caller_calls[key]
        mutect2 = calls.get("mutect2")
        mutserve = calls.get("mutserve")
        caller_names = sorted(calls)
        af_values = [calls[name]["af"] for name in caller_names if calls[name]["af"] is not None]

        coverage_depth = coverage.get((chrom, pos))
        consensus_af = sum(af_values) / len(af_values) if af_values else None

        depth_values = [calls[name]["dp"] for name in caller_names if calls[name]["dp"] is not None]
        if depth_values:
            consensus_depth = round(sum(depth_values) / len(depth_values))
        else:
            consensus_depth = coverage_depth

        alt_depth_values = [
            calls[name]["alt_depth"] for name in caller_names if calls[name]["alt_depth"] is not None
        ]
        if alt_depth_values:
            consensus_alt_depth = round(sum(alt_depth_values) / len(alt_depth_values))
        elif consensus_af is not None and consensus_depth is not None:
            consensus_alt_depth = round(consensus_af * consensus_depth)
        else:
            consensus_alt_depth = None

        ci_low, ci_high = wilson_interval(consensus_alt_depth, consensus_depth)
        tier, reason = confidence_tier(
            caller_names,
            calls,
            consensus_depth,
            consensus_alt_depth,
            af_values,
            args,
        )

        all_rows.append(
            {
                "sample_id": args.sample_id,
                "chrom": chrom,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "callers": ",".join(caller_names),
                "coverage_depth": coverage_depth if coverage_depth is not None else "",
                "consensus_depth": consensus_depth if consensus_depth is not None else "",
                "consensus_alt_depth": consensus_alt_depth if consensus_alt_depth is not None else "",
                "consensus_af_fraction": fmt_float(consensus_af),
                "consensus_heteroplasmy_pct_0_100": fmt_percent(consensus_af),
                "heteroplasmy_ci95_low_pct_0_100": fmt_percent(ci_low),
                "heteroplasmy_ci95_high_pct_0_100": fmt_percent(ci_high),
                "confidence_tier": tier,
                "confidence_reason": reason,
                "mutect2_filter": mutect2["filter"] if mutect2 else "",
                "mutect2_af_fraction": fmt_float(mutect2["af"]) if mutect2 else "",
                "mutect2_heteroplasmy_pct_0_100": fmt_percent(mutect2["af"]) if mutect2 else "",
                "mutect2_depth": mutect2["dp"] if mutect2 and mutect2["dp"] is not None else "",
                "mutect2_ref_depth": mutect2["ref_depth"] if mutect2 and mutect2["ref_depth"] is not None else "",
                "mutect2_alt_depth": mutect2["alt_depth"] if mutect2 and mutect2["alt_depth"] is not None else "",
                "mutserve_filter": mutserve["filter"] if mutserve else "",
                "mutserve_af_fraction": fmt_float(mutserve["af"]) if mutserve else "",
                "mutserve_heteroplasmy_pct_0_100": fmt_percent(mutserve["af"]) if mutserve else "",
                "mutserve_depth": mutserve["dp"] if mutserve and mutserve["dp"] is not None else "",
                "mutserve_ref_depth": mutserve["ref_depth"] if mutserve and mutserve["ref_depth"] is not None else "",
                "mutserve_alt_depth": mutserve["alt_depth"] if mutserve and mutserve["alt_depth"] is not None else "",
            }
        )

        for caller_name in caller_names:
            call = calls[caller_name]
            caller_row = {
                "sample_id": args.sample_id,
                "caller": caller_name,
                "chrom": chrom,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "filter": call["filter"],
                "af_fraction": fmt_float(call["af"]),
                "heteroplasmy_pct_0_100": fmt_percent(call["af"]),
                "depth": call["dp"] if call["dp"] is not None else "",
                "ref_depth": call["ref_depth"] if call["ref_depth"] is not None else "",
                "alt_depth": call["alt_depth"] if call["alt_depth"] is not None else "",
                "coverage_depth": coverage_depth if coverage_depth is not None else "",
                "confidence_tier": tier,
                "confidence_reason": reason,
                "depth_source": call["dp_source"],
                "ref_depth_source": call["ref_depth_source"],
                "alt_depth_source": call["alt_depth_source"],
            }
            if caller_name == "mutect2":
                mutect2_rows.append(caller_row)
            if caller_name == "mutserve":
                mutserve_rows.append(caller_row)

    confidence_filtered_rows = [
        row for row in all_rows if row["confidence_tier"] in ("HIGH", "MEDIUM")
    ]

    write_table(args.output, all_fieldnames, all_rows)
    if args.mutect2_output:
        write_table(args.mutect2_output, caller_fieldnames, mutect2_rows)
    if args.mutserve_output:
        write_table(args.mutserve_output, caller_fieldnames, mutserve_rows)
    if args.confidence_filtered_output:
        write_table(args.confidence_filtered_output, all_fieldnames, confidence_filtered_rows)


if __name__ == "__main__":
    main()
