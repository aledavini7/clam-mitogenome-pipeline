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

            for alt_idx, alt in enumerate(alts):
                if alt in ("", ".", "*") or alt.startswith("<"):
                    continue

                ref_depth = as_int(pick(ad_values, 0))
                alt_depth = as_int(pick(ad_values, alt_idx + 1))
                af = as_float(pick(af_values, alt_idx))

                if af is None and dp and alt_depth is not None and dp > 0:
                    af = alt_depth / dp

                if dp is None and ad_values:
                    numeric_ad = [as_int(value) for value in ad_values]
                    numeric_ad = [value for value in numeric_ad if value is not None]
                    dp = sum(numeric_ad) if numeric_ad else None

                key = (chrom, int(pos_text), ref, alt)
                calls[key] = {
                    "caller": caller,
                    "filter": filt,
                    "qual": qual,
                    "dp": dp,
                    "ref_depth": ref_depth,
                    "alt_depth": alt_depth,
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


def confidence_tier(callers, mutect2, coverage_depth, alt_depth, af_values, args):
    reasons = []

    if coverage_depth is None:
        reasons.append("missing_coverage")
    elif coverage_depth < args.medium_depth:
        reasons.append(f"coverage<{args.medium_depth}")

    if alt_depth is None:
        reasons.append("missing_alt_depth")
    elif alt_depth < args.medium_alt_reads:
        reasons.append(f"alt_reads<{args.medium_alt_reads}")

    if mutect2 and not is_pass(mutect2["filter"]):
        reasons.append(f"mutect2_filter={mutect2['filter']}")

    if len(callers) == 2 and len(af_values) == 2:
        discordance = abs(af_values[0] - af_values[1])
        if discordance > args.max_af_discordance:
            reasons.append(f"caller_af_discordance>{args.max_af_discordance}")
    else:
        discordance = None

    if (
        len(callers) == 2
        and mutect2
        and is_pass(mutect2["filter"])
        and coverage_depth is not None
        and coverage_depth >= args.high_depth
        and alt_depth is not None
        and alt_depth >= args.high_alt_reads
        and (discordance is None or discordance <= args.max_af_discordance)
    ):
        return "HIGH", "both_callers_pass_depth_alt_support"

    if (
        coverage_depth is not None
        and coverage_depth >= args.medium_depth
        and alt_depth is not None
        and alt_depth >= args.medium_alt_reads
        and (not mutect2 or is_pass(mutect2["filter"]))
    ):
        return "MEDIUM", "single_or_partial_support_with_acceptable_depth"

    return "LOW", ",".join(reasons) if reasons else "limited_support"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--mutect2-vcf", required=True)
    parser.add_argument("--mutserve-vcf", required=True)
    parser.add_argument("--coverage", required=True)
    parser.add_argument("--output", required=True)
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

    fieldnames = [
        "sample_id",
        "chrom",
        "position",
        "ref",
        "alt",
        "callers",
        "coverage_depth",
        "consensus_af",
        "heteroplasmy_percent",
        "consensus_alt_depth",
        "heteroplasmy_ci95_low",
        "heteroplasmy_ci95_high",
        "confidence_tier",
        "confidence_reason",
        "mutect2_filter",
        "mutect2_af",
        "mutect2_dp",
        "mutect2_ref_depth",
        "mutect2_alt_depth",
        "mutserve_filter",
        "mutserve_af",
        "mutserve_dp",
        "mutserve_ref_depth",
        "mutserve_alt_depth",
    ]

    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for key in sorted(caller_calls, key=lambda item: (item[0], item[1], item[2], item[3])):
            chrom, pos, ref, alt = key
            calls = caller_calls[key]
            mutect2 = calls.get("mutect2")
            mutserve = calls.get("mutserve")
            caller_names = sorted(calls)
            af_values = [calls[name]["af"] for name in caller_names if calls[name]["af"] is not None]

            coverage_depth = coverage.get((chrom, pos))
            if coverage_depth is None:
                coverage_depth = next((calls[name]["dp"] for name in caller_names if calls[name]["dp"] is not None), None)

            consensus_af = sum(af_values) / len(af_values) if af_values else None

            alt_depth_values = [
                calls[name]["alt_depth"] for name in caller_names if calls[name]["alt_depth"] is not None
            ]
            if alt_depth_values:
                consensus_alt_depth = round(sum(alt_depth_values) / len(alt_depth_values))
            elif consensus_af is not None and coverage_depth is not None:
                consensus_alt_depth = round(consensus_af * coverage_depth)
            else:
                consensus_alt_depth = None

            ci_low, ci_high = wilson_interval(consensus_alt_depth, coverage_depth)
            tier, reason = confidence_tier(
                caller_names,
                mutect2,
                coverage_depth,
                consensus_alt_depth,
                af_values,
                args,
            )

            writer.writerow(
                {
                    "sample_id": args.sample_id,
                    "chrom": chrom,
                    "position": pos,
                    "ref": ref,
                    "alt": alt,
                    "callers": ",".join(caller_names),
                    "coverage_depth": coverage_depth if coverage_depth is not None else "",
                    "consensus_af": fmt_float(consensus_af),
                    "heteroplasmy_percent": fmt_percent(consensus_af),
                    "consensus_alt_depth": consensus_alt_depth if consensus_alt_depth is not None else "",
                    "heteroplasmy_ci95_low": fmt_percent(ci_low),
                    "heteroplasmy_ci95_high": fmt_percent(ci_high),
                    "confidence_tier": tier,
                    "confidence_reason": reason,
                    "mutect2_filter": mutect2["filter"] if mutect2 else "",
                    "mutect2_af": fmt_float(mutect2["af"]) if mutect2 else "",
                    "mutect2_dp": mutect2["dp"] if mutect2 and mutect2["dp"] is not None else "",
                    "mutect2_ref_depth": mutect2["ref_depth"] if mutect2 and mutect2["ref_depth"] is not None else "",
                    "mutect2_alt_depth": mutect2["alt_depth"] if mutect2 and mutect2["alt_depth"] is not None else "",
                    "mutserve_filter": mutserve["filter"] if mutserve else "",
                    "mutserve_af": fmt_float(mutserve["af"]) if mutserve else "",
                    "mutserve_dp": mutserve["dp"] if mutserve and mutserve["dp"] is not None else "",
                    "mutserve_ref_depth": mutserve["ref_depth"] if mutserve and mutserve["ref_depth"] is not None else "",
                    "mutserve_alt_depth": mutserve["alt_depth"] if mutserve and mutserve["alt_depth"] is not None else "",
                }
            )


if __name__ == "__main__":
    main()
