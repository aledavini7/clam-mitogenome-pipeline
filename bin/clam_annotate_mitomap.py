#!/usr/bin/env python3
import argparse
import csv
import re
from collections import defaultdict


POSITION_ALIASES = {
    "pos",
    "position",
    "nt_position",
    "nucleotide_position",
    "mtdna_position",
    "rcrs_position",
    "locus",
}

REF_ALIASES = {"ref", "reference", "reference_allele", "ref_allele", "rcrs_base"}
ALT_ALIASES = {
    "alt",
    "allele",
    "variant_allele",
    "alternate",
    "alternate_allele",
    "alt_allele",
    "mutant",
    "mutant_allele",
}
VARIANT_ALIASES = {
    "variant",
    "mutation",
    "polymorphism",
    "nt_change",
    "nucleotide_change",
    "base_change",
    "substitution",
}


def normalize_header(value):
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def detect_dialect(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        sample = handle.read(8192)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;")
    except csv.Error:
        class TabDialect(csv.excel_tab):
            delimiter = "\t"

        return TabDialect


def as_int(value):
    if value in (None, ""):
        return None
    match = re.search(r"\d+", str(value).replace(",", ""))
    if not match:
        return None
    return int(match.group(0))


def first_present(row, aliases):
    for key, value in row.items():
        if normalize_header(key) in aliases and value not in (None, ""):
            return value
    return None


def parse_variant_text(value):
    if value in (None, ""):
        return None, None, None

    text = str(value).strip()

    match = re.search(r"m\.?\s*(\d+)\s*([ACGT]+)\s*>\s*([ACGT]+)", text, re.IGNORECASE)
    if match:
        return int(match.group(1)), match.group(2).upper(), match.group(3).upper()

    match = re.search(r"\b([ACGT]+)\s*(\d+)\s*([ACGT]+)\b", text, re.IGNORECASE)
    if match:
        return int(match.group(2)), match.group(1).upper(), match.group(3).upper()

    match = re.search(r"\b(\d+)\s*([ACGT]+)\b", text, re.IGNORECASE)
    if match:
        return int(match.group(1)), None, match.group(2).upper()

    match = re.search(r"\b(\d+)\b", text)
    if match:
        return int(match.group(1)), None, None

    return None, None, None


def collapse_values(rows, columns):
    values = {}
    for column in columns:
        seen = []
        for row in rows:
            value = row.get(column, "")
            value = str(value).strip()
            if value and value not in seen:
                seen.append(value)
        values[column] = "||".join(seen)
    return values


def choose_columns(fieldnames, limit):
    selected = []
    for field in fieldnames:
        normalized = normalize_header(field)
        if normalized in POSITION_ALIASES | REF_ALIASES | ALT_ALIASES | VARIANT_ALIASES:
            continue
        selected.append(field)
        if len(selected) >= limit:
            break
    return selected


def mitomap_record_keys(row):
    pos = as_int(first_present(row, POSITION_ALIASES))
    ref = first_present(row, REF_ALIASES)
    alt = first_present(row, ALT_ALIASES)

    if ref:
        ref = str(ref).strip().upper()
    if alt:
        alt = str(alt).strip().upper()

    variant_text = first_present(row, VARIANT_ALIASES)
    parsed_pos, parsed_ref, parsed_alt = parse_variant_text(variant_text)
    pos = pos or parsed_pos
    ref = ref or parsed_ref
    alt = alt or parsed_alt

    if pos is None:
        return None, None

    position_key = pos
    allele_key = None
    if ref and alt and re.fullmatch(r"[ACGT]+|-", ref) and re.fullmatch(r"[ACGT]+|-", alt):
        allele_key = (pos, ref, alt)

    return position_key, allele_key


def read_mitomap(path, keep_columns):
    dialect = detect_dialect(path)
    by_position = defaultdict(list)
    by_allele = defaultdict(list)

    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            return [], by_position, by_allele

        columns = choose_columns(reader.fieldnames, keep_columns)
        for row in reader:
            position_key, allele_key = mitomap_record_keys(row)
            if position_key is None:
                continue
            by_position[position_key].append(row)
            if allele_key is not None:
                by_allele[allele_key].append(row)

    return columns, by_position, by_allele


def annotate_row(row, mitomap_columns, by_position, by_allele):
    pos = as_int(row.get("position"))
    ref = str(row.get("ref", "")).upper()
    alt = str(row.get("alt", "")).upper()
    exact_rows = by_allele.get((pos, ref, alt), []) if pos is not None else []
    position_rows = by_position.get(pos, []) if pos is not None else []

    if exact_rows:
        match_rows = exact_rows
        match_type = "exact"
    elif position_rows:
        match_rows = position_rows
        match_type = "position"
    else:
        match_rows = []
        match_type = ""

    annotated = dict(row)
    annotated["mitomap_match"] = "YES" if match_rows else "NO"
    annotated["mitomap_match_type"] = match_type
    annotated["mitomap_match_count"] = len(match_rows)

    collapsed = collapse_values(match_rows, mitomap_columns) if match_rows else {}
    for column in mitomap_columns:
        annotated[f"mitomap_{normalize_header(column)}"] = collapsed.get(column, "")

    return annotated


def annotate_table(input_path, output_path, mitomap_columns, by_position, by_allele):
    with open(input_path, "r", encoding="utf-8", newline="") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        rows = [annotate_row(row, mitomap_columns, by_position, by_allele) for row in reader]
        fieldnames = list(reader.fieldnames or [])

    extra_fields = ["mitomap_match", "mitomap_match_type", "mitomap_match_count"]
    extra_fields.extend(f"mitomap_{normalize_header(column)}" for column in mitomap_columns)

    with open(output_path, "w", encoding="utf-8", newline="") as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=fieldnames + extra_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mitomap", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--keep-columns", type=int, default=12)
    args = parser.parse_args()

    mitomap_columns, by_position, by_allele = read_mitomap(args.mitomap, args.keep_columns)
    annotate_table(args.input, args.output, mitomap_columns, by_position, by_allele)


if __name__ == "__main__":
    main()
