#!/usr/bin/env python3
import argparse
import csv


MAF_FIELDS = [
    "Hugo_Symbol",
    "Entrez_Gene_Id",
    "Center",
    "NCBI_Build",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "CLAM_Source_Table",
    "CLAM_Callers",
    "CLAM_Filter",
    "CLAM_Confidence",
    "CLAM_Confidence_Reason",
    "CLAM_AF_Fraction",
    "CLAM_Heteroplasmy_Pct_0_100",
    "CLAM_Coverage_Depth",
]


def as_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def normalize_chrom(chrom):
    if chrom in ("NC_012920.1", "chrM", "M"):
        return "MT"
    return chrom


def trim_variant(pos, ref, alt):
    start = int(pos)
    ref = ref or ""
    alt = alt or ""

    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        start += 1

    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    maf_ref = ref if ref else "-"
    maf_alt = alt if alt else "-"

    if maf_ref == "-":
        end = start
    else:
        end = start + len(ref) - 1

    return start, end, maf_ref, maf_alt


def variant_type(ref, alt):
    if ref == "-":
        return "INS"
    if alt == "-":
        return "DEL"
    if len(ref) == len(alt):
        if len(ref) == 1:
            return "SNP"
        if len(ref) == 2:
            return "DNP"
        if len(ref) == 3:
            return "TNP"
        return "ONP"
    if len(ref) < len(alt):
        return "INS"
    if len(ref) > len(alt):
        return "DEL"
    return "ONP"


def first_present(row, names):
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return value
    return ""


def row_evidence(row, source_table):
    if "caller" in row:
        caller = row.get("caller", source_table)
        return {
            "callers": caller,
            "filter": row.get("filter", ""),
            "confidence": row.get("confidence_tier", ""),
            "confidence_reason": row.get("confidence_reason", ""),
            "af": row.get("af_fraction", ""),
            "heteroplasmy": row.get("heteroplasmy_pct_0_100", ""),
            "depth": row.get("depth", ""),
            "ref_count": row.get("ref_depth", ""),
            "alt_count": row.get("alt_depth", ""),
            "coverage": row.get("coverage_depth", ""),
        }

    depth = row.get("consensus_depth", "")
    alt_count = row.get("consensus_alt_depth", "")
    ref_count = ""
    depth_int = as_int(depth)
    alt_int = as_int(alt_count)
    if depth_int is not None and alt_int is not None:
        ref_count = max(depth_int - alt_int, 0)

    return {
        "callers": row.get("callers", source_table),
        "filter": first_present(row, ["mutect2_filter", "mutserve_filter"]),
        "confidence": row.get("confidence_tier", ""),
        "confidence_reason": row.get("confidence_reason", ""),
        "af": row.get("consensus_af_fraction", ""),
        "heteroplasmy": row.get("consensus_heteroplasmy_pct_0_100", ""),
        "depth": depth,
        "ref_count": ref_count,
        "alt_count": alt_count,
        "coverage": row.get("coverage_depth", ""),
    }


def to_maf_row(row, args):
    start, end, ref, alt = trim_variant(row["position"], row["ref"], row["alt"])
    evidence = row_evidence(row, args.source_table)

    return {
        "Hugo_Symbol": args.hugo_symbol,
        "Entrez_Gene_Id": "0",
        "Center": args.center,
        "NCBI_Build": args.ncbi_build,
        "Chromosome": normalize_chrom(row["chrom"]),
        "Start_Position": start,
        "End_Position": end,
        "Strand": "+",
        "Variant_Classification": args.variant_classification,
        "Variant_Type": variant_type(ref, alt),
        "Reference_Allele": ref,
        "Tumor_Seq_Allele1": ref,
        "Tumor_Seq_Allele2": alt,
        "Tumor_Sample_Barcode": row.get("sample_id") or args.sample_id,
        "t_depth": evidence["depth"],
        "t_ref_count": evidence["ref_count"],
        "t_alt_count": evidence["alt_count"],
        "CLAM_Source_Table": args.source_table,
        "CLAM_Callers": evidence["callers"],
        "CLAM_Filter": evidence["filter"],
        "CLAM_Confidence": evidence["confidence"],
        "CLAM_Confidence_Reason": evidence["confidence_reason"],
        "CLAM_AF_Fraction": evidence["af"],
        "CLAM_Heteroplasmy_Pct_0_100": evidence["heteroplasmy"],
        "CLAM_Coverage_Depth": evidence["coverage"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--source-table", required=True)
    parser.add_argument("--center", default="CLAM")
    parser.add_argument("--ncbi-build", default="rCRS")
    parser.add_argument("--hugo-symbol", default="MTDNA")
    parser.add_argument("--variant-classification", default="RNA")
    args = parser.parse_args()

    with open(args.input, "r", encoding="utf-8", newline="") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        rows = [to_maf_row(row, args) for row in reader]

    with open(args.output, "w", encoding="utf-8", newline="") as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=MAF_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
