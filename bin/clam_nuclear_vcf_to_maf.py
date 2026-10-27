#!/usr/bin/env python3
"""Convert a VEP-annotated VCF into a compact, MAF-compatible TSV."""

import argparse
import csv
import re
import sys


CONSEQUENCE_CLASS = {
    "transcript_ablation": "Splice_Site",
    "splice_acceptor_variant": "Splice_Site",
    "splice_donor_variant": "Splice_Site",
    "stop_gained": "Nonsense_Mutation",
    "frameshift_variant": "Frame_Shift_Del",
    "stop_lost": "Nonstop_Mutation",
    "start_lost": "Translation_Start_Site",
    "inframe_insertion": "In_Frame_Ins",
    "inframe_deletion": "In_Frame_Del",
    "missense_variant": "Missense_Mutation",
    "protein_altering_variant": "Missense_Mutation",
    "splice_region_variant": "Splice_Region",
    "synonymous_variant": "Silent",
    "5_prime_UTR_variant": "5'UTR",
    "3_prime_UTR_variant": "3'UTR",
    "non_coding_transcript_exon_variant": "RNA",
    "intron_variant": "Intron",
    "upstream_gene_variant": "5'Flank",
    "downstream_gene_variant": "3'Flank",
    "intergenic_variant": "IGR",
}
SEVERITY = list(CONSEQUENCE_CLASS)


def info_values(text):
    result = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def parse_sample(format_text, sample_text):
    keys = format_text.split(":")
    values = sample_text.split(":")
    return dict(zip(keys, values))


def choose_consequence(entries, alt):
    candidates = []
    for entry in entries:
        fields = dict(zip(entry[0], entry[1].split("|")))
        allele = fields.get("Allele", "")
        if allele and allele != alt and not (alt.startswith("<") and allele == alt):
            continue
        terms = fields.get("Consequence", "").split("&")
        rank = min((SEVERITY.index(term) for term in terms if term in SEVERITY), default=len(SEVERITY))
        canonical = fields.get("MANE_SELECT") or fields.get("CANONICAL") == "YES"
        candidates.append((not bool(canonical), rank, fields))
    return sorted(candidates, key=lambda item: (item[0], item[1]))[0][2] if candidates else {}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    columns = [
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
        "Start_Position", "End_Position", "Strand", "Variant_Classification",
        "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
        "Tumor_Sample_Barcode", "HGVSc", "HGVSp", "Transcript_ID", "Gene",
        "gnomAD_AF", "COSMIC_ID", "FILTER", "t_depth", "t_alt_count", "t_ref_count",
    ]
    csq_keys = []
    sample_column = None
    with open(args.input, encoding="utf-8") as source, open(args.output, "w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line in source:
            if line.startswith("##INFO=<ID=CSQ,"):
                match = re.search(r'Format: ([^\"]+)', line)
                if match:
                    csq_keys = match.group(1).split("|")
                continue
            if line.startswith("#CHROM\t"):
                sample_names = line.rstrip("\n").split("\t")[9:]
                if args.sample_id in sample_names:
                    sample_column = 9 + sample_names.index(args.sample_id)
                elif sample_names:
                    sample_column = 9
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt, _, filt, info = fields[:8]
            if filt not in ("PASS", "."):
                continue
            info_map = info_values(info)
            csq_entries = [(csq_keys, value) for value in info_map.get("CSQ", "").split(",") if value]
            consequence = choose_consequence(csq_entries, alt)
            terms = consequence.get("Consequence", "").split("&")
            classification = next((CONSEQUENCE_CLASS[t] for t in SEVERITY if t in terms), "Targeted_Region")
            sample_data = parse_sample(fields[8], fields[sample_column]) if sample_column is not None and len(fields) > sample_column else {}
            ad = sample_data.get("AD", "").split(",")
            ref_count = ad[0] if ad else ""
            alt_count = ad[1] if len(ad) > 1 else ""
            depth = sample_data.get("DP", "")
            if not depth and ref_count.isdigit() and alt_count.isdigit():
                depth = str(int(ref_count) + int(alt_count))
            variant_type = "SNP" if len(ref) == len(alt) == 1 else ("INS" if len(alt) > len(ref) else "DEL")
            row = {
                "Hugo_Symbol": consequence.get("SYMBOL", ""), "Entrez_Gene_Id": "",
                "Center": "CLAM", "NCBI_Build": "GRCh38", "Chromosome": chrom.removeprefix("chr"),
                "Start_Position": pos, "End_Position": str(int(pos) + len(ref) - 1), "Strand": "+",
                "Variant_Classification": classification, "Variant_Type": variant_type,
                "Reference_Allele": ref, "Tumor_Seq_Allele1": ref, "Tumor_Seq_Allele2": alt,
                "Tumor_Sample_Barcode": args.sample_id, "HGVSc": consequence.get("HGVSc", ""),
                "HGVSp": consequence.get("HGVSp", ""), "Transcript_ID": consequence.get("Feature", ""),
                "Gene": consequence.get("Gene", ""), "gnomAD_AF": consequence.get("gnomAD_AF", ""),
                "COSMIC_ID": consequence.get("COSMIC", ""), "FILTER": filt,
                "t_depth": depth, "t_alt_count": alt_count, "t_ref_count": ref_count,
            }
            writer.writerow(row)
    if not csq_keys:
        print("Warning: VCF has no VEP CSQ header; MAF will contain basic variant columns only", file=sys.stderr)


if __name__ == "__main__":
    main()
