#!/usr/bin/env bash
set -euo pipefail

# Run from the Nextflow launch directory after copying this script and
# clam_build_nuclear_report.py into that directory.
LAUNCHDIR="$(pwd -P)"
SAMPLE_ID="${1:-S69834_S6_100K_1_S17}"
CALLING_MODE="${2:-tumor_only}"
COMMON_AF_THRESHOLD="${COMMON_AF_THRESHOLD:-0.01}"
TOP_VARIANTS="${TOP_VARIANTS:-100}"
SAMPLE_DIR="${LAUNCHDIR}/results/${SAMPLE_ID}/nuclear"
ANNOTATION_DIR="${SAMPLE_DIR}/annotation"
CNV_DIR="${SAMPLE_DIR}/cnv/cnvkit_out"
REPORT_DIR="${SAMPLE_DIR}/report"
REPORT_BUILDER="$(dirname "$(realpath "$0")")/clam_build_nuclear_report.py"
VCF="${ANNOTATION_DIR}/${SAMPLE_ID}.annotated.vcf"
MAF="${ANNOTATION_DIR}/${SAMPLE_ID}.maf"
OUTPUT="${REPORT_DIR}/${SAMPLE_ID}_nuclear_report_preview.html"

for path in "$REPORT_BUILDER" "$VCF" "$MAF" "$CNV_DIR"; do
    if [[ ! -e "$path" ]]; then
        echo "Required input not found: $path" >&2
        echo "Expected report builder and outputs beneath: $LAUNCHDIR" >&2
        exit 1
    fi
done

mkdir -p "$REPORT_DIR"
python3 "$REPORT_BUILDER" \
    --sample-id "$SAMPLE_ID" \
    --calling-mode "$CALLING_MODE" \
    --vcf "$VCF" \
    --maf "$MAF" \
    --cnv-dir "$CNV_DIR" \
    --common-af-threshold "$COMMON_AF_THRESHOLD" \
    --top-variants "$TOP_VARIANTS" \
    --output "$OUTPUT"

echo "Preview report written to: $OUTPUT"
echo "Open it in a browser from the cluster or copy it back to your computer."
