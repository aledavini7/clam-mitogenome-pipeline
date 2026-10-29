#!/usr/bin/env bash
set -euo pipefail

# Rebuild the mitochondrial HTML report from an existing CLAM result folder.
# No alignment or variant calling is performed.
RESULT_DIR="$(cd "${1:-results_optimal_test/test_with_NUMT_correction/S69834_S6_100K_1_S17}" && pwd -P)"
SAMPLE_ID="${2:-$(basename "$RESULT_DIR")}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPORT_BUILDER="$SCRIPT_DIR/clam_build_report.py"
ANNOTATION="$RESULT_DIR/annotation"
MITOMAP="$ANNOTATION/mitomap"
VARIANTS="$RESULT_DIR/variants"
OUT="$RESULT_DIR/report"

required=(
  "$REPORT_BUILDER"
  "$RESULT_DIR/qc/${SAMPLE_ID}_bam_qc.tsv"
  "$RESULT_DIR/coverage/${SAMPLE_ID}_before.coverage"
  "$RESULT_DIR/coverage/${SAMPLE_ID}_after.coverage"
  "$VARIANTS/${SAMPLE_ID}_mutect2.vcf"
  "$VARIANTS/${SAMPLE_ID}_mutect2_filtered.vcf"
  "$VARIANTS/${SAMPLE_ID}_mutserve.vcf"
  "$VARIANTS/${SAMPLE_ID}_mutserve.txt"
  "$RESULT_DIR/haplogroups/${SAMPLE_ID}_mutect2_haplogroups.txt"
  "$RESULT_DIR/haplogroups/${SAMPLE_ID}_mutserve_haplogroups.txt"
  "$ANNOTATION/${SAMPLE_ID}_mutect2_variants.tsv"
  "$ANNOTATION/${SAMPLE_ID}_mutserve_variants.tsv"
  "$ANNOTATION/${SAMPLE_ID}_mitochondrial_variants_all.tsv"
  "$ANNOTATION/${SAMPLE_ID}_mitochondrial_variants_confidence_filtered.tsv"
  "$MITOMAP/${SAMPLE_ID}_mitochondrial_variants_all_mitomap.tsv"
  "$MITOMAP/${SAMPLE_ID}_mitochondrial_variants_confidence_filtered_mitomap.tsv"
)
for file in "${required[@]}"; do
  [[ -f "$file" ]] || { echo "Required report input missing: $file" >&2; exit 1; }
done

mkdir -p "$OUT"
python3 "$REPORT_BUILDER" \
  --sample-id "$SAMPLE_ID" \
  --analysis-mode wgs_numt_correction \
  --input-type fastq \
  --reference-set GRCh38+rCRS \
  --run-numt-correction true \
  --mt-length 16569 \
  --mt-contig NC_012920.1 \
  --bam-qc "$RESULT_DIR/qc/${SAMPLE_ID}_bam_qc.tsv" \
  --coverage-before "$RESULT_DIR/coverage/${SAMPLE_ID}_before.coverage" \
  --coverage-final "$RESULT_DIR/coverage/${SAMPLE_ID}_after.coverage" \
  --mutect2-vcf "$VARIANTS/${SAMPLE_ID}_mutect2.vcf" \
  --filtered-mutect2-vcf "$VARIANTS/${SAMPLE_ID}_mutect2_filtered.vcf" \
  --mutserve-vcf "$VARIANTS/${SAMPLE_ID}_mutserve.vcf" \
  --mutserve-summary "$VARIANTS/${SAMPLE_ID}_mutserve.txt" \
  --mutect2-haplogroups "$RESULT_DIR/haplogroups/${SAMPLE_ID}_mutect2_haplogroups.txt" \
  --mutserve-haplogroups "$RESULT_DIR/haplogroups/${SAMPLE_ID}_mutserve_haplogroups.txt" \
  --mutect2-tsv "$ANNOTATION/${SAMPLE_ID}_mutect2_variants.tsv" \
  --mutserve-tsv "$ANNOTATION/${SAMPLE_ID}_mutserve_variants.tsv" \
  --all-tsv "$ANNOTATION/${SAMPLE_ID}_mitochondrial_variants_all.tsv" \
  --confidence-filtered-tsv "$ANNOTATION/${SAMPLE_ID}_mitochondrial_variants_confidence_filtered.tsv" \
  --has-mitomap true \
  --mitomap-all-tsv "$MITOMAP/${SAMPLE_ID}_mitochondrial_variants_all_mitomap.tsv" \
  --mitomap-confidence-tsv "$MITOMAP/${SAMPLE_ID}_mitochondrial_variants_confidence_filtered_mitomap.tsv" \
  --has-numt-exploration false \
  --numt-summary-tsv "$ANNOTATION/${SAMPLE_ID}_mitochondrial_variants_all.tsv" \
  --numt-contig-tsv "$ANNOTATION/${SAMPLE_ID}_mitochondrial_variants_all.tsv" \
  --numt-window-tsv "$ANNOTATION/${SAMPLE_ID}_mitochondrial_variants_all.tsv" \
  --output-html "$OUT/${SAMPLE_ID}_CLAM_report_preview.html" \
  --output-qc "$OUT/${SAMPLE_ID}_CLAM_qc_summary_preview.tsv"

echo "Report-only preview written: $OUT/${SAMPLE_ID}_CLAM_report_preview.html"
