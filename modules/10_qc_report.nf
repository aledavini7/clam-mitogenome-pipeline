process bam_qc {

    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), val(bam_stages), path(bams)

    output:
    tuple val(sample_id), path("${sample_id}_bam_qc.tsv")

    script:
    """
    set -euo pipefail

    printf "sample_id\\tstage\\tbam\\ttotal_alignments\\tmapped_alignments\\tprimary_alignments\\tproperly_paired_alignments\\tduplicate_alignments\\tmt_contig_alignments\\n" > ${sample_id}_bam_qc.tsv

    IFS=',' read -ra stages <<< "$bam_stages"
    files=( $bams )

    for idx in "\${!files[@]}"; do
        bam="\${files[\$idx]}"
        stage="\${stages[\$idx]}"

        total=\$(samtools view -c "\$bam")
        mapped=\$(samtools view -c -F 4 "\$bam")
        primary=\$(samtools view -c -F 2304 "\$bam")
        paired=\$(samtools view -c -f 2 "\$bam")
        duplicates=\$(samtools view -c -f 1024 "\$bam")
        mt_contig=\$(samtools view -c "\$bam" "$params.mt_contig" 2>/dev/null || printf "%s" "\$mapped")

        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
            "$sample_id" "\$stage" "\$bam" "\$total" "\$mapped" "\$primary" "\$paired" "\$duplicates" "\$mt_contig" >> ${sample_id}_bam_qc.tsv
    done
    """

}

process no_mitomap_report_table {

    input:
    tuple val(sample_id),
        path(mutect2_tsv),
        path(mutserve_tsv),
        path(all_tsv),
        path(confidence_filtered_tsv)

    output:
    tuple val(sample_id), val(false), path("${sample_id}_mitomap_not_requested.tsv")

    script:
    """
    printf "mitomap_match\\tmitomap_match_type\\n" > ${sample_id}_mitomap_not_requested.tsv
    """

}

process no_numt_exploration_report_table {

    input:
    tuple val(sample_id),
        path(mutect2_tsv),
        path(mutserve_tsv),
        path(all_tsv),
        path(confidence_filtered_tsv)

    output:
    tuple val(sample_id),
        val(false),
        path("${sample_id}_numt_exploration_not_requested.tsv"),
        path("${sample_id}_numt_exploration_contigs_not_requested.tsv"),
        path("${sample_id}_numt_exploration_windows_not_requested.tsv")

    script:
    """
    printf "sample_id\\tpre_mitoscape_unique_read_ids\\tpost_mitoscape_unique_read_ids\\tretained_read_ids\\tmitoscape_rejected_read_ids\\tmitoscape_rejected_fraction_pct\\trejected_mt_alignments\\trejected_nuclear_alignments\\trejected_nuclear_primary_alignments\\trejected_nuclear_proper_pair_alignments\\trejected_nuclear_high_mapq_alignments\\tmin_mapq_threshold\\tcluster_window_bp\\tinterpretation\\n" > ${sample_id}_numt_exploration_not_requested.tsv
    printf "sample_id\\tnuclear_contig\\talignments\\thigh_mapq_alignments\\tmin_pos\\tmax_pos\\tmean_mapq\\n" > ${sample_id}_numt_exploration_contigs_not_requested.tsv
    printf "sample_id\\tnuclear_contig\\twindow_start\\twindow_end\\talignments\\tunique_read_ids\\thigh_mapq_alignments\\tmean_mapq\\n" > ${sample_id}_numt_exploration_windows_not_requested.tsv
    """

}

process build_report {

    publishDir "${params.outdir}/${sample_id}/report", mode: 'copy'

    input:
    tuple val(sample_id),
        path(bam_qc_tsv),
        path(coverage_before),
        path(coverage_final),
        path(mutect2_vcf),
        path(filtered_mutect2_vcf),
        path(mutserve_vcf),
        path(mutserve_txt),
        path(mutect2_haplogroups),
        path(mutserve_haplogroups),
        path(mutect2_tsv),
        path(mutserve_tsv),
        path(all_tsv),
        path(confidence_filtered_tsv),
        val(has_mitomap),
        path(mitomap_confidence_tsv),
        val(has_numt_exploration),
        path(numt_summary_tsv),
        path(numt_contig_tsv),
        path(numt_window_tsv)

    output:
    tuple val(sample_id),
        path("${sample_id}_CLAM_report.html"),
        path("${sample_id}_CLAM_qc_summary.tsv")

    script:
    """
    python3 ${projectDir}/bin/clam_build_report.py \\
        --sample-id "$sample_id" \\
        --analysis-mode "$params.analysis_mode" \\
        --input-type "$params.input_type" \\
        --reference-set "$params.reference_set" \\
        --run-numt-correction "$params.run_numt_correction" \\
        --mt-length "$params.report_mt_length" \\
        --mt-contig "$params.mt_contig" \\
        --bam-qc "$bam_qc_tsv" \\
        --coverage-before "$coverage_before" \\
        --coverage-final "$coverage_final" \\
        --mutect2-vcf "$mutect2_vcf" \\
        --filtered-mutect2-vcf "$filtered_mutect2_vcf" \\
        --mutserve-vcf "$mutserve_vcf" \\
        --mutserve-summary "$mutserve_txt" \\
        --mutect2-haplogroups "$mutect2_haplogroups" \\
        --mutserve-haplogroups "$mutserve_haplogroups" \\
        --mutect2-tsv "$mutect2_tsv" \\
        --mutserve-tsv "$mutserve_tsv" \\
        --all-tsv "$all_tsv" \\
        --confidence-filtered-tsv "$confidence_filtered_tsv" \\
        --has-mitomap "$has_mitomap" \\
        --mitomap-confidence-tsv "$mitomap_confidence_tsv" \\
        --has-numt-exploration "$has_numt_exploration" \\
        --numt-summary-tsv "$numt_summary_tsv" \\
        --numt-contig-tsv "$numt_contig_tsv" \\
        --numt-window-tsv "$numt_window_tsv" \\
        --output-html "${sample_id}_CLAM_report.html" \\
        --output-qc "${sample_id}_CLAM_qc_summary.tsv"
    """

}
