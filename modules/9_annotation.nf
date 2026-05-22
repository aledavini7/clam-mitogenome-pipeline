process merge_variant_calls {

    publishDir "${params.outdir}/${sample_id}/annotation", mode: 'copy'

    input:
    tuple val(sample_id), path(mutect2_vcf), path(mutserve_vcf), path(coverage_file)

    output:
    tuple val(sample_id),
        path("${sample_id}_mutect2_variants.tsv"),
        path("${sample_id}_mutserve_variants.tsv"),
        path("${sample_id}_mitochondrial_variants_all.tsv"),
        path("${sample_id}_mitochondrial_variants_confidence_filtered.tsv")

    script:
    """
    python3 ${projectDir}/bin/clam_merge_variants.py \
        --sample-id $sample_id \
        --mutect2-vcf $mutect2_vcf \
        --mutserve-vcf $mutserve_vcf \
        --coverage $coverage_file \
        --output ${sample_id}_mitochondrial_variants_all.tsv \
        --mutect2-output ${sample_id}_mutect2_variants.tsv \
        --mutserve-output ${sample_id}_mutserve_variants.tsv \
        --confidence-filtered-output ${sample_id}_mitochondrial_variants_confidence_filtered.tsv \
        --high-depth $params.annotation_high_depth \
        --medium-depth $params.annotation_medium_depth \
        --high-alt-reads $params.annotation_high_alt_reads \
        --medium-alt-reads $params.annotation_medium_alt_reads \
        --max-af-discordance $params.annotation_max_af_discordance
    """

}
