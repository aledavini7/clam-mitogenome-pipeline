process merge_variant_calls {

    publishDir "${params.outdir}/annotation/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(mutect2_vcf), path(mutserve_vcf), path(coverage_file)

    output:
    tuple val(sample_id), path("${sample_id}_mitochondrial_variants.tsv")

    script:
    """
    python3 ${projectDir}/bin/clam_merge_variants.py \
        --sample-id $sample_id \
        --mutect2-vcf $mutect2_vcf \
        --mutserve-vcf $mutserve_vcf \
        --coverage $coverage_file \
        --output ${sample_id}_mitochondrial_variants.tsv \
        --high-depth $params.annotation_high_depth \
        --medium-depth $params.annotation_medium_depth \
        --high-alt-reads $params.annotation_high_alt_reads \
        --medium-alt-reads $params.annotation_medium_alt_reads \
        --max-af-discordance $params.annotation_max_af_discordance
    """

}
