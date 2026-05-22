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

process annotation_mafs {

    publishDir "${params.outdir}/${sample_id}/annotation/mafs", mode: 'copy'

    input:
    tuple val(sample_id),
        path(mutect2_tsv),
        path(mutserve_tsv),
        path(all_tsv),
        path(confidence_filtered_tsv)

    output:
    tuple val(sample_id),
        path("${sample_id}_mutect2_variants.maf"),
        path("${sample_id}_mutserve_variants.maf"),
        path("${sample_id}_mitochondrial_variants_all.maf"),
        path("${sample_id}_mitochondrial_variants_confidence_filtered.maf")

    script:
    """
    python3 ${projectDir}/bin/clam_tsv_to_maf.py \
        --sample-id $sample_id \
        --input $mutect2_tsv \
        --output ${sample_id}_mutect2_variants.maf \
        --source-table mutect2 \
        --center $params.maf_center \
        --ncbi-build $params.maf_ncbi_build

    python3 ${projectDir}/bin/clam_tsv_to_maf.py \
        --sample-id $sample_id \
        --input $mutserve_tsv \
        --output ${sample_id}_mutserve_variants.maf \
        --source-table mutserve \
        --center $params.maf_center \
        --ncbi-build $params.maf_ncbi_build

    python3 ${projectDir}/bin/clam_tsv_to_maf.py \
        --sample-id $sample_id \
        --input $all_tsv \
        --output ${sample_id}_mitochondrial_variants_all.maf \
        --source-table consensus_all \
        --center $params.maf_center \
        --ncbi-build $params.maf_ncbi_build

    python3 ${projectDir}/bin/clam_tsv_to_maf.py \
        --sample-id $sample_id \
        --input $confidence_filtered_tsv \
        --output ${sample_id}_mitochondrial_variants_confidence_filtered.maf \
        --source-table consensus_confidence_filtered \
        --center $params.maf_center \
        --ncbi-build $params.maf_ncbi_build
    """

}

process mitomap_annotation {

    publishDir "${params.outdir}/${sample_id}/annotation/mitomap", mode: 'copy'

    input:
    tuple val(sample_id),
        path(mutect2_tsv),
        path(mutserve_tsv),
        path(all_tsv),
        path(confidence_filtered_tsv),
        path(mitomap_table)

    output:
    tuple val(sample_id),
        path("${sample_id}_mutect2_variants_mitomap.tsv"),
        path("${sample_id}_mutserve_variants_mitomap.tsv"),
        path("${sample_id}_mitochondrial_variants_all_mitomap.tsv"),
        path("${sample_id}_mitochondrial_variants_confidence_filtered_mitomap.tsv")

    script:
    """
    python3 ${projectDir}/bin/clam_annotate_mitomap.py \
        --mitomap $mitomap_table \
        --input $mutect2_tsv \
        --output ${sample_id}_mutect2_variants_mitomap.tsv \
        --keep-columns $params.mitomap_keep_columns

    python3 ${projectDir}/bin/clam_annotate_mitomap.py \
        --mitomap $mitomap_table \
        --input $mutserve_tsv \
        --output ${sample_id}_mutserve_variants_mitomap.tsv \
        --keep-columns $params.mitomap_keep_columns

    python3 ${projectDir}/bin/clam_annotate_mitomap.py \
        --mitomap $mitomap_table \
        --input $all_tsv \
        --output ${sample_id}_mitochondrial_variants_all_mitomap.tsv \
        --keep-columns $params.mitomap_keep_columns

    python3 ${projectDir}/bin/clam_annotate_mitomap.py \
        --mitomap $mitomap_table \
        --input $confidence_filtered_tsv \
        --output ${sample_id}_mitochondrial_variants_confidence_filtered_mitomap.tsv \
        --keep-columns $params.mitomap_keep_columns
    """

}
