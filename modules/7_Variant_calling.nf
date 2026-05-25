// Sixth module: Variant_calling.nf

process mutect2 {

    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_rg_bams), path(index_bai)

    output:
    tuple val(sample_id), path('*.vcf')
    tuple val(sample_id), path('*.vcf.idx')
    tuple val(sample_id), path('*.vcf.stats')

    script:
    """
    gatk Mutect2 --input $sorted_rg_bams --reference $params.fasta --mitochondria-mode true -L $params.mt_contig --output ${sample_id}_mutect2.vcf
    """

}

process filter_mutect2 {

    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs), path(vcf_idx), path(vcf_stats)

    output:
    tuple val(sample_id), path('*_mutect2_filtered.vcf')
    tuple val(sample_id), path('*_mutect2_filtered.vcf.idx')

    script:
    """
    gatk FilterMutectCalls \
        --mitochondria-mode true \
        --reference $params.fasta \
        --variant $vcfs \
        --stats $vcf_stats \
        --output ${sample_id}_mutect2_filtered.vcf
    """

}

process bgzip {

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path('*.vcf.gz')

    script:
    """
    bgzip -k $vcfs
    """

}


process mutserve {

    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_rg_bams), path(index_bai)

    output:
    tuple val(sample_id), path('*_mutserve.vcf')
    tuple val(sample_id), path('*_mutserve.txt')

    script:
    """
    $params.mutserve call --reference $params.mutserve_fasta --contig-name $params.mt_contig --output ${sample_id}_mutserve.vcf $sorted_rg_bams
    """

}


process mutserve_annotate {

    //publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path('*_mutserve_annotated.txt')

    script:
    """
    $params.mutserve annotate --input $vcfs --annotation $params.annotation --output ${sample_id}_mutserve_annotated.txt
    """

}


// These are steps needed for WXS data, in which we sometimes want to use the vcf coming from the MT alignment without NUMTs correction steps.

process mutserve_annotate_before {

    //publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'
    publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path('*_before_mutserve_annotated.txt')

    script:
    """
    $params.mutserve annotate --input $vcfs --annotation $params.annotation --output ${sample_id}_before_mutserve_annotated.txt
    """

}

process mutserve_annotate_after {

    publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path('*_after_mutserve_annotated.txt')

    script:
    """
    $params.mutserve annotate --input $vcfs --annotation $params.annotation --output ${sample_id}_after_mutserve_annotated.txt
    """

}
