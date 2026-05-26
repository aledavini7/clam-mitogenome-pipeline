

process index_0 {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bams)

    output:
    tuple val(sample_id), path('*.bam.bai')

    script:
    """
    samtools index $sorted_bams
    """

}

process index_1 {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bams)

    output:
    tuple val(sample_id), path('*.bam.bai')
    //tuple val(sample_id), path('*_sorted.bam.bai')

    script:
    """
    samtools index $sorted_bams
    """

}


process mutect2_0 {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_rg_bams), path(index_bai)

    output:
    tuple val(sample_id), path('*.vcf')
    tuple val(sample_id), path('*.vcf.idx')
    tuple val(sample_id), path('*.vcf.stats')

    script:
    """
    gatk Mutect2 --input $sorted_rg_bams --reference $params.fasta --mitochondria-mode true -L $params.mt_contig --output ${sample_id}_before_mutect2.vcf
    """

}

process mutect2_1 {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_rg_bams), path(index_bai)

    output:
    tuple val(sample_id), path('*.vcf')
    tuple val(sample_id), path('*.vcf.idx')
    tuple val(sample_id), path('*.vcf.stats')

    script:
    """
    gatk Mutect2 --input $sorted_rg_bams --reference $params.fasta --mitochondria-mode true -L $params.mt_contig --output ${sample_id}_after_mutect2.vcf
    """

}


process mutserve_0 {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_rg_bams), path(index_bai)

    output:
    tuple val(sample_id), path('*_mutserve.vcf')
    tuple val(sample_id), path('*_mutserve.txt')

    script:
    """
    $params.mutserve call --reference $params.mutserve_fasta --contig-name $params.mt_contig --output ${sample_id}_before_mutserve.vcf $sorted_rg_bams
    """

}


process mutserve_1 {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_rg_bams), path(index_bai)

    output:
    tuple val(sample_id), path('*_mutserve.vcf')
    tuple val(sample_id), path('*_mutserve.txt')

    script:
    """
    $params.mutserve call --reference $params.mutserve_fasta --contig-name $params.mt_contig --output ${sample_id}_after_mutserve.vcf $sorted_rg_bams
    """

}










