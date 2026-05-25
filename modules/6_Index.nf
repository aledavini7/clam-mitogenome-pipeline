// Fifth module: Index.nf

process index {

    publishDir "${params.outdir}/${sample_id}/bams", mode: 'copy'

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
