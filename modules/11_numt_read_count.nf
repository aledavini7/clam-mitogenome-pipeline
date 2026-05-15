

process index_nuclear_numt {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(numt_nuclear_bams)

    output:
    tuple val(sample_id), path('*_numt_nuc_sorted.bam.bai')

    script:
    """
    samtools index $numt_nuclear_bams
    """

}



process read_count {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(numt_nuclear_bams), path(index_bai)

    output:
    tuple val(sample_id), path('*_reads.txt')
    
    script:
    """
    touch ${sample_id}_reads.txt 

    chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

    for chromosome in "\${chromosomes[@]}"; do
        result=\$(samtools view -c -F 260 $numt_nuclear_bams "\$chromosome")
        echo "\$result" >> ${sample_id}_reads.txt 
    done
    """
}















