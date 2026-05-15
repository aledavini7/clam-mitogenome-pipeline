// Fourth module: Coverage.nf

process sorting {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(md_bams) // output of 'get_the_MD'
    tuple val(sample_id), path(mito_bams) // output of 'mitoscape'

    output:
    tuple val(sample_id), path('*_sorted_before.bam') // sorted before mitoscape
    tuple val(sample_id), path('*_sorted_after.bam') // sorted after mitoscape

    script:
    """
    samtools sort $md_bams > ${sample_id}_sorted_before.bam
    samtools sort $mito_bams > ${sample_id}_sorted_after.bam
    """

}

process coverage {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(before_bams)
    tuple val(sample_id), path(after_bams)

    output:
    tuple val(sample_id), path('*_before.coverage') // coverage before mitoscape
    tuple val(sample_id), path('*_after.coverage') // coverage after mitoscape
   

    script:
    """
    samtools depth $before_bams > ${sample_id}_before.coverage
    samtools depth $after_bams > ${sample_id}_after.coverage
    """

}


// Sorting and Coverage single file (if I have mitochondrial sequencing data, I have only information about mtDNA, so I perform just single steps)

process sorting_single {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bams) // output of 'Alignment_MT'

    output:
    tuple val(sample_id), path('*_sorted.bam') // sorted

    script:
    """
    samtools sort $bams > ${sample_id}_sorted.bam
    """

}

process coverage_single {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bams)

    output:
    tuple val(sample_id), path('*.coverage') // coverage 
   

    script:
    """
    samtools depth $sorted_bams > ${sample_id}.coverage
    """

}