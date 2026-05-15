// Prepare data to obtain the reads which will be used for the pipeline.
// This step is needed when we start from a bam (such as for public data obtained with dbGaP or TCGA). In this case we need to sort the bam and convert it in fastq (R1, R2)
// We start from a bam that have benn previously aligned to the nuclear human reference genome. 
// Since we want to re-align it to the mitochodnrial genome only, we have to restart from the raw fastq file and realing it.

process sorting_1 {

    input:
    tuple val(sample_id), path(bams)

    output:
    tuple val(sample_id), path('*_sorted.bam')

    script:
    """
    samtools sort -n $bams > ${sample_id}_sorted.bam
    """

}


process bam2fqgz {

	input:
	tuple val(sample_id), path(bams)

	output:
	tuple val(sample_id), path('*.fastq.gz')

	script:
	"""
	samtools fastq \
	    -0 /dev/null -1 ${sample_id}_R1.fastq.gz -2 ${sample_id}_R2.fastq.gz \
	    -f 0x3 \
	    -N \
	    $bams 
	"""

}
