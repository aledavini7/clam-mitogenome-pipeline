// Second module: NUC_aignment.nf

process bam2fq {

	input:
	tuple val(sample_id), path(bams)

	output:
	tuple val(sample_id), path('*.fastq')

	script:
	"""
	samtools fastq \
	    -0 /dev/null -1 ${sample_id}_R1.fastq -2 ${sample_id}_R2.fastq \
	    -f 0x3 \
	    -N \
	    $bams 
	"""

}


process alignment_Nuc {

	publishDir "${params.outdir}/${sample_id}/bams", mode: 'copy'

	input:
	tuple val(sample_id), path(fastqs)

	output:
    tuple val(sample_id), path('*_NT.bam')

	script:
	"""
	gsnap \
        -D $params.genomedir \
        -d $params.nuc_gsnap_db \
        --input-buffer-size=10000 \
        --nthreads=$task.cpus \
        --npaths=1 \
        --print-snps \
        -A sam \
		--read-group-id=$sample_id \
		--read-group-name=MT_$sample_id \
        --read-group-library=bar \
        --read-group-platform=illumina \
        ${fastqs[0]} ${fastqs[1]} | \
        samtools view -f 0x3 -u | \
        samtools sort -o ${sample_id}_NT.bam -m 3G -O BAM
	"""

}
