// Third module: Mitoscape.nf

process mitoscape {

	publishDir "${params.outdir}/${sample_id}/bams", mode: 'copy'

	input:
	tuple val(sample_id), path(mt_bams), path(nuc_bams), path(bams)

	output:
    tuple val(sample_id), path('*_MTDNA.bam')

	script:
	"""
    java -Xmx20G -jar ${params.classifier} \
	--threads $task.cpus \
	--prob 0.5 \
	--ld ${params.mitomap} \
	--numt ${params.numt} \
	--classifier ${params.model} \
	--prefix ${sample_id} \
	--out ${sample_id}_MTDNA.bam
	"""

}