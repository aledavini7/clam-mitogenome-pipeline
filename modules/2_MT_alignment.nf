// First module: MT_alignment
// First step. We take the new reads (fastq R1 and R2 obtained from bam2fq) and align them only to the mitochondrial chromosome.
// I use gsnap as suggested in the github page of MitoScape (https://github.com/larryns/MitoScape/blob/master/README.md)
// MitoScape was tested on gsnap, which is able to handle the circular mtDNA chromosome, but in theory any aligner should do.

process alignment_MT {

	publishDir "${params.outdir}/${sample_id}/bams", mode: 'copy'

        input:
	tuple val(sample_id), path(reads)
	
        output:
	tuple val(sample_id), path('*_MT.bam')

        script:
	"""
	gsnap \
        -D $params.datadir \
        -d $params.mt_gsnap_db \
        --input-buffer-size=25000 \
        --gunzip \
        --nthreads=$task.cpus \
        --npaths=1 \
        -A sam \
        --read-group-id=$sample_id \
	--read-group-name=MT_$sample_id \
        --read-group-library=bar \
        --read-group-platform=illumina \
        ${reads[0]} ${reads[1]} | \
        samtools view -f 0x3 -u -| \
        samtools sort -n -o ${sample_id}_MT.bam -m 3G -O BAM
	"""	

}

process get_the_MD {

	publishDir "${params.outdir}/${sample_id}/bams", mode: 'copy'

	input:
	tuple val(sample_id), path(bams)

	output:
	tuple val(sample_id), path('*_MT_MD.bam')

	script:
	"""
	samtools calmd -e --output-fmt BAM $bams $params.fasta > ${sample_id}_MT_MD.bam
	"""
	
}
