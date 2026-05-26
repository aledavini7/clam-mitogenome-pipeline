
process sort_and_view {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(first_mito_bams) // output of 'get_the_MD'
    tuple val(sample_id), path(second_mito_bams) // output of 'mitoscape'

    output:
    tuple val(sample_id), path('*_MT_MD.txt')
    tuple val(sample_id), path('*_MTDNA.txt')

    script:
    """
    samtools view $first_mito_bams | cut -f 1 | sort > ${sample_id}_MT_MD.txt
    samtools view $second_mito_bams | cut -f 1 | sort > ${sample_id}_MTDNA.txt
    """

}


process remove_duplicates {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(first_mito_reads_id) 
    tuple val(sample_id), path(second_mito_reads_id)

    output:
    tuple val(sample_id), path('*_numts.txt')

    script:
    """
    comm -23 $first_mito_reads_id $second_mito_reads_id > ${sample_id}_numts.txt
    """

}

process subtraction {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(numts_reads_id), path(mito_bams) 
    //tuple val(sample_id), path(mito_bams)

    output:
    tuple val(sample_id), path('*_numts.bam')

    script:
    """
    samtools view -N $numts_reads_id -o ${sample_id}_numts.bam $mito_bams
    """

}

process alignment_numt_Nuc {

	publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

	input:
	tuple val(sample_id), path(fastqs)

	output:
    tuple val(sample_id), path('*_numt_NT.bam')

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
        samtools sort -o ${sample_id}_numt_NT.bam -m 3G -O BAM
	"""

}

process sorting_numts {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(numt_bams)

    output:
    tuple val(sample_id), path('*_numt_sorted.bam')

    script:
    """
    samtools sort $numt_bams > ${sample_id}_numt_sorted.bam
    """

}

process sorting_nuclear_numts {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(numt_nuc_bams)

    output:
    tuple val(sample_id), path('*_numt_nuc_sorted.bam')

    script:
    """
    samtools sort $numt_nuc_bams > ${sample_id}_numt_nuc_sorted.bam
    """

}


process index_numt {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bams)

    output:
    tuple val(sample_id), path('*_numt_sorted.bam.bai')

    script:
    """
    samtools index $sorted_bams
    """

}


process coverage_numts {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(numt_bams)

    output:
    tuple val(sample_id), path('*_numt.coverage')
   

    script:
    """
    samtools depth $numt_bams > ${sample_id}_numt.coverage
    """

}

process coverage_nuclear_numts {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(numt_nuc_bams)

    output:
    tuple val(sample_id), path('*_numt_nuc.coverage') 
   

    script:
    """
    samtools depth $numt_nuc_bams > ${sample_id}_numt_nuc.coverage
    """

}


process mutect2_numts {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(numt_bams), path(index_bai)

    output:
    tuple val(sample_id), path('*.vcf')
    tuple val(sample_id), path('*.vcf.idx')
    tuple val(sample_id), path('*.vcf.stats')

    script:
    """
    gatk Mutect2 --input $numt_bams --reference $params.fasta --mitochondria-mode true -L $params.mt_contig --output ${sample_id}_numts_mutect2.vcf
    """

}

process mutserve_numts {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(numt_bams), path(index_bai)

    output:
    tuple val(sample_id), path('*_mutserve.vcf')
    tuple val(sample_id), path('*_mutserve.txt')

    script:
    """
    $params.mutserve call --reference $params.mutserve_fasta --contig-name $params.mt_contig --output ${sample_id}_mutserve.vcf $numt_bams
    """

}


process haplogrep_numts {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(gz_vcfs)

    output:
    tuple val(sample_id), path('*_mutect2_numts_haplogroups.txt')

    script:
    """
    $params.haplogrep classify --tree $params.haplogrep_tree --in $gz_vcfs --out ${sample_id}_mutect2_numts_haplogroups.txt
    """

}

process haplogrep1_numts {

    publishDir "${params.outdir}/${sample_id}/numts", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path('*_mutserve_numts_haplogroups.txt')

    script:
    """
    $params.haplogrep classify --tree $params.haplogrep_tree --in $vcfs --out ${sample_id}_mutserve_numts_haplogroups.txt
    """

}










