// 1st module of the annotation process. I take all the vcfs obtained from the previous workflow and I convert them in maf files.

process vcf2maf_mutect2 {

    //publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    conda '/hpcnfs/home/ieo5898/miniconda3/envs/vcf2maf'

    input:
    tuple val(sample_id), path(vcfs_mutect2)

    output:
    tuple val(sample_id), path('*_mutect2.maf')

    script:
    """
    perl $params.vcf2maf --input-vcf $vcfs_mutect2 --output-maf ${sample_id}_mutect2.maf --ref-fasta $params.mutect2_fasta --vep-path $params.vep --vep-data $params.cache --ncbi-build GRCh38 --tumor-id $sample_id --vep-overwrite
    """

}


process vcf2maf_mutserve {

    //publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    conda '/hpcnfs/home/ieo5898/miniconda3/envs/vcf2maf'

    input:
    tuple val(sample_id), path(vcfs_mutserve)

    output:
    tuple val(sample_id), path('*_mutserve.maf')

    script:
    """
    perl $params.vcf2maf --input-vcf $vcfs_mutserve --output-maf ${sample_id}_mutserve.maf --ref-fasta $params.mutserve_fasta --vep-path $params.vep --vep-data $params.cache --ncbi-build GRCh38 --tumor-id $sample_id --vep-overwrite
    """

}



// These are steps needed for WXS data, in which we sometimes want to use the vcf coming from the MT alignment without NUMTs correction steps.
process vcf2maf_mutect2_before {

    publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'

    conda '/hpcnfs/home/ieo5898/miniconda3/envs/vcf2maf'

    input:
    tuple val(sample_id), path(vcfs_mutect2)

    output:
    tuple val(sample_id), path('*_before_mutect2.maf')

    script:
    """
    perl $params.vcf2maf --input-vcf $vcfs_mutect2 --output-maf ${sample_id}_before_mutect2.maf --ref-fasta $params.mutect2_fasta --vep-path $params.vep --vep-data $params.cache --ncbi-build GRCh38 --tumor-id $sample_id --vep-overwrite
    """

}

process vcf2maf_mutect2_after {

    publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'

    conda '/hpcnfs/home/ieo5898/miniconda3/envs/vcf2maf'

    input:
    tuple val(sample_id), path(vcfs_mutect2)

    output:
    tuple val(sample_id), path('*_after_mutect2.maf')

    script:
    """
    perl $params.vcf2maf --input-vcf $vcfs_mutect2 --output-maf ${sample_id}_after_mutect2.maf --ref-fasta $params.mutect2_fasta --vep-path $params.vep --vep-data $params.cache --ncbi-build GRCh38 --tumor-id $sample_id --vep-overwrite
    """

}

process vcf2maf_mutserve_before {

    //publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'
    publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'

    conda '/hpcnfs/home/ieo5898/miniconda3/envs/vcf2maf'

    input:
    tuple val(sample_id), path(vcfs_mutserve)

    output:
    tuple val(sample_id), path('*_before_mutserve.maf')

    script:
    """
    perl $params.vcf2maf --input-vcf $vcfs_mutserve --output-maf ${sample_id}_before_mutserve.maf --ref-fasta $params.mutserve_fasta --vep-path $params.vep --vep-data $params.cache --ncbi-build GRCh38 --tumor-id $sample_id --vep-overwrite
    """

}

process vcf2maf_mutserve_after {

    publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'

    conda '/hpcnfs/home/ieo5898/miniconda3/envs/vcf2maf'

    input:
    tuple val(sample_id), path(vcfs_mutserve)

    output:
    tuple val(sample_id), path('*_after_mutserve.maf')

    script:
    """
    perl $params.vcf2maf --input-vcf $vcfs_mutserve --output-maf ${sample_id}_after_mutserve.maf --ref-fasta $params.mutserve_fasta --vep-path $params.vep --vep-data $params.cache --ncbi-build GRCh38 --tumor-id $sample_id --vep-overwrite
    """

}
