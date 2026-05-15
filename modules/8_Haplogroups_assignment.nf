// Seventh module: Haplogroups_assignment.nf

process haplogrep {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(gz_vcfs)

    output:
    tuple val(sample_id), path('*_haplogroups.txt')

    script:
    """
    $params.haplogrep classify --tree $params.haplogrep_tree --in $gz_vcfs --out ${sample_id}_mutect2_haplogroups.txt
    """

}

process haplogrep1 {

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path('*_mutserve_haplogroups.txt')

    script:
    """
    $params.haplogrep classify --tree $params.haplogrep_tree --in $vcfs --out ${sample_id}_mutserve_haplogroups.txt
    """

}
