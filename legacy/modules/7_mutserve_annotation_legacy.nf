// Legacy mutserve annotation processes removed from the active CLAM workflow.
//
// These relied on local mutserve annotation resources and were superseded by the
// current CLAM annotation tables, MITOMAP overlay, MAF generation, and HTML
// report. They are kept here only as historical reference.

process mutserve_annotate {

    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path('*_mutserve_annotated.txt')

    script:
    """
    $params.mutserve annotate --input $vcfs --annotation $params.annotation --output ${sample_id}_mutserve_annotated.txt
    """

}

process mutserve_annotate_before {

    publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path('*_before_mutserve_annotated.txt')

    script:
    """
    $params.mutserve annotate --input $vcfs --annotation $params.annotation --output ${sample_id}_before_mutserve_annotated.txt
    """

}

process mutserve_annotate_after {

    publishDir "${params.outdir}/analysis/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path('*_after_mutserve_annotated.txt')

    script:
    """
    $params.mutserve annotate --input $vcfs --annotation $params.annotation --output ${sample_id}_after_mutserve_annotated.txt
    """

}
