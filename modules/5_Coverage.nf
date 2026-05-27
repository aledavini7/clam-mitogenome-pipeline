// Fourth module: Coverage.nf

process sorting {

    publishDir "${params.outdir}/${sample_id}/bams", mode: 'copy'

    input:
    tuple val(sample_id), path(md_bams), path(mito_bams)

    output:
    tuple val(sample_id), path('*_sorted_before.bam') // sorted before mitoscape
    tuple val(sample_id), path('*_sorted_after.bam') // sorted after mitoscape

    script:
    """
    set -euo pipefail

    case "\$(basename "$md_bams")" in
        ${sample_id}_*) ;;
        *) echo "ERROR: before-MitoScape BAM '$md_bams' does not match sample '$sample_id'" >&2; exit 1 ;;
    esac

    case "\$(basename "$mito_bams")" in
        ${sample_id}_*) ;;
        *) echo "ERROR: after-MitoScape BAM '$mito_bams' does not match sample '$sample_id'" >&2; exit 1 ;;
    esac

    samtools sort $md_bams > ${sample_id}_sorted_before.bam
    samtools sort $mito_bams > ${sample_id}_sorted_after.bam
    """

}

process coverage {

    publishDir "${params.outdir}/${sample_id}/coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(before_bams), path(after_bams)

    output:
    tuple val(sample_id), path('*_before.coverage') // coverage before mitoscape
    tuple val(sample_id), path('*_after.coverage') // coverage after mitoscape
   

    script:
    """
    set -euo pipefail

    case "\$(basename "$before_bams")" in
        ${sample_id}_*) ;;
        *) echo "ERROR: before-coverage BAM '$before_bams' does not match sample '$sample_id'" >&2; exit 1 ;;
    esac

    case "\$(basename "$after_bams")" in
        ${sample_id}_*) ;;
        *) echo "ERROR: after-coverage BAM '$after_bams' does not match sample '$sample_id'" >&2; exit 1 ;;
    esac

    samtools depth $before_bams > ${sample_id}_before.coverage
    samtools depth $after_bams > ${sample_id}_after.coverage
    """

}


// Sorting and Coverage single file (if I have mitochondrial sequencing data, I have only information about mtDNA, so I perform just single steps)

process sorting_single {

    publishDir "${params.outdir}/${sample_id}/bams", mode: 'copy'

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

    publishDir "${params.outdir}/${sample_id}/coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bams)

    output:
    tuple val(sample_id), path('*.coverage') // coverage 
   

    script:
    """
    samtools depth $sorted_bams > ${sample_id}.coverage
    """

}
