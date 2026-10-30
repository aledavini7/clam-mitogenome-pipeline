process nuclear_manta {
    tag "$sample_id ($calling_mode)"
    publishDir "${params.outdir}/${sample_id}/nuclear/sv", mode: 'copy'
    container params.manta_container

    input:
    tuple val(sample_id), val(calling_mode), path(bams), path(bais), path(reference), path(reference_fai)

    output:
    tuple val(sample_id), path('manta_out')

    script:
    def configuration = calling_mode == 'paired' ? "--tumorBam ${bams[0]} --normalBam ${bams[1]}" : calling_mode == 'germline' ? "--bam ${bams[0]}" : "--tumorBam ${bams[0]}"
    def result_vcf = calling_mode == 'paired' ? 'somaticSV.vcf.gz' : calling_mode == 'germline' ? 'diploidSV.vcf.gz' : 'tumorSV.vcf.gz'

    """
    set -euo pipefail

    # Restrict the workflow to primary chromosomes. Manta documents this as a
    # way to avoid excess work on GRCh38 decoy and unplaced contigs. The exome
    # mode remains enabled, while the call regions stay broad enough to include
    # breakpoints outside the capture BED itself.
    awk -F '\\t' '
        {
            name = \$1
            sub(/^chr/, "", name)
            if (name ~ /^([1-9]|1[0-9]|2[0-2]|X|Y)\$/) print \$1 "\\t0\\t" \$2
        }
    ' $reference_fai > manta_call_regions.bed
    bgzip -c manta_call_regions.bed > manta_call_regions.bed.gz
    tabix -f -p bed manta_call_regions.bed.gz

    configManta.py \\
        ${configuration} \\
        --referenceFasta $reference \\
        --runDir manta_run \\
        --callRegions manta_call_regions.bed.gz \\
        --exome
    manta_run/runWorkflow.py -j $task.cpus

    mkdir -p manta_out
    cp -a manta_run/results/variants/. manta_out/
    test -s manta_out/${result_vcf}
    """
}
