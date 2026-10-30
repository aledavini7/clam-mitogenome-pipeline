process nuclear_reference_index {
    tag 'GRCh38'
    container params.nuclear_tool_container

    input:
    path reference

    output:
    tuple path(reference), path("${reference.name}.fai"), path('nuclear_index.*')

    script:
    """
    samtools faidx $reference
    awk '/^>/{sub(/[[:space:]].*/, "")} {print}' $reference > nuclear_reference.bwa.fa
    bwa-mem2 index -p nuclear_index nuclear_reference.bwa.fa
    """
}

process nuclear_alignment {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/nuclear/nuclear_bams", mode: 'copy'
    container params.nuclear_tool_container

    input:
    tuple val(sample_id), path(reads)
    tuple path(reference), path(reference_fai), path(index_files)

    output:
    tuple val(sample_id), path("${sample_id}.nuclear.bam"), path("${sample_id}.nuclear.bam.bai")

    script:
    """
    samtools dict -o nuclear_reference.dict $reference
    set -o pipefail
    bwa-mem2 mem -t $task.cpus \
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \
        nuclear_index ${reads[0]} ${reads[1]} | \
        awk -v dict=nuclear_reference.dict '
            BEGIN {
                while ((getline line < dict) > 0) {
                    if (substr(line, 1, 3) == "@SQ") sq[++n] = line
                }
                close(dict)
            }
            /^@SQ/ { next }
            /^@HD/ {
                print
                for (i = 1; i <= n; i++) print sq[i]
                inserted = 1
                next
            }
            /^@/ {
                if (!inserted) {
                    for (i = 1; i <= n; i++) print sq[i]
                    inserted = 1
                }
                print
                next
            }
            {
                if (!inserted) {
                    for (i = 1; i <= n; i++) print sq[i]
                    inserted = 1
                }
                print
            }
            END {
                if (!inserted) {
                    for (i = 1; i <= n; i++) print sq[i]
                }
                if (n == 0) {
                    print "ERROR: samtools dict produced no @SQ records" > "/dev/stderr"
                    exit 2
                }
            }
        ' | \
        samtools sort -n -@ $task.cpus -o ${sample_id}.name.bam -
    samtools fixmate -m ${sample_id}.name.bam ${sample_id}.fixmate.bam
    samtools sort -@ $task.cpus -o ${sample_id}.coord.bam ${sample_id}.fixmate.bam
    samtools markdup -@ $task.cpus ${sample_id}.coord.bam ${sample_id}.nuclear.bam
    samtools index ${sample_id}.nuclear.bam
    samtools quickcheck -v ${sample_id}.nuclear.bam
    """
}

process nuclear_prepare_bam {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/nuclear/nuclear_bams", mode: 'copy'
    container params.nuclear_tool_container

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("${sample_id}.nuclear.bam"), path("${sample_id}.nuclear.bam.bai")

    script:
    """
    samtools quickcheck -v $input_bam
    samtools sort -n -@ $task.cpus -o ${sample_id}.name.bam $input_bam
    samtools fixmate -m ${sample_id}.name.bam ${sample_id}.fixmate.bam
    samtools sort -@ $task.cpus -o ${sample_id}.coord.bam ${sample_id}.fixmate.bam
    samtools markdup -@ $task.cpus ${sample_id}.coord.bam ${sample_id}.nuclear.bam
    samtools index ${sample_id}.nuclear.bam
    samtools quickcheck -v ${sample_id}.nuclear.bam
    """
}

process nuclear_tumor_only_call {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/nuclear/nuclear_variants", mode: 'copy'
    container params.mutect2_container

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference), path(reference_fai), path(target_bed)
    path germline_resource
    path pon

    output:
    tuple val(sample_id), path("${sample_id}.somatic.unfiltered.vcf"), path("${sample_id}.somatic.filtered.vcf")

    script:
    def germline_arg = germline_resource ? "--germline-resource ${germline_resource[0]}" : ''
    def pon_arg = pon ? "--panel-of-normals ${pon[0]}" : ''
    """
    gatk CreateSequenceDictionary -R $reference -O ${reference.baseName}.dict
    gatk Mutect2 \
        -R $reference \
        -I $bam \
        -tumor $sample_id \
        -L $target_bed \
        $germline_arg $pon_arg \
        -O ${sample_id}.somatic.unfiltered.vcf
    gatk FilterMutectCalls \
        -R $reference \
        -V ${sample_id}.somatic.unfiltered.vcf \
        -O ${sample_id}.somatic.filtered.vcf
    """
}

process nuclear_paired_call {
    tag "$pair_id"
    publishDir "${params.outdir}/${tumor_id}/nuclear/nuclear_variants", mode: 'copy'
    container params.mutect2_container

    input:
    tuple val(pair_id), val(tumor_id), path(tumor_bam), path(tumor_bai), val(normal_id), path(normal_bam), path(normal_bai), path(reference), path(reference_fai), path(target_bed)
    path germline_resource
    path pon

    output:
    tuple val(tumor_id), path("${tumor_id}.somatic.unfiltered.vcf"), path("${tumor_id}.somatic.filtered.vcf")

    script:
    def germline_arg = germline_resource ? "--germline-resource ${germline_resource[0]}" : ''
    def pon_arg = pon ? "--panel-of-normals ${pon[0]}" : ''
    """
    gatk CreateSequenceDictionary -R $reference -O ${reference.baseName}.dict
    gatk Mutect2 \
        -R $reference \
        -I $tumor_bam -tumor $tumor_id \
        -I $normal_bam -normal $normal_id \
        -L $target_bed \
        $germline_arg $pon_arg \
        -O ${tumor_id}.somatic.unfiltered.vcf
    gatk FilterMutectCalls \
        -R $reference \
        -V ${tumor_id}.somatic.unfiltered.vcf \
        -O ${tumor_id}.somatic.filtered.vcf
    """
}

process nuclear_germline_call {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/nuclear/nuclear_variants", mode: 'copy'
    container params.mutect2_container

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference), path(reference_fai), path(target_bed)

    output:
    tuple val(sample_id), path("${sample_id}.germline.g.vcf.gz"), path("${sample_id}.germline.filtered.vcf")

    script:
    """
    gatk CreateSequenceDictionary -R $reference -O ${reference.baseName}.dict
    gatk HaplotypeCaller \
        -R $reference -I $bam -L $target_bed \
        -ERC GVCF -O ${sample_id}.germline.g.vcf.gz
    gatk GenotypeGVCFs \
        -R $reference -V ${sample_id}.germline.g.vcf.gz \
        -O ${sample_id}.germline.raw.vcf
    gatk VariantFiltration \
        -R $reference -V ${sample_id}.germline.raw.vcf \
        --filter-name SNP_HARD_FILTER \
        --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
        --filter-name INDEL_HARD_FILTER \
        --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' \
        -O ${sample_id}.germline.filtered.vcf
    """
}

process nuclear_annotate {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/nuclear/nuclear_annotation", mode: 'copy'
    container params.nuclear_tool_container

    input:
    tuple val(sample_id), path(raw_vcf), path(reference), path(reference_fai)
    path vep_cache
    path gnomad_vcf
    path cosmic_vcf

    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf"), path("${sample_id}.maf")

    script:
    def custom = []
    if (gnomad_vcf) custom << "--custom file=${gnomad_vcf[0]},short_name=gnomAD,format=vcf,type=exact,coords=0,fields=AF"
    if (cosmic_vcf) custom << "--custom file=${cosmic_vcf[0]},short_name=COSMIC,format=vcf,type=exact,coords=0,fields=ID"
    """
    vep \
        --input_file $raw_vcf \
        --output_file ${sample_id}.annotated.vcf \
        --format vcf --vcf --offline --cache --dir_cache $vep_cache \
        --species homo_sapiens --assembly GRCh38 --cache_version 112 --fasta $reference \
        --symbol --canonical --mane --hgvs --numbers --fork $task.cpus \
        --force_overwrite \
        ${custom.join(' ')}
    python3 ${projectDir}/bin/clam_nuclear_vcf_to_maf.py \
        --sample-id $sample_id --input ${sample_id}.annotated.vcf --output ${sample_id}.maf
    """
}

process nuclear_cnvkit_flat {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/nuclear/nuclear_cnv", mode: 'copy'
    container params.nuclear_tool_container

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference), path(reference_fai), path(target_bed)

    output:
    tuple val(sample_id), path('cnvkit_out')

    script:
    """
    mkdir -p cnvkit_out
    cnvkit.py batch $bam -n \
        --targets $target_bed --fasta $reference \
        --output-reference cnvkit_out/flat_reference.cnn \
        --output-dir cnvkit_out --processes $task.cpus
    for cns in cnvkit_out/*.cns; do
        cnvkit.py call \$cns -o \${cns%.cns}.call.cns
        cnvkit.py export seg \${cns%.cns}.call.cns -o \${cns%.cns}.seg
    done
    for cnr in cnvkit_out/*.cnr; do
        [[ -f "\$cnr" ]] || continue
        cns="\${cnr%.cnr}.call.cns"
        [[ -s "\$cns" ]] || continue
        plot="\${cnr%.cnr}.scatter.png"
        cnvkit.py scatter "\$cnr" -s "\$cns" -o "\$plot" --title "${sample_id} genome-wide log2 copy ratio"
    done
    """
}

process nuclear_cnvkit_paired {
    tag "$pair_id"
    publishDir "${params.outdir}/${tumor_id}/nuclear/nuclear_cnv", mode: 'copy'
    container params.nuclear_tool_container

    input:
    tuple val(pair_id), val(tumor_id), path(tumor_bam), path(tumor_bai), val(normal_id), path(normal_bam), path(normal_bai), path(reference), path(reference_fai), path(target_bed)

    output:
    tuple val(tumor_id), path('cnvkit_out')

    script:
    """
    mkdir -p cnvkit_out
    cnvkit.py batch $tumor_bam --normal $normal_bam \
        --targets $target_bed --fasta $reference \
        --output-reference cnvkit_out/paired_reference.cnn \
        --output-dir cnvkit_out --processes $task.cpus
    for cns in cnvkit_out/*.cns; do
        cnvkit.py call \$cns -o \${cns%.cns}.call.cns
        cnvkit.py export seg \${cns%.cns}.call.cns -o \${cns%.cns}.seg
    done
    for cnr in cnvkit_out/*.cnr; do
        [[ -f "\$cnr" ]] || continue
        cns="\${cnr%.cnr}.call.cns"
        [[ -s "\$cns" ]] || continue
        plot="\${cnr%.cnr}.scatter.png"
        cnvkit.py scatter "\$cnr" -s "\$cns" -o "\$plot" --title "${tumor_id} genome-wide log2 copy ratio"
    done
    """
}

process nuclear_report {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/nuclear/nuclear_report", mode: 'copy'
    container params.annotation_container

    input:
    tuple val(sample_id), val(calling_mode), path(annotated_vcf), path(maf), path(cnvkit_out)

    output:
    tuple val(sample_id), path("${sample_id}_nuclear_report.html")

    script:
    """
    python3 ${projectDir}/bin/clam_build_nuclear_report.py \
        --sample-id $sample_id --calling-mode $calling_mode --vcf $annotated_vcf --maf $maf \
        --cnv-dir $cnvkit_out --common-af-threshold ${params.nuclear_common_af_threshold} \
        --top-variants ${params.nuclear_report_top_variants} \
        --output ${sample_id}_nuclear_report.html
    """
}
