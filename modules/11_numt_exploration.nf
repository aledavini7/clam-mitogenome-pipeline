process numt_like_exploration {

    publishDir "${params.outdir}/${sample_id}/numt_exploration", mode: 'copy'

    input:
    tuple val(sample_id), path(before_bam), path(after_bam), path(nuclear_bam)

    output:
    tuple val(sample_id),
        path("${sample_id}_mitoscape_rejected_read_ids.txt"),
        path("${sample_id}_mitoscape_rejected_mt.bam"),
        path("${sample_id}_mitoscape_rejected_mt.bam.bai"),
        path("${sample_id}_mitoscape_rejected_nuclear.bam"),
        path("${sample_id}_mitoscape_rejected_nuclear.bam.bai"),
        path("${sample_id}_numt_like_summary.tsv"),
        path("${sample_id}_numt_like_nuclear_alignments.tsv"),
        path("${sample_id}_numt_like_contig_summary.tsv"),
        path("${sample_id}_numt_like_window_summary.tsv")

    script:
    """
    set -euo pipefail
    export LC_ALL=C

    min_mapq=$params.numt_exploration_min_mapq
    window_bp=$params.numt_exploration_cluster_window

    samtools view "$before_bam" | cut -f 1 | awk '{ sub(/\\/[12]\$/, "", \$0); print }' | sort -u > pre_mitoscape.ids
    samtools view "$after_bam" | cut -f 1 | awk '{ sub(/\\/[12]\$/, "", \$0); print }' | sort -u > post_mitoscape.ids

    comm -23 pre_mitoscape.ids post_mitoscape.ids > ${sample_id}_mitoscape_rejected_read_ids.txt
    awk '{
        print
        print \$0 "/1"
        print \$0 "/2"
    }' ${sample_id}_mitoscape_rejected_read_ids.txt | sort -u > rejected_read_query.ids

    retained_reads=\$(comm -12 pre_mitoscape.ids post_mitoscape.ids | wc -l | awk '{print \$1}')
    pre_reads=\$(wc -l < pre_mitoscape.ids | awk '{print \$1}')
    post_reads=\$(wc -l < post_mitoscape.ids | awk '{print \$1}')
    rejected_reads=\$(wc -l < ${sample_id}_mitoscape_rejected_read_ids.txt | awk '{print \$1}')
    rejected_fraction_pct=\$(awk -v rejected="\$rejected_reads" -v pre="\$pre_reads" 'BEGIN { if (pre > 0) printf "%.4f", rejected / pre * 100; else printf "0.0000" }')

    make_empty_bam() {
        template="\$1"
        output="\$2"
        samtools view -H "\$template" | samtools view -b -o "\$output" -
        samtools index "\$output"
    }

    if [ "\$rejected_reads" -gt 0 ]; then
        samtools view -h -N rejected_read_query.ids -b "$before_bam" | \\
            samtools sort -@ $task.cpus -o ${sample_id}_mitoscape_rejected_mt.bam
        samtools view -h -N rejected_read_query.ids -b "$nuclear_bam" | \\
            samtools sort -@ $task.cpus -o ${sample_id}_mitoscape_rejected_nuclear.bam
        samtools index ${sample_id}_mitoscape_rejected_mt.bam
        samtools index ${sample_id}_mitoscape_rejected_nuclear.bam
    else
        make_empty_bam "$before_bam" ${sample_id}_mitoscape_rejected_mt.bam
        make_empty_bam "$nuclear_bam" ${sample_id}_mitoscape_rejected_nuclear.bam
    fi

    mt_alignments=\$(samtools view -c ${sample_id}_mitoscape_rejected_mt.bam)
    nuclear_alignments=\$(samtools view -c ${sample_id}_mitoscape_rejected_nuclear.bam)
    nuclear_primary=\$(samtools view -c -F 2304 ${sample_id}_mitoscape_rejected_nuclear.bam)
    nuclear_proper_pairs=\$(samtools view -c -f 2 ${sample_id}_mitoscape_rejected_nuclear.bam)
    nuclear_high_mapq=\$(samtools view -c -q "\$min_mapq" ${sample_id}_mitoscape_rejected_nuclear.bam)

    if [ "\$rejected_reads" -eq 0 ]; then
        interpretation="no_mitoscape_rejected_reads"
    elif [ "\$nuclear_high_mapq" -gt 0 ]; then
        interpretation="numt_like_signal_review_recommended"
    else
        interpretation="mitoscape_rejected_low_nuclear_mapping_support"
    fi

    printf "sample_id\\tpre_mitoscape_unique_read_ids\\tpost_mitoscape_unique_read_ids\\tretained_read_ids\\tmitoscape_rejected_read_ids\\tmitoscape_rejected_fraction_pct\\trejected_mt_alignments\\trejected_nuclear_alignments\\trejected_nuclear_primary_alignments\\trejected_nuclear_proper_pair_alignments\\trejected_nuclear_high_mapq_alignments\\tmin_mapq_threshold\\tcluster_window_bp\\tinterpretation\\n" > ${sample_id}_numt_like_summary.tsv
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "$sample_id" "\$pre_reads" "\$post_reads" "\$retained_reads" "\$rejected_reads" "\$rejected_fraction_pct" \\
        "\$mt_alignments" "\$nuclear_alignments" "\$nuclear_primary" "\$nuclear_proper_pairs" "\$nuclear_high_mapq" \\
        "\$min_mapq" "\$window_bp" "\$interpretation" >> ${sample_id}_numt_like_summary.tsv

    samtools view ${sample_id}_mitoscape_rejected_nuclear.bam | \\
        gawk -v sample="$sample_id" -v min_mapq="\$min_mapq" 'BEGIN {
            OFS="\\t"
            print "sample_id","read_id","flag","nuclear_contig","position","mapq","cigar","mate_contig","mate_position","template_length","is_primary","is_high_mapq"
        }
        {
            is_primary = (and(\$2, 256) == 0 && and(\$2, 2048) == 0 ? "YES" : "NO")
            is_high_mapq = (\$5 >= min_mapq ? "YES" : "NO")
            print sample,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,is_primary,is_high_mapq
        }' > ${sample_id}_numt_like_nuclear_alignments.tsv

    samtools view ${sample_id}_mitoscape_rejected_nuclear.bam | \\
        awk -v sample="$sample_id" -v min_mapq="\$min_mapq" 'BEGIN { OFS="\\t" }
        {
            contig = \$3
            if (contig == "*" || contig == "") next
            count[contig]++
            sum_mapq[contig] += \$5
            if (\$5 >= min_mapq) high_mapq[contig]++
            if (!(contig in min_pos) || \$4 < min_pos[contig]) min_pos[contig] = \$4
            if (!(contig in max_pos) || \$4 > max_pos[contig]) max_pos[contig] = \$4
        }
        END {
            print "sample_id","nuclear_contig","alignments","high_mapq_alignments","min_pos","max_pos","mean_mapq"
            for (contig in count) {
                printf "%s\\t%s\\t%d\\t%d\\t%d\\t%d\\t%.2f\\n", sample, contig, count[contig], high_mapq[contig] + 0, min_pos[contig], max_pos[contig], sum_mapq[contig] / count[contig]
            }
        }' > ${sample_id}_numt_like_contig_summary.tsv

    samtools view ${sample_id}_mitoscape_rejected_nuclear.bam | \\
        awk -v sample="$sample_id" -v min_mapq="\$min_mapq" -v window_bp="\$window_bp" 'BEGIN { OFS="\\t" }
        {
            contig = \$3
            pos = \$4
            if (contig == "*" || pos <= 0) next
            start = int((pos - 1) / window_bp) * window_bp + 1
            end = start + window_bp - 1
            key = contig SUBSEP start SUBSEP end
            count[key]++
            sum_mapq[key] += \$5
            if (\$5 >= min_mapq) high_mapq[key]++
            read_id = \$1
            sub(/\\/[12]\$/, "", read_id)
            seen[key SUBSEP read_id] = 1
        }
        END {
            for (item in seen) {
                split(item, parts, SUBSEP)
                key = parts[1] SUBSEP parts[2] SUBSEP parts[3]
                unique_reads[key]++
            }
            print "sample_id","nuclear_contig","window_start","window_end","alignments","unique_read_ids","high_mapq_alignments","mean_mapq"
            for (key in count) {
                split(key, parts, SUBSEP)
                printf "%s\\t%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%.2f\\n", sample, parts[1], parts[2], parts[3], count[key], unique_reads[key] + 0, high_mapq[key] + 0, sum_mapq[key] / count[key]
            }
        }' > ${sample_id}_numt_like_window_summary.tsv
    """

}
