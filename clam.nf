/*
 * pipeline input parameters
 */

// Example:
// nextflow run clam.nf --input 'data/*_{R1,R2}.fastq.gz' --input_type fastq
// nextflow run clam.nf --input 'data/*.bam' --input_type bam --analysis_mode wgs_numt_correction

nextflow.enable.dsl = 2

def normalizeAnalysisMode(value) {
    def requested = value?.toString()?.trim()
    if (!requested) {
        return null
    }

    def aliases = [
        'wgs_numt_correction': 'wgs_numt_correction',
        'wgs-numt-correction': 'wgs_numt_correction',
        'wgs'                : 'wgs_numt_correction',
        'numt'               : 'wgs_numt_correction',
        'full'               : 'wgs_numt_correction',
        'grch38'             : 'wgs_numt_correction',
        'hg38'               : 'wgs_numt_correction',
        'rcrs_only'          : 'rcrs_only',
        'rcrs-only'          : 'rcrs_only',
        'mitochondrial_only' : 'rcrs_only',
        'mitochondrial-only' : 'rcrs_only',
        'short'              : 'rcrs_only',
        'rcrs'               : 'rcrs_only'
    ]

    return aliases[requested.toLowerCase()]
}

def sanitizeSampleId(value) {
    def sample_id = value?.toString()?.trim()
    if (!sample_id) {
        return null
    }

    return sample_id
        .replaceAll(/[^A-Za-z0-9_.-]+/, '_')
        .replaceAll(/^_+|_+$/, '')
}

def integerParamError(name, value, minValue) {
    try {
        def parsed = Integer.parseInt(value?.toString())
        if (parsed < minValue) {
            return "--${name} must be >= ${minValue}; received '${value}'"
        }
        return null
    } catch (Exception ignored) {
        return "--${name} must be an integer; received '${value}'"
    }
}

def numberParamError(name, value, minValue, maxValue) {
    try {
        def parsed = Double.parseDouble(value?.toString())
        if (parsed < minValue || parsed > maxValue) {
            return "--${name} must be between ${minValue} and ${maxValue}; received '${value}'"
        }
        return null
    } catch (Exception ignored) {
        return "--${name} must be numeric; received '${value}'"
    }
}

include { sorting_1; bam2fqgz } from './modules/1_data_preparation.nf'
include { alignment_MT; get_the_MD; get_the_MD_mitoscape } from './modules/2_MT_alignment.nf'
include { bam2fq; alignment_Nuc } from './modules/3_NUC_alignment.nf'
include { mitoscape } from './modules/4_Mitoscape.nf'
include { sorting; coverage; sorting_single; coverage_single } from './modules/5_Coverage.nf'
include { index } from './modules/6_Index.nf'
include { mutect2; filter_mutect2; bgzip; mutserve } from './modules/7_Variant_calling.nf'
include { haplogrep; haplogrep1 } from './modules/8_Haplogroups_assignment.nf'
include { merge_variant_calls; annotation_mafs; mitomap_annotation } from './modules/9_annotation.nf'
include { bam_qc; no_mitomap_report_table; no_numt_exploration_report_table; build_report } from './modules/10_qc_report.nf'
include { numt_like_exploration } from './modules/11_numt_exploration.nf'

workflow {

    input_pattern = params.input?.toString()?.trim()
    if (!input_pattern) {
        error "Missing required parameter: --input"
    }

    input_type = params.input_type?.toString()?.trim()?.toLowerCase()
    if (!input_type) {
        error "Missing required parameter: --input_type. Supported values: fastq, bam"
    }

    if (!['fastq', 'bam'].contains(input_type)) {
        error "Unsupported --input_type '${params.input_type}'. Supported values: fastq, bam"
    }

    if (!params.outdir?.toString()?.trim()) {
        error "Missing required parameter: --outdir"
    }

    analysis_mode = normalizeAnalysisMode(params.analysis_mode ?: params.genome)
    if (!analysis_mode) {
        error "Unsupported --analysis_mode '${params.analysis_mode ?: params.genome}'. Supported values: wgs_numt_correction, rcrs_only"
    }

    if (!params.analysis_modes || !params.analysis_modes[analysis_mode]) {
        error "Missing configuration for analysis mode '${analysis_mode}'. Check conf/genomes.config"
    }

    run_numt_correction = params.run_numt_correction != null && params.run_numt_correction.toString().toBoolean()
    run_numt_exploration = params.run_numt_exploration != null && params.run_numt_exploration.toString().toBoolean()
    reference_set = params.reference_set ?: (run_numt_correction ? 'GRCh38+rCRS' : 'rCRS')

    if (run_numt_exploration && !run_numt_correction) {
        error "--run_numt_exploration is only available with --analysis_mode wgs_numt_correction"
    }

    if (!params.datadir || !params.mt_gsnap_db || !params.fasta || !params.mutserve_fasta || !params.mt_contig) {
        error "Incomplete mitochondrial reference configuration for analysis mode '${analysis_mode}'"
    }

    if (run_numt_correction && (!params.genomedir || !params.nuc_gsnap_db || !params.numt)) {
        error "Incomplete nuclear/NUMTs reference configuration for analysis mode '${analysis_mode}'"
    }

    [
        ['annotation_high_depth', params.annotation_high_depth, 1],
        ['annotation_medium_depth', params.annotation_medium_depth, 1],
        ['annotation_high_alt_reads', params.annotation_high_alt_reads, 1],
        ['annotation_medium_alt_reads', params.annotation_medium_alt_reads, 1],
        ['mitomap_keep_columns', params.mitomap_keep_columns, 0],
        ['report_mt_length', params.report_mt_length, 1],
        ['numt_exploration_min_mapq', params.numt_exploration_min_mapq, 0],
        ['numt_exploration_cluster_window', params.numt_exploration_cluster_window, 1],
    ].each { name, value, minValue ->
        def message = integerParamError(name, value, minValue)
        if (message) {
            error message
        }
    }

    def af_message = numberParamError('annotation_max_af_discordance', params.annotation_max_af_discordance, 0.0, 1.0)
    if (af_message) {
        error af_message
    }

    if ((params.annotation_high_depth as Integer) < (params.annotation_medium_depth as Integer)) {
        error "--annotation_high_depth must be >= --annotation_medium_depth"
    }

    if ((params.annotation_high_alt_reads as Integer) < (params.annotation_medium_alt_reads as Integer)) {
        error "--annotation_high_alt_reads must be >= --annotation_medium_alt_reads"
    }

    mitomap_variant_table = params.mitomap_variant_table?.toString()?.trim()

    log.info """\
        C L A M - N F   P I P E L I N E
        ===============================
        input           : ${input_pattern}
        input_type      : ${input_type}
        analysis_mode   : ${analysis_mode}
        reference_set   : ${reference_set}
        outdir          : ${params.outdir}
        mt_gsnap_dir    : ${params.datadir}
        mt_gsnap_db     : ${params.mt_gsnap_db}
        mt_fasta        : ${params.fasta}
        mt_contig       : ${params.mt_contig}
        numt_correction : ${run_numt_correction}
        nuc_gsnap_dir   : ${params.genomedir ?: 'not used'}
        nuc_gsnap_db    : ${params.nuc_gsnap_db ?: 'not used'}
        numt_regions    : ${params.numt ?: 'not used'}
        numt_explore    : ${run_numt_exploration}
        mitomap_table   : ${mitomap_variant_table ?: 'not used'}
        """
        .stripIndent()

    if (input_type == 'fastq') {
        reads_ch = Channel
            .fromFilePairs(input_pattern, size: 2, checkIfExists: true)
            .ifEmpty { error "No FASTQ pairs matched --input '${input_pattern}'. Use a paired-end pattern such as '*_{R1,R2}.fastq.gz'." }
            .map { sample_id, reads ->
                def clean_id = sanitizeSampleId(sample_id)
                if (!clean_id) {
                    error "Could not derive a valid sample ID from FASTQ pair '${sample_id}'"
                }
                [clean_id, reads]
            }
    }

    if (input_type == 'bam') {
        bam_ch = Channel
            .fromPath(input_pattern, checkIfExists: true)
            .ifEmpty { error "No BAM files matched --input '${input_pattern}'." }
            .map { bam ->
                def clean_id = sanitizeSampleId(bam.simpleName)
                if (!clean_id) {
                    error "Could not derive a valid sample ID from BAM '${bam}'"
                }
                [clean_id, bam]
            }

        sort_bam_ch = sorting_1(bam_ch)
        reads_ch = bam2fqgz(sort_bam_ch)
    }

    bams_ch = alignment_MT(reads_ch)
    md_ch = get_the_MD(bams_ch)

    if (run_numt_correction) {
        md_mitoscape_ch = get_the_MD_mitoscape(bams_ch)
        fq_ch = bam2fq(bams_ch)
        nuc_ch = alignment_Nuc(fq_ch)

        md_mitoscape_ch
            .combine(nuc_ch, by: 0)
            .set { mito_ch }

        mito_ch
            .combine(bams_ch, by: 0)
            .set { final_mito_ch }

        mt_ch = mitoscape(final_mito_ch)

        sort_ch = sorting(md_ch, mt_ch)
        cov_ch = coverage(sort_ch)
        final_coverage_ch = cov_ch[1]
        ind_ch = index(sort_ch[1])

        sort_ch[0]
            .combine(sort_ch[1], by: 0)
            .map { sample_id, before_bam, after_bam -> [sample_id, 'before,after', [before_bam, after_bam]] }
            .set { report_bam_input_ch }

        cov_ch[0]
            .combine(cov_ch[1], by: 0)
            .set { report_coverage_ch }

        sort_ch[1]
            .combine(ind_ch, by: 0)
            .set { var_ch }

        if (run_numt_exploration) {
            sort_ch[0]
                .combine(sort_ch[1], by: 0)
                .combine(nuc_ch, by: 0)
                .set { numt_exploration_input_ch }

            numt_exploration_raw_ch = numt_like_exploration(numt_exploration_input_ch)
        }
    } else {
        sort_single_ch = sorting_single(md_ch)
        cov_ch = coverage_single(sort_single_ch)
        final_coverage_ch = cov_ch
        ind_ch = index(sort_single_ch)

        sort_single_ch
            .map { sample_id, sorted_bam -> [sample_id, 'final', [sorted_bam]] }
            .set { report_bam_input_ch }

        cov_ch
            .map { sample_id, coverage_file -> [sample_id, coverage_file, coverage_file] }
            .set { report_coverage_ch }

        sort_single_ch
            .combine(ind_ch, by: 0)
            .set { var_ch }
    }

    raw_mutect2_ch = mutect2(var_ch)
    vcf_2_ch = mutserve(var_ch)

    raw_mutect2_ch[0]
        .combine(raw_mutect2_ch[1], by: 0)
        .combine(raw_mutect2_ch[2], by: 0)
        .set { mutect2_filter_input_ch }

    filtered_mutect2_ch = filter_mutect2(mutect2_filter_input_ch)

    filtered_mutect2_ch[0]
        .combine(vcf_2_ch[0], by: 0)
        .combine(final_coverage_ch, by: 0)
        .set { annotation_input_ch }

    variant_summary_ch = merge_variant_calls(annotation_input_ch)
    variant_maf_ch = annotation_mafs(variant_summary_ch)

    if (mitomap_variant_table) {
        mitomap_table_ch = Channel.value(file(mitomap_variant_table, checkIfExists: true))

        variant_summary_ch
            .combine(mitomap_table_ch)
            .set { mitomap_input_ch }

        mitomap_annotation_ch = mitomap_annotation(mitomap_input_ch)
        mitomap_annotation_ch
            .map { sample_id, mutect2_mitomap_tsv, mutserve_mitomap_tsv, all_mitomap_tsv, confidence_mitomap_tsv -> [sample_id, true, confidence_mitomap_tsv] }
            .set { mitomap_report_ch }
    } else {
        mitomap_report_ch = no_mitomap_report_table(variant_summary_ch)
    }

    if (run_numt_exploration) {
        numt_exploration_raw_ch
            .map { sample_id, rejected_ids, rejected_mt_bam, rejected_mt_bai, rejected_nuclear_bam, rejected_nuclear_bai, summary_tsv, alignments_tsv, contig_tsv, window_tsv ->
                [sample_id, true, summary_tsv, contig_tsv, window_tsv]
            }
            .set { numt_report_ch }
    } else {
        numt_report_ch = no_numt_exploration_report_table(variant_summary_ch)
    }
    
    gz_ch = bgzip(filtered_mutect2_ch[0])
    
    hplg_ch = haplogrep(gz_ch)
    hplg_2_ch = haplogrep1(vcf_2_ch[0])

    bam_qc_ch = bam_qc(report_bam_input_ch)

    bam_qc_ch
        .combine(report_coverage_ch, by: 0)
        .combine(raw_mutect2_ch[0], by: 0)
        .combine(filtered_mutect2_ch[0], by: 0)
        .combine(vcf_2_ch[0], by: 0)
        .combine(vcf_2_ch[1], by: 0)
        .combine(hplg_ch, by: 0)
        .combine(hplg_2_ch, by: 0)
        .combine(variant_summary_ch, by: 0)
        .combine(mitomap_report_ch, by: 0)
        .combine(numt_report_ch, by: 0)
        .set { report_input_ch }

    report_ch = build_report(report_input_ch)


}
