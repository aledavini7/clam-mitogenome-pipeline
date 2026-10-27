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
        'with_numt_correction'    : 'wgs_numt_correction',
        'without_numt_correction' : 'rcrs_only',
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
include { nuclear_reference_index; nuclear_alignment; nuclear_prepare_bam; nuclear_tumor_only_call; nuclear_paired_call; nuclear_germline_call; nuclear_annotate; nuclear_cnvkit_flat; nuclear_cnvkit_paired; nuclear_report } from './modules/12_nuclear_wes.nf'

process merge_fastq_lanes {
    tag "$sample_id"

    input:
    tuple val(sample_id), val(lanes), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz")

    script:
    def r1 = lanes.collect { lane -> reads.find { it.name == lane.r1 } }
    def r2 = lanes.collect { lane -> reads.find { it.name == lane.r2 } }
    """
    cat ${r1.join(' ')} > ${sample_id}_R1.fastq.gz
    cat ${r2.join(' ')} > ${sample_id}_R2.fastq.gz
    """
}

def fastqLaneKey(file) {
    def matcher = file.name =~ /^(.*)_L(\d{3})_(R[12])_(\d{3})\.fastq\.gz$/
    if (matcher.matches()) return [sample: matcher[0][1], lane: matcher[0][2], mate: matcher[0][3]]
    matcher = file.name =~ /^(.*)_(R[12])\.fastq\.gz$/
    if (matcher.matches()) return [sample: matcher[0][1], lane: '000', mate: matcher[0][2]]
    return null
}

def normalizeFastqPairs(files, sample) {
    def lanes = [:].withDefault { [:].withDefault { [] } }
    files.each { file ->
        def key = fastqLaneKey(file)
        if (key && key.sample == sample) lanes[key.lane][key.mate] << file
    }
    def sorted = lanes.keySet().sort()
    def lane_names = sorted.collect { lane ->
        if (lanes[lane].R1.size() != 1 || lanes[lane].R2.size() != 1)
            error "Sample '${sample}' lane '${lane}' must have exactly one R1 and one R2 FASTQ."
        [r1: lanes[lane].R1[0].name, r2: lanes[lane].R2[0].name]
    }
    tuple(lane_names, lanes.values().collectMany { it.values().flatten() })
}

workflow {

    fastq_input_pattern = params.input?.toString()?.trim()
    // Seqera may preserve shell-style escapes around glob characters in the text field.
    fastq_input_pattern = fastq_input_pattern?.replace('\\*', '*')?.replace('\\?', '?')
    if (!fastq_input_pattern) {
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

    analysis_target = (params.analysis_target ?: 'mtDNA').toString().trim()
    if (!(analysis_target in ['mtDNA', 'nDNA', 'both'])) {
        error "Unsupported --analysis_target '${analysis_target}'. Supported values: mtDNA, nDNA, both"
    }
    run_mt_analysis = analysis_target in ['mtDNA', 'both']
    run_nuclear_analysis = analysis_target in ['nDNA', 'both']

    nuclear_calling_mode = (params.nuclear_calling_mode ?: 'tumor_only').toString().trim().toLowerCase()
    if (!(nuclear_calling_mode in ['tumor_only', 'paired', 'germline'])) {
        error "Unsupported --nuclear_calling_mode '${nuclear_calling_mode}'. Supported values: tumor_only, paired, germline"
    }

    if (run_nuclear_analysis) {
        if (!params.nuclear_reference?.toString()?.trim() || !params.nuclear_target_bed?.toString()?.trim()) {
            error "nDNA analysis requires --nuclear_reference and --nuclear_target_bed"
        }
        if (!file(params.nuclear_reference).exists()) {
            error "nDNA reference FASTA does not exist: ${params.nuclear_reference}"
        }
        if (!file(params.nuclear_target_bed).exists()) {
            error "nDNA target BED does not exist: ${params.nuclear_target_bed}"
        }
        run_nuclear_annotation = params.nuclear_annotate == null || params.nuclear_annotate.toString().toBoolean()
        if (run_nuclear_annotation && !params.nuclear_vep_cache?.toString()?.trim()) {
            error "nDNA variant annotation requires --nuclear_vep_cache pointing to an offline Ensembl VEP cache"
        }
        if (run_nuclear_annotation && !file(params.nuclear_vep_cache).exists()) {
            error "nDNA VEP cache directory does not exist: ${params.nuclear_vep_cache}"
        }
        [
            ['nuclear_germline_resource', params.nuclear_germline_resource, 'nuclear_germline_resource_tbi', params.nuclear_germline_resource_tbi],
            ['nuclear_pon', params.nuclear_pon, 'nuclear_pon_tbi', params.nuclear_pon_tbi],
            ['nuclear_gnomad_vcf', params.nuclear_gnomad_vcf, 'nuclear_gnomad_tbi', params.nuclear_gnomad_tbi],
            ['nuclear_cosmic_vcf', params.nuclear_cosmic_vcf, 'nuclear_cosmic_tbi', params.nuclear_cosmic_tbi]
        ].each { name, path_value, index_name, index_value ->
            if (path_value?.toString()?.trim() && !index_value?.toString()?.trim()) {
                error "${name} requires its matching ${index_name} path"
            }
            if (!path_value?.toString()?.trim() && index_value?.toString()?.trim()) {
                error "${index_name} cannot be supplied without ${name}"
            }
            [[name, path_value], [index_name, index_value]].each { resource_name, resource_path ->
                if (resource_path?.toString()?.trim() && !file(resource_path).exists()) {
                    error "${resource_name} does not exist: ${resource_path}"
                }
            }
        }
        if (nuclear_calling_mode == 'paired' && !params.nuclear_sample_manifest?.toString()?.trim()) {
            error "Paired nDNA calling requires --nuclear_sample_manifest with sample_id,role,pair_id columns"
        }
        if (nuclear_calling_mode == 'paired' && !file(params.nuclear_sample_manifest).exists()) {
            error "nDNA sample manifest does not exist: ${params.nuclear_sample_manifest}"
        }
    }

    analysis_mode = normalizeAnalysisMode(params.analysis_mode ?: params.genome)
    if (run_mt_analysis && !analysis_mode) {
        error "Unsupported --analysis_mode '${params.analysis_mode ?: params.genome}'. Supported values: wgs_numt_correction, rcrs_only"
    }

    if (run_mt_analysis && (!params.analysis_modes || !params.analysis_modes[analysis_mode])) {
        error "Missing configuration for analysis mode '${analysis_mode}'. Check conf/genomes.config"
    }

    run_numt_correction = params.run_numt_correction != null && params.run_numt_correction.toString().toBoolean()
    run_numt_exploration = params.run_numt_exploration != null && params.run_numt_exploration.toString().toBoolean()
    reference_set = params.reference_set ?: (run_numt_correction ? 'GRCh38+rCRS' : 'rCRS')

    if (run_mt_analysis && run_numt_exploration && !run_numt_correction) {
        error "--run_numt_exploration is only available with --analysis_mode wgs_numt_correction"
    }

    if (run_mt_analysis && (!params.datadir || !params.mt_gsnap_db || !params.fasta || !params.mutserve_fasta || !params.mt_contig)) {
        error "Incomplete mitochondrial reference configuration for analysis mode '${analysis_mode}'"
    }

    if (run_mt_analysis && run_numt_correction && (!params.genomedir || !params.nuc_gsnap_db || !params.numt)) {
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

    if (run_mt_analysis && (params.annotation_high_depth as Integer) < (params.annotation_medium_depth as Integer)) {
        error "--annotation_high_depth must be >= --annotation_medium_depth"
    }

    if (run_mt_analysis && (params.annotation_high_alt_reads as Integer) < (params.annotation_medium_alt_reads as Integer)) {
        error "--annotation_high_alt_reads must be >= --annotation_medium_alt_reads"
    }

    mitomap_variant_table = params.mitomap_variant_table?.toString()?.trim()

    log.info """\
        C L A M - N F   P I P E L I N E
        ===============================
        input           : ${fastq_input_pattern}
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
        fastq_file_ch = Channel.fromPath(fastq_input_pattern, checkIfExists: true).collect()
        fastq_groups_ch = fastq_file_ch.flatMap { files ->
            def groups = [:].withDefault { [] }
            files.each { file ->
                def lane_key = fastqLaneKey(file)
                if (lane_key) groups[lane_key.sample] << file
            }
            groups.collect { sample, sample_files ->
                def clean_id = sanitizeSampleId(sample)
                if (!clean_id) error "Could not derive a valid sample ID from FASTQ files '${sample}'"
                def pair_data = normalizeFastqPairs(sample_files, sample)
                tuple(clean_id, pair_data[0], pair_data[1])
            }
        }
        reads_ch = merge_fastq_lanes(fastq_groups_ch)
            .map { sample_id, r1, r2 -> [sample_id, [r1, r2]] }
            .ifEmpty { error "No paired-end R1/R2 FASTQ files matched --input '${fastq_input_pattern}'." }
    }

    if (input_type == 'bam') {
        bam_ch = Channel
            .fromPath(fastq_input_pattern, checkIfExists: true)
            .ifEmpty { error "No BAM files matched --input '${fastq_input_pattern}'." }
            .map { bam ->
                def clean_id = sanitizeSampleId(bam.simpleName)
                if (!clean_id) {
                    error "Could not derive a valid sample ID from BAM '${bam}'"
                }
                [clean_id, bam]
            }

        if (run_mt_analysis) {
            sort_bam_ch = sorting_1(bam_ch)
            reads_ch = bam2fqgz(sort_bam_ch)
        }
    }

    if (run_nuclear_analysis) {
        nuclear_reference_input_ch = Channel.value(file(params.nuclear_reference, checkIfExists: true))
        nuclear_reference_index_ch = nuclear_reference_index(nuclear_reference_input_ch)
        nuclear_reference_ch = nuclear_reference_index_ch.map { reference, reference_fai, index_files -> tuple(reference, reference_fai) }
        nuclear_bed_ch = Channel.value(file(params.nuclear_target_bed, checkIfExists: true))

        if (input_type == 'fastq') {
            nuclear_bams_ch = nuclear_alignment(reads_ch, nuclear_reference_index_ch)
        } else {
            nuclear_bams_ch = nuclear_prepare_bam(bam_ch)
        }

        nuclear_call_input_ch = nuclear_bams_ch
            .combine(nuclear_reference_ch)
            .combine(nuclear_bed_ch)
            .map { sample_id, bam, bai, reference, reference_fai, target_bed -> tuple(sample_id, bam, bai, reference, reference_fai, target_bed) }
        nuclear_germline_resource_ch = Channel.value(params.nuclear_germline_resource ? [file(params.nuclear_germline_resource, checkIfExists: true), file(params.nuclear_germline_resource_tbi, checkIfExists: true)] : [])
        nuclear_pon_ch = Channel.value(params.nuclear_pon ? [file(params.nuclear_pon, checkIfExists: true), file(params.nuclear_pon_tbi, checkIfExists: true)] : [])

        if (nuclear_calling_mode == 'germline') {
            nuclear_variants_ch = nuclear_germline_call(nuclear_call_input_ch)
            nuclear_cnv_input_ch = nuclear_bams_ch
                .combine(nuclear_reference_ch)
                .combine(nuclear_bed_ch)
                .map { sample_id, bam, bai, reference, reference_fai, target_bed -> tuple(sample_id, bam, bai, reference, target_bed) }
            nuclear_cnv_ch = nuclear_cnvkit_flat(nuclear_cnv_input_ch)
        } else if (nuclear_calling_mode == 'paired') {
            manifest_rows = file(params.nuclear_sample_manifest).text.readLines()
                .findAll { it.trim() && !it.trim().startsWith('#') }
            if (manifest_rows.size() < 2) {
                error "--nuclear_sample_manifest must contain a header and at least one sample pair"
            }
            manifest_header = manifest_rows[0].split('\t', -1).collect { it.trim().toLowerCase() }
            required_manifest_columns = ['sample_id', 'role', 'pair_id']
            if (!required_manifest_columns.every { manifest_header.contains(it) }) {
                error "--nuclear_sample_manifest must be tab-delimited with sample_id, role, and pair_id columns"
            }
            manifest = manifest_rows.drop(1).collect { line ->
                def values = line.split('\t', -1).collect { it.trim() }
                [sample_id: values[manifest_header.indexOf('sample_id')], role: values[manifest_header.indexOf('role')].toLowerCase(), pair_id: values[manifest_header.indexOf('pair_id')]]
            }
            if (manifest.any { !it.sample_id || !it.pair_id || !(it.role in ['tumor', 'normal']) }) {
                error "Each manifest row needs a sample_id, pair_id, and role of tumor or normal"
            }
            pair_records = manifest.groupBy { it.pair_id }.collect { pair_id, rows ->
                def tumors = rows.findAll { it.role == 'tumor' }
                def normals = rows.findAll { it.role == 'normal' }
                if (tumors.size() != 1 || normals.size() != 1) {
                    error "Paired manifest group '${pair_id}' must contain exactly one tumor and one normal"
                }
                [pair_id: pair_id, tumor_id: tumors[0].sample_id, normal_id: normals[0].sample_id]
            }
            nuclear_pair_bams_ch = nuclear_bams_ch.collect().flatMap { rows ->
                def bam_by_sample = rows.collectEntries { row -> [(row[0]): [bam: row[1], bai: row[2]] ] }
                pair_records.collect { pair ->
                    if (!bam_by_sample[pair.tumor_id] || !bam_by_sample[pair.normal_id]) {
                        error "Manifest pair '${pair.pair_id}' references samples absent from the input: ${pair.tumor_id}, ${pair.normal_id}"
                    }
                    tuple(pair.pair_id, pair.tumor_id, bam_by_sample[pair.tumor_id].bam, bam_by_sample[pair.tumor_id].bai,
                        pair.normal_id, bam_by_sample[pair.normal_id].bam, bam_by_sample[pair.normal_id].bai)
                }
            }
            nuclear_pair_call_input_ch = nuclear_pair_bams_ch
                .combine(nuclear_reference_ch)
                .combine(nuclear_bed_ch)
                .map { pair_id, tumor_id, tumor_bam, tumor_bai, normal_id, normal_bam, normal_bai, reference, reference_fai, target_bed ->
                    tuple(pair_id, tumor_id, tumor_bam, tumor_bai, normal_id, normal_bam, normal_bai, reference, reference_fai, target_bed)
                }
            nuclear_variants_ch = nuclear_paired_call(nuclear_pair_call_input_ch, nuclear_germline_resource_ch, nuclear_pon_ch)
            nuclear_cnv_input_ch = nuclear_pair_bams_ch
                .combine(nuclear_reference_ch)
                .combine(nuclear_bed_ch)
                .map { pair_id, tumor_id, tumor_bam, tumor_bai, normal_id, normal_bam, normal_bai, reference, reference_fai, target_bed ->
                    tuple(pair_id, tumor_id, tumor_bam, tumor_bai, normal_id, normal_bam, normal_bai, reference, reference_fai, target_bed)
                }
            nuclear_cnv_ch = nuclear_cnvkit_paired(nuclear_cnv_input_ch)
        } else {
            nuclear_variants_ch = nuclear_tumor_only_call(nuclear_call_input_ch, nuclear_germline_resource_ch, nuclear_pon_ch)
            nuclear_cnv_input_ch = nuclear_bams_ch
                .combine(nuclear_reference_ch)
                .combine(nuclear_bed_ch)
                .map { sample_id, bam, bai, reference, reference_fai, target_bed -> tuple(sample_id, bam, bai, reference, reference_fai, target_bed) }
            nuclear_cnv_ch = nuclear_cnvkit_flat(nuclear_cnv_input_ch)
        }

        if (run_nuclear_annotation) {
            nuclear_annotated_ch = nuclear_variants_ch
                .combine(nuclear_reference_ch)
                .map { sample_id, unfiltered_vcf, filtered_vcf, reference, reference_fai -> tuple(sample_id, filtered_vcf, reference, reference_fai) }
            nuclear_vep_cache_ch = Channel.value(file(params.nuclear_vep_cache, checkIfExists: true))
            nuclear_gnomad_vcf_ch = Channel.value(params.nuclear_gnomad_vcf ? [file(params.nuclear_gnomad_vcf, checkIfExists: true), file(params.nuclear_gnomad_tbi, checkIfExists: true)] : [])
            nuclear_cosmic_vcf_ch = Channel.value(params.nuclear_cosmic_vcf ? [file(params.nuclear_cosmic_vcf, checkIfExists: true), file(params.nuclear_cosmic_tbi, checkIfExists: true)] : [])
            nuclear_annotated_ch = nuclear_annotate(nuclear_annotated_ch, nuclear_vep_cache_ch, nuclear_gnomad_vcf_ch, nuclear_cosmic_vcf_ch)
            nuclear_annotated_ch
                .join(nuclear_cnv_ch, by: 0)
                .map { sample_id, annotated_vcf, maf, cnvkit_out -> tuple(sample_id, nuclear_calling_mode, annotated_vcf, maf, cnvkit_out) }
                .set { nuclear_report_input_ch }
            nuclear_report(nuclear_report_input_ch)
        } else {
            log.info 'nDNA annotation disabled: VEP and MAF generation will be skipped.'
        }
    }

    if (run_mt_analysis) {
    bams_ch = alignment_MT(reads_ch)
    md_ch = get_the_MD(bams_ch)

    if (run_numt_correction) {
        md_mitoscape_ch = get_the_MD_mitoscape(bams_ch)
        fq_ch = bam2fq(bams_ch)
        nuc_ch = alignment_Nuc(fq_ch)

        md_mitoscape_ch
            .join(nuc_ch, by: 0)
            .set { mito_ch }

        mito_ch
            .join(bams_ch, by: 0)
            .set { final_mito_ch }

        mt_ch = mitoscape(final_mito_ch)

        md_ch
            .join(mt_ch, by: 0)
            .set { sorting_input_ch }

        sort_ch = sorting(sorting_input_ch)

        sort_ch[0]
            .join(sort_ch[1], by: 0)
            .set { sorted_pair_ch }

        cov_ch = coverage(sorted_pair_ch)
        final_coverage_ch = cov_ch[1]
        ind_ch = index(sort_ch[1])

        sorted_pair_ch
            .map { sample_id, before_bam, after_bam -> [sample_id, 'before,after', [before_bam, after_bam]] }
            .set { report_bam_input_ch }

        cov_ch[0]
            .join(cov_ch[1], by: 0)
            .set { report_coverage_ch }

        sort_ch[1]
            .join(ind_ch, by: 0)
            .set { var_ch }

        if (run_numt_exploration) {
            sorted_pair_ch
                .join(nuc_ch, by: 0)
                .set { numt_exploration_input_ch }

            numt_exploration_raw_ch = numt_like_exploration(numt_exploration_input_ch)
        }
    } else {
        sort_single_ch = sorting_single(md_ch)
        cov_ch = coverage_single(sort_single_ch)
        final_coverage_ch = cov_ch[1]
        ind_ch = index(sort_single_ch)

        sort_single_ch
            .map { sample_id, sorted_bam -> [sample_id, 'final', [sorted_bam]] }
            .set { report_bam_input_ch }

        cov_ch[0]
            .join(cov_ch[1], by: 0)
            .map { sample_id, coverage_before, coverage_final -> [sample_id, coverage_before, coverage_final] }
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


}
