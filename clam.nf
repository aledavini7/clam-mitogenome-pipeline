/*
 * pipeline input parameters
 */

// Example:
// nextflow run clam.nf --input 'data/*_{R1,R2}.fastq.gz' --input_type fastq
// nextflow run clam.nf --input 'data/*.bam' --input_type bam --genome GRCh38

nextflow.enable.dsl = 2

def normalizeGenomeName(value) {
    def requested = value?.toString()
    if (!requested) {
        return null
    }

    def aliases = [
        'grch38': 'GRCh38',
        'hg38'  : 'GRCh38',
        'grch37': 'GRCh37',
        'hg19'  : 'GRCh37',
        'rcrs'  : 'rCRS'
    ]

    return aliases[requested.toLowerCase()]
}

include { sorting_1; bam2fqgz } from './modules/1_data_preparation.nf'
include { alignment_MT; get_the_MD } from './modules/2_MT_alignment.nf'
include { bam2fq; alignment_Nuc } from './modules/3_NUC_alignment.nf'
include { mitoscape } from './modules/4_Mitoscape.nf'
include { sorting; coverage; sorting_single; coverage_single } from './modules/5_Coverage.nf'
include { index } from './modules/6_Index.nf'
include { mutect2; filter_mutect2; bgzip; mutserve } from './modules/7_Variant_calling.nf'
include { haplogrep; haplogrep1 } from './modules/8_Haplogroups_assignment.nf'
include { merge_variant_calls } from './modules/9_annotation.nf'

workflow {

    if (!params.input) {
        error "Missing required parameter: --input"
    }

    input_type = params.input_type?.toString()?.toLowerCase()
    if (!input_type) {
        error "Missing required parameter: --input_type. Supported values: fastq, bam"
    }

    if (!['fastq', 'bam'].contains(input_type)) {
        error "Unsupported --input_type '${params.input_type}'. Supported values: fastq, bam"
    }

    genome = normalizeGenomeName(params.genome)
    if (!genome) {
        error "Unsupported --genome '${params.genome}'. Supported values: GRCh38, GRCh37, rCRS"
    }

    if (!params.genomes || !params.genomes[genome]) {
        error "Missing configuration for genome '${genome}'. Check conf/genomes.config"
    }

    run_numt_correction = params.run_numt_correction != null && params.run_numt_correction.toString().toBoolean()

    if (!params.datadir || !params.mt_gsnap_db || !params.fasta || !params.mutserve_fasta || !params.mt_contig) {
        error "Incomplete mitochondrial reference configuration for genome '${genome}'"
    }

    if (run_numt_correction && (!params.genomedir || !params.nuc_gsnap_db || !params.numt)) {
        error "Incomplete nuclear/NUMTs reference configuration for genome '${genome}'"
    }

    log.info """\
        C L A M - N F   P I P E L I N E
        ===============================
        input           : ${params.input}
        input_type      : ${input_type}
        genome          : ${genome}
        outdir          : ${params.outdir}
        mt_gsnap_dir    : ${params.datadir}
        mt_gsnap_db     : ${params.mt_gsnap_db}
        mt_fasta        : ${params.fasta}
        mt_contig       : ${params.mt_contig}
        numt_correction : ${run_numt_correction}
        nuc_gsnap_dir   : ${params.genomedir ?: 'not used'}
        nuc_gsnap_db    : ${params.nuc_gsnap_db ?: 'not used'}
        numt_regions    : ${params.numt ?: 'not used'}
        """
        .stripIndent()

    if (input_type == 'fastq') {
        reads_ch = Channel
            .fromFilePairs(params.input, size: 2, checkIfExists: true)
            .map { sample_id, reads -> [sample_id.toString(), reads] }
    }

    if (input_type == 'bam') {
        bam_ch = Channel
            .fromPath(params.input, checkIfExists: true)
            .map { bam -> [bam.simpleName, bam] }

        sort_bam_ch = sorting_1(bam_ch)
        reads_ch = bam2fqgz(sort_bam_ch)
    }

    bams_ch = alignment_MT(reads_ch)
    md_ch = get_the_MD(bams_ch)

    if (run_numt_correction) {
        fq_ch = bam2fq(bams_ch)
        nuc_ch = alignment_Nuc(fq_ch)

        md_ch
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

        sort_ch[1]
            .combine(ind_ch, by: 0)
            .set { var_ch }
    } else {
        sort_single_ch = sorting_single(md_ch)
        cov_ch = coverage_single(sort_single_ch)
        final_coverage_ch = cov_ch
        ind_ch = index(sort_single_ch)

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
    
    gz_ch = bgzip(filtered_mutect2_ch[0])
    
    hplg_ch = haplogrep(gz_ch)
    hplg_2_ch = haplogrep1(vcf_2_ch[0])


}
