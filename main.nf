#!/usr/bin/env nextflow
/*
 * bcl-qc: Illumina BCL demux, Dragen alignment, and QC summary for clinical NGS.
 * Nextflow DSL2 port of ucladx/bcl-qc bclqc.py (github.com/ucladx/bcl-qc).
 *
 * Usage:
 *   nextflow run . --run_dir <path> --sampleinfo <file.tsv> [--steps demux,align,qc]
 */

nextflow.enable.dsl = 2

// --- Parameter defaults (override on CLI or via -params-file) ---------------
params.run_dir     = null
params.sampleinfo  = null
params.fastqs_dir  = '/staging/hot/reads'
params.bams_dir    = '/mnt/pns/bams'
params.steps       = 'demux,align,qc'

// Index flavor -> Dragen panel. Any flavor absent from this map is skipped
// (e.g. I8N2_N10 exome samples are not aligned by this pipeline).
params.flavor_to_panel = [
    'I10'      : 'PCP',
    'U5N2_I10' : 'HEME',
]

// --- Helpers ---------------------------------------------------------------
def dragenBin() {
    def bin = params.dragen_versions[params.dragen_version]
    if( !bin || !file(bin).exists() )
        error "Dragen ${params.dragen_version} not found at '${bin}'. Update conf/panels.config or pass --dragen_version."
    return bin
}

def flavorFromSamplesheet(String ss) {
    // "SampleSheet_I10.csv" -> "I10"  (matches bclqc.py get_index)
    return file(ss).name.replaceFirst(/^SampleSheet_/, '').replaceFirst(/\.csv$/, '')
}

// Heme subpanels inherit heme_comp and override target_intervals
// (matches bclqc.py QCSumInfo.get_config). YAML is small, load-on-call is fine.
def resolveQcsumCfg(String panel_key) {
    def y = new org.yaml.snakeyaml.Yaml().load(file("${projectDir}/config/qcsum_config.yaml").text) as Map
    if( !y.containsKey(panel_key) )
        error "Panel '${panel_key}' not in config/qcsum_config.yaml"
    if( panel_key.startsWith('heme_') && panel_key != 'heme_comp' ) {
        def cfg = new LinkedHashMap(y.heme_comp)
        def sub = y[panel_key] ?: [:]
        if( sub.target_intervals ) cfg.target_intervals = sub.target_intervals
        return cfg
    }
    return y[panel_key]
}

// ============================================================================
// PROCESSES
// ============================================================================

process BCL_CONVERT {
    tag "${flavor}"
    publishDir "${params.fastqs_dir}/${run_name}", mode: 'copy'

    input:
        tuple val(flavor), path(samplesheet)
        val   run_dir_abs
        val   run_name

    output:
        tuple val(flavor), path("${flavor}/Reports/fastq_list.csv"), emit: fastq_list
        path "${flavor}",                                            emit: outdir

    script:
        def dragen = dragenBin()
        """
        ${dragen} \\
            --bcl-conversion-only true \\
            --bcl-use-hw false \\
            --bcl-only-matched-reads true \\
            --bcl-input-directory ${run_dir_abs} \\
            --sample-sheet ${samplesheet} \\
            --output-directory ${flavor}
        """

    stub:
        """
        mkdir -p ${flavor}/Reports
        echo "RGID,RGSM,RGLB,Lane,Read1File,Read2File" > ${flavor}/Reports/fastq_list.csv
        echo "stub.1,STUB_SAMPLE,stub,1,stub_R1.fastq.gz,stub_R2.fastq.gz" >> ${flavor}/Reports/fastq_list.csv
        """
}

process DRAGEN_ALIGN {
    tag "${sample_id}/${panel}"
    publishDir "${params.bams_dir}/${run_name}/${sample_id}", mode: 'copy'

    input:
        tuple val(sample_id), val(panel), path(fastq_list)
        val   run_name

    output:
        tuple val(sample_id), val(panel), path("${sample_id}.cram"), emit: cram
        path "*",                                                     emit: all

    script:
        def dragen   = dragenBin()
        def cfg      = params.panels[panel]
        def bed      = cfg.bed
        def ref      = cfg.ref
        def heme_ext = panel == 'HEME' ?
            '--umi-enable true --umi-source qname --umi-correction-scheme random ' +
            '--umi-min-supporting-reads 1 ' +
            "--umi-metrics-interval-file ${bed} " +
            '--vc-enable-umi-germline true --vc-enable-high-sensitivity-mode true' :
            '--enable-duplicate-marking true'
        """
        ${dragen} \\
            --intermediate-results-dir ${params.intermediate_results_dir} \\
            --enable-map-align true --enable-map-align-output true \\
            --output-format CRAM --generate-sa-tags true --enable-sort true \\
            --soft-read-trimmers polyg,quality --trim-min-quality 2 \\
            --ref-dir ${ref} \\
            --qc-coverage-tag-1 target_bed --qc-coverage-region-1 ${bed} \\
            --qc-coverage-reports-1 cov_report --qc-coverage-ignore-overlaps true \\
            --enable-variant-caller true --vc-combine-phased-variants-distance 6 \\
            --vc-emit-ref-confidence GVCF --enable-hla true \\
            --fastq-list ${fastq_list} --fastq-list-sample-id ${sample_id} \\
            --output-directory . --output-file-prefix ${sample_id} \\
            ${heme_ext}
        """

    stub:
        """
        touch ${sample_id}.cram ${sample_id}.cram.crai ${sample_id}.hard-filtered.vcf.gz
        """
}

process SAVE_OCC_PF_PLOT {
    tag "${run_name}"
    publishDir "${params.bams_dir}/${run_name}", mode: 'copy'

    input:
        val run_dir_abs
        val run_name

    output:
        path "occ_pf_lane_mqc.jpg"

    script:
        """
        save_occ_pf_plot.py --run-dir ${run_dir_abs} --output-dir .
        """

    stub:
        "touch occ_pf_lane_mqc.jpg"
}

process PICARD_QCSUM {
    tag "${sample_id}/${panel_key}"
    publishDir "${params.bams_dir}/${run_name}/${sample_id}", mode: 'copy'

    input:
        // cfg is a resolved map from qcsum_config.yaml (see resolveQcsumCfg in workflow)
        tuple val(sample_id), val(panel_key), path(bam), val(cfg)
        val   run_name

    output:
        path "${sample_id}.hsm.txt",   emit: hsm
        path "${sample_id}.qcsum.txt", emit: qcsum

    script:
        """
        picard CollectHsMetrics \\
            I=${bam} \\
            O=${sample_id}.hsm.txt \\
            R=${params.picard_ref} \\
            BAIT_INTERVALS=${cfg.bait_intervals} \\
            TARGET_INTERVALS=${cfg.target_intervals}

        perl ${projectDir}/bin/qcsum_metrics.pl \\
            --prefix ${sample_id} \\
            --qcfolder . \\
            --pipeline_version ${cfg.pipeline_version} \\
            --platform ${cfg.platform} \\
            --pass_min_align_pct ${cfg.pass_min_align_pct} \\
            --fail_min_align_pct ${cfg.fail_min_align_pct} \\
            --covered ${cfg.covered} \\
            --pass_min_roi_pct ${cfg.pass_min_roi_pct} \\
            --fail_min_roi_pct ${cfg.fail_min_roi_pct} \\
            --pass_min_avgcov ${cfg.pass_min_avgcov} \\
            --fail_min_avgcov ${cfg.fail_min_avgcov} \\
            --pass_min_reads ${cfg.pass_min_reads} \\
            --fail_min_reads ${cfg.fail_min_reads} \\
            --capture ${cfg.capture} \\
            --capture_version ${cfg.capture_version}
        """

    stub:
        """
        touch ${sample_id}.hsm.txt
        echo "Sample" > ${sample_id}.qcsum.txt
        echo "${sample_id},NovaSeq6000,vSTUB,PASS,PASS,0,0,stub,0,0,0,0,0,0,0,0,0,0,0,0" >> ${sample_id}.qcsum.txt
        """
}

process AGGREGATE_QCSUM {
    tag "${run_name}"
    publishDir "${params.bams_dir}/${run_name}", mode: 'copy'

    input:
        path qcsum_files
        val  run_name

    output:
        path "qcsum_mqc.csv"

    script:
        """
        aggregate_qcsum.py --out qcsum_mqc.csv ${qcsum_files}
        """

    stub:
        "cp ${projectDir}/assets/qcsum_mqc_header.txt qcsum_mqc.csv"
}

process MULTIQC {
    tag "${run_name}"
    publishDir "${params.bams_dir}/${run_name}", mode: 'copy'

    input:
        val run_name
        path extras  // qcsum_mqc.csv, occ_pf_lane_mqc.jpg, etc. (staged alongside bams/fastqs)

    output:
        path "multiqc_report.html"
        path "multiqc_data"

    script:
        """
        # Drop wgs coverage reports; irrelevant for targeted panels (matches bclqc.py)
        find ${params.bams_dir}/${run_name} -name "*.wgs_*.csv" -delete || true
        multiqc --force --config ${projectDir}/config/multiqc_config.yaml --outdir . \\
            ${params.bams_dir}/${run_name} \\
            ${params.fastqs_dir}/${run_name} \\
            .
        """

    stub:
        """
        touch multiqc_report.html
        mkdir -p multiqc_data && touch multiqc_data/.keep
        """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow {
    if( !params.run_dir ) error "Missing --run_dir"
    def run_dir_abs = file(params.run_dir).toAbsolutePath().toString()
    def run_name    = file(params.run_dir).name
    def steps       = params.steps.tokenize(',')

    // ---------------- demux ------------------------------------------------
    def fastq_list_ch
    if( 'demux' in steps ) {
        def samplesheets = Channel.fromPath("${run_dir_abs}/SampleSheet_*.csv")
            .ifEmpty { error "No SampleSheet_*.csv in ${run_dir_abs}" }
            .map { ss -> tuple(flavorFromSamplesheet(ss.toString()), ss) }

        fastq_list_ch = BCL_CONVERT(samplesheets, run_dir_abs, run_name).fastq_list
    } else {
        // Resume or re-enter: reuse previously published fastq_lists.
        fastq_list_ch = Channel.fromPath("${params.fastqs_dir}/${run_name}/*/Reports/fastq_list.csv")
            .map { fl -> tuple(file(fl).parent.parent.name, fl) }
    }

    // ---------------- align ------------------------------------------------
    if( 'align' in steps ) {
        def align_ch = fastq_list_ch
            .filter { flavor, _fl -> params.flavor_to_panel.containsKey(flavor) }
            .flatMap { flavor, fl ->
                def panel = params.flavor_to_panel[flavor]
                def rows  = file(fl).splitCsv(header: true)
                rows.collect { tuple(it.RGSM, panel, fl) }.unique { it[0] }
            }
        DRAGEN_ALIGN(align_ch, run_name)
    }

    // ---------------- qc ---------------------------------------------------
    if( 'qc' in steps ) {
        if( !params.sampleinfo ) error "Missing --sampleinfo (required for qc step)"

        SAVE_OCC_PF_PLOT(run_dir_abs, run_name)

        def sampleinfo_to_qcsum = [
            'Comprehensive Heme Panel'          : 'heme_comp',
            'MPN Screen Panel (JAK2, CALR, MPL)': 'heme_mpn',
            'JAK2'                              : 'heme_jak2',
            'Peripheral Blood Lymphoma Panel'   : 'heme_pblp',
            'Tumor'                             : 'pcp_tumor',
            'Normal'                            : 'pcp_normal',
        ]
        def qcsum_in = Channel.fromPath(params.sampleinfo)
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def raw   = row.Panel ?: row.Tumor
                def panel = sampleinfo_to_qcsum[raw] ?: raw
                tuple(row.Samples, panel, file(row['BAM Path']), resolveQcsumCfg(panel))
            }

        def picard = PICARD_QCSUM(qcsum_in, run_name)
        def agg    = AGGREGATE_QCSUM(picard.qcsum.collect(), run_name)

        MULTIQC(run_name, agg.mix(SAVE_OCC_PF_PLOT.out).collect())
    }
}
