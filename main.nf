#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input         = null
params.outdir        = 'results'
params.threads       = 4
params.star_index    = null
params.salmon_index  = null
params.skip_trimming = false
params.skip_qc       = false

if (!params.input) {
    error "❌ Please specify --input samplesheet.csv"
}

include { FASTQC as FASTQC_RAW     } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { FASTP                     } from './modules/fastp'
include { STAR_ALIGN                } from './modules/star'
include { SALMON_QUANT              } from './modules/salmon'
include { MULTIQC                   } from './modules/multiqc'

workflow {
    
    ch_input = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            
            if (row.containsKey('fastq_2') && row.fastq_2 && row.fastq_2 != '') {
                meta.single_end = false
                [meta, [file(row.fastq_1), file(row.fastq_2)]]
            } else {
                meta.single_end = true
                [meta, file(row.fastq_1)]
            }
        }

    if (!params.skip_qc) {
        FASTQC_RAW(ch_input)
        ch_raw_fastqc = FASTQC_RAW.out.zip
    }

    if (!params.skip_trimming) {
        FASTP(ch_input)
        ch_trimmed = FASTP.out.reads
        ch_fastp_json = FASTP.out.json
    } else {
        ch_trimmed = ch_input
    }

    if (!params.skip_qc && !params.skip_trimming) {
        FASTQC_TRIMMED(ch_trimmed)
        ch_trimmed_fastqc = FASTQC_TRIMMED.out.zip
    }

    if (params.star_index) {
        STAR_ALIGN(
            ch_trimmed,
            Channel.fromPath("${params.star_index}/*").collect()
        )
        ch_bam = STAR_ALIGN.out.bam
        ch_star_log = STAR_ALIGN.out.log
    }

    if (params.salmon_index) {
        SALMON_QUANT(
            ch_trimmed,
            file(params.salmon_index)
        )
        ch_salmon_results = SALMON_QUANT.out.results
    }

    ch_multiqc_files = Channel.empty()
    if (!params.skip_qc) {
        if (!params.skip_trimming) {
            ch_multiqc_files = ch_multiqc_files.mix(
                FASTQC_RAW.out.zip.map { meta, file -> file }.collect().ifEmpty([]),
                FASTP.out.json.map { meta, file -> file }.collect().ifEmpty([]),
                FASTQC_TRIMMED.out.zip.map { meta, file -> file }.collect().ifEmpty([])
            )
        }
        
        MULTIQC(ch_multiqc_files.collect())
    }
}

workflow.onComplete {
    println """
    ====================================
    Pipeline completed!
    ====================================
    Status:   ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Results:  ${params.outdir}
    Duration: ${workflow.duration}
    ====================================
    """.stripIndent()
}