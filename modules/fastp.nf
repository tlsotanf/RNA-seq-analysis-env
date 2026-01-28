process FASTP {
    
    tag "$meta.id"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.json")            , emit: json
    tuple val(meta), path("*.html")            , emit: html

    script:
    def prefix = meta.id
    
    // Single-end 체크 (reads가 1개만 있으면 single-end)
    if (reads instanceof List && reads.size() == 1 || meta.single_end) {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${prefix}.trimmed.fastq.gz \\
            --thread ${params.threads} \\
            --cut_front \\
            --cut_tail \\
            --cut_mean_quality 20 \\
            --length_required 30 \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${prefix}_1.trimmed.fastq.gz \\
            --out2 ${prefix}_2.trimmed.fastq.gz \\
            --thread ${params.threads} \\
            --cut_front \\
            --cut_tail \\
            --cut_mean_quality 20 \\
            --length_required 30 \\
            --detect_adapter_for_pe \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html
        """
    }
}
