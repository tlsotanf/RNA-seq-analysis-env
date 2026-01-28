process MULTIQC {
    
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"        , emit: data

    script:
    """
    multiqc . -f
    """
}
