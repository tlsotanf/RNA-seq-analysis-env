process SALMON_QUANT {
    
    tag "$meta.id"
    publishDir "${params.outdir}/salmon", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(reads)
    path salmon_index

    output:
    tuple val(meta), path("${meta.id}"), emit: results, optional: true

    script:
    def prefix = meta.id
    def read_files = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def orphans = meta.single_end ? "" : "--allowOrphans"
    """
    salmon quant \\
        -i ${salmon_index} \\
        -l A \\
        ${read_files} \\
        -p ${params.threads} \\
        --validateMappings \\
        --minAssignedFrags 1 \\
        ${orphans} \\
        -o ${prefix}
    """
}