process STAR_ALIGN {
    
    tag "$meta.id"
    publishDir "${params.outdir}/star", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path star_index

    output:
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("*.Log.final.out")                , emit: log
    tuple val(meta), path("*.SJ.out.tab")                   , emit: junction

    script:
    def prefix = meta.id
    """
    STAR \\
        --genomeDir ${star_index} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --readFilesCommand zcat \\
        --runThreadN ${params.threads} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --quantMode TranscriptomeSAM \\
        --outFileNamePrefix ${prefix}.
    """
}
