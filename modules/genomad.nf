process GENOMAD {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/genomad:1.11.0--pyhdfd78af_0'

    publishDir path: { "${params.outdir}/identify/genomad/${meta.id}" }, mode: 'copy'
    // Flat copy of the gene-level TSV for DRAM-v's --genomad_genes glob.
    // Drops the ${meta.id}_summary/ prefix so the file lands directly under
    // dramv_input/genomad_genes/<sample>_virus_genes.tsv.
    publishDir path: { "${params.outdir}/dramv_input/genomad_genes" },
               mode: 'copy', pattern: "**/*_virus_genes.tsv",
               saveAs: { fn -> file(fn).name }

    input:
    tuple val(meta), path(contigs)
    path  genomad_db

    output:
    tuple val(meta), path("${meta.id}_summary/${meta.id}_virus.fna"),                   emit: virus_fasta
    tuple val(meta), path("${meta.id}_summary/${meta.id}_virus_summary.tsv"),           emit: virus_summary
    tuple val(meta), path("${meta.id}_summary/${meta.id}_virus_genes.tsv"),             emit: virus_genes,    optional: true
    tuple val(meta), path("${meta.id}_summary/${meta.id}_virus_proteins.faa"),          emit: virus_proteins, optional: true
    tuple val(meta), path("${meta.id}_summary/${meta.id}_taxonomy.tsv"),                emit: virus_taxonomy, optional: true
    path "versions.yml",                                                                emit: versions

    script:
    def args = task.ext.args ?: '--cleanup'
    """
    set -euo pipefail

    # Stage input as <sample>.fna so geNomad's auto-naming uses our sample id.
    # geNomad accepts gzipped input; staging unzipped keeps file basenames predictable.
    if [[ "${contigs}" == *.gz ]]; then
        zcat ${contigs} > ${meta.id}.fna
    else
        cp ${contigs} ${meta.id}.fna
    fi

    genomad end-to-end \\
        ${meta.id}.fna \\
        . \\
        ${genomad_db} \\
        --threads ${task.cpus} \\
        ${args}

    rm -f ${meta.id}.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(genomad --version 2>&1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
