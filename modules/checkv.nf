process CHECKV {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/checkv:1.0.3--pyhdfd78af_0'

    publishDir path: { "${params.outdir}/identify/checkv/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(virus_fasta)
    path  checkv_db

    output:
    tuple val(meta), path("checkv/quality_summary.tsv"), emit: quality
    tuple val(meta), path("checkv/viruses.fna"),         emit: viruses
    tuple val(meta), path("checkv/proviruses.fna"),      emit: proviruses
    tuple val(meta), path("checkv/contamination.tsv"),   emit: contamination,  optional: true
    tuple val(meta), path("checkv/completeness.tsv"),    emit: completeness,   optional: true
    path "versions.yml",                                 emit: versions

    script:
    """
    set -euo pipefail

    checkv end_to_end \\
        ${virus_fasta} \\
        checkv \\
        -d ${checkv_db} \\
        -t ${task.cpus}

    # CheckV always emits viruses.fna; touch proviruses.fna if the run produced none
    # so downstream channels never see a missing file.
    [ -f checkv/proviruses.fna ] || touch checkv/proviruses.fna
    [ -f checkv/viruses.fna ]    || touch checkv/viruses.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv --version 2>&1 | sed 's/^.*version //')
    END_VERSIONS
    """
}
