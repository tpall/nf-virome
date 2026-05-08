process EXTRACT_REPS {
    label 'process_low'

    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    publishDir "${params.outdir}/cluster", mode: 'copy'

    input:
    path catalog_fasta
    path reps_ids

    output:
    path "votu_catalog.fa", emit: catalog

    script:
    """
    set -euo pipefail
    seqkit grep -f ${reps_ids} ${catalog_fasta} > votu_catalog.fa
    """
}
