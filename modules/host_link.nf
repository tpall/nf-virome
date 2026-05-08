process HOST_LINK {
    label 'process_low'

    container 'python:3.12'

    publishDir "${params.outdir}/host_couple", mode: 'copy'

    input:
    path blast_hits
    path spacer_mapping

    output:
    path "votu_host_links.tsv",   emit: hits
    path "votu_host_summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    host_link.py \\
        --blast       ${blast_hits} \\
        --mapping     ${spacer_mapping} \\
        --out-hits    votu_host_links.tsv \\
        --out-summary votu_host_summary.tsv
    """
}
