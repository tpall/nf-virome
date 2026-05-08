process VOTU_CLUSTER {
    label 'process_low'

    // Stdlib-only script. Use the full python image (Debian-based) so `ps`
    // from procps is available — Nextflow needs it for task metrics, and the
    // slim variant omits it. Singularity auto-prefixes `docker://`.
    container 'python:3.12'

    publishDir "${params.outdir}/cluster", mode: 'copy'

    input:
    path catalog_fasta
    path skani_ani
    val  min_ani
    val  min_af

    output:
    path "votu_clusters.tsv", emit: clusters
    path "votu_reps.txt",     emit: reps_ids

    script:
    """
    set -euo pipefail
    votu_cluster.py \\
        --fasta ${catalog_fasta} \\
        --skani ${skani_ani} \\
        --ani ${min_ani} \\
        --af ${min_af} \\
        --out-clusters votu_clusters.tsv \\
        --out-reps votu_reps.txt
    """
}
