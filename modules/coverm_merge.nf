process COVERM_MERGE {
    label 'process_low'

    container 'python:3.12'

    publishDir "${params.outdir}/quantify", mode: 'copy'

    input:
    // Stage files with their original names (EV01_coverage.tsv, ...).
    // A `stageAs: "input/*_coverage.tsv"` would substitute `*` with a numeric
    // index when given multiple files, which would strip the sample id.
    path coverages

    output:
    path "votu_count_matrix.tsv", emit: count
    path "votu_tpm_matrix.tsv",   emit: tpm
    path "votu_long.tsv",         emit: long

    script:
    """
    set -euo pipefail
    coverm_merge.py \\
        --indir . \\
        --out-count votu_count_matrix.tsv \\
        --out-tpm   votu_tpm_matrix.tsv \\
        --out-long  votu_long.tsv
    """
}
