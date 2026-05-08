process CONCAT_FASTAS {
    label 'process_low'

    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    publishDir "${params.outdir}/cluster", mode: 'copy'

    input:
    path fastas, stageAs: "input/*"

    output:
    path "all_filtered.fna", emit: combined

    script:
    """
    set -euo pipefail
    : > all_filtered.fna
    # Each input is named <sample>_filtered.fna; prefix every contig header with
    # the sample id so the concatenated catalog has globally unique sequence ids.
    for f in input/*_filtered.fna; do
        sample=\$(basename "\$f" _filtered.fna)
        if [ -s "\$f" ]; then
            awk -v s="\$sample" '/^>/ {sub(/^>/, ">"s"_"); print; next} {print}' "\$f" >> all_filtered.fna
        fi
    done
    """
}
