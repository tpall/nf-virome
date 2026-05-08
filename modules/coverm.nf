process COVERM {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/coverm:0.7.0--hcb7b614_4'

    publishDir "${params.outdir}/quantify/per_sample", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path catalog
    val  min_covered_fraction   // 0..1
    val  min_read_pid           // 0..1

    output:
    tuple val(meta), path("${meta.id}_coverage.tsv"), emit: coverage
    path "versions.yml",                              emit: versions

    script:
    def reads_arg = meta.single_end ? "--single ${reads}" : "--coupled ${reads}"
    """
    set -euo pipefail

    # CoverM rejects --min-covered-fraction > 0 when 'count' is among the methods
    # (the count estimator does not natively compute covered_fraction). We instead
    # collect raw values and apply the threshold in awk below — this keeps counts,
    # TPM, and covered_fraction all derived from the same alignments and lets the
    # threshold be retuned later without re-mapping.
    coverm contig \\
        --reference ${catalog} \\
        ${reads_arg} \\
        -m count tpm covered_fraction \\
        --min-read-percent-identity ${min_read_pid} \\
        --output-format dense \\
        -t ${task.cpus} \\
        -o coverm_raw.tsv

    # CoverM column layout (with -m count tpm covered_fraction):
    #   contig, "<reads> Read Count", "<reads> TPM", "<reads> Covered Fraction"
    # Zero out count and TPM where covered_fraction < threshold, so vOTUs that
    # only catch a stray spurious read aren't called "present".
    awk -F'\\t' -v thr=${min_covered_fraction} 'BEGIN{OFS="\\t"}
        NR==1 { print "vOTU","count","tpm","covered_fraction"; next }
        {
            cf = \$4 + 0
            if (cf < thr) { print \$1, 0, 0, cf } else { print \$1, \$2, \$3, cf }
        }
    ' coverm_raw.tsv > ${meta.id}_coverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version 2>&1 | head -1 | sed 's/^coverm //')
    END_VERSIONS
    """
}
