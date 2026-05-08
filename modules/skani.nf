process SKANI_TRIANGLE {
    label 'process_high'

    container 'quay.io/biocontainers/skani:0.3.1--ha6fb395_0'

    publishDir "${params.outdir}/cluster", mode: 'copy'

    input:
    path catalog_fasta

    output:
    path "skani_ani.tsv",  emit: ani
    path "versions.yml",   emit: versions

    script:
    """
    set -euo pipefail

    # skani triangle treats each input FILE as one genome. To get pairwise
    # contig-level ANI we split the catalog into one file per contig and
    # pass the list to skani.
    mkdir -p contigs
    awk 'BEGIN{out=""} /^>/ { if (out!="") close(out); name=substr(\$1,2); out="contigs/"name".fa" } { print > out }' \\
        ${catalog_fasta}

    # Avoid argv overflow with thousands of files: hand a list file to skani.
    find contigs -name '*.fa' > input_list.txt

    skani triangle \\
        -l input_list.txt \\
        --sparse \\
        --min-af 50 \\
        -t ${task.cpus} \\
        -o skani_ani.tsv

    rm -rf contigs input_list.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version 2>&1 | sed 's/^.*skani //')
    END_VERSIONS
    """
}
