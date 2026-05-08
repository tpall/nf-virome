process MINCED {
    tag "$bin_name"
    label 'process_low'

    container 'quay.io/biocontainers/minced:0.3.2--0'

    input:
    path bin_fa

    output:
    tuple val(bin_name), path("${bin_name}_spacers.fa"), emit: spacers
    tuple val(bin_name), path("${bin_name}.crisprs"),    emit: arrays,  optional: true
    path "versions.yml",                                 emit: versions

    script:
    // Strip both the .gz (if any) and the {.fa, .fasta, .fna} suffix.
    bin_name = bin_fa.name.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, "")
    """
    set -euo pipefail

    # minced needs an uncompressed FASTA. Decompress on the fly if needed.
    if [[ "${bin_fa}" == *.gz ]]; then
        zcat ${bin_fa} > ${bin_name}.fa
        INPUT=${bin_name}.fa
    else
        INPUT=${bin_fa}
    fi

    # minced 0.3.2 CLI: `minced [options] file.fa [outputFile]` — one positional
    # output max. With `-spacers`, minced derives a spacer filename from the
    # input basename: <basename>_spacers.fa. Since INPUT is ${bin_name}.fa,
    # the spacer file is ${bin_name}_spacers.fa — already our canonical name.
    # Discard the stdout summary; OK to fail (some bins have no CRISPR arrays).
    minced -spacers \$INPUT > /dev/null 2>&1 || true

    # Ensure the canonical file exists even when minced wrote nothing.
    [ -f ${bin_name}_spacers.fa ] || : > ${bin_name}_spacers.fa

    # Clean up the staged uncompressed FASTA if we made one.
    [ "\$INPUT" != "${bin_fa}" ] && rm -f "\$INPUT" || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minced: \$(minced --version 2>&1 | sed -n 's/^minced \\(.*\\)/\\1/p' | head -1)
    END_VERSIONS
    """
}
