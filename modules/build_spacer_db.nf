process BUILD_SPACER_DB {
    label 'process_low'

    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    publishDir "${params.outdir}/host_couple", mode: 'copy'

    input:
    path spacer_files       // one *_spacers.fa per bin

    output:
    path "all_spacers.fa",     emit: db
    path "spacer_to_bin.tsv",  emit: mapping
    path "bin_spacer_stats.tsv", emit: stats

    script:
    """
    set -euo pipefail

    : > all_spacers.fa
    : > spacer_to_bin.tsv
    printf "bin\\tn_spacers\\n" > bin_spacer_stats.tsv

    # Each input is named <bin>_spacers.fa. Prefix every spacer header with the
    # bin id so the downstream BLAST hit's qseqid encodes the source bin.
    # Also accumulate spacer->bin mapping and per-bin counts.
    #
    # NOTE: don't use `\$(grep -c '^>' f || echo 0)` — when grep finds no
    # matches it both prints "0" and exits 1, so the `|| echo 0` *also* runs
    # and the substitution captures "0\\n0", inflating the stats file.
    for f in *_spacers.fa; do
        # Skip the output file we're building if it already exists in cwd
        [ "\$f" = "all_spacers.fa" ] && continue
        bin=\$(basename "\$f" _spacers.fa)
        n=\$(grep -c '^>' "\$f" 2>/dev/null) || n=0
        printf "%s\\t%s\\n" "\$bin" "\$n" >> bin_spacer_stats.tsv
        if [ "\$n" -gt 0 ]; then
            seqkit replace -p '^(\\S+).*' -r "\${bin}__spacer_{nr}" "\$f" >> all_spacers.fa
            seqkit fx2tab -n -i "\$f" | awk -v bin="\$bin" 'BEGIN{n=0} {n++; printf "%s__spacer_%d\\t%s\\n", bin, n, bin}' >> spacer_to_bin.tsv
        fi
    done

    n_total=\$(grep -c '^>' all_spacers.fa || echo 0)
    echo "Built spacer DB: \${n_total} spacers from \$(awk 'NR>1 {n++} END{print n+0}' bin_spacer_stats.tsv) bins" >&2
    """
}
