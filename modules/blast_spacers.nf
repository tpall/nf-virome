process BLAST_SPACERS {
    label 'process_medium'

    container 'quay.io/biocontainers/blast:2.17.0--h66d330f_0'

    publishDir "${params.outdir}/host_couple", mode: 'copy'

    input:
    path spacers
    path catalog

    output:
    path "spacer_blast_filtered.tsv", emit: hits
    path "spacer_blast_raw.tsv",      emit: raw
    path "versions.yml",              emit: versions

    script:
    """
    set -euo pipefail

    # vOTU catalog as the BLAST DB
    makeblastdb -in ${catalog} -dbtype nucl -out catalog_db -parse_seqids

    # blastn-short for ~30 bp spacer queries. -word_size 7 catches near-perfect
    # matches; we will further filter for >=95% identity and <=1 mismatch on
    # the full spacer length.
    blastn \\
        -task blastn-short \\
        -query ${spacers} \\
        -db catalog_db \\
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \\
        -perc_identity 90 \\
        -word_size 7 \\
        -dust no \\
        -num_threads ${task.cpus} \\
        -out spacer_blast_raw.tsv

    # Strict CRISPR-spacer host-prediction filter (Edwards et al. 2016 style):
    #   >=95% identity, <=1 mismatch, alignment covers >=95% of the spacer length
    awk -F'\\t' 'BEGIN{OFS="\\t"} \$3+0>=95 && \$5+0<=1 && \$13>0 && (\$4/\$13)>=0.95' \\
        spacer_blast_raw.tsv > spacer_blast_filtered.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | head -1 | sed 's/^blastn: //')
    END_VERSIONS
    """
}
