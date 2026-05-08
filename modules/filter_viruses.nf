process FILTER_VIRUSES {
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    publishDir path: { "${params.outdir}/identify/filtered/${meta.id}" }, mode: 'copy'
    // Flat copy for downstream DRAM-v consumption: every sample's filtered fasta
    // ends up in one directory so a single --input_fasta arg works without nested
    // globs. See README "AMG annotation" section.
    publishDir path: { "${params.outdir}/dramv_input/fastas" },
               mode: 'copy', pattern: "*_filtered.fna"

    input:
    tuple val(meta), path(viruses), path(proviruses), path(quality)
    val   min_length
    val   keep_quality

    output:
    tuple val(meta), path("${meta.id}_filtered.fna"),  emit: filtered
    tuple val(meta), path("${meta.id}_kept_ids.txt"),  emit: kept_ids
    tuple val(meta), path("${meta.id}_filter_stats.tsv"), emit: stats

    script:
    // keep_quality is a comma-separated string; build alternation regex
    def quality_re = keep_quality.tokenize(',').collect{ it.trim() }.join('|')
    """
    set -euo pipefail

    # CheckV viruses.fna = trimmed-input viral contigs; proviruses.fna = excised provirus regions
    cat ${viruses} ${proviruses} > all.fna

    # Keep contigs whose checkv_quality matches AND contig_length >= min_length.
    # Header columns (CheckV >=1.0): contig_id, contig_length, provirus, proviral_length,
    #   gene_count, viral_genes, host_genes, checkv_quality, miuvig_quality, completeness,
    #   completeness_method, contamination, kmer_freq, warnings
    awk -F'\\t' -v re="^(${quality_re})\$" -v minlen=${min_length} '
        NR==1 {
            for (i=1; i<=NF; i++) {
                if (\$i == "contig_id")      cid=i;
                if (\$i == "contig_length")  cl=i;
                if (\$i == "checkv_quality") cq=i;
            }
            next
        }
        \$cq ~ re && \$cl+0 >= minlen { print \$cid }
    ' ${quality} > ${meta.id}_kept_ids.txt

    if [ -s ${meta.id}_kept_ids.txt ]; then
        seqkit grep -f ${meta.id}_kept_ids.txt all.fna > ${meta.id}_filtered.fna
    else
        : > ${meta.id}_filtered.fna
    fi

    # Per-sample stats line
    n_in=\$(grep -c '^>' all.fna || true)
    n_out=\$(grep -c '^>' ${meta.id}_filtered.fna || true)
    printf "sample\\tn_input\\tn_kept\\tmin_length\\tquality_kept\\n%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${meta.id}" "\$n_in" "\$n_out" "${min_length}" "${keep_quality}" > ${meta.id}_filter_stats.tsv
    """
}
