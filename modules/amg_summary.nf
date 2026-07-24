process AMG_SUMMARY {
    label 'process_low'

    container 'python:3.12'

    publishDir "${params.outdir}/annotate/summary", mode: 'copy'

    input:
    path annotations         // ANNOTATE/annotations_with_flags.tsv
    path summary             // SUMMARIZE/summarized_genomes.tsv
    path quant_long          // QUANTIFY/votu_long.tsv

    output:
    path "votu_amg.tsv",                  emit: votu_amg
    path "votu_amg_genes.tsv",            emit: votu_amg_genes
    path "sample_amg_abundance.tsv",      emit: sample_amg
    path "sample_amg_category_tpm.tsv",   emit: sample_category
    path "amg_stats.tsv",                 emit: amg_stats

    script:
    """
    set -euo pipefail
    amg_summary.py \\
        --annotations         ${annotations} \\
        --summary             ${summary} \\
        --quant-long          ${quant_long} \\
        --out-votu-amg        votu_amg.tsv \\
        --out-votu-amg-genes  votu_amg_genes.tsv \\
        --out-sample-amg      sample_amg_abundance.tsv \\
        --out-sample-category sample_amg_category_tpm.tsv \\
        --out-amg-stats       amg_stats.tsv
    """
}
