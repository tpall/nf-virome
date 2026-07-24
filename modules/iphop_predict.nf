process IPHOP_PREDICT {
    label 'process_high'
    label 'process_long'

    container 'quay.io/biocontainers/iphop:1.4.2--pyhdfd78af_0'

    publishDir "${params.outdir}/host_iphop", mode: 'copy'

    input:
    path catalog        // cluster/votu_catalog.fa (vOTU representatives)
    path iphop_db       // iPHoP database dir (stock, or add_to_db-augmented with study MAGs)

    output:
    path "iphop_host_genus.csv",        emit: genus
    path "iphop_host_genome.csv",       emit: genome
    path "Detailed_output_by_tool.csv", emit: detailed, optional: true
    path "versions.yml",                emit: versions

    script:
    // iPHoP integrates CRISPR + blast homology + genome k-mer/codon-usage
    // (Random Forest) into one calibrated host call with a confidence score,
    // giving broad coverage beyond the exact-spacer HOST_COUPLE stage. Default
    // --min_score 90 → filenames carry an _m90 suffix; normalize by glob.
    """
    set -euo pipefail
    iphop predict \\
        --fa_file ${catalog} \\
        --db_dir ${iphop_db} \\
        --out_dir iphop_out \\
        --num_threads ${task.cpus}

    cp iphop_out/Host_prediction_to_genus_m*.csv  iphop_host_genus.csv
    cp iphop_out/Host_prediction_to_genome_m*.csv iphop_host_genome.csv
    if [ -f iphop_out/Detailed_output_by_tool.csv ]; then
        cp iphop_out/Detailed_output_by_tool.csv .
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iphop: 1.4.2
    END_VERSIONS
    """
}
