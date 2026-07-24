process PHATYP {
    label 'process_medium'

    container 'quay.io/biocontainers/phabox:2.1.13--pyhdfd78af_1'

    publishDir "${params.outdir}/lifestyle", mode: 'copy'

    input:
    path catalog        // cluster/votu_catalog.fa (vOTU representatives)
    path phabox_db      // PhaBOX2 database directory (--dbdir)

    output:
    path "votu_lifestyle.tsv", emit: lifestyle
    path "versions.yml",       emit: versions

    script:
    // PhaBOX2 PhaTYP task: temperate/virulent prediction per contig. Verified
    // end-to-end (PhaBOX2 2.1.13, phabox_db_v2_2, 2026-07-24): writes
    // <outpth>/final_prediction/phatyp_prediction.tsv with columns
    // Accession (= vOTU rep id), Length, TYPE (virulent|temperate), PhaTYPScore.
    // PhaTYP's default --len 3000 drops nothing here (catalog is >= min_length 5000).
    """
    set -euo pipefail
    phabox2 --task phatyp \\
        --contigs ${catalog} \\
        --outpth phabox_out \\
        --dbdir ${phabox_db} \\
        --threads ${task.cpus}

    cp phabox_out/final_prediction/phatyp_prediction.tsv votu_lifestyle.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phabox: 2.1.13
    END_VERSIONS
    """
}
