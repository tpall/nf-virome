include { CONCAT_FASTAS  } from '../modules/concat_fastas'
include { SKANI_TRIANGLE } from '../modules/skani'
include { VOTU_CLUSTER   } from '../modules/votu_cluster'
include { EXTRACT_REPS   } from '../modules/extract_reps'

workflow CLUSTER {
    take:
    filtered_ch     // channel: [meta, filtered.fna] from IDENTIFY
    min_ani         // val  (percent)
    min_af          // val  (percent)

    main:
    // Drop meta and gather every per-sample filtered FASTA into one channel
    fastas = filtered_ch.map { meta, fa -> fa }.collect()

    CONCAT_FASTAS( fastas )
    SKANI_TRIANGLE( CONCAT_FASTAS.out.combined )
    VOTU_CLUSTER(
        CONCAT_FASTAS.out.combined,
        SKANI_TRIANGLE.out.ani,
        min_ani,
        min_af
    )
    EXTRACT_REPS( CONCAT_FASTAS.out.combined, VOTU_CLUSTER.out.reps_ids )

    emit:
    catalog  = EXTRACT_REPS.out.catalog
    clusters = VOTU_CLUSTER.out.clusters
    reps_ids = VOTU_CLUSTER.out.reps_ids
    ani      = SKANI_TRIANGLE.out.ani
    versions = SKANI_TRIANGLE.out.versions
}
