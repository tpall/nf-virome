include { IPHOP_PREDICT } from '../modules/iphop_predict'

/*
 * Stage 4b (HOST_IPHOP): integrated host prediction for the vOTU catalog.
 * Complements the exact-spacer HOST_COUPLE stage — iPHoP fuses CRISPR, blast
 * homology, and genome k-mer/codon-usage similarity into a calibrated call with
 * a confidence score, so it reaches the many vOTUs that carry no matching spacer.
 * Point --iphop_db at the stock db, or an add_to_db-augmented db that includes
 * the cohort's own MAGs (bin/build_iphop_db.sh) so hosts land on the study taxa.
 */
workflow HOST_IPHOP {
    take:
    catalog             // path: cluster/votu_catalog.fa
    iphop_db            // path: iPHoP database dir

    main:
    IPHOP_PREDICT(catalog, iphop_db)

    emit:
    genus    = IPHOP_PREDICT.out.genus
    genome   = IPHOP_PREDICT.out.genome
    versions = IPHOP_PREDICT.out.versions
}
