include { COVERM       } from '../modules/coverm'
include { COVERM_MERGE } from '../modules/coverm_merge'

workflow QUANTIFY {
    take:
    reads_ch              // channel: [meta, reads]  (reads = [r1] or [r1,r2])
    catalog               // path: votu_catalog.fa
    min_covered_fraction  // val (0..1)
    min_read_pid          // val (0..1)

    main:
    COVERM( reads_ch, catalog, min_covered_fraction, min_read_pid )

    coverages = COVERM.out.coverage.map { meta, f -> f }.collect()

    COVERM_MERGE( coverages )

    emit:
    count_matrix = COVERM_MERGE.out.count
    tpm_matrix   = COVERM_MERGE.out.tpm
    long_table   = COVERM_MERGE.out.long
    versions     = COVERM.out.versions
}
