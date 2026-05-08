include { GENOMAD }        from '../modules/genomad'
include { CHECKV }         from '../modules/checkv'
include { FILTER_VIRUSES } from '../modules/filter_viruses'

workflow IDENTIFY {
    take:
    contigs_ch        // channel: [meta, contigs.fa.gz]
    genomad_db        // path
    checkv_db         // path
    min_length        // val
    keep_quality      // val (comma-separated string)

    main:
    GENOMAD( contigs_ch, genomad_db )
    CHECKV(  GENOMAD.out.virus_fasta, checkv_db )

    filter_input = CHECKV.out.viruses
        .join(CHECKV.out.proviruses, by: 0)
        .join(CHECKV.out.quality,    by: 0)

    FILTER_VIRUSES( filter_input, min_length, keep_quality )

    emit:
    filtered = FILTER_VIRUSES.out.filtered
    kept_ids = FILTER_VIRUSES.out.kept_ids
    stats    = FILTER_VIRUSES.out.stats
    quality  = CHECKV.out.quality
    summary  = GENOMAD.out.virus_summary
    taxonomy = GENOMAD.out.virus_taxonomy
    versions = GENOMAD.out.versions.mix(CHECKV.out.versions)
}
