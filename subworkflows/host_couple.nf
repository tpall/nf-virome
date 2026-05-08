include { MINCED          } from '../modules/minced'
include { BUILD_SPACER_DB } from '../modules/build_spacer_db'
include { BLAST_SPACERS   } from '../modules/blast_spacers'
include { HOST_LINK       } from '../modules/host_link'

workflow HOST_COUPLE {
    take:
    bins_ch     // channel of bin FASTA files (one per element)
    catalog     // path: votu_catalog.fa

    main:
    MINCED( bins_ch )

    // Gather all per-bin spacer files (drop bin-name tag for the global DB).
    spacer_files = MINCED.out.spacers
        .map { name, f -> f }
        .collect()

    BUILD_SPACER_DB( spacer_files )
    BLAST_SPACERS  ( BUILD_SPACER_DB.out.db, catalog )
    HOST_LINK      ( BLAST_SPACERS.out.hits, BUILD_SPACER_DB.out.mapping )

    emit:
    spacer_db        = BUILD_SPACER_DB.out.db
    spacer_to_bin    = BUILD_SPACER_DB.out.mapping
    bin_spacer_stats = BUILD_SPACER_DB.out.stats
    blast_hits       = BLAST_SPACERS.out.hits
    host_links       = HOST_LINK.out.hits
    host_summary     = HOST_LINK.out.summary
    versions         = MINCED.out.versions.mix( BLAST_SPACERS.out.versions )
}
