include { AMG_SUMMARY } from '../modules/amg_summary'

/*
 * Stage 4-summarize: join external DRAM-v annotations to per-sample CoverM
 * quantification. AMG annotation runs as a separate `nextflow run tpall/DRAM
 * -r dev ...` against the vOTU catalog (see README); this subworkflow gathers
 * its outputs into vOTU- and sample-level AMG tables.
 *
 * If neither DRAM-v output table is present, the subworkflow logs a notice
 * and emits no channels rather than failing — Stage 4 is optional and the
 * rest of the pipeline must still run.
 */
workflow SUMMARIZE {
    take:
    quant_long              // path: votu_long.tsv from QUANTIFY
    annotations_path        // string or null: ANNOTATE/annotations_with_flags.tsv
    summary_path            // string or null: SUMMARIZE/summarized_genomes.tsv

    main:
    def ann_file = annotations_path ? file(annotations_path) : null
    def sum_file = summary_path     ? file(summary_path)     : null

    if (ann_file?.exists() && sum_file?.exists()) {
        AMG_SUMMARY(ann_file, sum_file, quant_long)
        votu_amg        = AMG_SUMMARY.out.votu_amg
        votu_amg_genes  = AMG_SUMMARY.out.votu_amg_genes
        sample_amg      = AMG_SUMMARY.out.sample_amg
        sample_category = AMG_SUMMARY.out.sample_category
        amg_stats       = AMG_SUMMARY.out.amg_stats
    }
    else {
        log.warn "SUMMARIZE: skipping AMG summary — DRAM-v outputs not found " +
                 "(annotations=${annotations_path ?: 'unset'}, " +
                 "summary=${summary_path ?: 'unset'}). " +
                 "Run tpall/DRAM dev against \${outdir}/cluster/votu_catalog.fa, " +
                 "then re-run with -resume."
        votu_amg        = Channel.empty()
        votu_amg_genes  = Channel.empty()
        sample_amg      = Channel.empty()
        sample_category = Channel.empty()
        amg_stats       = Channel.empty()
    }

    emit:
    votu_amg
    votu_amg_genes
    sample_amg
    sample_category
    amg_stats
}
