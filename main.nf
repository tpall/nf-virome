#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { IDENTIFY     } from './subworkflows/identify'
include { CLUSTER      } from './subworkflows/cluster'
include { QUANTIFY     } from './subworkflows/quantify'
include { HOST_COUPLE  } from './subworkflows/host_couple'
include { SUMMARIZE    } from './subworkflows/summarize'
include { LIFESTYLE    } from './subworkflows/lifestyle'

workflow {

    if (!params.samplesheet) {
        error "Missing --samplesheet (CSV with: sample, contigs, reads_1, reads_2, bins_dir, single_end)"
    }

    // Per-sample manifest channel.
    // Required columns: sample, contigs, reads_1.
    // Optional: reads_2 (paired-end), bins_dir (host coupling), single_end ("true"/"false").
    samplesheet_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id         : row.sample,
                single_end : Boolean.parseBoolean((row.single_end ?: 'false').toString().trim())
            ]
            [ meta, row ]
        }

    contigs_ch = samplesheet_ch.map { meta, row ->
        [ meta, file(row.contigs, checkIfExists: true) ]
    }

    reads_ch = samplesheet_ch.map { meta, row ->
        def reads = meta.single_end ?
            [ file(row.reads_1, checkIfExists: true) ] :
            [ file(row.reads_1, checkIfExists: true),
              file(row.reads_2, checkIfExists: true) ]
        [ meta, reads ]
    }

    // Stage 1: per-sample geNomad + CheckV + length/quality filter.
    IDENTIFY(
        contigs_ch,
        file(params.genomad_db, checkIfExists: true),
        file(params.checkv_db,  checkIfExists: true),
        params.min_length,
        params.keep_quality
    )

    // Concatenated per-sample identification stats.
    IDENTIFY.out.stats
        .map { meta, f -> f }
        .collectFile(
            name: 'identify_stats.tsv',
            storeDir: "${params.outdir}/identify",
            keepHeader: true,
            skip: 1
        )

    // Stage 2: cluster filtered viruses across all samples into a vOTU catalog.
    CLUSTER(
        IDENTIFY.out.filtered,
        params.votu_min_ani,
        params.votu_min_af
    )

    // Stage 3: quantify by mapping reads back to the catalog.
    QUANTIFY(
        reads_ch,
        CLUSTER.out.catalog,
        params.coverm_min_covered_fraction,
        params.coverm_min_read_pid
    )

    // Stage 4: host coupling via CRISPR-spacers from existing MAGs.
    // Each unique bins_dir from the samplesheet contributes its bin FASTAs.
    // Bins are typically gzipped — minced will decompress on the fly.
    // Skip the stage entirely if --skip_host_couple is set, or if no
    // bins_dir column exists / all rows omit it.
    if (!params.skip_host_couple) {
        bins_ch = samplesheet_ch
            .map { meta, row -> row.bins_dir }
            .filter { it != null && it.toString().trim() != '' }
            .unique()
            .flatMap { d ->
                def dir = file(d, checkIfExists: true)
                dir.list()
                   .findAll { it ==~ /.+\.(fa|fasta|fna)(\.gz)?$/ }
                   .collect { dir.resolve(it) }
            }

        HOST_COUPLE(
            bins_ch,
            CLUSTER.out.catalog
        )
    }

    // Stage 5 (LIFESTYLE): temperate/virulent prediction on the vOTU catalog
    // representatives with PhaTYP (PhaBOX2). Optional — runs only when the
    // --phabox_db is set. Replaces the earlier R-side provirus-suffix proxy.
    if (params.phabox_db) {
        LIFESTYLE(
            CLUSTER.out.catalog,
            file(params.phabox_db, checkIfExists: true)
        )
    }

    // Stage 6 (SUMMARIZE): join external DRAM-v (tpall/DRAM -r dev --use_dramv)
    // AMG annotations to the CoverM quantification. DRAM-v runs separately
    // against ${outdir}/cluster/votu_catalog.fa; this stage skips gracefully
    // (warning, no failure) until its two output tables exist, so a first pass
    // runs discovery and a `-resume` pass after DRAM-v fires the AMG summary.
    def dramv_ann = params.dramv2_annotations ?:
        (params.dramv2_outdir ? "${params.dramv2_outdir}/ANNOTATE/annotations_with_flags.tsv" : null)
    def dramv_sum = params.dramv2_summary ?:
        (params.dramv2_outdir ? "${params.dramv2_outdir}/SUMMARIZE/summarized_genomes.tsv" : null)

    SUMMARIZE(
        QUANTIFY.out.long_table,
        dramv_ann,
        dramv_sum
    )
}
