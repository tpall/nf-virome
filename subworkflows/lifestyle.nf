include { PHATYP } from '../modules/phatyp'

/*
 * Stage 5 (LIFESTYLE): predict temperate vs virulent lifestyle for the vOTU
 * catalog representatives with PhaTYP (PhaBOX2). Replaces the earlier R-side
 * `|provirus_` CheckV-suffix proxy with a real per-vOTU call. Kept as a
 * subworkflow so a second caller (e.g. BACPHLIP as a cross-check) can be added
 * without touching main.nf.
 */
workflow LIFESTYLE {
    take:
    catalog             // path: cluster/votu_catalog.fa
    phabox_db           // path: PhaBOX2 database dir

    main:
    PHATYP(catalog, phabox_db)

    emit:
    lifestyle = PHATYP.out.lifestyle
    versions  = PHATYP.out.versions
}
