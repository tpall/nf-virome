#!/usr/bin/env bash
#
# Stage 6 half 1: run DRAM-v (tpall/DRAM dev) on a finished nf-virome vOTU
# catalog to produce the AMG annotation tables that nf-virome's SUMMARIZE
# subworkflow consumes. This closes the loop the two-pass design leaves open:
#
#   1. bin/run_pipeline.sh  ...            # discovery; SUMMARIZE no-ops
#   2. bin/run_dramv.sh <nf-virome-outdir> # THIS script; writes <outdir>/dramv
#   3. DRAMV2_OUTDIR=<outdir>/dramv bin/run_pipeline.sh ...   # -resume fires SUMMARIZE
#
# Why it's safe to run now (tpall/DRAM issue #13): PR #12's global `-T 10`
# hmmsearch prefilter inflates raw VOG / V-flag counts ~2x, but the AMG endpoint
# is UNCHANGED (auxiliary_score is geNomad/VirSorter-derived; the AMG gate keys
# off amg_flags + auxiliary_score, not VOG). A validation run measured 294 -> 294
# AMG calls across `-E` vs `-T`. SUMMARIZE consumes exactly that AMG endpoint, so
# the open per-DB fix (issue #13) does not affect our output.
#
# Inputs come straight from the nf-virome outdir:
#   <outdir>/cluster/votu_catalog.fa                     (the catalog)
#   <outdir>/dramv_input/genomad_genes/*_virus_genes.tsv (per-sample geNomad
#       genes; DRAM-v's --genomad_genes matches them to catalog contigs by their
#       sample-prefixed ids, so no per-catalog re-run of geNomad is needed)
#
# Flags mirror the issue-#13 validation run (catalog mode, DRAM-v).
#
# Environment (DB config + repo are site-specific — tpall/DRAM owns its DB setup):
#   DRAM_REPO     nextflow project for tpall/DRAM (default: tpall/DRAM, pulled by NF)
#   DRAM_REF      git ref (default: dev — carries PR #12)
#   DRAM_PROFILE  nextflow profile(s) (default: singularity)
#   DRAM_CONFIG   optional -c config supplying DRAM DB paths (e.g. ~/dram_dbs.config)
#
# Usage:
#   bin/run_dramv.sh <nf-virome-outdir>
# Example:
#   DRAM_CONFIG=~/dram_dbs.config bin/run_dramv.sh results_aging620
#
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <nf-virome-outdir>"
    exit 1
fi

OUTDIR="$1"
CATALOG="${OUTDIR}/cluster/votu_catalog.fa"
GENOMAD_GENES="${OUTDIR}/dramv_input/genomad_genes/*_virus_genes.tsv"
DRAMV_OUT="${OUTDIR}/dramv"

if [[ ! -f "${CATALOG}" ]]; then
    echo "Catalog not found: ${CATALOG}"
    echo "Run discovery first (bin/run_pipeline.sh) so CLUSTER emits the catalog."
    exit 1
fi

cfg=()
[[ -n "${DRAM_CONFIG:-}" ]] && cfg+=(-c "${DRAM_CONFIG}")

nextflow run "${DRAM_REPO:-tpall/DRAM}" -r "${DRAM_REF:-dev}" \
    --input_fasta   "${CATALOG}" \
    --genomad_genes "${GENOMAD_GENES}" \
    --use_dramv --use_kofam --use_merops \
    --max_auxiliary_score 3 \
    --call --annotate --summarize \
    --outdir "${DRAMV_OUT}" \
    -profile "${DRAM_PROFILE:-singularity}" \
    "${cfg[@]}"

echo
echo "DRAM-v done. Tables for SUMMARIZE:"
echo "  ${DRAMV_OUT}/ANNOTATE/annotations_with_flags.tsv"
echo "  ${DRAMV_OUT}/SUMMARIZE/summarized_genomes.tsv"
echo "Now re-run discovery with DRAMV2_OUTDIR=${DRAMV_OUT} to fire SUMMARIZE."
