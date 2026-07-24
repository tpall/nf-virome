#!/usr/bin/env bash
#
# Launch nf-virome against a samplesheet, parking outputs under a named --outdir
# with -with-report / -with-trace / -with-timeline alongside. Wrap a tmux/screen
# session around this for multi-day runs so detaching doesn't kill the orchestrator.
#
# Prerequisites
# -------------
# - On HPC, with `module load singularity/3.8.5` and `conda activate java` (nextflow).
# - Set the required DB paths in the environment (or edit the flags below):
#     export GENOMAD_DB=/path/to/genomad_db
#     export CHECKV_DB=/path/to/checkv-db-v1.5
#   Optional:
#     export PHABOX_DB=/path/to/phabox_db   # enables Stage 5 LIFESTYLE (PhaTYP)
#     export DRAMV2_OUTDIR=/path/to/dramv    # enables Stage 6 SUMMARIZE (AMG)
# - Build the samplesheet first (see assets/build_samplesheet.R):
#     # 620-sample aging cohort (8 projects, one shared vOTU catalog).
#     # Run this ON THE HPC so --check-files can flag single-end samples.
#     Rscript assets/build_samplesheet.R assets/samplesheet_aging620.csv \
#       --samples assets/aging_samples_620.txt --eluring ../eluring --check-files
#
# Usage
# -----
#   bin/run_pipeline.sh <samplesheet.csv> <outdir>
#
# Example
# -------
#   tmux new -s virome
#   bin/run_pipeline.sh assets/samplesheet_aging620.csv results_aging620
#   # detach: Ctrl+b d   reattach: tmux attach -t virome
#
# Two-pass AMG summary
# --------------------
# DRAM-v (tpall/DRAM -r dev --use_dramv) runs SEPARATELY against
# ${outdir}/cluster/votu_catalog.fa. Stage 6 SUMMARIZE skips gracefully until
# its tables exist, so:
#   1. First pass: run discovery (this script; SUMMARIZE no-ops).
#   2. Run DRAM-v on the catalog.
#   3. Re-run this exact command with DRAMV2_OUTDIR set; Nextflow -resume reuses
#      the cached discovery work and only fires SUMMARIZE.
#
# Resume after interruption: re-run the same command (Nextflow auto-resumes from
# the cached work dir) or pin a session id with `-resume <uuid>`.
#
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <samplesheet.csv> <outdir>"
    exit 1
fi

SAMPLESHEET="$1"
OUTDIR="$2"

if [[ ! -f "${SAMPLESHEET}" ]]; then
    echo "Samplesheet not found: ${SAMPLESHEET}"
    echo "Build it first with assets/build_samplesheet.R — see the header of this script."
    exit 1
fi

mkdir -p "${OUTDIR}"

# Optional stage databases: only pass the flag when the env var is set, so the
# graceful-skip gating in main.nf stays intact when they're absent.
extra=()
[[ -n "${PHABOX_DB:-}"      ]] && extra+=(--phabox_db      "${PHABOX_DB}")
[[ -n "${DRAMV2_OUTDIR:-}"  ]] && extra+=(--dramv2_outdir  "${DRAMV2_OUTDIR}")

nextflow run . \
    -profile slurm,singularity \
    --samplesheet "${SAMPLESHEET}" \
    --outdir      "${OUTDIR}" \
    --genomad_db  "${GENOMAD_DB:?set GENOMAD_DB in the environment}" \
    --checkv_db   "${CHECKV_DB:?set CHECKV_DB in the environment}" \
    "${extra[@]}" \
    -with-report   "${OUTDIR}/report.html" \
    -with-trace    "${OUTDIR}/trace.txt" \
    -with-timeline "${OUTDIR}/timeline.html"
