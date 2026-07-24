#!/usr/bin/env bash
#
# Augment the stock iPHoP database with the study's own MAGs so HOST_IPHOP
# predicts hosts among the cohort's bacteria (the 620-cohort winners/losers),
# not just generic GTDB references. Runs `iphop add_to_db` once; point
# --iphop_db at the resulting OUT_DB for the pipeline run.
#
# Run this ON THE HPC — it reads the refined bins and per-sample GTDB-Tk
# summaries straight from the assembly outputs (paths in eluring data_paths.csv).
#
# Inputs it stages (for the 620 aging samples, or any --samples list):
#   MAGs      : <HPC_ROOT>/<folder>/results/binrefine/final_bins/*.fa[.gz]
#   taxonomy  : <HPC_ROOT>/<bac120-col>  (per-sample GTDB-Tk bac120 summary)
# iPHoP requires the fasta basenames to match the `user_genome` column of the
# merged GTDB-Tk summary; the assembly pipeline already names them consistently.
#
# Environment:
#   BASE_DB   stock iPHoP db dir (downloaded once; --db_dir)         [required]
#   OUT_DB    output augmented db dir (--out_dir)                    [required]
#   HPC_ROOT  projects root (default /gpfs/helios/home/taavi74/Projects)
#   THREADS   threads for add_to_db (default 8)
#
# Usage:
#   BASE_DB=~/db/iphop_db OUT_DB=~/db/iphop_db_plus620 \
#     bin/build_iphop_db.sh assets/aging_samples_620.txt ../eluring
#
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <samples.txt> [eluring_root]"
    exit 1
fi
SAMPLES="$1"
ELURING="${2:-../eluring}"
DATA_PATHS="${ELURING}/data/data_paths.csv"
HPC_ROOT="${HPC_ROOT:-/gpfs/helios/home/taavi74/Projects}"
THREADS="${THREADS:-8}"
: "${BASE_DB:?set BASE_DB to the stock iPHoP db dir}"
: "${OUT_DB:?set OUT_DB to the output augmented db dir}"

[[ -f "${SAMPLES}" ]]    || { echo "samples list not found: ${SAMPLES}"; exit 1; }
[[ -f "${DATA_PATHS}" ]] || { echo "data_paths.csv not found: ${DATA_PATHS} (pass eluring_root)"; exit 1; }

STAGE="$(mktemp -d)"
MAGS="${STAGE}/mags"; GTDB="${STAGE}/gtdb"
mkdir -p "${MAGS}" "${GTDB}"
trap 'rm -rf "${STAGE}"' EXIT

# data_paths.csv columns: 1 sample, 2 project, 3 folder, 4 covstats,
# 5 contig_to_bin, 6 bac120, 7 ar53. Select the requested samples.
awk -F, -v OFS='\t' 'NR==FNR{want[$1]=1;next} FNR==1{next} ($1 in want){print $3,$6,$7}' \
    "${SAMPLES}" "${DATA_PATHS}" > "${STAGE}/sel.tsv"
n_sel=$(wc -l < "${STAGE}/sel.tsv")
echo "Selected ${n_sel} samples from ${DATA_PATHS}"

# Stage MAGs: one final_bins dir per project (folder), symlink/decompress once.
awk '{print $1}' "${STAGE}/sel.tsv" | sort -u | while read -r folder; do
    bins="${HPC_ROOT}/${folder}/results/binrefine/final_bins"
    [[ -d "${bins}" ]] || { echo "  WARN: no final_bins for ${folder}"; continue; }
    for f in "${bins}"/*.fa "${bins}"/*.fna; do
        [[ -e "${f}" ]] || continue
        ln -sf "${f}" "${MAGS}/"
    done
    for f in "${bins}"/*.fa.gz "${bins}"/*.fna.gz; do
        [[ -e "${f}" ]] || continue
        b=$(basename "${f}" .gz)
        [[ -e "${MAGS}/${b}" ]] || zcat "${f}" > "${MAGS}/${b}"
    done
done
echo "Staged $(ls "${MAGS}" | wc -l) MAG fastas"

# Merge per-sample GTDB-Tk bac120 (and ar53 if present) summaries, one header.
merge_summaries() {  # $1 = column index in sel.tsv, $2 = output name
    local col="$1" out="${GTDB}/$2" first=1
    awk -v c="${col}" '{print $c}' "${STAGE}/sel.tsv" | while read -r rel; do
        [[ "${rel}" == "NA" || -z "${rel}" ]] && continue
        local p="${HPC_ROOT}/${rel}"
        [[ -f "${p}" ]] || continue
        if [[ ${first} -eq 1 ]]; then cat "${p}" > "${out}"; first=0
        else tail -n +2 "${p}" >> "${out}"; fi
    done
    [[ -f "${out}" ]] && echo "  ${2}: $(( $(wc -l < "${out}") - 1 )) genomes"
}
merge_summaries 2 gtdbtk.bac120.summary.tsv
merge_summaries 3 gtdbtk.ar53.summary.tsv

echo "Running iphop add_to_db → ${OUT_DB}"
iphop add_to_db \
    --fna_dir "${MAGS}" \
    --gtdb_dir "${GTDB}" \
    --out_dir "${OUT_DB}" \
    --db_dir "${BASE_DB}" \
    --num_threads "${THREADS}"

echo "Done. Run the pipeline with --iphop_db ${OUT_DB}"
