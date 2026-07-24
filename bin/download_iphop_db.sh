#!/usr/bin/env bash
#
# Download + extract the iPHoP base database on the HPC. Resumable (wget -c) and
# retry-tolerant, because the only official host (NERSC portal) has intermittent
# outages. Fetches the current default db (latest_rw = Jun_2025_pub_rw at time of
# writing). Point --iphop_db (or build_iphop_db.sh BASE_DB) at the extracted dir.
#
# Run detached so it survives logout and keeps retrying if NERSC is down:
#   screen -dmS iphop_dl bin/download_iphop_db.sh /gpfs/space/projects/preterm/databases/iphop
#   screen -r iphop_dl        # reattach
#   tail -f <dest>/iphop_download.log
#
# Idempotent: re-running resumes a partial tarball or no-ops if already complete.
#
set -uo pipefail

DEST="${1:-/gpfs/space/projects/preterm/databases/iphop}"
URL="https://portal.nersc.gov/cfs/m342/iphop/db/iPHoP.latest_rw.tar.gz"
LOG="${DEST}/iphop_download.log"
MAX_HOURS="${MAX_HOURS:-24}"      # give up after this many hours of retrying
RETRY_SLEEP="${RETRY_SLEEP:-1200}" # 20 min between retries while the portal is down

mkdir -p "${DEST}"
cd "${DEST}"
echo "$(date) === start; dest=${DEST} url=${URL}" | tee -a "${LOG}"

deadline=$(( $(date +%s) + MAX_HOURS * 3600 ))
until wget -c --timeout=60 --tries=3 "${URL}" >> "${LOG}" 2>&1; do
    if [ "$(date +%s)" -ge "${deadline}" ]; then
        echo "$(date) gave up after ${MAX_HOURS}h (NERSC portal still unreachable)" | tee -a "${LOG}"
        exit 1
    fi
    echo "$(date) download attempt failed (NERSC portal likely down); retry in $((RETRY_SLEEP/60)) min" | tee -a "${LOG}"
    sleep "${RETRY_SLEEP}"
done

TARBALL="${DEST}/iPHoP.latest_rw.tar.gz"
echo "$(date) download complete ($(du -h "${TARBALL}" | cut -f1)); extracting" | tee -a "${LOG}"
tar -C "${DEST}" -zxvf "${TARBALL}" >> "${LOG}" 2>&1

DBDIR=$(tar -tzf "${TARBALL}" 2>/dev/null | head -1 | cut -d/ -f1)
echo "$(date) === DONE. iphop_db = ${DEST}/${DBDIR}" | tee -a "${LOG}"
echo "Point the pipeline at:  --iphop_db ${DEST}/${DBDIR}"
echo "(or build_iphop_db.sh with BASE_DB=${DEST}/${DBDIR})"
