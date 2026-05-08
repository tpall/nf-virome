#!/usr/bin/env python3
"""
Build per-(vOTU, bin) host links from CRISPR-spacer BLAST hits.

Inputs
------
--blast    : filtered BLAST TSV with columns
             qseqid sseqid pident length mismatch gapopen qstart qend
             sstart send evalue bitscore qlen slen
             (qseqid format: <bin>__spacer_<N>; sseqid is the vOTU id)
--mapping  : spacer_to_bin.tsv emitted by BUILD_SPACER_DB (spacer_id \\t bin)

Outputs
-------
--out-hits         : long-format `vOTU \\t bin \\t spacer \\t pident \\t mismatch`
--out-summary      : per-vOTU summary
                     `vOTU \\t n_bins \\t bins (semicolon-joined)`
                     plus n_spacer_hits.

Bin taxonomy is intentionally not joined here — that lives in the project's
GTDB outputs and is loaded in the R analysis layer alongside the bacterial
phyloseq data.
"""
import argparse
import csv
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--blast",       required=True)
    p.add_argument("--mapping",     required=True)
    p.add_argument("--out-hits",    required=True)
    p.add_argument("--out-summary", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    # Load spacer -> bin mapping. (We could parse the bin from qseqid via the
    # __spacer_ prefix, but the explicit mapping is more robust against header
    # collisions.)
    spacer_bin = {}
    with open(args.mapping) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                spacer_bin[parts[0]] = parts[1]

    # Walk BLAST hits.
    per_vOTU_bins = defaultdict(set)        # vOTU -> set of bins
    per_vOTU_hits = defaultdict(int)        # vOTU -> spacer-hit count
    long_rows = []

    cols = ("qseqid sseqid pident length mismatch gapopen qstart qend "
            "sstart send evalue bitscore qlen slen").split()

    with open(args.blast) as f:
        for line in f:
            row = line.rstrip("\n").split("\t")
            if len(row) < len(cols):
                continue
            r = dict(zip(cols, row))
            spacer = r["qseqid"]
            votu   = r["sseqid"]
            # Fall back to the __spacer_ prefix if the spacer is missing from
            # the mapping (shouldn't happen, but defensive)
            bin_id = spacer_bin.get(spacer, spacer.split("__spacer_")[0])

            per_vOTU_bins[votu].add(bin_id)
            per_vOTU_hits[votu] += 1
            long_rows.append((votu, bin_id, spacer, r["pident"], r["mismatch"]))

    with open(args.out_hits, "w") as f:
        f.write("vOTU\tbin\tspacer\tpident\tmismatch\n")
        for r in long_rows:
            f.write("\t".join(r) + "\n")

    with open(args.out_summary, "w") as f:
        f.write("vOTU\tn_bins\tn_spacer_hits\tbins\n")
        for votu in sorted(per_vOTU_bins):
            bins = sorted(per_vOTU_bins[votu])
            f.write(f"{votu}\t{len(bins)}\t{per_vOTU_hits[votu]}\t{';'.join(bins)}\n")

    sys.stderr.write(
        f"Linked {len(per_vOTU_bins)} vOTUs to "
        f"{sum(len(v) for v in per_vOTU_bins.values())} (vOTU, bin) pairs "
        f"from {len(long_rows)} spacer hits\n"
    )


if __name__ == "__main__":
    main()
