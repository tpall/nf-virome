#!/usr/bin/env python3
"""
Cluster viral contigs at MIUViG thresholds (95% ANI, 85% AF on the shorter sequence)
from skani triangle --sparse output.

Method: greedy centroid clustering. Sequences sorted by length descending; the
longest unassigned sequence becomes a representative and absorbs any unassigned
neighbours that meet both thresholds. Equivalent to the standard "aniclust"
procedure used in IMG/VR and CheckV.

Inputs
------
--fasta : multi-FASTA whose headers are the sequence IDs to cluster.
--skani : skani triangle output (sparse mode) with the columns:
          Ref_file Query_file ANI Align_fraction_ref Align_fraction_query Ref_name Query_name
--ani   : minimum ANI (percent), default 95.0.
--af    : minimum aligned fraction of the shorter sequence (percent), default 85.0.

Outputs
-------
--out-clusters : TSV with header `representative\tmember`.
--out-reps     : one representative ID per line.
"""
import argparse
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta",        required=True)
    p.add_argument("--skani",        required=True)
    p.add_argument("--ani",          type=float, default=95.0)
    p.add_argument("--af",           type=float, default=85.0)
    p.add_argument("--out-clusters", required=True)
    p.add_argument("--out-reps",     required=True)
    return p.parse_args()


def read_lengths(fa_path):
    lens = {}
    cur, n = None, 0
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    lens[cur] = n
                cur = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
        if cur is not None:
            lens[cur] = n
    return lens


def main():
    args = parse_args()
    lens = read_lengths(args.fasta)

    # Build undirected adjacency from skani sparse output.
    # AF on the shorter sequence = max(af_ref, af_query): the shorter sequence's
    # alignment fraction is whichever side has the higher AF (full coverage of
    # the shorter implies higher AF for that side).
    adj = defaultdict(set)
    with open(args.skani) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        # Required columns
        for col in ("ANI", "Align_fraction_ref", "Align_fraction_query",
                    "Ref_name", "Query_name"):
            if col not in idx:
                sys.exit(f"skani output missing required column: {col}")

        for line in f:
            parts = line.rstrip("\n").split("\t")
            ani = float(parts[idx["ANI"]])
            af  = max(float(parts[idx["Align_fraction_ref"]]),
                      float(parts[idx["Align_fraction_query"]]))
            r   = parts[idx["Ref_name"]]
            q   = parts[idx["Query_name"]]
            if r == q:
                continue
            if ani >= args.ani and af >= args.af:
                adj[r].add(q)
                adj[q].add(r)

    # Greedy: longest-first becomes representative; absorbs unassigned neighbours
    sorted_seqs = sorted(lens, key=lens.get, reverse=True)
    assigned = {}
    reps = []
    for s in sorted_seqs:
        if s in assigned:
            continue
        assigned[s] = s
        reps.append(s)
        for nb in adj.get(s, ()):
            if nb not in assigned:
                assigned[nb] = s

    with open(args.out_clusters, "w") as f:
        f.write("representative\tmember\n")
        for member, rep in sorted(assigned.items()):
            f.write(f"{rep}\t{member}\n")

    with open(args.out_reps, "w") as f:
        for r in reps:
            f.write(f"{r}\n")

    sys.stderr.write(
        f"Clustered {len(assigned)} sequences into {len(reps)} vOTUs "
        f"(ANI≥{args.ani}, AF≥{args.af})\n"
    )


if __name__ == "__main__":
    main()
