#!/usr/bin/env python3
"""
Merge per-sample CoverM TSVs into wide count and TPM matrices keyed on vOTU.

Each input file is named `<sample>_coverage.tsv` and has the canonical layout
emitted by modules/coverm.nf:

    vOTU    count   tpm   covered_fraction

Outputs
-------
--out-count : wide TSV (rows=vOTU, cols=samples) of read counts
--out-tpm   : wide TSV (rows=vOTU, cols=samples) of TPM values
--out-long  : long TSV `vOTU \t sample \t count \t tpm \t covered_fraction`
              for diagnostic / sparse-data use.
"""
import argparse
import csv
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--indir",     required=True, help="directory of <sample>_coverage.tsv files")
    p.add_argument("--out-count", required=True)
    p.add_argument("--out-tpm",   required=True)
    p.add_argument("--out-long",  required=True)
    return p.parse_args()


def main():
    args = parse_args()
    indir = Path(args.indir)
    files = sorted(indir.glob("*_coverage.tsv"))
    if not files:
        sys.exit(f"No *_coverage.tsv files in {indir}")

    samples = []
    counts, tpms = {}, {}      # vOTU -> {sample -> value}
    long_rows = []

    for f in files:
        sample = f.name.replace("_coverage.tsv", "")
        samples.append(sample)
        with f.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                votu = row["vOTU"]
                c = int(float(row["count"]))
                t = float(row["tpm"])
                cf = float(row["covered_fraction"])
                counts.setdefault(votu, {})[sample] = c
                tpms.setdefault(votu, {})[sample]   = t
                long_rows.append((votu, sample, c, t, cf))

    votus = sorted(counts.keys())

    def write_wide(path, table, fmt):
        with open(path, "w") as out:
            out.write("vOTU\t" + "\t".join(samples) + "\n")
            for v in votus:
                row = [fmt(table[v].get(s, 0)) for s in samples]
                out.write(v + "\t" + "\t".join(row) + "\n")

    write_wide(args.out_count, counts, lambda x: str(int(x)))
    write_wide(args.out_tpm,   tpms,   lambda x: f"{x:.6g}")

    with open(args.out_long, "w") as out:
        out.write("vOTU\tsample\tcount\ttpm\tcovered_fraction\n")
        for v, s, c, t, cf in long_rows:
            out.write(f"{v}\t{s}\t{c}\t{t:.6g}\t{cf:.6g}\n")

    sys.stderr.write(
        f"Merged {len(samples)} samples × {len(votus)} vOTUs "
        f"({len(long_rows)} long-format rows)\n"
    )


if __name__ == "__main__":
    main()
