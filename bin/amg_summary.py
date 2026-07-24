#!/usr/bin/env python3
"""
Join DRAM-v annotations to the per-sample CoverM long table to produce
per-vOTU and per-sample AMG summaries.

Inputs
------
--annotations  ANNOTATE/annotations_with_flags.tsv
                  gene-level: query_id, scaffold, amg_flags, is_transposon,
                  kofam_id, pfam_id, dbcan_id, merops_id, ...
--summary      SUMMARIZE/summarized_genomes.tsv
                  wide AMG distillate: 6 metadata cols
                  (gene_id, gene_description, topic_ecosystem, category,
                  subcategory, pathway) followed by one column per vOTU
                  containing the gene-call count of that KO in that vOTU.
--quant-long   QUANTIFY/votu_long.tsv
                  vOTU x sample long table (count, tpm, covered_fraction).

Outputs (all TSVs)
------------------
--out-votu-amg
    Long: vOTU x KO from `summarized_genomes.tsv` where count > 0.
    Cols: vOTU, gene_id, n_genes, gene_description, topic_ecosystem,
          category, subcategory, pathway.
--out-votu-amg-genes
    Long: gene-call rows from `annotations_with_flags.tsv` where
    amg_flags is non-empty and is_transposon is false. The DRAM-v
    'scaffold' is the vOTU representative.
    Cols: vOTU, query_id, gene_number, rank, amg_flags,
          kofam_id, kofam_EC, pfam_id, dbcan_id, merops_id.
--out-sample-amg
    Long: sample x KO. sum_count and sum_tpm aggregate
        Sigma_{vOTU} (gene_count_in_vOTU * vOTU_quant_in_sample)
    Cols: sample, gene_id, gene_description, category, subcategory,
          pathway, sum_count, sum_tpm.
--out-sample-category
    Wide: rows=sample, cols=category, values=sum_tpm aggregated
    across all KOs in that category and all vOTUs.
--out-amg-stats
    One row per sample: n_amg_kos_detected, n_amg_categories,
    n_amg_votus_present (count > 0 in that sample).
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annotations",         required=True)
    p.add_argument("--summary",             required=True)
    p.add_argument("--quant-long",          required=True)
    p.add_argument("--out-votu-amg",        required=True)
    p.add_argument("--out-votu-amg-genes",  required=True)
    p.add_argument("--out-sample-amg",      required=True)
    p.add_argument("--out-sample-category", required=True)
    p.add_argument("--out-amg-stats",       required=True)
    return p.parse_args()


def load_quant_long(path):
    """Return {vOTU: {sample: (count, tpm)}} and the ordered sample list."""
    quant = defaultdict(dict)
    samples_seen = []
    seen = set()
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            s = row["sample"]
            if s not in seen:
                seen.add(s)
                samples_seen.append(s)
            quant[row["vOTU"]][s] = (
                int(float(row["count"])),
                float(row["tpm"]),
            )
    return quant, sorted(samples_seen)


def load_summary(path):
    """
    Parse summarized_genomes.tsv.

    Returns
    -------
    ko_meta : dict gene_id -> dict of metadata
    votu_ko : dict vOTU -> dict gene_id -> int count (count > 0 only)
    """
    ko_meta = {}
    votu_ko = defaultdict(dict)
    META_COLS = ("gene_id", "gene_description", "topic_ecosystem",
                 "category", "subcategory", "pathway")
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # First six columns must be metadata in this order (DRAM-v convention).
        if tuple(header[:6]) != META_COLS:
            sys.exit(f"Unexpected header in {path}: first 6 cols were {header[:6]}, "
                     f"expected {META_COLS}")
        votu_cols = header[6:]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            ko = parts[0]
            ko_meta[ko] = {
                "gene_description": parts[1],
                "topic_ecosystem":  parts[2],
                "category":         parts[3],
                "subcategory":      parts[4],
                "pathway":          parts[5],
            }
            for v, cell in zip(votu_cols, parts[6:]):
                if not cell:
                    continue
                try:
                    n = int(float(cell))
                except ValueError:
                    continue
                if n > 0:
                    votu_ko[v][ko] = n
    return ko_meta, votu_ko


def write_votu_amg(path, ko_meta, votu_ko):
    cols = ["vOTU", "gene_id", "n_genes", "gene_description",
            "topic_ecosystem", "category", "subcategory", "pathway"]
    with open(path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for v in sorted(votu_ko):
            for ko in sorted(votu_ko[v]):
                m = ko_meta[ko]
                out.write("\t".join([
                    v, ko, str(votu_ko[v][ko]),
                    m["gene_description"], m["topic_ecosystem"],
                    m["category"], m["subcategory"], m["pathway"],
                ]) + "\n")


def write_votu_amg_genes(path, annotations_path):
    """Stream the gene-level annotations and emit AMG-flagged, non-transposon rows."""
    out_cols = ["vOTU", "query_id", "gene_number", "rank", "amg_flags",
                "kofam_id", "kofam_EC", "pfam_id", "dbcan_id", "merops_id"]
    n_in = n_out = 0
    with open(annotations_path) as fh, open(path, "w") as out:
        reader = csv.DictReader(fh, delimiter="\t")
        out.write("\t".join(out_cols) + "\n")
        for row in reader:
            n_in += 1
            if not row.get("amg_flags"):
                continue
            if (row.get("is_transposon") or "").lower() == "true":
                continue
            out.write("\t".join([
                row["scaffold"], row["query_id"], row.get("gene_number", ""),
                row.get("rank", ""), row["amg_flags"],
                row.get("kofam_id", ""), row.get("kofam_EC", ""),
                row.get("pfam_id", ""), row.get("dbcan_id", ""),
                row.get("merops_id", ""),
            ]) + "\n")
            n_out += 1
    return n_in, n_out


def aggregate_sample_amg(ko_meta, votu_ko, quant):
    """
    For each (sample, KO), sum gene_count_in_vOTU * vOTU_quant_in_sample
    across vOTUs that contain that KO and are detected in that sample.

    Returns dict (sample, ko) -> [sum_count, sum_tpm].
    """
    agg = defaultdict(lambda: [0, 0.0])
    for v, kos in votu_ko.items():
        sample_quant = quant.get(v)
        if not sample_quant:
            continue
        for ko, n_genes in kos.items():
            for sample, (c, t) in sample_quant.items():
                if c == 0 and t == 0.0:
                    continue
                key = (sample, ko)
                agg[key][0] += n_genes * c
                agg[key][1] += n_genes * t
    return agg


def write_sample_amg(path, agg, ko_meta):
    cols = ["sample", "gene_id", "gene_description", "category",
            "subcategory", "pathway", "sum_count", "sum_tpm"]
    with open(path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for (s, ko) in sorted(agg):
            m = ko_meta[ko]
            sc, st = agg[(s, ko)]
            out.write("\t".join([
                s, ko, m["gene_description"], m["category"],
                m["subcategory"], m["pathway"],
                str(sc), f"{st:.6g}",
            ]) + "\n")


def write_sample_category(path, agg, ko_meta, samples):
    """Wide: rows=sample, cols=category, values=sum_tpm aggregated over KOs."""
    by_sc = defaultdict(float)   # (sample, category) -> sum_tpm
    cats = set()
    for (s, ko), (_sc, st) in agg.items():
        cat = ko_meta[ko]["category"] or "(uncategorized)"
        by_sc[(s, cat)] += st
        cats.add(cat)
    cats = sorted(cats)
    with open(path, "w") as out:
        out.write("sample\t" + "\t".join(cats) + "\n")
        for s in samples:
            row = [f"{by_sc[(s, c)]:.6g}" for c in cats]
            out.write(s + "\t" + "\t".join(row) + "\n")


def write_amg_stats(path, agg, ko_meta, votu_ko, quant, samples):
    """Per-sample: n_amg_kos_detected, n_amg_categories, n_amg_votus_present."""
    kos_per_sample = defaultdict(set)
    cats_per_sample = defaultdict(set)
    for (s, ko), (sc, _st) in agg.items():
        if sc <= 0:
            continue
        kos_per_sample[s].add(ko)
        cat = ko_meta[ko]["category"] or "(uncategorized)"
        cats_per_sample[s].add(cat)

    amg_votus = set(votu_ko)
    votus_per_sample = defaultdict(int)
    for v in amg_votus:
        for s, (c, _t) in quant.get(v, {}).items():
            if c > 0:
                votus_per_sample[s] += 1

    cols = ["sample", "n_amg_kos_detected", "n_amg_categories",
            "n_amg_votus_present"]
    with open(path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for s in samples:
            out.write("\t".join([
                s,
                str(len(kos_per_sample[s])),
                str(len(cats_per_sample[s])),
                str(votus_per_sample[s]),
            ]) + "\n")


def main():
    args = parse_args()

    for label, path in [("--annotations", args.annotations),
                        ("--summary",     args.summary),
                        ("--quant-long",  args.quant_long)]:
        if not Path(path).is_file():
            sys.exit(f"{label} not found: {path}")

    quant, samples = load_quant_long(args.quant_long)
    ko_meta, votu_ko = load_summary(args.summary)

    write_votu_amg(args.out_votu_amg, ko_meta, votu_ko)
    n_in, n_out = write_votu_amg_genes(args.out_votu_amg_genes, args.annotations)

    agg = aggregate_sample_amg(ko_meta, votu_ko, quant)
    write_sample_amg(args.out_sample_amg, agg, ko_meta)
    write_sample_category(args.out_sample_category, agg, ko_meta, samples)
    write_amg_stats(args.out_amg_stats, agg, ko_meta, votu_ko, quant, samples)

    sys.stderr.write(
        f"AMG summary: {len(votu_ko)} AMG-bearing vOTUs x {len(ko_meta)} KOs; "
        f"{n_out}/{n_in} flagged-non-transposon gene calls; "
        f"{len(samples)} samples; {len(agg)} (sample, KO) cells.\n"
    )


if __name__ == "__main__":
    main()
