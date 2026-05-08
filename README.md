# nf-virome

Generic Nextflow workflow for viral identification, clustering, quantification, and host coupling from metagenomic assemblies.

## What it does

| Stage | Tool | Output |
| --- | --- | --- |
| **1. IDENTIFY** | geNomad → CheckV → length/quality filter | per-sample filtered viral fasta + per-gene TSV |
| **2. CLUSTER** | skani → leiden vOTU clustering | `votu_catalog.fa` + cluster table |
| **3. QUANTIFY** | CoverM mean-read mapping | per-sample × vOTU coverage / RPKM / TPM matrices + long table |
| **4. HOST_COUPLE** | minced → blast vs vOTUs | CRISPR-spacer-based host predictions |

Stage 4 runs only if the samplesheet's `bins_dir` column is populated for at least one sample.

## Quick start

```bash
nextflow run tpall/nf-virome \
  --samplesheet samples.csv \
  --outdir results \
  --genomad_db /path/to/genomad-db \
  --checkv_db  /path/to/checkv-db-v1.5 \
  -profile singularity,slurm
```

Sample sheet schema (CSV header required):

| Column | Required? | Notes |
| --- | --- | --- |
| `sample` | yes | unique id, used as the per-sample namespace |
| `contigs` | yes | path to the assembled contigs (`.fa`/`.fa.gz`/etc.) |
| `reads_1` | yes | forward reads (or single-end reads if `single_end=true`) |
| `reads_2` | conditional | reverse reads — required when `single_end=false`, ignored otherwise |
| `bins_dir` | optional | directory of bin FASTAs for HOST_COUPLE; leave blank to skip |
| `single_end` | optional | `true`/`false` (default `false`) |

See `assets/samplesheet_example.csv` for a minimal example.

## Outputs

```
results/
├── identify/
│   ├── identify_stats.tsv               (per-sample CheckV count summary)
│   ├── filtered/<sample>/<sample>_filtered.fna
│   ├── genomad/<sample>/<sample>_summary/
│   └── checkv/<sample>/
├── dramv_input/                         (flat handoff for downstream DRAM-v)
│   ├── fastas/<sample>_filtered.fna
│   └── genomad_genes/<sample>_virus_genes.tsv
├── cluster/
│   ├── votu_catalog.fa
│   └── votu_clusters.tsv
├── quantify/
│   ├── per_sample/<sample>.coverm.tsv
│   ├── votu_relab.tsv  votu_rpkm.tsv  votu_tpm.tsv
│   └── votu_long.tsv
└── host_couple/
    ├── spacer_db/, blast_hits.tsv
    └── host_summary.tsv
```

## Downstream: DRAM-v AMG annotation

`dramv_input/` is published flat (one file per sample, no nested dirs) so [tpall/DRAM](https://github.com/tpall/DRAM) Phase 2 can be invoked with simple globs:

```bash
nextflow run tpall/DRAM -r dev \
  --input_fasta results/dramv_input/fastas \
  --fasta_fmt "*.fna" \
  --genomad_genes "results/dramv_input/genomad_genes/*.tsv" \
  --use_dramv --call --annotate --summarize \
  -profile singularity
```

Or, for the production catalog-mode launch, run on `results/cluster/votu_catalog.fa` directly. The per-sample mode is mainly useful for DRAM-v development / phase testing where gene-id alignment matters.

## Profiles

- `standard` — local executor, bring your own resources
- `slurm` — HPC SLURM via `conf/slurm.config`; per-process labels `process_low/medium/high/long`
- `singularity` — enables singularity, autoMounts, common cache dir. Add HPC bind paths via `singularity.runOptions` in your environment config

Combine: `-profile singularity,slurm`.

## Defaults

| Param | Default | Source |
| --- | --- | --- |
| `--min_length` | 5000 | bp |
| `--keep_quality` | `Medium-quality,High-quality,Complete` | CheckV |
| `--votu_min_ani` | 95.0 | MIUViG |
| `--votu_min_af` | 85.0 | MIUViG (shorter sequence) |
| `--coverm_min_covered_fraction` | 0.70 | Nayfach 2021 / IMG/VR |
| `--coverm_min_read_pid` | 0.95 | Nayfach 2021 / IMG/VR |

## Citing

If you use nf-virome, please cite the underlying tools (geNomad, CheckV, skani, CoverM, minced, blast) as listed in the `versions.yml` files emitted next to each output.

## Origin

Extracted from [tpall/eluring-virome](https://github.com/tpall/eluring-virome) v0.1.0. The eluring-virome repo retains the cohort-specific AMG summary stage that consumes nf-virome's outputs.
