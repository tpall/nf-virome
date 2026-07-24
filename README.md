# nf-virome

Generic Nextflow workflow for viral identification, clustering, quantification, host coupling, lifestyle prediction, and AMG summary from metagenomic assemblies.

## What it does

| Stage | Tool | Output |
| --- | --- | --- |
| **1. IDENTIFY** | geNomad → CheckV → length/quality filter | per-sample filtered viral fasta + per-gene TSV |
| **2. CLUSTER** | skani → greedy centroid ANI clustering (aniclust-equivalent) | `votu_catalog.fa` + cluster table |
| **3. QUANTIFY** | CoverM read mapping | per-sample × vOTU count / TPM matrices + long table |
| **4. HOST_COUPLE** | minced → blast vs vOTUs | exact CRISPR-spacer host links |
| **4b. HOST_IPHOP** | iPHoP (CRISPR + blast + k-mer/codon RF) | calibrated host genus/genome calls |
| **5. LIFESTYLE** | PhaTYP (PhaBOX2) | per-vOTU temperate/virulent call |
| **6. SUMMARIZE** | join external DRAM-v AMG tables | vOTU / sample AMG summaries |

Stage 4 runs only if the samplesheet's `bins_dir` column is populated for at least one sample (skip with `--skip_host_couple`). Stage 4b runs only if `--iphop_db` is set. Stage 5 runs only if `--phabox_db` is set. Stage 6 runs only once the external DRAM-v tables exist (see below); it skips gracefully otherwise. Stages 4 and 4b are complementary: HOST_COUPLE gives high-confidence exact-spacer links against the study MAGs, HOST_IPHOP gives broad, calibrated coverage for the majority of vOTUs that carry no matching spacer.

## Quick start

Raw invocation:

```bash
nextflow run tpall/nf-virome \
  --samplesheet samples.csv \
  --outdir results \
  --genomad_db /path/to/genomad-db \
  --checkv_db  /path/to/checkv-db-v1.5 \
  -profile singularity,slurm
```

For the eluring cohorts, build the samplesheet from the shared manifest and launch via the wrapper (parks report/trace/timeline alongside the outputs, handles the optional-stage flags):

```bash
# 620-sample aging cohort (8 projects, one shared vOTU catalog).
# Run on the HPC so --check-files can flag single-end samples.
Rscript assets/build_samplesheet.R assets/samplesheet_aging620.csv \
  --samples assets/aging_samples_620.txt --eluring ../eluring --check-files

export GENOMAD_DB=/path/to/genomad-db CHECKV_DB=/path/to/checkv-db-v1.5
tmux new -s virome
bin/run_pipeline.sh assets/samplesheet_aging620.csv results_aging620
```

The 620-sample id list (`assets/aging_samples_620.txt`) is the aging-cohort membership from `shared/output/ps_aging.rds` (age ≥ 12), so the virome catalog is built over exactly the samples in the bacterial paper.

Sample sheet schema (CSV header required):

| Column | Required? | Notes |
| --- | --- | --- |
| `sample` | yes | unique id, used as the per-sample namespace |
| `contigs` | yes | path to the assembled contigs (`.fa`/`.fa.gz`/etc.) |
| `reads_1` | yes | forward reads (or single-end reads if `single_end=true`) |
| `reads_2` | conditional | reverse reads, required when `single_end=false`, ignored otherwise |
| `bins_dir` | optional | directory of bin FASTAs for HOST_COUPLE; leave blank to skip |
| `single_end` | optional | `true`/`false` (default `false`) |
| `project`, `cohort` | optional | provenance only; ignored by the pipeline |

`assets/build_samplesheet.R` emits this schema from `eluring/data/data_paths.csv`; see `assets/samplesheet_example.csv` for a minimal hand-written example.

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
│   ├── all_filtered.fna                 (sample-prefixed concatenation)
│   ├── votu_catalog.fa                  (representative sequences)
│   ├── votu_clusters.tsv                (representative → member)
│   └── votu_reps.txt
├── quantify/
│   ├── per_sample/<sample>_coverage.tsv
│   ├── votu_count_matrix.tsv  votu_tpm_matrix.tsv
│   └── votu_long.tsv
├── host_couple/
│   ├── all_spacers.fa  spacer_to_bin.tsv  bin_spacer_stats.tsv
│   ├── votu_host_links.tsv
│   └── votu_host_summary.tsv
├── host_iphop/                          (Stage 4b, if --iphop_db)
│   ├── iphop_host_genus.csv  iphop_host_genome.csv
│   └── Detailed_output_by_tool.csv
├── lifestyle/
│   └── votu_lifestyle.tsv               (PhaTYP temperate/virulent)
├── dramv/                               (external DRAM-v run, via bin/run_dramv.sh)
│   └── ANNOTATE/collected_trnas.tsv     (per-vOTU tRNA counts; free genome feature)
└── annotate/summary/                    (Stage 6, once DRAM-v has run)
    ├── votu_amg.tsv  votu_amg_genes.tsv
    ├── sample_amg_abundance.tsv  sample_amg_category_tpm.tsv
    └── amg_stats.tsv
```

## Host prediction (Stages 4 + 4b)

Two complementary callers. `HOST_COUPLE` (Stage 4) matches CRISPR spacers mined
from the study MAGs against the catalog: exact and high-confidence, but only for
the minority of vOTUs a spacer covers. `HOST_IPHOP` (Stage 4b) runs iPHoP, which
fuses CRISPR, blast homology, and genome k-mer/codon-usage similarity into a
calibrated genus/genome call with a confidence score, reaching the rest.

To land iPHoP's calls on the cohort's own bacteria instead of generic GTDB
references, augment the stock db with the study MAGs once (reads the refined bins
and per-sample GTDB-Tk summaries from the assembly outputs), then point
`--iphop_db` at the result:

```bash
BASE_DB=~/db/iphop_db OUT_DB=~/db/iphop_db_plus620 \
  bin/build_iphop_db.sh assets/aging_samples_620.txt ../eluring
```

## DRAM-v AMG annotation (Stage 6, two passes)

DRAM-v ([tpall/DRAM](https://github.com/tpall/DRAM) dev) runs SEPARATELY against
the finished catalog; Stage 6 only joins its tables to the quantification. Its
inputs are already published by discovery: the catalog (`cluster/votu_catalog.fa`)
and the per-sample geNomad genes (`dramv_input/genomad_genes/*_virus_genes.tsv`,
which DRAM-v matches to catalog contigs by their sample-prefixed ids). Run it via
the launcher, which uses the validated catalog-mode flags:

```bash
DRAM_CONFIG=~/dram_dbs.config bin/run_dramv.sh results_aging620
```

Then re-run discovery with `DRAMV2_OUTDIR` set; Nextflow `-resume` reuses the
cached work and only fires SUMMARIZE:

```bash
DRAMV2_OUTDIR=results_aging620/dramv bin/run_pipeline.sh \
  assets/samplesheet_aging620.csv results_aging620
```

(Equivalently, pass `--dramv2_outdir`, or point `--dramv2_annotations` /
`--dramv2_summary` at the two tables directly.)

Safe to run despite tpall/DRAM issue #13: the global `-T` prefilter inflates raw
VOG / V-flag counts but leaves the AMG endpoint unchanged (auxiliary_score is
geNomad/VirSorter-derived; a validation run measured 294 → 294 AMG calls). Stage 6
consumes only that AMG endpoint.

## Profiles

- `standard` — local executor, bring your own resources
- `slurm` — HPC SLURM via `conf/slurm.config`; per-process labels `process_low/medium/high/long`; heavy per-sample stages submit as one job array (`--array_size`)
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
| `--iphop_db` | null | iPHoP DB dir; enables Stage 4b HOST_IPHOP |
| `--phabox_db` | null | PhaBOX2 DB dir; enables Stage 5 LIFESTYLE |
| `--dramv2_outdir` | null | external DRAM-v run dir; enables Stage 6 SUMMARIZE |

## Citing

If you use nf-virome, please cite the underlying tools (geNomad, CheckV, skani, CoverM, minced, blast, PhaBOX/PhaTYP, DRAM-v) as listed in the `versions.yml` files emitted next to each output.

## Origin

Consolidates [tpall/eluring-virome](https://github.com/tpall/eluring-virome) v0.1.0: the generic discovery core was extracted here, and the cohort-specific AMG summary stage (SUMMARIZE) was folded back in as Stage 6. nf-virome is now the single canonical pipeline; eluring-virome is retired. DRAM-v itself remains external (tpall/DRAM `-r dev --use_dramv`).
