# CLAM Mitochondrial Genome Pipeline

CLAM is a Nextflow pipeline for mitochondrial genome analysis from human WGS-like
data and a separate, optional nuclear WES analysis branch. The mitochondrial
workflow remains the default and retains its existing NUMT-aware correction and
variant-calling behavior.

This repository is a modernization of an older local CLAM workflow. The current
goal is to make the pipeline reproducible, containerized, cluster-friendly, and
ready for future Seqera Launchpad support.

## Pipeline Logic

The core idea is simple: regardless of whether the starting point is FASTQ or an
already mapped BAM, CLAM rebuilds a consistent mitochondrial analysis space.

```mermaid
flowchart TD
    A[FASTQ pairs] --> C[Align reads to rCRS with GSNAP]
    B[Mapped BAM] --> B1[Sort BAM]
    B1 --> B2[Convert BAM to FASTQ]
    B2 --> C
    C --> D[Add MD tags with samtools calmd]
    D --> E{Analysis mode}
    E -->|with_NUMT_correction| F[Convert candidate MT BAM to FASTQ]
    F --> G[Re-align candidate reads to nuclear-only genome]
    G --> H[MitoScape classification]
    H --> I[Final NUMT-corrected MT BAM]
    E -->|without_NUMT_correction| J[Use rCRS MT BAM directly]
    I --> K[Sort, index, coverage]
    J --> K
    K --> L[Mutect2 mitochondrial calling]
    K --> M[mutserve calling]
    L --> L1[FilterMutectCalls mitochondrial filtering]
    L1 --> N[HaploGrep3 haplogroup assignment]
    M --> O[HaploGrep3 haplogroup assignment]
    L1 --> P[Unified heteroplasmy table]
    M --> P
```

## Why Remap?

Public datasets often arrive in different states: raw FASTQ, mapped BAMs, or
alignments produced against different human references. CLAM keeps the current
core workflow conservative by accepting:

- FASTQ pairs, which are mapped directly into the CLAM mitochondrial workflow.
- Mapped BAMs, which are first converted back to FASTQ and then remapped.

This makes downstream mitochondrial analyses more comparable because every
sample passes through the same mapping and filtering logic.

CLAM intentionally creates two MD-tagged mitochondrial BAM representations when
NUMT correction is enabled. The normal `samtools calmd` BAM keeps real read
bases for Mutect2 and mutserve. A separate MitoScape-only `samtools calmd -e`
BAM uses `=` for reference-matching bases, preserving compatibility with
MitoScape's MD parser without exposing that representation to variant callers.

## NUMT-Aware Strategy

Nuclear mitochondrial DNA segments, or NUMTs, can mimic mitochondrial signal in
WGS data. In `with_NUMT_correction` mode, CLAM uses the classic MitoScape logic:

1. Align reads to the mitochondrial reference.
2. Extract candidate mitochondrial reads.
3. Re-align those candidate reads to a nuclear-only human genome.
4. Use MitoScape, known NUMT regions, linkage disequilibrium data, and the
   bundled classifier model to retain likely true mitochondrial reads.

For `--analysis_mode without_NUMT_correction`, the pipeline intentionally skips the nuclear
remapping and MitoScape branch. This mode is useful for mitochondrial-only
analyses, WES-style tests, and quick checks, but it is not NUMT-corrected.

CLAM can also run an optional exploratory layer with
`--run_numt_exploration true`. This does not call novel NUMT insertions.
Instead, it subtracts the MitoScape-retained read IDs from the pre-MitoScape
candidate mitochondrial read set, extracts the MitoScape-rejected reads, and
summarizes where those reads align in the nuclear remapping BAM. The resulting
tables should be interpreted as NUMT-like or MitoScape-rejected evidence that
may guide follow-up review, not as confirmed breakpoint-level NUMT events.

## Supported Inputs

| Input type | Parameter | Expected pattern | Notes |
| --- | --- | --- | --- |
| FASTQ pairs | `--input_type fastq` | `data/*.fastq.gz` | Accepts one R1/R2 pair per sample or lane-named files such as `*_L001_R1_001.fastq.gz`; lane files are concatenated by mate in lane order. |
| BAM | `--input_type bam` | `data/*.bam` | BAMs are sorted, converted to FASTQ, and remapped inside CLAM. |

Example runs:

```bash
nextflow run clam.nf \
  --input 'data/*.fastq.gz' \
  --input_type fastq \
  --analysis_mode wgs_numt_correction
```

```bash
nextflow run clam.nf \
  --input 'data/*.bam' \
  --input_type bam \
  --analysis_mode wgs_numt_correction
```

## Analysis Modes

`--analysis_target` selects which analysis branch runs:

| Target | Behavior |
| --- | --- |
| `mtDNA` (default) | Existing mitochondrial workflow only. |
| `nDNA` | Nuclear WES workflow only, using BWA-MEM2 for FASTQ input or preparing the supplied BAM. |
| `both` | Runs both branches from the supplied FASTQ/BAM inputs. |

For `nDNA` and `both`, provide `--nuclear_reference` (GRCh38 FASTA) and
`--nuclear_target_bed` (capture BED path). By default, VEP annotation and MAF
generation are enabled and require `--nuclear_vep_cache` (offline GRCh38 VEP
cache). For a temporary alignment/calling/CNV test before the cache is ready,
set `--nuclear_annotate false`; this skips VEP, MAF generation, and the nuclear
HTML report. Re-enable it for annotated outputs. The BED remains an external input so it can be selected as a
Seqera file parameter; it is not vendored in this repository. Confirm that the
BED describes the exact enrichment design used. Twist Exome 2.0 alone and
Exome 2.0 plus Comprehensive Exome Spike-in have different target BEDs.

Choose one nuclear calling mode with `--nuclear_calling_mode`:

| Mode | Caller | Use |
| --- | --- | --- |
| `tumor_only` | Mutect2 | Somatic candidate calls when only tumor samples are supplied. |
| `paired` | Mutect2 | Somatic calls with matched tumor-normal samples and a tab-separated manifest. |
| `germline` | HaplotypeCaller | Germline-style variant calls from each sample. Variants in T cells may also be acquired somatic mutations and are not proven constitutional by this mode. |

The paired manifest requires the columns `sample_id`, `role`, and `pair_id`,
with exactly one `tumor` and one `normal` row per pair. Tumor-only and germline
modes also run CNVkit with a flat reference; paired mode uses the matched normal
for CNVkit. Optional Mutect2 germline-resource and panel-of-normals VCF paths,
plus optional gnomAD and COSMIC VCF paths, can be supplied separately with
their matching `.tbi` indexes. COSMIC requires an appropriately licensed
resource. The nDNA outputs include filtered
VCFs, VEP-annotated VCFs, MAF-compatible tables, CNVkit results, and a per-sample
HTML summary. Structural-variant calling is not included in this initial WES
branch.

Example nuclear-only run (replace every example path with an accessible path):

```bash
nextflow run clam.nf \
  --input '/path/to/fastqs/*.fastq.gz' \
  --input_type fastq \
  --analysis_target nDNA \
  --nuclear_reference /path/to/GRCh38.primary_assembly.genome.fa \
  --nuclear_target_bed /path/to/twist_capture_targets_hg38.bed \
  --nuclear_vep_cache /path/to/vep_cache \
  --nuclear_annotate true \
  --nuclear_calling_mode germline
```

The reference/index and BED must use compatible GRCh38 contig names and
coordinates. For BAM input, the BAM must already be aligned to the supplied
reference build.

| Mode | NUMT correction | Internal references | Intended use |
| --- | --- | --- | --- |
| `with_NUMT_correction` | Yes | rCRS plus GRCh38 nuclear-only remapping | Runs nuclear remapping and MitoScape filtering. |
| `without_NUMT_correction` | No | rCRS only | Skips nuclear remapping and MitoScape filtering. |

Large nuclear reference genomes and nuclear-only GSNAP indexes are not bundled
in the containers. The small rCRS GSNAP index is bundled in `clam-core` to avoid
GSNAP index-version mismatches during mitochondrial alignment. The normalized
rCRS FASTA used by mutserve is also indexed inside the image, so mutserve does
not need to write index files into the read-only container filesystem.
Nuclear-only GSNAP indexes for NUMT correction are configured in
`conf/genomes.config` and currently point to IEO cluster reference paths.
The `wgs_numt_correction` mode always uses a GRCh38 nuclear-only index built
with `clam-core:0.1.2` / GSNAP 2025-07-31 from
`GRCh38.primary_assembly.nuclear_only.fa`, derived from GENCODE release 49
`GRCh38.primary_assembly.genome.fa.gz`, and stored at:

```text
/hpcnfs/scratch/ED/genome/clam_refs/gsnap_GRCh38_gencode49
```

## Containers

The pipeline is designed to run with Singularity/Apptainer on the cluster while
using Docker/OCI images as the source:

| Image | Purpose |
| --- | --- |
| `ghcr.io/aledavini7/clam-core:0.1.2` | GSNAP/GMAP, samtools, htslib/bgzip, MitoScape, mutserve, HaploGrep3, a current-format rCRS GSNAP index, indexed rCRS FASTA resources, and bundled small CLAM resources. |
| `ghcr.io/aledavini7/clam-mitoscape:0.1.0` | Java 8 runtime for the Spark-based MitoScape NUMT classification step. |
| `ghcr.io/aledavini7/clam-mutect2:0.1.0` | GATK4/Mutect2 runtime. |
| `ghcr.io/aledavini7/clam-annotation:0.1.1` | Lightweight Python runtime for mitochondrial variant summary and heteroplasmy tables. |
| `ghcr.io/aledavini7/clam-nuclear:0.1.0` | BWA-MEM2, samtools, CNVkit, Ensembl VEP 112, and Python for nuclear WES. |

Both images are built for `linux/amd64`.
MitoScape runs in its own Java 8 container because its bundled Spark/Scala
runtime is not fully compatible with the Java 17 runtime used by `clam-core`.

Container build definitions live in:

- `containers/clam-core/`
- `containers/mitoscape/`
- `containers/mutect2/`
- `containers/annotation/`
- `containers/nuclear/`

The GitHub Actions workflow in `.github/workflows/build-clam-core.yml` builds
and publishes the container images to GHCR.

## Repository Layout

Active Nextflow modules live in `modules/` and are the only module files imported
by `clam.nf`. Older exploratory CLAM modules have been moved to
`legacy/modules/` so the current workflow stays easier to inspect while the old
implementation remains traceable.

## Seqera Launchpad

The repository includes `nextflow_schema.json` in the root directory. Seqera
uses this file to build the Launchpad parameter form for `input`, `input_type`,
`analysis_mode`, and `outdir`.

Suggested Launchpad settings for the first test runs:

| Field | Value |
| --- | --- |
| Pipeline | `https://github.com/aledavini7/clam-mitogenome-pipeline` |
| Revision | `main` or a pinned commit SHA |
| Main script | inferred from `manifest.mainScript = 'clam.nf'` |
| Config profile | leave empty for the current IEO SLURM/Singularity defaults |
| Work directory | an accessible cluster work path, for example `/hpcscratch/ieo/ieo5898/clam-mitogenome-pipeline-work` |

Minimal run parameters:

```json
{
  "input": "/path/to/data/*_{R1,R2}.fastq.gz",
  "input_type": "fastq",
  "analysis_mode": "wgs_numt_correction",
  "outdir": "results"
}
```

For BAM input:

```json
{
  "input": "/path/to/data/*.bam",
  "input_type": "bam",
  "analysis_mode": "wgs_numt_correction",
  "outdir": "results"
}
```

## Main Outputs

For each sample, CLAM writes results under:

```text
results/<sample_id>/
```

Core output categories include:

- BAMs under `results/<sample_id>/bams/`: rCRS alignment BAMs,
  MD-tagged BAMs, MitoScape-corrected BAMs, final coordinate-sorted BAMs, and
  BAM index files
- coverage files under `results/<sample_id>/coverage/`
- raw and FilterMutectCalls-filtered Mutect2 VCFs under
  `results/<sample_id>/variants/`
- mutserve VCFs and summaries under `results/<sample_id>/variants/`
- HaploGrep3 haplogroup reports under `results/<sample_id>/haplogroups/`
- first-pass QC summaries under `results/<sample_id>/qc/`
- a self-contained HTML report under `results/<sample_id>/report/`, combining
  coverage, BAM QC, variant caller summaries, confidence tiers, haplogroups, and
  optional MITOMAP matches
- annotation TSVs under `results/<sample_id>/annotation/`: caller-specific
  Mutect2 and mutserve tables, an all-calls consensus table, and a
  confidence-filtered consensus table retaining `HIGH` and `MEDIUM` calls with
  final coverage, heteroplasmy estimates, Wilson confidence intervals, and
  confidence tiers
- lightweight maftools-compatible MAF files under
  `results/<sample_id>/annotation/mafs/`, generated from the same annotation
  TSVs. These first-pass MAFs use conservative mitochondrial placeholders
  (`Hugo_Symbol=MTDNA`, `Variant_Classification=RNA`) until richer MITOMAP/VEP
  annotation is added.
- optional MITOMAP-annotated TSVs under
  `results/<sample_id>/annotation/mitomap/` when `--mitomap_variant_table` is
  provided. CLAM expects a local MITOMAP-derived TSV/CSV file and performs exact
  allele matching when possible, falling back to position-level matching.
- optional NUMT-like exploration outputs under
  `results/<sample_id>/numt_exploration/` when
  `--run_numt_exploration true` is used in `wgs_numt_correction` mode. These
  include MitoScape-rejected read IDs, extracted rejected-read BAMs, nuclear
  alignment summaries, and window-level nuclear remapping summaries. These are
  exploratory outputs and should not be reported as confirmed NUMT insertions
  without independent breakpoint-level support.

MITOMAP is maintained as a human mitochondrial genome database and reports
published data on human mtDNA variation. Because MITOMAP resources and terms of
use should be managed explicitly by the user/institute, CLAM does not bundle a
MITOMAP table. Use `--mitomap_variant_table /path/to/mitomap_table.tsv` to add
this annotation layer.

## Current Development Status

Implemented in the modernized core:

- DSL2 Nextflow workflow in `clam.nf`
- FASTQ and BAM input modes
- analysis mode selection with `wgs_numt_correction` and `rcrs_only`
- containerized CLAM core runtime
- containerized Mutect2 runtime
- containerized MitoScape Java 8 runtime
- FilterMutectCalls mitochondrial filtering
- first-pass unified heteroplasmy summary table
- lightweight maftools-compatible MAF generation from CLAM annotation TSVs
- optional MITOMAP table annotation for CLAM annotation TSVs
- optional MitoScape-rejected / NUMT-like read exploration for WGS mode
- SLURM/Singularity-oriented configuration
- `nextflow_schema.json` for Seqera Launchpad

Still planned:

- polish the HTML report with richer plots and run-level JSON summaries
- add a second, stricter NUMT candidate layer based on discordant/split-read
  breakpoint clustering once validation criteria are finalized
- add richer VEP/vcf2maf-style consequence annotation once the cache/resource
  strategy is finalized
- add test profiles with small synthetic fixtures
- add Python and nf-test checks for merge, annotation, and reporting logic
- add nf-core-style metadata, documentation, and CI checks

## Development Notes

This repository is being intentionally modernized step by step. The preferred
pattern is:

1. Update one logical piece of the workflow.
2. Run a Nextflow preview or a small test.
3. Commit the change.
4. Move to the next module.

That keeps the old CLAM behavior traceable while making the new implementation
easier to review, run, and eventually launch through Seqera.
