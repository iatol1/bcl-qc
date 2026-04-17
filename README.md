# bcl-qc

Nextflow DSL2 pipeline for demultiplexing, aligning, and QC-summarising
Illumina clinical NGS runs on a Dragen server. A direct port of
[ucladx/bcl-qc](https://github.com/ucladx/bcl-qc) `bclqc.py`, same tools and
same outputs — just wrapped in Nextflow for resumability, provenance, and
container-pinned reproducibility.

## What it does

```
run_dir/SampleSheet_*.csv ──► BCL_CONVERT (Dragen BCLConvert)
                                   │
                                   ▼
               fastq_list.csv ──► DRAGEN_ALIGN  (panel: PCP or HEME)
                                        │
                                        ▼
sampleinfo.tsv  ──► PICARD_QCSUM (Picard CollectHsMetrics + qcsum_metrics.pl)
                                        │
                                        ▼
run_dir/InterOp ──► SAVE_OCC_PF_PLOT    AGGREGATE_QCSUM ──► qcsum_mqc.csv
                                                │
                                                ▼
                                             MULTIQC (multiqc_report.html)
```

- Index flavor → panel mapping (from SampleSheet filename): `I10` → `PCP`,
  `U5N2_I10` → `HEME`. Exome flavor `I8N2_N10` is skipped for alignment,
  matching upstream. Override via `params.flavor_to_panel`.
- Dragen binary is resolved through `params.dragen_versions[params.dragen_version]`,
  so multiple Dragen installs (4.2/4.3/4.4) can coexist without clobbering `PATH`.

## Quick start

```bash
nextflow run . \
    --run_dir     /mnt/pns/runs/240401_A01234_0001_BXYZ \
    --sampleinfo  /mnt/pns/sampleinfo/240401.tsv \
    --steps       demux,align,qc \
    --dragen_version 4.3 \
    -profile      singularity \
    -resume
```

Required arguments:

| flag              | meaning                                                         |
|-------------------|-----------------------------------------------------------------|
| `--run_dir`       | Illumina run directory (contains `SampleSheet_*.csv`, InterOp/) |
| `--sampleinfo`    | TSV with columns `Samples`, `Panel` or `Tumor`, `BAM Path` (qc step only) |

Optional:

| flag                | default                  | meaning                              |
|---------------------|--------------------------|--------------------------------------|
| `--fastqs_dir`      | `/staging/hot/reads`     | parent dir for FASTQ outputs         |
| `--bams_dir`        | `/mnt/pns/bams`          | parent dir for CRAM/QC outputs       |
| `--steps`           | `demux,align,qc`         | subset of stages to run              |
| `--dragen_version`  | `4.3`                    | key into `params.dragen_versions`    |

## Profiles

- `-profile conda` — reuses upstream's `config/conda_env.yaml` (picard 3.3.0,
  multiqc 1.23, interop 1.3.2, …). Dragen is expected on the host.
- `-profile docker` / `-profile singularity` — per-process biocontainer images
  pinned in `conf/containers.config`. Dragen still runs on the host.
- `-profile test` — stub run; no Dragen required. Combine with `-stub`.

## Provenance

Every run writes `<bams_dir>/pipeline_info/{report,timeline,dag,trace}` for
audit. Work dirs are content-hashed; `-resume` is safe.

## Verification

**Offline unit test** (no Dragen / Picard / MultiQC required):

```bash
bash tests/test_qcsum_metrics.sh
```

Runs `bin/qcsum_metrics.pl` against a frozen Picard HsMetrics fixture and
diffs against a golden `*.qcsum.txt`. Because `qcsum_metrics.pl` is copied
verbatim from ucladx/bcl-qc, a byte-identical result proves per-sample QC
parity with upstream.

**Demux parity** (operator, on Dragen host with demo BaseSpace data):

```bash
# 1. Upstream reference
git clone https://github.com/ucladx/bcl-qc /tmp/upstream
python /tmp/upstream/bclqc.py --run-dir demo_run --steps demux \
       --fastqs-dir /tmp/ref/reads --bams-dir /tmp/ref/bams

# 2. This pipeline
nextflow run . --run_dir demo_run --steps demux \
       --fastqs_dir /tmp/nf/reads --bams_dir /tmp/nf/bams -profile conda

# 3. Diff
diff <(md5sum /tmp/ref/reads/demo_run/I8/**/*.fastq.gz | sort) \
     <(md5sum /tmp/nf/reads/demo_run/I8/**/*.fastq.gz | sort)
```

Same command pattern for `--steps qc`: `qcsum_mqc.csv` should be byte-identical
since aggregation preserves upstream ordering (sorted by sample id).

## Layout

```
main.nf                    DSL2 workflow
nextflow.config            profiles, process defaults, provenance
conf/panels.config         Dragen paths + panel refs/BEDs
conf/containers.config     biocontainer pins
config/qcsum_config.yaml   QC thresholds per panel (verbatim from upstream)
config/multiqc_config.yaml MultiQC module order (verbatim from upstream)
bin/save_occ_pf_plot.py    InterOp → occupancy scatter (from bclqc.py)
bin/qcsum_metrics.pl       Picard HsMetrics → qcsum.txt (verbatim upstream)
bin/aggregate_qcsum.py     Concat per-sample qcsum.txt → qcsum_mqc.csv
demo_run/SampleSheet_I8.csv Demo samplesheet (upstream)
tests/                     Offline qcsum parity test + fixtures
```

## Upstream

Based on `ucladx/bcl-qc@f5b4381` (latest as of 2026-04-17).
