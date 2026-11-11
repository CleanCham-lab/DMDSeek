# DMDSeek — DMD exon deletion/duplication caller for PacBio HiFi BAMs

**DMDSeek** calls exon-level deletions and duplications in *DMD* (GRCh38, MANE NM_004006.3, 79 exons)
from one PacBio HiFi BAM, either whole-genome (WGS) or Twist targeted capture. It compares the
sample's depth with a fixed reference panel and confirms boundaries with breakpoint reads (Sawfish),
then writes an HTML report with an exon-level (MLPA-like) chart.

Research use only. Confirm reportable findings with an orthogonal method (MLPA, qPCR) and review the
reads in IGV.

## Before you start (public release)

The public bundle (`dmdseek-<version>-linux-x86_64-public.tar.gz`) contains the program only.
Two things have to be added before the first run; with both in place DMDSeek behaves exactly like
the version validated by the authors.

**1. Install Sawfish** (breakpoint caller, not bundled because it has its own licence)

```bash
conda install -c bioconda sawfish      # tested with sawfish 2.2.1
```

DMDSeek finds `sawfish` on PATH (or give `--sawfish /path/to/sawfish`). Without Sawfish the result is
depth-only: no tier A calls and less precise boundaries.

**2. Build a reference panel from your own cohort**

Calls are made by comparing each sample with a reference panel, so the panel must come from data like
yours (same sequencing chemistry; for targeted data, the same capture design and protocol).

| data | minimum | recommended |
|---|---|---|
| WGS | 8 males **and** 8 females | ≥ 20 per sex |
| targeted capture | 8 samples per capture batch | ≥ 20 per batch |

Samples carrying a DMD deletion/duplication may be included (near-zero bins are masked automatically;
known events can be masked with `--wgs-events` / `--target-events`).

List your cohort in a **sample sheet** - CSV, TSV or Excel (.xlsx, first sheet is read). Templates:
`dmdseek/templates/sample_sheet_template.csv` / `.xlsx`.

| column | | |
|---|---|---|
| `sample_id` | required | your sample name |
| `bam` | required | indexed GRCh38 BAM (`.bai`/`.csi` next to it): absolute path, path relative to the sheet, or a Windows path (`E:\data\S1.bam` becomes `/mnt/e/data/S1.bam` on WSL) |
| `sex` | optional | `M` / `F`; empty = inferred from the reads (recommended for WGS) |
| `batch` | optional | capture batch / pool for targeted data; one panel is built per batch |

Header names are matched loosely ("Sample ID", "BAM path", "Gender" also work).

```bash
./dmdseek/dmdseek build-reference --wgs-bams wgs_cohort.xlsx    -o dmdseek/reference/dmd_reference.npz   # WGS only
./dmdseek/dmdseek build-reference --target-bams tgt_cohort.csv  -o dmdseek/reference/dmd_reference.npz   # targeted only
./dmdseek/dmdseek build-reference --wgs-bams wgs_cohort.xlsx --target-bams tgt_cohort.csv \
    -o dmdseek/reference/dmd_reference.npz   # both in one panel file
```

Saved as `dmdseek/reference/dmd_reference.npz` it is used automatically; elsewhere, pass it with
`--reference my_panel.npz`.

**3. Check it on your own samples** with known MLPA results (see *Validation*).

## Quick start

```bash
tar -xzf dmdseek-1.0.0-linux-x86_64.tar.gz
./dmdseek/dmdseek call sample.bam --ref GRCh38.fa -o results/
# several BAMs: one sub-folder per sample + results/dmdseek_summary.tsv
./dmdseek/dmdseek call a.bam b.bam c.bam --ref GRCh38.fa -o results/ -t 8
# or a sample sheet (CSV/TSV/XLSX with columns sample_id, bam) so that sample names are yours
./dmdseek/dmdseek call --list samples.xlsx --ref GRCh38.fa -o results/ -t 8
# WGS or targeted is detected from the reads; -wgs / -target force it
./dmdseek/dmdseek call sample.bam --ref GRCh38.fa -o results/ -target
./dmdseek/dmdseek -help            # all commands;  ./dmdseek/dmdseek call -help  for one command
```

* Input: coordinate-sorted, indexed BAM (or CRAM with `--ref`) aligned to GRCh38 with `chr` names.
  WGS vs targeted and the sex are detected from the reads (override with `--mode`, `--sex`).
* `--ref` (GRCh38 FASTA with `.fai`) enables the breakpoint step. Without it the result is
  depth-only (tier B at best).
* Runs on Linux x86_64 / WSL2 with no other installation, apart from Sawfish and the reference
  panel described above.

## Output (per sample)

| file | content |
|---|---|
| `<sample>.dmdseek.report.html` | result, exon chart, gene and zoom plots, events, QC |
| `<sample>.dmdseek.events.tsv` | all events with tier, exons, reading frame, ratio, breakpoints |
| `<sample>.dmdseek.exons.tsv` | 79 exons: depth, ratio, z, state, review flag |
| `<sample>.dmdseek.summary.json` | machine-readable summary |
| `<sample>.dmdseek.bins.tsv.gz` | 200-bp depth ratios |

Tiers: **A** exonic, clean copy-number fit and a concordant breakpoint junction; **B** exonic, needs
review (no junction, complex junction or depth-only); **C** intronic with junction; **D** weak (table
only). *Exon-level review flags* mark single exons with unusual depth that are not called (e.g. small
heterozygous deletions in women, probe-site variants seen by MLPA).

## How it works

1. Depth per 200-bp bin and per exon with pysam (MAPQ ≥ 20; unmapped/secondary/QC-fail/duplicate
   reads excluded; deleted bases not counted) — the same code builds the reference.
2. Normalisation: WGS by the median depth of ±1 Mb flanks; targeted by the on-target depth level after
   dividing by the capture profile of the best-matching reference batch.
3. Ratio to the reference median of the same sex (WGS) or batch (targeted); sex-aware HMM
   (male CN 0/1/2, female CN 0/1/2/3, plus half steps for mosaicism).
4. Sawfish on a chrX:29–35.5 Mb regional BAM; DEL/DUP records and BND pairs define exact boundaries
   whenever the depth over the junction interval agrees.
5. Exons, reading frame (exon-length rule), tier, recurrent intronic variants from the reference cohort.

## Reference panel

`reference/dmd_reference.npz` holds only aggregated per-bin statistics (median, MAD, capture
profile) per sex / batch, the exon coordinates and recurrent intronic events. It contains no sample
identifiers and no per-sample depth.

A panel is valid only for data like the cohort it was built from (for targeted data: the same
capture design and protocol). `dmdseek info` shows the panels in use. How to build one, and how many
samples it needs, is described in *Before you start*.

Carriers in the reference cohort are fine: their event bins are masked (`--*-events`) and near-zero
bins are removed automatically.

## Validation

```bash
# truth sheet (CSV/TSV/XLSX): sample_id, type (DEL/DUP/NONE, or EXCLUDE), exons (e.g. 45~47)
./dmdseek/dmdseek validate --truth truth.xlsx --summary results/dmdseek_summary.tsv -o validation.tsv
```

Calling the samples that built the reference is optimistic. `tools/crossval_target.sh` gives an
honest estimate: K folds with families kept together; each fold is called with a reference built
without it (`tools/make_lists.py` turns a cohort manifest into the BAM list).

```bash
python3 tools/make_lists.py --manifest manifest.tsv --bam-col mapped_bam -o tgt_bams.tsv
EXON_BED=DMD_MANE_79exons.bed FASTA=GRCh38.fa TARGET_EVENTS=events.tsv \
    bash tools/crossval_target.sh ./dmdseek/dmdseek tgt_bams.tsv truth.tsv cv_out 5
```

## Building the executable

```bash
bash build_dmdseek.sh --exon-bed DMD_MANE_79exons.bed --sawfish $(which sawfish)   # bundle (Python, pysam, matplotlib inside)
build/dist/dmdseek/dmdseek build-reference ... -o build/dist/dmdseek/reference/dmd_reference.npz
bash build_dmdseek.sh --package                            # smoke test + tar.gz
```

Build on the oldest Linux you want to support (glibc compatibility). `--logo FILE` (png/svg) puts an
institutional logo at the top of every report; it is not part of this repository. Terminal colours
follow the terminal (`NO_COLOR=1` turns them off, `DMDSEEK_PLAIN=1` also drops icons).

## Authors

Developed at the Molecular Diagnostics Laboratory, Department of Laboratory Medicine, Seoul National
University Hospital.

* Seung Won Chae — developer
* Moon-woo Seong, MD, PhD — principal investigator

Affiliations: Department of Laboratory Medicine, Seoul National University Hospital; Seoul National
University Cancer Research Institute.

## License

DMDSeek is released under the MIT License (see `LICENSE`); citation information is in `CITATION.cff`.
The reference panel shipped in the release bundle is a separate data file built from patient samples
(aggregated statistics only).

## Third-party software

Bundled: Python, NumPy (BSD), matplotlib (PSF-based), pysam/htslib (MIT), Sawfish (Pacific
Biosciences; see `bin/sawfish-license`). Licence texts of the bundled Python packages are in `licenses/`. Sawfish: Saunders et al., *Bioinformatics* 2025.
