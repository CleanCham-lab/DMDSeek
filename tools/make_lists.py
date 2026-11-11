#!/usr/bin/env python3
"""Make a dmdseek BAM list (sample_id, bam, sex, batch, family) from a cohort manifest TSV and check
that every BAM and its index exist.  Files are only read, never modified.

  python3 make_lists.py --manifest DMD_target_manifest.tsv --bam-col mapped_bam -o tgt_bams.tsv
  python3 make_lists.py --manifest CJH_manifest.tsv --bam-col bam -o wgs_bams.tsv
"""
import argparse
import csv
import os
import sys

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--manifest', required=True)
ap.add_argument('--bam-col', default='bam', help='column holding the BAM path (targeted: mapped_bam)')
ap.add_argument('--exclude', nargs='*', default=[], help='text files with sample IDs to leave out')
ap.add_argument('-o', '--out', required=True)
a = ap.parse_args()

excl = set()
for p in a.exclude:
    excl |= {l.split()[0] for l in open(p) if l.strip() and not l.startswith('#')}
rows, bad = [], []
with open(a.manifest, encoding='utf8') as f:
    for r in csv.DictReader((l for l in f if not l.startswith('##')), delimiter='\t'):
        sid = r['sample_id'].strip()
        if sid in excl:
            continue
        bam = (r.get(a.bam_col) or '').strip()
        if not bam or bam == '.':
            bad.append((sid, f'no value in column {a.bam_col}'))
            continue
        if not os.path.isfile(bam):
            bad.append((sid, f'missing {bam}'))
            continue
        if not any(os.path.isfile(bam + s) for s in ('.bai', '.csi')) and \
                not os.path.isfile(bam[:-4] + '.bai'):
            bad.append((sid, f'no index for {bam}'))
            continue
        sex = (r.get('sex') or '').strip().upper()[:1]
        rows.append([sid, bam, sex if sex in ('M', 'F') else '', (r.get('batch') or '').strip(),
                     (r.get('family') or '').strip() or sid])
with open(a.out, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['sample_id', 'bam', 'sex', 'batch', 'family'])
    w.writerows(rows)
print(f'[lists] {len(rows)} BAMs -> {a.out}  (M={sum(r[2] == "M" for r in rows)}, F={sum(r[2] == "F" for r in rows)}, '
      f'no sex={sum(r[2] == "" for r in rows)}; excluded {len(excl)})')
for sid, why in bad:
    print(f'[lists] SKIPPED {sid}: {why}', file=sys.stderr)
sys.exit(1 if bad and not rows else 0)
