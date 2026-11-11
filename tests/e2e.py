#!/usr/bin/env python3
"""End-to-end test: simulated cohorts -> PoN files (dmd_pon / dmd_target_prep formats) ->
`dmdseek build-reference` -> `dmdseek call` on held-out BAMs (with a fake sawfish) -> checks."""
import csv
import gzip
import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import sim_bams as sb                                   # noqa: E402
from dmdseek import bamio                               # noqa: E402
from dmdseek.reference import crude_scale, mode_level   # noqa: E402

W = Path(sys.argv[1])
EXON_BED = Path(sys.argv[2])
W.mkdir(parents=True, exist_ok=True)
PROF = sb.profiles()
ex = sorted((int(l.split()[1]), int(l.split()[2]), int(l.split()[3].split('_')[1])) for l in open(EXON_BED))
EX = {n: (a, b) for a, b, n in ex}


def span(first, last):
    """genomic interval covering exons first..last with intronic breakpoints (minus strand)."""
    lo = (EX[last][0] + EX[last + 1][1]) // 2 if last < 79 else EX[79][0] - 3000
    hi = (EX[first - 1][0] + EX[first][1]) // 2 if first > 1 else EX[1][1] + 3000
    return lo, hi


def depth_arrays(bam_path, starts, ends):
    bam = bamio.open_bam(bam_path)
    w0, w1 = int(min(starts[0], min(a for a, b, n in ex))), int(max(ends[-1], max(b for a, b, n in ex)))
    cov, _ = bamio.coverage(bam, w0, w1)
    D = bamio.interval_means(cov, w0, starts, ends)
    E = bamio.interval_means(cov, w0, [a for a, b, n in ex], [b for a, b, n in ex])
    st = bamio.index_stats(bam)
    return D, E, st


def events_tsv(path, rows):
    with open(path, 'w') as f:
        f.write('sample_id\ttier\ttype\tstart\tend\texons\tpartial_exons\n')
        for sid, typ, a, b in rows:
            f.write(f'{sid}\tA_CONFIRMED\t{typ}\t{a}\t{b}\tex1\t.\n')


# ------------------------------------------------------------------ WGS reference cohort
rng = np.random.default_rng(11)
ws, we = np.arange(sb.WGS_WIN[0], sb.WGS_WIN[1] - 199, 200), None
we = ws + 200
cohort = W / 'cohort_wgs'
cohort.mkdir(exist_ok=True)
ids, sexes, Ds, Es, evrows = [], [], [], [], []
for i in range(20):
    sex = 'M' if i % 2 == 0 else 'F'
    sid = f'W{i:02d}'
    evs = []
    if i == 3:
        lo, hi = span(45, 50); evs = [('DEL', lo, hi, 1)]; evrows.append((sid, 'DEL', lo, hi))
    if i == 8:
        lo, hi = span(3, 7); evs = [('DEL', lo, hi, 1)]; evrows.append((sid, 'DEL', lo, hi))
    lo, hi = 32_969_276, 32_970_876        # common intronic deletion in 4 samples
    if i in (2, 5, 9, 14):
        evs = evs + [('DEL', lo, hi, 1)]; evrows.append((sid, 'DEL', lo, hi))
    p = cohort / f'{sid}.bam'
    if not p.exists():
        sb.simulate(p, 'wgs', sex, rng.uniform(22, 34), evs, seed=100 + i, prof=PROF)
    D, E, _ = depth_arrays(p, ws, we)
    ids.append(sid); sexes.append(sex); Ds.append(D); Es.append(E)
D = np.vstack(Ds); E = np.vstack(Es)
flank = (we <= sb.DMD_START) | (ws >= sb.DMD_END)
fl = np.median(D[:, flank], axis=1)
np.savez_compressed(W / 'CJH_DMD_PoN.npz', sample_ids=np.array(ids), sex=np.array(sexes), start=ws, end=we,
                    raw_depth=D, normalized_depth=D / fl[:, None], flank_depth=fl)
np.savez_compressed(W / 'CJH_DMD_PoN_exons.npz', sample_ids=np.array(ids), exons=np.array([f'exon_{n:02d}' for a, b, n in ex]),
                    starts=np.array([a for a, b, n in ex]), ends=np.array([b for a, b, n in ex]),
                    raw_depth=E, normalized_depth=E / fl[:, None])
events_tsv(W / 'wgs_events.tsv', evrows)
# intronic C-tier rows for the recurrent-event annotation
with open(W / 'wgs_events.tsv', 'a') as f:
    for sid, typ, a, b in evrows:
        if a == 32_969_276:
            f.write(f'{sid}\tC_INTRONIC\t{typ}\t{a}\t{b}\t.\t.\n')

# ------------------------------------------------------------------ targeted reference cohort
ts = np.arange(sb.TGT_WIN[0], sb.TGT_WIN[1], 200)
te = np.minimum(ts + 200, sb.TGT_WIN[1])
cohort = W / 'cohort_tgt'
cohort.mkdir(exist_ok=True)
ids, sexes, batches, Ds, Es, offs, evrows = [], [], [], [], [], [], []
hot = [(45, 47), (45, 50), (51, 51), (48, 50), (45, 47), (52, 52)]
for i in range(24):
    b = 'B1' if i < 12 else 'B2'
    sex = 'M' if i % 3 else 'F'
    sid = f'T{i:02d}'
    evs = []
    if i % 2 == 0:                                       # positive-rich cohort
        f_, l_ = hot[(i // 2) % len(hot)]
        lo, hi = span(f_, l_)
        evs = [('DEL', lo, hi, 1)]
        evrows.append((sid, 'DEL', lo, hi))
    p = cohort / f'{sid}.bam'
    if not p.exists():
        sb.simulate(p, 'targeted', sex, rng.uniform(25, 80), evs, seed=300 + i, bg_reads=60000, batch=b, prof=PROF)
    Dd, Ee, st = depth_arrays(p, ts, te)
    off = sum(st.values()) - st.get('chrX', 0) - st.get('chrY', 0)
    ids.append(sid); sexes.append(sex); batches.append(b); Ds.append(Dd); Es.append(Ee); offs.append(off)
D = np.vstack(Ds); E = np.vstack(Es); batches = np.array(batches)
s0 = np.array([crude_scale(d) for d in D])
N0 = D / s0[:, None]
scale = np.zeros(len(ids))
for b in ('B1', 'B2'):
    m = batches == b
    prof = np.nanmedian(N0[m], axis=0)
    for i in np.flatnonzero(m):
        ok = (prof > 0.2) & (N0[i] > 0.15)
        scale[i] = s0[i] * mode_level((N0[i] / prof)[ok])
np.savez_compressed(W / 'DMD_target_PoN.npz', sample_ids=np.array(ids), sex=np.array(sexes), batch=batches,
                    family=np.array(ids), prep_flag=np.array(['PASS'] * len(ids)), start=ts, end=te,
                    raw_depth=D, normalized_depth=D / scale[:, None], flank_depth=scale)
np.savez_compressed(W / 'DMD_target_PoN_exons.npz', sample_ids=np.array(ids),
                    exons=np.array([f'exon_{n:02d}' for a, b, n in ex]),
                    starts=np.array([a for a, b, n in ex]), ends=np.array([b for a, b, n in ex]),
                    raw_depth=E, normalized_depth=E / scale[:, None])
with open(W / 'DMD_target_input_QC.tsv', 'w') as f:
    f.write('sample_id\tofftarget_reads\tscale\n')
    for sid, o, sc in zip(ids, offs, scale):
        f.write(f'{sid}\t{o}\t{sc}\n')
events_tsv(W / 'tgt_events.tsv', evrows)
print('[E2E] cohorts simulated')

# ------------------------------------------------------------------ test samples (held out)
tests = {
    'wgs_M_del45_50': ('wgs', 'M', 30, [('DEL',) + span(45, 50) + (1,)], {}),
    'wgs_F_het_del45_47': ('wgs', 'F', 30, [('DEL',) + span(45, 47) + (1,)], {}),
    'wgs_M_normal_common_intronic': ('wgs', 'M', 28, [('DEL', 32_969_276, 32_970_876, 1)], {}),
    'wgs_F_dup2_7': ('wgs', 'F', 30, [('DUP',) + span(2, 7) + (1,)], {}),
    'tgt_M_del3_7': ('targeted', 'M', 50, [('DEL',) + span(3, 7) + (1,)], {'batch': 'B1'}),
    'tgt_F_het_del51_54': ('targeted', 'F', 60, [('DEL',) + span(51, 54) + (1,)], {'batch': 'B2'}),
    'tgt_F_dup1_29': ('targeted', 'F', 60, [('DUP',) + span(1, 29) + (1,)], {'batch': 'B2'}),
    'tgt_M_normal': ('targeted', 'M', 45, [], {'batch': 'B1'}),
    'tgt_M_wholegene': ('targeted', 'M', 45, [('DEL', sb.TGT_WIN[0] - 40000, sb.TGT_WIN[1] + 40000, 1)], {'batch': 'B1'}),
    'tgt_F_small_het57_nojunction': ('targeted', 'F', 60, [('DEL', (EX[57][0] + EX[57][1]) // 2 - 2000,
                                                            (EX[57][0] + EX[57][1]) // 2 + 2000, 1)], {'batch': 'B1'}),
}
tdir = W / 'tests'
tdir.mkdir(exist_ok=True)
fake = W / 'fake_vcf'
fake.mkdir(exist_ok=True)
for k, (name, (mode, sex, dep, evs, kw)) in enumerate(tests.items()):
    p = tdir / f'{name}.bam'
    if not p.exists():
        sb.simulate(p, mode, sex, dep, evs, seed=900 + k, prof=PROF,
                    bg_reads=(2_000_000 if mode == 'wgs' else 70000), **kw)
    # fake Sawfish output: large events as BND pairs, small as DEL/DUP records; none for the 'nojunction' case
    with gzip.open(fake / f'{name}.vcf.gz', 'wt') as f:
        f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n')
        for j, (t, a, b, n) in enumerate(evs):
            if 'nojunction' in name or 'wholegene' in name:
                continue
            if b - a > 100_000:
                alo = f'N[chrX:{b}[' if t == 'DEL' else f']chrX:{b}]N'
                ahi = f']chrX:{a}]N' if t == 'DEL' else f'N[chrX:{a}['
                f.write(f'chrX\t{a}\tsawfish:0:{j}:0:0:0\tN\t{alo}\t999\tPASS\tSVTYPE=BND\tGT\t{"1" if sex == "M" else "0/1"}\n')
                f.write(f'chrX\t{b}\tsawfish:0:{j}:0:0:1\tN\t{ahi}\t999\tPASS\tSVTYPE=BND\tGT\t{"1" if sex == "M" else "0/1"}\n')
            else:
                f.write(f'chrX\t{a}\tsawfish:0:{j}:0:0\tN\t<{t}>\t999\tPASS\tSVTYPE={t};END={b}\tGT\t'
                        f'{"1" if sex == "M" else "0/1"}\n')
sawfish = W / 'fake_sawfish'
sawfish.write_text(f'''#!/usr/bin/env bash
set -e
if [[ "$1" == --version ]]; then echo "sawfish 2.0.4 (fake for tests)"; exit 0; fi
cmd=$1; shift
while [[ $# -gt 0 ]]; do case $1 in --output-dir) out=$2; shift 2;; --sample) smp=$2; shift 2;; --bam) bam=$2; shift 2;; *) shift;; esac; done
mkdir -p "$out"
if [[ $cmd == discover ]]; then samtools_free=1; echo "$bam" > "$out/bam.txt"; exit 0; fi
name=$(basename "$(dirname "$smp")"); name=${{name%.sawfish}}
cp "{fake}/$name.vcf.gz" "$out/genotyped.sv.vcf.gz"
''')
sawfish.chmod(sawfish.stat().st_mode | stat.S_IEXEC)
(W / 'GRCh38_fake.fa').write_text('>chrX\nN\n')
print('[E2E] test BAMs ready')
json.dump({k: [list(map(str, e)) for e in v[3]] for k, v in tests.items()}, open(W / 'tests_truth.json', 'w'))
