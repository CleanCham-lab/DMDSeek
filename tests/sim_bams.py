#!/usr/bin/env python3
"""Synthetic PacBio-like BAMs (GRCh38 header) for end-to-end tests of dmdseek.

Reads are placed per 200-bp bin with Poisson starts proportional to a shared depth profile x copy
number; reads crossing a deletion are written as primary + supplementary with SA tags.
"""
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pysam

GRCH38 = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
          'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
          'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
          'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
          'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
CONTIGS = list(GRCH38)
TID = {c: i for i, c in enumerate(CONTIGS)}
HEADER = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': c, 'LN': GRCH38[c]} for c in CONTIGS],
          'RG': [{'ID': 'sim', 'SM': 'sim', 'PL': 'PACBIO'}]}
DMD_START, DMD_END = 31_119_221, 33_339_388
WGS_WIN = (DMD_START - 1_000_000, DMD_END + 1_000_000)
TGT_WIN = (31_097_676, 33_339_609)


def ar1(rng, n, rho, sd):
    e = rng.normal(0, sd * np.sqrt(1 - rho ** 2), n)
    x = np.empty(n)
    x[0] = rng.normal(0, sd)
    for t in range(1, n):
        x[t] = rho * x[t - 1] + e[t]
    return x


def profiles(seed=5):
    rng = np.random.default_rng(seed)
    nw = (WGS_WIN[1] - WGS_WIN[0]) // 200 + 1
    nt = (TGT_WIN[1] - TGT_WIN[0]) // 200 + 1
    wgs = np.exp(ar1(rng, nw, 0.98, 0.12))
    tgt = np.exp(ar1(rng, nt, 0.98, 0.35))
    tgt[rng.random(nt) < 0.02] *= 0.05                         # capture dropouts
    batch = {b: np.exp(ar1(rng, nt, 0.98, 0.12)) for b in ('B1', 'B2')}
    return wgs, tgt, batch


def _rec(name, chrom, pos, length, flag=0, sa=None, mapq=60):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    a.reference_id = TID[chrom]
    a.reference_start = int(pos)
    a.mapping_quality = mapq
    a.cigartuples = [(0, int(length))]
    a.query_sequence = None
    tags = [('RG', 'sim')]
    if sa:
        tags.append(('SA', sa))
    a.set_tags(tags)
    return a


def simulate(path, mode, sex, depth, events=(), seed=1, bg_reads=0, batch='B1', read_len=None, prof=None):
    """events: [(type, start, end, copies_affected)] e.g. ('DEL', a, b, 1) removes one copy."""
    rng = np.random.default_rng(seed)
    wgs_p, tgt_p, bmod = prof or profiles()
    win = WGS_WIN if mode == 'wgs' else TGT_WIN
    shape = wgs_p if mode == 'wgs' else tgt_p * bmod[batch]
    rl_mean = read_len or (15000 if mode == 'wgs' else 8000)
    copies = 1 if sex == 'M' else 2
    recs = []
    k = 0
    pad = 30000
    lo_w, hi_w = win[0] - pad, win[1] + pad

    def lam_at(centers):
        idx = np.clip((centers - win[0]) // 200, 0, len(shape) - 1)
        lam = depth / copies * shape[idx] * 200 / rl_mean
        if mode == 'targeted':
            lam = np.where((centers < win[0]) | (centers > win[1]), lam * 0.02, lam)
        return lam

    def emit(name, pieces):
        pieces = [(a, b) for a, b in pieces if b - a >= 500]
        if not pieces:
            return
        total = sum(b - a for a, b in pieces)
        offs, o = [], 0
        for a, b in pieces:
            offs.append(o)
            o += b - a
        def cig(i):
            a, b = pieces[i]
            pre, post = offs[i], total - offs[i] - (b - a)
            return (f'{pre}S' if pre else '') + f'{b - a}M' + (f'{post}S' if post else '')
        for i, (a, b) in enumerate(pieces):
            sa = ''.join(f'chrX,{pieces[j][0] + 1},+,{cig(j)},60,0;' for j in range(len(pieces)) if j != i) or None
            recs.append(_rec(name, 'chrX', a, b - a, 0 if i == 0 else 2048, sa))

    for hap in range(copies):
        dels = sorted((a, b) for t, a, b, n in events if t == 'DEL' and hap < n)
        dups = sorted((a, b) for t, a, b, n in events if t == 'DUP' and hap < n)
        # reference-coordinate copy of the haplotype, deleted intervals skipped
        nb = (hi_w - lo_w) // 200
        centers = lo_w + 200 * np.arange(nb) + 100
        lam = lam_at(centers)
        for a, b in dels:
            lam[(centers >= a) & (centers < b)] = 0
        for j in np.flatnonzero(rng.poisson(lam)):
            s = int(lo_w + 200 * j + rng.integers(0, 200))
            ln = int(max(1000, rng.lognormal(np.log(rl_mean), 0.3)))
            pieces, cur, rem = [], s, ln
            for a, b in dels:
                if a <= cur:
                    continue
                stop = min(cur + rem, a)
                pieces.append((cur, stop))
                rem -= stop - cur
                cur = b
                if rem <= 0:
                    break
            if rem > 0:
                pieces.append((cur, cur + rem))
            k += 1
            emit(f'r{seed}_{k}', pieces)
        # extra tandem copies: reads crossing the b -> a junction are split
        for a, b in dups:
            c2 = a + 200 * np.arange((b - a) // 200) + 100
            for j in np.flatnonzero(rng.poisson(lam_at(c2))):
                s = int(a + 200 * j + rng.integers(0, 200))
                ln = int(max(1000, rng.lognormal(np.log(rl_mean), 0.3)))
                pieces = [(s, min(b, s + ln))]
                if s + ln > b:
                    pieces.append((a, min(b, a + (s + ln - b))))
                k += 1
                emit(f'd{seed}_{k}', pieces)
    if bg_reads:
        lens = np.array([GRCH38[c] for c in CONTIGS], dtype=float)
        w = lens.copy()
        w[TID['chrX']] *= copies / 2
        w[TID['chrY']] *= (0.4 / 2) if sex == 'M' else 0.0
        w /= w.sum()
        counts = rng.multinomial(bg_reads, w)
        if sex == 'F':
            counts[TID['chrY']] += int(bg_reads * 2e-5)          # mismapped reads
        for c, n in zip(CONTIGS, counts):
            pos = np.sort(rng.integers(0, GRCH38[c] - 20000, n))
            for q in pos:
                if c == 'chrX' and win[0] - pad <= q <= win[1] + pad:
                    continue
                k += 1
                recs.append(_rec(f'b{seed}_{k}', c, q, 1000))
    recs.sort(key=lambda r: (r.reference_id, r.reference_start))
    path = Path(path)
    tmp = str(path) + '.tmp.bam'
    with pysam.AlignmentFile(tmp, 'wb', header=HEADER) as out:
        for r in recs:
            out.write(r)
    Path(tmp).rename(path)
    pysam.index(str(path))
    return len(recs)


if __name__ == '__main__':
    print(simulate(sys.argv[1], 'wgs', 'M', 30, [('DEL', 31_800_000, 32_000_000, 1)], bg_reads=1000))
