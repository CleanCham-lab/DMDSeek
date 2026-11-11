"""Reference panel: build (from cohort PoN files) and use (normalise one sample against it).

The reference file holds only aggregated per-bin statistics (median, MAD, capture profile) per sex /
panel, the exon coordinates and recurrent intronic events - no sample IDs and no per-sample depth.
"""
import csv
import datetime
import json
import sys
import math
import os
from collections import defaultdict
from pathlib import Path

import numpy as np

from . import __version__
from .constants import (CHROM, DMD_END, DMD_START, MIN_MAPQ, REFERENCE_FORMAT, TILE)
from .model import robust_mad

MIN_PANEL_N = 8          # samples needed for a panel (sex / batch)
MIN_BIN_N = 4            # unmasked samples needed at one bin (carrier-rich bins keep fewer)


# ---------------------------------------------------------------- helpers shared by build and call
def mode_level(R):
    """Most populated depth level (smoothed histogram peak of log2 ratio), refined by the median of
    bins within +/-20%. A DEL/DUP covering 30-49% of the gene cannot pull it (unlike the median)."""
    R = R[np.isfinite(R) & (R > 0)]
    if R.size < 100:
        return float(np.median(R)) if R.size else float('nan')
    L = np.log2(R)
    edges = np.arange(-3.0, 2.0 + 1e-9, 0.01)
    h, _ = np.histogram(L, bins=edges)
    k = np.exp(-0.5 * (np.arange(-24, 25) / 8.0) ** 2)
    m = (edges[:-1] + 0.005)[int(np.argmax(np.convolve(h, k / k.sum(), mode='same')))]
    sel = np.abs(L - m) < 0.263
    return float(2 ** np.median(L[sel])) if sel.sum() >= 100 else float(np.median(R))


def crude_scale(D):
    hi = np.nanpercentile(D, 90)
    return float(np.nanmedian(D[D > 0.1 * hi])) if hi > 0 else float('nan')


def load_exon_bed(path):
    ex = []
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                c, a, b, name = line.rstrip('\n').split('\t')[:4]
                num = int(''.join(ch for ch in name if ch.isdigit()))
                ex.append((num, int(a), int(b), name))
    ex.sort(key=lambda t: t[1])
    if len(ex) != 79:
        raise ValueError(f'{path}: {len(ex)} exons, expected 79 (MANE NM_004006.3)')
    by = {n: (a, b) for n, a, b, _ in ex}
    intron44 = by[44][0] - by[45][1]          # minus strand: exon 44 lies above exon 45
    if intron44 < 150_000:
        print(f'[WARN] {path}: intron 44 is {intron44:,} bp (GRCh38 MANE: ~248 kb) - '
              'is this the real NM_004006.3 GRCh38 exon BED?', file=sys.stderr)
    return ex


# ---------------------------------------------------------------- build
def _event_masks(paths, ids, starts, ends, tiers=('A', 'B', 'C')):
    """sample index -> boolean bin mask of called events (from dmd_deldup_caller events TSVs)."""
    pos = {s: i for i, s in enumerate(ids)}
    mask = np.zeros((len(ids), len(starts)), dtype=bool)
    n = 0
    for p in paths or []:
        with open(p, encoding='utf8') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                if r['sample_id'] in pos and r['tier'][:1] in tiers:
                    m = (ends > int(r['start'])) & (starts < int(r['end']))
                    mask[pos[r['sample_id']]] |= m
                    n += 1
    return mask, n


def _cn0_mask(norm, groups):
    """Near-zero values (< 0.3 x group 90th percentile): carriers in positive-rich cohorts."""
    mask = np.zeros(norm.shape, dtype=bool)
    for g in np.unique(groups):
        gm = groups == g
        with np.errstate(all='ignore'):
            p90 = np.nanpercentile(norm[gm], 90, axis=0)
            mask[gm] = (norm[gm] / np.maximum(p90, 1e-3) < 0.3) & (p90 > 0.15)
    return mask


def _stats(mat, mask):
    """Per-bin median / MAD of unmasked values. Bins where carriers were masked keep fewer values:
    with MIN_BIN_N..MIN_PANEL_N-1 values the MAD is floored at the panel's typical relative MAD
    (a MAD from a handful of values is unreliable); below MIN_BIN_N the bin is unusable."""
    m = np.where(mask, np.nan, mat)
    with np.errstate(all='ignore'):
        n = np.sum(np.isfinite(m), axis=0)
        med = np.nanmedian(m, axis=0)
        mad = robust_mad(m)
        rel = mad / med
        good = (n >= MIN_PANEL_N) & np.isfinite(rel) & (med > 0)
        typical = np.nanmedian(rel[good]) if good.any() else np.nanmedian(rel)
        few = (n < MIN_PANEL_N) & (n >= MIN_BIN_N)
        if np.isfinite(typical):
            mad[few] = np.fmax(mad[few], typical * med[few])
    med[n < MIN_BIN_N] = np.nan
    mad[n < MIN_BIN_N] = np.nan
    return med.astype(np.float32), mad.astype(np.float32), n.astype(np.int16)


def _exon_mask(bin_mask, b_starts, b_ends, e_starts, e_ends):
    out = np.zeros((bin_mask.shape[0], len(e_starts)), dtype=bool)
    for j, (a, b) in enumerate(zip(e_starts, e_ends)):
        sel = (b_ends > a) & (b_starts < b)
        out[:, j] = bin_mask[:, sel].any(axis=1)
    return out


def _recurrent_intronic(paths, cohort, min_samples=3, tol=600):
    ev = []
    for p in paths or []:
        with open(p, encoding='utf8') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                if r.get('exons', '.') in ('.', '') and r.get('partial_exons', '.') in ('.', ''):
                    ev.append((r['type'], int(r['start']), int(r['end']), r['sample_id']))
    clusters = []
    for typ, a, b, sid in sorted(ev, key=lambda t: (t[0], t[1])):
        for c in clusters:
            if c['type'] == typ and abs(c['a'] - a) <= tol and abs(c['b'] - b) <= tol:
                c['s'].add(sid)
                c['as'].append(a)
                c['bs'].append(b)
                break
        else:
            clusters.append(dict(type=typ, a=a, b=b, s={sid}, **{'as': [a], 'bs': [b]}))
    return [(c['type'], int(np.median(c['as'])), int(np.median(c['bs'])), len(c['s']), cohort)
            for c in clusters if len(c['s']) >= min_samples]


def _tiles(a, b):
    st = np.arange(a, b, TILE, dtype=np.int64)
    return st, np.minimum(st + TILE, b)


def _read_list(path):
    """Sample sheet (CSV / TSV / XLSX): sample_id, bam [, sex, batch]; all BAMs must exist and be indexed."""
    from .tables import read_sample_sheet
    rows = read_sample_sheet(path)
    bad = []
    for r in rows:
        b = r['bam']
        if not os.path.isfile(b):
            bad.append(f"{r['sample_id']}: not found {b}")
        elif not any(os.path.isfile(x) for x in (b + '.bai', b + '.csi', b[:-4] + '.bai', b + '.crai')):
            bad.append(f"{r['sample_id']}: no index (.bai/.csi) for {b}")
    if bad:
        raise ValueError(f'{path}: {len(bad)} of {len(rows)} BAMs unusable:\n  ' + '\n  '.join(bad[:20])
                         + ('\n  ...' if len(bad) > 20 else ''))
    return rows


def _bam_depths(rows, mode, st, en, e_st, e_en, threads, log):
    """Depth per bin / exon and index stats for each listed BAM (same code as `dmdseek call`)."""
    from concurrent.futures import ThreadPoolExecutor
    from . import bamio

    w0, w1 = int(min(st[0], e_st.min())), int(max(en[-1], e_en.max()))

    def one(r):
        try:
            return _one(r)
        except Exception as exc:
            raise ValueError(f"{r['sample_id']} ({r['bam']}): {exc}") from exc

    def _one(r):
        bam = bamio.open_bam(r['bam'])
        stats = bamio.index_stats(bam)
        cov, _ = bamio.coverage(bam, w0, w1)
        D = bamio.interval_means(cov, w0, st, en)
        E = bamio.interval_means(cov, w0, e_st, e_en)
        sex = (r.get('sex') or '').strip().upper()[:1]
        if sex not in ('M', 'F'):
            sex = bamio.sex_from_stats(bam, stats, mode)[0]
        bam.close()
        return D, E, stats, sex

    res = []
    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        for k, out in enumerate(ex.map(one, rows), 1):
            res.append(out)
            if k % 10 == 0 or k == len(rows):
                log(f'  depth {k}/{len(rows)}')
    return res


def _cohort_wgs(a, out, exons, log):
    """-> ids, sex, st, en, norm (bin depth / flank median), enorm (exon depth / flank median)."""
    if a.wgs_bams:
        rows = _read_list(a.wgs_bams)
        st, en = _tiles(DMD_START - 1_000_000, DMD_END + 1_000_000)
        res = _bam_depths(rows, 'wgs', st, en, out['exon_start'], out['exon_end'], a.threads, log)
        flank = (en <= DMD_START) | (st >= DMD_END)
        D = np.vstack([x[0] for x in res])
        E = np.vstack([x[1] for x in res])
        fl = np.median(D[:, flank], axis=1)
        return (np.array([r['sample_id'] for r in rows]), np.array([x[3] for x in res]), st, en,
                D / fl[:, None], E / fl[:, None], 'pysam')
    p = np.load(a.wgs_pon)
    ids = p['sample_ids'].astype(str)
    enorm = np.full((len(ids), 79), np.nan)
    if a.wgs_pon_exons:
        e = np.load(a.wgs_pon_exons)
        epos = {s: i for i, s in enumerate(e['sample_ids'].astype(str))}
        order = np.argsort(e['starts'])
        if not np.array_equal(e['starts'][order], out['exon_start']):
            raise ValueError('WGS exon PoN coordinates differ from --exon-bed')
        for i, s in enumerate(ids):
            if s in epos:
                enorm[i] = e['normalized_depth'][epos[s]][order]
    return (ids, p['sex'].astype(str), p['start'].astype(np.int64), p['end'].astype(np.int64),
            p['normalized_depth'].astype(np.float64), enorm, 'PoN (mosdepth)')


def _cohort_tgt(a, out, exons, log):
    """-> ids, sex, batch, flags, st, en, D, E(raw), scale, offtarget(dict)."""
    if a.target_bams:
        rows = _read_list(a.target_bams)
        w = [int(x) for x in a.target_window.replace(',', '').split('-')]
        st, en = _tiles(w[0], w[1])
        res = _bam_depths(rows, 'targeted', st, en, out['exon_start'], out['exon_end'], a.threads, log)
        ids = np.array([r['sample_id'] for r in rows])
        batch = np.array([(r.get('batch') or 'B1').strip() for r in rows])
        D = np.vstack([x[0] for x in res])
        E = np.vstack([x[1] for x in res])
        off = {}
        for sid, (_, _, stats, _) in zip(ids, res):
            off[sid] = sum(stats.values()) - stats.get(CHROM, 0) - stats.get('chrY', 0)
        s0 = np.array([crude_scale(d) for d in D])
        scale = np.full(len(ids), np.nan)
        for g in np.unique(batch):
            m = (batch == g) & np.isfinite(s0)
            with np.errstate(all='ignore'):
                prof = np.nanmedian(D[m] / s0[m, None], axis=0)
            for i in np.flatnonzero(m):
                with np.errstate(all='ignore'):
                    R = D[i] / s0[i] / prof
                ok = np.isfinite(R) & (prof > 0.2) & (R > 0.15)
                scale[i] = s0[i] * mode_level(R[ok])
        sex = np.array([x[3] for x in res])
        # chrY-undetermined: X dosage from on-target scale per off-target reads vs batch males
        for g in np.unique(batch):
            m = batch == g
            cap = np.array([scale[i] / (off[ids[i]] / 1e6) if off[ids[i]] else np.nan for i in range(len(ids))])
            ref_m = np.nanmedian(cap[m & (sex == 'M')]) if (m & (sex == 'M')).any() else np.nan
            for i in np.flatnonzero(m & (sex == 'U')):
                sex[i] = 'F' if cap[i] / ref_m > 1.6 else ('M' if cap[i] / ref_m < 1.3 else 'U')
        flags = np.array(['PASS' if np.isfinite(scale[i]) else 'NO_ONTARGET' for i in range(len(ids))])
        return ids, sex, batch, flags, st, en, D, E, scale, off, 'pysam'
    p = np.load(a.target_pon)
    ids = p['sample_ids'].astype(str)
    E = np.full((len(ids), 79), np.nan)
    if a.target_pon_exons:
        e = np.load(a.target_pon_exons)
        epos = {s: i for i, s in enumerate(e['sample_ids'].astype(str))}
        order = np.argsort(e['starts'])
        if not np.array_equal(e['starts'][order], out['exon_start']):
            raise ValueError('targeted exon PoN coordinates differ from --exon-bed')
        for i, s in enumerate(ids):
            if s in epos:
                E[i] = e['raw_depth'][epos[s]][order]
    off = {}
    if a.target_qc:
        with open(a.target_qc, encoding='utf8') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                if r.get('offtarget_reads') not in (None, '', 'NA'):
                    off[r['sample_id']] = float(r['offtarget_reads'])
    flags = p['prep_flag'].astype(str) if 'prep_flag' in p.files else np.array(['PASS'] * len(ids))
    return (ids, p['sex'].astype(str), p['batch'].astype(str), flags, p['start'].astype(np.int64),
            p['end'].astype(np.int64), p['raw_depth'].astype(np.float64), E,
            p['flank_depth'].astype(np.float64), off, 'PoN (mosdepth)')


def build_reference(a):
    exons = load_exon_bed(a.exon_bed)
    excl = set()
    for p in a.exclude_samples or []:
        excl |= {l.strip() for l in open(p) if l.strip()}
    out = dict(exon_num=np.array([e[0] for e in exons], dtype=np.int16),
               exon_start=np.array([e[1] for e in exons], dtype=np.int64),
               exon_end=np.array([e[2] for e in exons], dtype=np.int64))
    import hashlib
    bed_md5 = hashlib.md5(open(a.exon_bed, 'rb').read()).hexdigest()
    meta = dict(format=REFERENCE_FORMAT, tool_version=__version__, genome='GRCh38', chrom=CHROM,
                exon_bed_md5=bed_md5,
                dmd_start=DMD_START, dmd_end=DMD_END, tile=TILE, min_mapq=MIN_MAPQ,
                created=datetime.date.today().isoformat(), label=a.label, panels={})
    recurrent = []
    log = print

    if a.wgs_pon or a.wgs_bams:
        ids, sex, st, en, norm, enorm, src = _cohort_wgs(a, out, exons, log)
        keep = np.array([s not in excl for s in ids]) & np.isin(sex, ['M', 'F'])
        evm, nev = _event_masks(a.wgs_events, ids, st, en)
        mask = evm | _cn0_mask(norm, sex)
        emask = _exon_mask(mask, st, en, out['exon_start'], out['exon_end'])
        out['wgs_start'], out['wgs_end'] = st, en
        cnt = {}
        for sx in ('M', 'F'):
            sel = keep & (sex == sx)
            out[f'wgs_{sx}_med'], out[f'wgs_{sx}_mad'], out[f'wgs_{sx}_n'] = _stats(norm[sel], mask[sel])
            emed, emad, _ = _stats(enorm[sel], emask[sel])
            out[f'wgs_{sx}_exon_med'], out[f'wgs_{sx}_exon_mad'] = emed, emad
            cnt[sx] = int(sel.sum())
        nU = int((~np.isin(sex, ['M', 'F'])).sum())
        if nU:
            log(f'[WARN] WGS: {nU} sample(s) with undetermined sex left out - add a "sex" column to the list')
        for sx in ('M', 'F'):
            if cnt[sx] < MIN_PANEL_N:
                log(f'[WARN] WGS {sx} panel has {cnt[sx]} samples (< {MIN_PANEL_N}): {sx} WGS samples cannot be called')
        meta['panels']['wgs'] = dict(kind='wgs', n_M=cnt['M'], n_F=cnt['F'], n_bins=len(st), depth_source=src,
                                     normalisation='median depth of bins outside DMD (+/-1 Mb flanks)')
        recurrent += _recurrent_intronic(a.wgs_events, 'wgs')
        log(f'[WGS] samples M={cnt["M"]} F={cnt["F"]} (excluded {int((~keep).sum())}); depth from {src}; '
            f'event masks={nev}; masked {(mask[keep].mean() if keep.any() else 0.0):.2%} of values')

    if a.target_pon or a.target_bams:
        ids, sex, batch, flags, st, en, D, E, scale, off, src = _cohort_tgt(a, out, exons, log)
        keep = (np.array([s not in excl for s in ids]) & np.isfinite(scale) & np.isin(sex, ['M', 'F'])
                & np.array(['NO_ONTARGET' not in f and 'CHECK_X_DOSAGE' not in f for f in flags]))
        norm = D / scale[:, None]
        evm, nev = _event_masks(a.target_events, ids, st, en)
        mask = evm | _cn0_mask(norm, batch)
        emask = _exon_mask(mask, st, en, out['exon_start'], out['exon_end'])
        out['tgt_start'], out['tgt_end'] = st, en
        groups = sorted(set(batch[keep]))
        names = []
        for k, g in enumerate(groups + ['ALL'], 1):
            sel = keep & ((batch == g) if g != 'ALL' else True)
            if sel.sum() < MIN_PANEL_N or (g == 'ALL' and len(groups) == 1):
                continue
            name = f'T{k}' if g != 'ALL' else 'T_pooled'
            names.append(name)
            s0 = np.array([crude_scale(D[i]) for i in np.flatnonzero(sel)])
            with np.errstate(all='ignore'):
                prof = np.nanmedian(D[sel] / s0[:, None], axis=0)
            med, mad, n = _stats(norm[sel], mask[sel])
            out[f'{name}_prof'], out[f'{name}_med'], out[f'{name}_mad'], out[f'{name}_n'] = \
                prof.astype(np.float32), med, mad, n
            emed, emad, _ = _stats(E[sel] / scale[sel][:, None], emask[sel])
            out[f'{name}_exon_med'], out[f'{name}_exon_mad'] = emed, emad
            caps = defaultdict(list)
            for i in np.flatnonzero(sel):
                if off.get(ids[i]):
                    caps[sex[i]].append(scale[i] / (off[ids[i]] / 1e6))
            meta['panels'][name] = dict(
                kind='targeted', n=int(sel.sum()), n_M=int((sel & (sex == 'M')).sum()),
                n_F=int((sel & (sex == 'F')).sum()), depth_source=src,
                cap_M=float(np.median(caps['M'])) if caps['M'] else None,
                cap_F=float(np.median(caps['F'])) if caps['F'] else None,
                normalisation='on-target depth level (mode) after dividing by the panel capture profile')
        log(f'[TARGETED] panels {names} from {int(keep.sum())} samples (excluded {int((~keep).sum())}); '
            f'depth from {src}; event masks={nev}; masked {(mask[keep].mean() if keep.any() else 0.0):.2%} of values')
        recurrent += _recurrent_intronic(a.target_events, 'targeted')

    if recurrent:
        out['recurrent'] = np.array([f'{t}|{x}|{y}|{n}|{c}' for t, x, y, n, c in recurrent])
    out['meta'] = np.array(json.dumps(meta))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).resolve().parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, **out)
    log(f'[REFERENCE] {a.out}  panels: {", ".join(meta["panels"])}  recurrent intronic events: {len(recurrent)}')
    log('[REFERENCE] contains aggregated statistics only (no sample IDs, no per-sample depth)')


# ---------------------------------------------------------------- use
class Reference:
    def __init__(self, path):
        z = np.load(path, allow_pickle=False)
        self.path = str(path)
        self.z = {k: z[k] for k in z.files}
        self.meta = json.loads(str(self.z['meta']))
        if self.meta.get('format') != REFERENCE_FORMAT:
            raise ValueError(f'{path}: reference format {self.meta.get("format")} (tool expects {REFERENCE_FORMAT})')
        self.exons = [dict(num=int(n), start=int(a), end=int(b), name=f'exon_{int(n):02d}')
                      for n, a, b in zip(self.z['exon_num'], self.z['exon_start'], self.z['exon_end'])]
        self.recurrent = []
        for s in self.z.get('recurrent', []):
            t, x, y, n, c = str(s).split('|')
            self.recurrent.append(dict(type=t, start=int(x), end=int(y), n=int(n), cohort=c))

    def panels(self, kind):
        return [k for k, v in self.meta['panels'].items() if v['kind'] == kind]

    def has(self, kind):
        return bool(self.panels(kind))

    # ---- WGS
    def wgs_bins(self):
        return self.z['wgs_start'], self.z['wgs_end']

    def wgs_normalise(self, D, E, sex):
        st, en = self.wgs_bins()
        flank = (en <= DMD_START) | (st >= DMD_END)
        scale = float(np.median(D[flank]))
        if not scale > 0:
            scale = float('nan')
        g = sex if sex in ('M', 'F') else 'F'
        if self.meta['panels']['wgs'].get(f'n_{g}', 0) < MIN_PANEL_N:
            raise ValueError(f'reference has no usable {g} WGS panel (n={self.meta["panels"]["wgs"].get(f"n_{g}")})')
        return dict(panel='wgs', scale=scale, norm=D / scale, exon_norm=E / scale,
                    med=self.z[f'wgs_{g}_med'].astype(float), mad=self.z[f'wgs_{g}_mad'].astype(float),
                    exon_med=self.z.get(f'wgs_{g}_exon_med', np.full(79, np.nan)).astype(float),
                    exon_mad=self.z.get(f'wgs_{g}_exon_mad', np.full(79, np.nan)).astype(float),
                    flags=[], panel_corr=float('nan'))

    # ---- targeted
    def tgt_bins(self):
        return self.z['tgt_start'], self.z['tgt_end']

    def tgt_normalise(self, D, E, sex, offtarget_reads, force_panel=None):
        s0 = crude_scale(D)
        N0 = D / s0
        best, best_r = None, -2.0
        cands = [force_panel] if force_panel else [p for p in self.panels('targeted')]
        for name in cands:
            prof = self.z[f'{name}_prof'].astype(float)
            ok = np.isfinite(prof) & (prof > 0.2) & (N0 > 0.15)
            rr = float(np.corrcoef(np.log(N0[ok]), np.log(prof[ok]))[0, 1]) if ok.sum() > 100 else -1
            # prefer a batch panel over the pooled one unless clearly better
            if name == 'T_pooled':
                rr -= 0.01
            if rr > best_r:
                best, best_r = name, rr
        prof = self.z[f'{best}_prof'].astype(float)
        ok = np.isfinite(prof) & (prof > 0.2)
        with np.errstate(all='ignore'):
            R = N0 / prof
        s2 = mode_level(R[ok & (R > 0.15)])
        scale = s0 * s2
        flags = []
        pm = self.meta['panels'][best]
        x_dos = float('nan')
        if offtarget_reads and pm.get('cap_M'):
            cap = scale / (offtarget_reads / 1e6)
            x_dos = cap / pm['cap_M']
            exp = pm['cap_M'] if sex == 'M' else (pm.get('cap_F') or 2 * pm['cap_M'])
            rel = math.log2(cap / exp) if exp else 0.0
            if abs(rel) > 0.7:
                Rn = R / s2
                second = float(np.mean(Rn[ok] > 1.6) if rel < 0 else np.mean(Rn[ok] < 0.65))
                if second >= 0.10 and abs(rel - math.copysign(1.0, rel)) < 0.35:
                    scale *= 2 ** -math.copysign(1.0, rel)
                    flags.append(f'RESCALED_x{2 ** -math.copysign(1.0, rel):g}(large_event_shifted_level)')
                else:
                    flags.append(f'CHECK_X_DOSAGE(on/off-target={2 ** rel:.2f}x_expected; one X? capture?)')
        return dict(panel=best, panel_corr=best_r, scale=scale, norm=D / scale, exon_norm=E / scale,
                    med=self.z[f'{best}_med'].astype(float), mad=self.z[f'{best}_mad'].astype(float),
                    exon_med=self.z.get(f'{best}_exon_med', np.full(79, np.nan)).astype(float),
                    exon_mad=self.z.get(f'{best}_exon_mad', np.full(79, np.nan)).astype(float),
                    flags=flags, x_dosage=x_dos)
