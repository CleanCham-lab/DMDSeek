"""Copy-number segmentation of one sample against a stored reference (sex-aware HMM)."""
import math
from collections import Counter
from dataclasses import dataclass

import numpy as np

from .constants import DMD_START, DMD_END


@dataclass
class Params:
    p_switch: float = 1e-6          # HMM state change probability per bin
    outlier_eps: float = 1e-3       # uniform outlier mixture weight
    max_rel_mad: float = 0.35       # bins with unstable reference (MAD/median) are masked
    min_bins: int = 5               # min valid 200-bp bins per event (1 kb)
    min_abs_z: float = 3.0
    min_llr: float = 15.0
    fit_tol: float = 0.35           # CLEAN if |ratio - expected| <= fit_tol * |expected - 1|
    min_null_ratio: float = 2.5     # long low-amplitude segment must beat same-length windows by this
    min_effect_sd: float = 4.0      # depth-only exonic event -> tier B if effect >= this (or null_ratio)
    min_bins_depth_only: int = 10
    min_sv_bp: int = 500
    min_ro: float = 0.5
    bp_tol: int = 2000
    exon_flag_z: float = 3.0        # exon review flag (not a call): |z| vs reference AND vs own exon noise


def robust_mad(a, axis=0):
    med = np.nanmedian(a, axis=axis)
    return 1.4826 * np.nanmedian(np.abs(a - np.expand_dims(med, axis=axis)), axis=axis)


def ratio_to_reference(norm, ref_med, ref_mad, p):
    """Depth ratio, z and relative reference spread for one normalised sample."""
    with np.errstate(all='ignore'):
        r = norm / np.maximum(ref_med, 1e-4)
        z = (norm - ref_med) / np.maximum(ref_mad, 0.07 * ref_med)
        rel_mad = ref_mad / np.maximum(ref_med, 1e-3)
    valid = (np.isfinite(norm) & np.isfinite(ref_med) & np.isfinite(ref_mad)
             & (ref_med > 0.15) & (rel_mad < p.max_rel_mad))
    return r, z, np.maximum(np.nan_to_num(rel_mad, nan=1.0), 0.03), valid


def noise_inflation(r, sd_rel, valid):
    """k: how much noisier this sample is than the reference spread; sd: robust SD of ratio."""
    m = valid & (np.abs(r - 1) < 0.35)
    if m.sum() < 200:
        return 1.0, float('nan')
    resid = (r[m] - 1) / sd_rel[m]
    k = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    sd = 1.4826 * np.median(np.abs(r[m] - np.median(r[m])))
    return float(min(max(k, 1.0), 6.0)), float(sd)


def state_means(sex):
    """Expected depth ratio per state; males get half steps (mosaic / mis-sexed female carrier)."""
    if sex == 'M':
        return np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0]), 1
    return np.array([0.0, 0.5, 1.0, 1.5, 2.0]), 2


def emissions(r, sd_rel, k, mu, eps):
    sd = np.sqrt((k * sd_rel[:, None]) ** 2 * np.maximum(mu[None, :], 0.15) + 0.03 ** 2)
    ll = -0.5 * ((r[:, None] - mu[None, :]) / sd) ** 2 - np.log(sd) - 0.5 * math.log(2 * math.pi)
    return np.logaddexp(math.log(1 - eps) + ll, math.log(eps) - math.log(4.0))


def viterbi(E, base_state, p_switch):
    n, S = E.shape
    lt = np.full((S, S), math.log(p_switch / (S - 1)))
    np.fill_diagonal(lt, math.log(1 - p_switch))
    V = np.full(S, math.log(p_switch / (S - 1)))
    V[base_state] = math.log(1 - p_switch)
    V = V + E[0]
    bp = np.empty((n, S), dtype=np.int8)
    for t in range(1, n):
        M = V[:, None] + lt
        bp[t] = M.argmax(axis=0)
        V = M.max(axis=0) + E[t]
    path = np.empty(n, dtype=np.int8)
    path[-1] = int(V.argmax())
    for t in range(n - 1, 0, -1):
        path[t - 1] = bp[t, path[t]]
    return path


def segment(r, z, sd_rel, valid, starts, ends, sex, p, junctions=()):
    """HMM segmentation + gap merge + junction-guided refinement.
    junctions: [(lo, hi, 'DEL'|'DUP')] from the breakpoint caller."""
    mu, base = state_means(sex)
    base_state = int(np.flatnonzero(mu == 1.0)[0])
    k, sample_sd = noise_inflation(r, sd_rel, valid)
    vidx = np.flatnonzero(valid)
    qc = dict(n_valid=int(valid.sum()), n_bins=len(valid), noise_k=k, sample_sd=sample_sd,
              gap_merges=0, junction_guided=0)
    if len(vidx) < 100:
        return [], qc
    E = emissions(r[vidx], sd_rel[vidx], k, mu, p.outlier_eps)
    rv = r[vidx]
    csum = np.concatenate([[0.0], np.cumsum(rv - 1.0)])

    def null_ratio(a, b):
        """|mean(r-1)| of the event vs the 99.5th pct of same-length windows elsewhere."""
        n = b - a
        sw = np.arange(0, len(rv) - n + 1)
        keep = (sw + n <= a) | (sw >= b)
        if keep.sum() < 20:
            return float('nan')
        means = np.abs(csum[sw[keep] + n] - csum[sw[keep]]) / n
        return abs(csum[b] - csum[a]) / n / max(float(np.percentile(means, 99.5)), 1e-3)

    path = viterbi(E, base_state, p.p_switch)
    direction = np.sign(mu[path] - 1.0).astype(int)
    segs, j = [], 0
    while j < len(path):
        if direction[j] == 0:
            j += 1
            continue
        d, jj = direction[j], j + 1
        while jj < len(path) and direction[jj] == d:
            jj += 1
        segs.append((j, jj, int(d)))
        j = jj

    merged = []
    for a, b, d in segs:
        if merged and merged[-1][2] == d and a - merged[-1][1] <= max(25, 0.01 * (b - merged[-1][0])):
            merged[-1] = (merged[-1][0], b, d)
            continue
        merged.append((a, b, d))
    qc['gap_merges'] = len(segs) - len(merged)
    segs = merged

    vst, ven = starts[vidx], ends[vidx]
    jg = set()
    for lo, hi, typ in sorted(junctions, key=lambda x: x[0] - x[1]):     # largest first
        d = -1 if typ == 'DEL' else 1
        a = int(np.searchsorted(ven, lo, side='right'))
        b = int(np.searchsorted(vst, hi, side='left'))
        if b - a < 3 or d * (float(np.median(r[vidx[a:b]])) - 1.0) < 0.25:
            continue
        if any(x[0] < b and x[1] > a for x in jg):
            continue
        segs = [x for x in segs if not (x[2] == d and x[0] < b and x[1] > a)]
        segs.append((a, b, d))
        jg.add((a, b))
    segs.sort()
    qc['junction_guided'] = len(jg)

    out = []
    for a, b, d in segs:
        bins = vidx[a:b]
        guided = (a, b) in jg
        if guided:
            cand = np.flatnonzero(np.sign(mu - 1.0) == d)
            state = int(cand[np.argmin(np.abs(mu[cand] - float(np.median(r[bins]))))])
        else:
            state = Counter(path[a:b].tolist()).most_common(1)[0][0]
        med_r = float(np.median(r[bins]))
        mu_s = float(mu[state])
        llr = float(np.sum(E[a:b, state] - E[a:b, base_state]))
        cn_state = mu_s * base
        integer_state = abs(cn_state - round(cn_state)) < 1e-6
        effect_sd = abs(med_r - 1.0) / float(k * np.median(sd_rel[bins]))
        n_valid = b - a
        med_z = float(np.nanmedian(z[bins]))
        st, en = int(starts[bins[0]]), int(ends[bins[-1]])
        clean = integer_state and abs(med_r - mu_s) <= p.fit_tol * abs(mu_s - 1)
        nr = null_ratio(a, b)
        keep = (n_valid >= p.min_bins and llr >= p.min_llr
                and (abs(med_z) >= p.min_abs_z or (np.isfinite(nr) and nr >= p.min_null_ratio)))
        if guided:
            keep = n_valid >= 3 and llr > 0
        if not keep:
            continue
        beyond = ('outside_DMD' if (en <= DMD_START or st >= DMD_END) else
                  '+'.join(x for x, c in (('3prime_flank', st < DMD_START), ('5prime_flank', en > DMD_END)) if c))
        out.append(dict(
            type='DEL' if d < 0 else 'DUP', start=st, end=en, size_bp=en - st,
            n_bins=n_valid, n_masked_bins=int(bins[-1] - bins[0] + 1) - n_valid,
            median_ratio=med_r, est_CN=max(0.0, base * med_r),
            hmm_CN=(int(round(cn_state)) if integer_state else cn_state), hmm_mu=mu_s, baseline_CN=base,
            cn_fit='CLEAN' if clean else 'INTERMEDIATE', effect_sd=effect_sd, null_ratio=nr,
            llr_score=llr, median_z=med_z, junction_guided=guided, beyond_dmd=beyond))
    return out, qc
