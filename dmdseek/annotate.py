"""Exon / reading-frame annotation, junction matching and evidence tiers."""
from .constants import CHROM, FIRST_LAST_EXONS


def recip_overlap(a0, a1, b0, b1):
    ov = min(a1, b1) - max(a0, b0)
    if ov <= 0:
        return 0.0
    return min(ov / max(a1 - a0, 1), ov / max(b1 - b0, 1))


def exon_overlap(st, en, exons):
    full, partial = [], []
    for e in exons:
        if e['start'] >= st and e['end'] <= en:
            full.append(e['num'])
        elif e['start'] < en and e['end'] > st:
            partial.append(e['num'])
    return full, partial


def exon_label(nums):
    if not nums:
        return '.'
    nums = sorted(nums)
    runs, s, q = [], nums[0], nums[0]
    for n in nums[1:]:
        if n == q + 1:
            q = n
            continue
        runs.append((s, q))
        s = q = n
    runs.append((s, q))
    return ','.join(f'ex{a}' if a == b else f'ex{a}-{b}' for a, b in runs)


def reading_frame(full, partial, exons):
    if not full and not partial:
        return 'NA_INTRONIC'
    if partial:
        return 'UNDETERMINED_PARTIAL_EXON'
    if len(full) == len(exons):
        return 'WHOLE_GENE'
    if set(full) & FIRST_LAST_EXONS:
        return 'NA_FIRST_OR_LAST_EXON'
    length = sum(e['end'] - e['start'] for e in exons if e['num'] in set(full))
    return 'IN_FRAME' if length % 3 == 0 else 'OUT_OF_FRAME'


def bnd_pair_kind(sv):
    """Intrachromosomal BND -> (DEL|DUP|INV, lo, hi) from the VCF ALT orientation."""
    alt = sv.get('alt', '') or ''
    mates = [q for c, q in sv.get('mates', []) if c == CHROM]
    if not mates or ('[' not in alt and ']' not in alt):
        return None
    mpos, pos = mates[0], sv['start']
    lo, hi, lower = min(pos, mpos), max(pos, mpos), pos < mpos
    if alt.endswith('['):
        kind = 'DEL' if lower else 'DUP'
    elif alt.startswith(']'):
        kind = 'DUP' if lower else 'DEL'
    else:
        kind = 'INV'
    return kind, lo, hi


def junction_intervals(svs, p):
    out, seen = [], set()
    for sv in svs:
        if sv['filt'] not in ('PASS', '.'):
            continue
        if sv['type'] in ('DEL', 'DUP'):
            lo, hi, typ = sv['start'], sv['end'], sv['type']
        elif sv['type'] == 'BND':
            pair = bnd_pair_kind(sv)
            if not pair or pair[0] not in ('DEL', 'DUP'):
                continue
            typ, lo, hi = pair
        else:
            continue
        key = (typ, lo // 50, hi // 50)
        if hi - lo >= p.min_sv_bp and key not in seen:
            seen.add(key)
            out.append((lo, hi, typ))
    return out


def sv_label(sv):
    return f"{sv['type']}:{sv['start']}-{sv['end']}:{sv['id']}"


def match_junctions(ev, svs, p):
    tol = max(p.bp_tol, 0.1 * ev['size_bp'])
    concordant, discordant, complex_hits = [], [], []
    for sv in svs:
        if sv['filt'] not in ('PASS', '.'):
            continue
        typ = sv['type']
        if typ in ('DEL', 'DUP'):
            if sv['end'] - sv['start'] < p.min_sv_bp:
                continue
            ro = recip_overlap(ev['start'], ev['end'], sv['start'], sv['end'])
            ends_ok = abs(sv['start'] - ev['start']) <= tol and abs(sv['end'] - ev['end']) <= tol
            if typ == ev['type'] and (ro >= p.min_ro or ends_ok):
                concordant.append(('BREAKPOINT' if ends_ok else 'OVERLAP', sv))
            elif typ != ev['type'] and ro >= p.min_ro:
                discordant.append(sv)
        elif typ == 'INS' and ev['type'] == 'DUP' and sv.get('ins_len', 0) >= p.min_sv_bp:
            near = ev['start'] - tol <= sv['start'] <= ev['end'] + tol
            if near and 0.5 <= sv['ins_len'] / max(ev['size_bp'], 1) <= 2.0:
                concordant.append(('TANDEM_DUP_AS_INS', sv))
        elif typ in ('BND', 'INV'):
            pair = bnd_pair_kind(sv) if typ == 'BND' else None
            if pair and pair[0] in ('DEL', 'DUP') and pair[2] - pair[1] >= p.min_sv_bp:
                kind, lo, hi = pair
                if kind == ev['type'] and abs(lo - ev['start']) <= tol and abs(hi - ev['end']) <= tol:
                    concordant.append(('BREAKPOINT_BND', dict(sv, start=lo, end=hi)))
                    continue
                if kind != ev['type'] and recip_overlap(ev['start'], ev['end'], lo, hi) >= p.min_ro:
                    discordant.append(dict(sv, type=f'BND:{kind}', start=lo, end=hi))
                    continue
            pts = [sv['start'], sv['end']] + [q for c, q in sv.get('mates', []) if c == CHROM]
            if any(abs(q - ev['start']) <= tol or abs(q - ev['end']) <= tol for q in pts):
                complex_hits.append(sv)
    used = {c[1]['id'].rsplit(':', 1)[0] for c in concordant if c[0] == 'BREAKPOINT_BND'}
    complex_hits = [c for c in complex_hits if c['id'].rsplit(':', 1)[0] not in used]
    rank = {'BREAKPOINT': 0, 'BREAKPOINT_BND': 0, 'OVERLAP': 1, 'TANDEM_DUP_AS_INS': 2}
    concordant.sort(key=lambda t: rank[t[0]])
    return (concordant[0] if concordant else None), concordant, discordant, complex_hits


def assign_tier(ev, p):
    exonic = ev['n_exons'] > 0
    junction = ev['junction_match'] != 'NONE'
    notes = []
    if ev.get('junction_guided'):
        notes.append('bounds_from_junction')
    if ev['beyond_dmd'] == 'outside_DMD':
        notes.append('outside_DMD')
    elif ev['beyond_dmd']:
        notes.append(f"extends_beyond_DMD({ev['beyond_dmd']})")
    if ev['cn_fit'] == 'INTERMEDIATE':
        notes.append('intermediate_CN(mosaic/contamination/sex?)')
    if ev['discordant_junctions'] != '.':
        notes.append('opposite_type_junction')
    if ev['complex_junctions'] != '.':
        notes.append('BND/INV_near_boundary(complex?)')
    if ev['type'] == 'DUP' and not junction:
        notes.append('tandem_status_unknown')
    # depth-only evidence: large effect or stands out from the sample's own same-length windows; an
    # intermediate copy number must also stand out (otherwise it looks like the sample's noise)
    distinct = ev['null_ratio'] >= p.min_null_ratio
    strong = ev['n_bins'] >= p.min_bins_depth_only and (ev['effect_sd'] >= p.min_effect_sd or distinct) \
        and (ev['cn_fit'] == 'CLEAN' or distinct)
    if exonic and junction and ev['cn_fit'] == 'CLEAN' and ev['discordant_junctions'] == '.':
        tier = 'A'
    elif exonic and (junction or strong or ev['complex_junctions'] != '.'):
        tier = 'B'
    elif exonic:
        tier = 'D'
        notes.append('weak_depth_only')
    elif junction:
        tier = 'C'
    else:
        tier = 'D'
    return tier, notes


TIER_TEXT = {'A': 'A (depth + breakpoint)', 'B': 'B (review)', 'C': 'C (intronic)', 'D': 'D (low)'}
