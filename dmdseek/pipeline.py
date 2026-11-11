"""One BAM -> DMD del/dup calls, tables, plots and HTML report."""
import csv
import gzip
import json
import math
from pathlib import Path

import numpy as np

from . import __version__, annotate, bamio, report
from . import sawfish as sf
from . import ui
from .annotate import (assign_tier, exon_label, exon_overlap, junction_intervals, match_junctions,
                       reading_frame, recip_overlap, sv_label)
from .model import ratio_to_reference, segment


def sample_name(bam_path):
    p = Path(bam_path)
    stem = p.name
    for suf in ('.bam', '.cram'):
        if stem.endswith(suf):
            stem = stem[:-len(suf)]
    return p.parent.name if stem in ('mapped', 'aligned', 'hifi_reads', 'reads') else stem


def _fmt(v):
    if isinstance(v, float):
        return 'NA' if not math.isfinite(v) else f'{v:.4g}'
    return str(v)


MODE_NAME = {'wgs': 'WGS (whole genome)', 'targeted': 'targeted capture'}


def call_sample(bam_path, ref, p, out_dir, sample=None, mode='auto', sex=None, fasta=None,
                sawfish_exe=None, threads=4, cram_ref=None, log=print, forced_mode=False):
    sample = sample or sample_name(bam_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    flags = []
    bam = bamio.open_bam(bam_path, cram_ref)
    stats = bamio.index_stats(bam)
    auto_mode, frac, total, in_dmd = bamio.detect_mode(bam, stats)
    if mode == 'auto' and auto_mode == 'targeted' and frac < 0.01:
        flags.append(f'NO_ONTARGET_READS(targeted assumed: {total:,} mapped reads, DMD fraction {frac:.5f})')
    mode = auto_mode if mode == 'auto' else mode
    if not ref.has(mode):
        raise bamio.InputError(f'reference {ref.path} has no {mode} panel')
    sex_call, sx = bamio.sex_from_stats(bam, stats, mode)
    sex_source = 'chrY + X dosage' if mode == 'wgs' else 'chrY (off-target)'
    if sx.get('flag'):
        flags.append(sx['flag'])

    starts, ends = ref.wgs_bins() if mode == 'wgs' else ref.tgt_bins()
    exons = ref.exons
    e_st = np.array([e['start'] for e in exons])
    e_en = np.array([e['end'] for e in exons])
    w0, w1 = int(min(starts[0], e_st.min())), int(max(ends[-1], e_en.max()))
    cov, nreads = bamio.coverage(bam, w0, w1)
    D = bamio.interval_means(cov, w0, starts, ends)
    E = bamio.interval_means(cov, w0, e_st, e_en)
    bam.close()
    how = ('set by user' if forced_mode else f'auto-detected: {frac:.2%} of mapped reads in DMD')
    if forced_mode and mode != auto_mode:
        flags.append(f'MODE_SET_BY_USER({mode}; reads look {auto_mode})')
    log(ui.item('dna', 'data type', f'{ui.bold(MODE_NAME[mode])}  {ui.dim("· " + how)}'))
    log(ui.item('chart', 'depth', f'{nreads:,} reads in the DMD window'))

    if mode == 'wgs':
        use_sex = sex or (sex_call if sex_call in ('M', 'F') else 'F')
        if not sex and sex_call == 'U':
            flags.append('SEX_UNDETERMINED(treated_as_F)')
        nz = ref.wgs_normalise(D, E, use_sex)
        x_dos = sx['x_dosage']
    else:
        off = sx.get('offtarget_reads')
        use_sex = sex or sex_call
        if use_sex not in ('M', 'F'):
            probe = ref.tgt_normalise(D, E, 'F', off)
            xd = probe.get('x_dosage', float('nan'))
            use_sex = 'F' if (math.isfinite(xd) and xd > 1.6) else ('M' if math.isfinite(xd) and xd < 1.3 else 'F')
            sex_source = f'X dosage {xd:.2f}x male (chrY inconclusive)'
            if not (math.isfinite(xd) and (xd > 1.6 or xd < 1.3)):
                flags.append('SEX_UNDETERMINED(treated_as_F)')
        nz = ref.tgt_normalise(D, E, use_sex, off)
        x_dos = nz.get('x_dosage', float('nan'))
    if sex:
        sex_source = 'user-specified'
    flags += nz['flags']
    log(ui.item('person', 'sex', f"{ui.bold({'M': 'male', 'F': 'female'}.get(use_sex, use_sex))}  {ui.dim('· from ' + sex_source)}"))

    # breakpoints
    svs, sawfish_status = [], 'not run (no --ref FASTA)'
    if fasta and sawfish_exe:
        try:
            svs = sf.run(bam_path, out_dir / f'{sample}.sawfish', fasta, sawfish_exe, threads, cram_ref, log)
            sawfish_status = f'Sawfish ({sf.sawfish_version(sawfish_exe) or "version ?"}): {len(svs)} records in DMD'
            log(ui.item('link', '', f"{len(svs)} SV record{'' if len(svs) == 1 else 's'} in DMD"))
        except Exception as exc:          # depth-only result is still useful
            sawfish_status = f'FAILED: {exc}'
            flags.append('SAWFISH_FAILED(depth_only)')
            log(ui.item('warn', 'breakpoints', ui.yellow(f'Sawfish failed - depth only ({exc})')))
    elif fasta and not sawfish_exe:
        sawfish_status = 'not run (sawfish executable not found)'
    junctions = junction_intervals(svs, p)

    r, z, sd_rel, valid = ratio_to_reference(nz['norm'], nz['med'], nz['mad'], p)
    segs, qc = segment(r, z, sd_rel, valid, starts, ends, use_sex, p, junctions)

    events = []
    for ev in segs:
        best, conc, disc, cplx = match_junctions(ev, svs, p)
        if best:
            kind, sv = best
            ev['junction_match'] = kind
            if kind != 'TANDEM_DUP_AS_INS':
                ev['bp_start'], ev['bp_end'] = sv['start'], sv['end']
                ev['bp_text'] = f"chrX:{sv['start']:,}-{sv['end']:,} ({kind}, GT {sv.get('gt', '.')})"
            else:
                ev['bp_text'] = f"tandem DUP as INS at chrX:{sv['start']:,}"
        else:
            ev['junction_match'] = 'NONE'
            ev['bp_text'] = 'none'
        ev['discordant_junctions'] = ';'.join(sv_label(s) for s in disc) or '.'
        ev['complex_junctions'] = ';'.join(sv_label(s) for s in cplx[:3]) or '.'
        a = ev.get('bp_start', ev['start'])
        b = ev.get('bp_end', ev['end'])
        full, partial = exon_overlap(a, b, exons)
        ev['exon_nums'] = sorted(full + partial)
        ev['exons'] = exon_label(full) if full else '.'
        ev['partial_exons'] = exon_label(partial) if partial else '.'
        ev['exons_long'] = (('exon ' + ev['exons'].replace('ex', '')) if full else '') + \
            ((' + partial ' + ev['partial_exons'].replace('ex', '')) if partial else '') or 'intronic'
        ev['n_exons'] = len(full) + len(partial)
        ev['frame'] = reading_frame(full, partial, exons)
        tier, notes = assign_tier(ev, p)
        for rc in ref.recurrent:
            if rc['type'] == ev['type'] and recip_overlap(ev['start'], ev['end'], rc['start'], rc['end']) >= 0.8:
                notes.append(f"seen_in_{rc['n']}_{rc['cohort']}_reference_samples(common_intronic_variant?)")
        ev['tier'] = tier
        ev['notes'] = ';'.join(notes) or '.'
        events.append(ev)
    events.sort(key=lambda e: ('ABCD'.index(e['tier']), e['start']))

    # exon-level table + review flags (single-exon het events can fall below the bin-level model)
    with np.errstate(all='ignore'):
        er = nz['exon_norm'] / np.maximum(nz['exon_med'], 1e-4)
        ez = (nz['exon_norm'] - nz['exon_med']) / np.maximum(nz['exon_mad'], 0.07 * nz['exon_med'])
    called = {n: e for e in events if e['tier'] in ('A', 'B') for n in e['exon_nums']}
    # sample-specific exon noise (WGS exon ratios are much noisier than capture ones)
    free = np.array([float(er[i]) for i, e in enumerate(exons)
                     if e['num'] not in called and np.isfinite(er[i]) and er[i] > 0.2])
    ex_sd = max(1.4826 * float(np.median(np.abs(free - np.median(free)))), 0.05) if free.size >= 20 else 0.2
    exon_table = []
    for i, e in enumerate(exons):
        rr, zz = float(er[i]), float(ez[i])
        zs = (rr - 1.0) / ex_sd if math.isfinite(rr) else float('nan')
        state = (f"{called[e['num']]['type']} [{called[e['num']]['tier']}]" if e['num'] in called else 'normal')
        flag = ''
        # both: unusual for this exon in the reference AND for this sample's own exon noise
        if (e['num'] not in called and math.isfinite(zs) and math.isfinite(zz)
                and abs(zs) >= p.exon_flag_z and abs(zz) >= p.exon_flag_z and np.sign(zs) == np.sign(zz)):
            if rr < 0.7:
                flag = 'low exon depth: small deletion? probe-site artefact? check IGV / MLPA'
            elif rr > 1.4:
                flag = 'high exon depth: small duplication? check IGV / MLPA'
        exon_table.append(dict(num=e['num'], start=e['start'], end=e['end'], depth=float(E[i]), ratio=rr, z=zz,
                               z_sample=zs, state=state, flag=flag))

    sv_only = []
    for lo, hi, typ in junctions:
        if any(e['type'] == typ and recip_overlap(lo, hi, e.get('bp_start', e['start']), e.get('bp_end', e['end'])) >= 0.5
               for e in events):
            continue
        full, partial = exon_overlap(lo, hi, exons)
        if full or partial:
            m = valid & (ends > lo) & (starts < hi)
            sv_only.append(dict(type=typ, start=lo, end=hi, exons=exon_label(full + partial),
                                depth_ratio=float(np.median(r[m])) if m.any() else float('nan')))

    ab = [e for e in events if e['tier'] in ('A', 'B')]
    final = '; '.join(f"{e['type']} {e['exons'] if e['exons'] != '.' else e['partial_exons']} "
                      f"{e['frame']} [{e['tier']}]" for e in ab) or 'NO_EXONIC_DELDUP'
    if nz['scale'] <= 0 or not math.isfinite(nz['scale']):
        final = 'NO_DEPTH (whole-gene deletion in a male, or capture failure) - check IGV / MLPA'
        flags.append('NO_DEPTH')
    if math.isfinite(qc['sample_sd']) and qc['noise_k'] >= 2.0:
        flags.append(f"NOISY(k={qc['noise_k']:.1f})")

    res = dict(sample=sample, bam=str(bam_path), mode=mode, mode_source=how, mode_fraction=frac, sex=use_sex, sex_source=sex_source,
               x_dosage=x_dos if x_dos is not None else float('nan'), y_fraction=sx['y_fraction'],
               panel=nz['panel'], panel_corr=nz['panel_corr'], scale=nz['scale'],
               n_valid=qc['n_valid'], n_bins=qc['n_bins'], noise_k=qc['noise_k'], sample_sd=qc['sample_sd'],
               sawfish_status=sawfish_status, flags=flags, final_call=final, events=events, exon_table=exon_table,
               sv_only=sv_only, junctions=junctions, exons=exons, bins_start=starts, bins_end=ends, ratio=r,
               valid=valid, reference_label=ref.meta.get('label', '.'), reference_created=ref.meta.get('created', '.'))
    write_outputs(res, out_dir)
    return res


EVENT_COLS = ['sample', 'tier', 'type', 'exons', 'partial_exons', 'frame', 'chrom', 'start_1based', 'end',
              'size_bp', 'median_ratio', 'est_CN', 'hmm_CN', 'baseline_CN', 'cn_fit', 'n_bins', 'llr_score',
              'median_z', 'effect_sd', 'null_ratio', 'junction_match', 'bp_start', 'bp_end',
              'discordant_junctions', 'complex_junctions', 'notes']


def write_outputs(res, out_dir):
    s = res['sample']
    with open(out_dir / f'{s}.dmdseek.events.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(EVENT_COLS)
        for e in res['events']:
            w.writerow([_fmt(x) for x in [s, e['tier'], e['type'], e['exons'], e['partial_exons'], e['frame'], 'chrX',
                                          e['start'] + 1, e['end'], e['size_bp'], e['median_ratio'], e['est_CN'],
                                          e['hmm_CN'], e['baseline_CN'], e['cn_fit'], e['n_bins'], e['llr_score'],
                                          e['median_z'], e['effect_sd'], e['null_ratio'], e['junction_match'],
                                          e.get('bp_start', '.'), e.get('bp_end', '.'), e['discordant_junctions'],
                                          e['complex_junctions'], e['notes']]])
    with open(out_dir / f'{s}.dmdseek.exons.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['sample', 'exon', 'start_1based', 'end', 'mean_depth', 'ratio', 'z_vs_reference',
                    'z_vs_sample_noise', 'state', 'review_flag'])
        for e in res['exon_table']:
            w.writerow([_fmt(x) for x in [s, e['num'], e['start'] + 1, e['end'], e['depth'], e['ratio'], e['z'],
                                          e['z_sample'], e['state'], e['flag'] or '.']])
    with gzip.open(out_dir / f'{s}.dmdseek.bins.tsv.gz', 'wt') as f:
        f.write('chrom\tstart_1based\tend\tratio\tusable\n')
        for a, b, rr, v in zip(res['bins_start'], res['bins_end'], res['ratio'], res['valid']):
            f.write(f'chrX\t{a + 1}\t{b}\t{_fmt(float(rr))}\t{int(v)}\n')
    summ = {k: res[k] for k in ('sample', 'bam', 'mode', 'mode_source', 'mode_fraction', 'sex', 'sex_source', 'x_dosage',
                                'y_fraction', 'panel', 'panel_corr', 'scale', 'n_valid', 'n_bins', 'noise_k',
                                'sample_sd', 'sawfish_status', 'flags', 'final_call', 'reference_label',
                                'reference_created')}
    summ['tool_version'] = __version__
    summ['events_AB'] = [{k: e[k] for k in ('tier', 'type', 'exons', 'partial_exons', 'frame', 'start', 'end',
                                            'median_ratio', 'junction_match')} for e in res['events']
                         if e['tier'] in ('A', 'B')]
    summ['exon_review_flags'] = [e['num'] for e in res['exon_table'] if e['flag']]
    with open(out_dir / f'{s}.dmdseek.summary.json', 'w') as f:
        json.dump(summ, f, indent=1, default=lambda o: None if isinstance(o, float) and not math.isfinite(o) else str(o))
    gene_png = out_dir / f'{s}.dmdseek.gene.png'
    exon_png = out_dir / f'{s}.dmdseek.exons.png'
    report.gene_plot(gene_png, res)
    report.exon_plot(exon_png, res)
    report.html_report(out_dir / f'{s}.dmdseek.report.html', res, gene_png, exon_png)
