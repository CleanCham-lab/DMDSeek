"""Plots (gene / zoom, MLPA-style exon chart) and a self-contained HTML report."""
import base64
import html
import io

import numpy as np

from . import __version__
from .branding import AFFILIATION, DEPARTMENT, INSTITUTION, LAB, PEOPLE
from .constants import DMD_END, DMD_START

COL_DEL, COL_DUP = '#d03b3b', '#2a78d6'
COL_POINT, COL_INK, COL_INK2, COL_GRID, COL_SURF = '#9a9893', '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _style(ax):
    ax.set_facecolor(COL_SURF)
    ax.tick_params(labelsize=7, colors=COL_INK2)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(COL_GRID)


def gene_plot(path, res):
    plt = _plt()
    st, en, r, valid = res['bins_start'], res['bins_end'], res['ratio'], res['valid']
    sex, exons, events = res['sex'], res['exons'], res['events']
    base = 1 if sex == 'M' else 2
    mid = (st + en) / 2 / 1e6
    shown = [e for e in events if e['tier'] in ('A', 'B')][:2] or [e for e in events if e['tier'] == 'C'][:1]
    fig, axes = plt.subplots(1 + len(shown), 1, figsize=(12, 3.0 * (1 + len(shown))), facecolor=COL_SURF,
                             squeeze=False)
    axes = axes[:, 0]

    def panel(ax, x0, x1, title):
        _style(ax)
        m = valid & (mid >= x0) & (mid <= x1)
        ax.scatter(mid[m], np.clip(r[m], 0, 2.6), s=2, color=COL_POINT, alpha=0.55, linewidths=0, rasterized=True)
        for c in range(0, 5):
            y = c / base
            if y <= 2.6:
                ax.axhline(y, color=COL_INK2 if c == base else COL_GRID, lw=0.8 if c == base else 0.6, zorder=0)
                ax.text(x1, y, f' CN{c}', va='center', ha='left', fontsize=7, color=COL_INK2)
        for e in events:
            a, b = e['start'] / 1e6, e['end'] / 1e6
            if e['tier'] in ('A', 'B', 'C') and b >= x0 and a <= x1:
                col = COL_DEL if e['type'] == 'DEL' else COL_DUP
                ax.axvspan(a, b, color=col, alpha=0.12, lw=0)
                ax.hlines(e['median_ratio'], a, b, color=col, lw=2.2)
        for lo, hi, typ in res['junctions']:
            for q in (lo, hi):
                if x0 <= q / 1e6 <= x1:
                    ax.axvline(q / 1e6, color=COL_DEL if typ == 'DEL' else COL_DUP, lw=0.8, ls=(0, (3, 2)))
        last = -1e9
        for e in sorted(exons, key=lambda e: e['start']):
            xm = (e['start'] + e['end']) / 2 / 1e6
            if x0 <= xm <= x1:
                ax.plot([xm, xm], [-0.17, -0.07], color=COL_INK2, lw=0.7)
                want = e['num'] % 5 == 0 or e['num'] in (1, 79) or (x1 - x0) < 0.5
                if want and xm - last >= 0.012 * (x1 - x0):
                    ax.text(xm, -0.24, str(e['num']), ha='center', va='top', fontsize=6, color=COL_INK2)
                    last = xm
        ax.set_xlim(x0, x1)
        ax.set_ylim(-0.42, 2.7)
        ax.set_ylabel('depth ratio\n(sample / reference)', fontsize=8, color=COL_INK2)
        ax.set_title(title, fontsize=9, color=COL_INK, loc='left')

    calls = '; '.join(f"{e['type']} {e['exons']} [{e['tier']}]" for e in events if e['tier'] in ('A', 'B')) \
        or 'no exonic DEL/DUP'
    panel(axes[0], DMD_START / 1e6 - 0.05, DMD_END / 1e6 + 0.05,
          f"{res['sample']} ({sex}, {res['mode']})  |  {calls}  — DMD on minus strand: exon 1 at right; "
          f"dashed = breakpoints")
    for k, e in enumerate(shown):
        pad = max(e['size_bp'] * 0.6, 20_000) / 1e6
        panel(axes[1 + k], e['start'] / 1e6 - pad, e['end'] / 1e6 + pad,
              f"zoom: {e['type']} {e['exons']}  chrX:{e['start'] + 1:,}-{e['end']:,}  ratio {e['median_ratio']:.2f}"
              f"  breakpoint: {e['junction_match']}  frame: {e['frame']}")
    axes[-1].set_xlabel('GRCh38 chrX (Mb)', fontsize=8, color=COL_INK2)
    fig.tight_layout()
    fig.savefig(path, dpi=120, facecolor=COL_SURF)
    plt.close(fig)


def exon_plot(path, res):
    """MLPA-like view: exon 1-79 ratio with expected copy-number levels and +/-30% cut-offs."""
    plt = _plt()
    ex = res['exon_table']
    x = np.array([e['num'] for e in ex])
    y = np.array([e['ratio'] if np.isfinite(e['ratio']) else np.nan for e in ex])
    sex = res['sex']
    fig, ax = plt.subplots(figsize=(12, 3.2), facecolor=COL_SURF)
    _style(ax)
    levels = (0, 1, 2) if sex == 'M' else (0, 0.5, 1, 1.5, 2)
    for lv in levels:
        ax.axhline(lv, color=COL_INK2 if lv == 1 else COL_GRID, lw=0.8 if lv == 1 else 0.6, zorder=0)
    for cut in (0.7, 1.3):
        ax.axhline(cut, color=COL_INK2, lw=0.6, ls=(0, (2, 3)), zorder=0)
    for e in res['events']:
        if e['tier'] in ('A', 'B') and e['exon_nums']:
            col = COL_DEL if e['type'] == 'DEL' else COL_DUP
            ax.axvspan(min(e['exon_nums']) - 0.5, max(e['exon_nums']) + 0.5, color=col, alpha=0.10, lw=0)
    col = np.where(y < 0.7, COL_DEL, np.where(y > 1.3, COL_DUP, COL_POINT))
    ax.scatter(x, np.clip(y, 0, 2.4), s=22, c=col, zorder=3, linewidths=0)
    ax.set_xlim(0.3, 79.7)
    ax.set_ylim(-0.08, 2.5)
    ax.set_xticks(range(1, 80, 2))
    ax.set_xticklabels([str(i) for i in range(1, 80, 2)], fontsize=6)
    ax.set_xlabel('DMD exon (MANE NM_004006.3)', fontsize=8, color=COL_INK2)
    ax.set_ylabel('exon depth ratio', fontsize=8, color=COL_INK2)
    exp = 'male: 0 = deleted, 2 = duplicated' if sex == 'M' else 'female: 0.5 = het deletion, 1.5 = het duplication'
    ax.set_title(f"{res['sample']} — exon-level ratio ({exp}; dotted = 0.7 / 1.3)", fontsize=9, color=COL_INK,
                 loc='left')
    fig.tight_layout()
    fig.savefig(path, dpi=120, facecolor=COL_SURF)
    plt.close(fig)


def _img(path):
    with open(path, 'rb') as f:
        return 'data:image/png;base64,' + base64.b64encode(f.read()).decode()


def _e(x):
    return html.escape(str(x))


CSS = """
:root{--ink:#0b0b0b;--ink2:#52514e;--line:#e4e3df;--surf:#fcfcfb;--card:#ffffff;--del:#d03b3b;--dup:#2a78d6;
--a:#0ca30c;--b:#c98500;--c:#6250d6}
*{box-sizing:border-box}body{margin:0;background:var(--surf);color:var(--ink);
font:14px/1.5 -apple-system,"Segoe UI","Noto Sans KR","Malgun Gothic",Roboto,sans-serif}
main{max-width:1180px;margin:0 auto;padding:24px 16px 48px}
h1{font-size:20px;margin:0 0 4px}h2{font-size:15px;margin:28px 0 8px}.sub{color:var(--ink2);font-size:12px}
.result{margin:16px 0;padding:14px 16px;border:1px solid var(--line);border-left:6px solid var(--ink2);
background:var(--card);border-radius:6px}.result.pos{border-left-color:var(--del)}.result.dup{border-left-color:var(--dup)}
.result .big{font-size:17px;font-weight:600}
table{border-collapse:collapse;width:100%;background:var(--card);font-size:12.5px}
th,td{border-bottom:1px solid var(--line);padding:5px 8px;text-align:left;vertical-align:top}
th{color:var(--ink2);font-weight:600;background:#f6f5f2}
.kv td{overflow-wrap:anywhere}.kv td:first-child{color:var(--ink2);width:30%;min-width:110px}
.result .sub,.result .big{overflow-wrap:anywhere}
.tier{display:inline-block;min-width:22px;text-align:center;border-radius:4px;color:#fff;font-weight:600;padding:0 5px}
.tA{background:var(--a)}.tB{background:var(--b)}.tC{background:var(--c)}.tD{background:#9a9893}
.flag{color:#a14d00}img{max-width:100%;border:1px solid var(--line);border-radius:4px;background:#fff}
details{margin-top:8px}summary{cursor:pointer;color:var(--ink2)}
.scroll{overflow-x:auto}footer{margin-top:36px;color:var(--ink2);font-size:11.5px}
.brand{display:flex;align-items:center;gap:16px;flex-wrap:wrap;padding:14px 0 12px;border-bottom:2px solid var(--ink);
margin-bottom:18px}.brand img{height:52px;width:auto;border:0;background:none}
.brand .org{flex:1 1 280px;min-width:0}.brand .inst{font-size:17px;font-weight:700;letter-spacing:.2px}
.brand .dept{font-size:13px;color:var(--ink2)}.brand .tool{text-align:right;font-size:12px;color:var(--ink2)}
.brand .tool b{display:block;font-size:20px;color:var(--ink);letter-spacing:.3px}
.credit{font-size:12px;color:var(--ink2);margin:-8px 0 18px;line-height:1.6}.credit b{color:var(--ink);font-weight:600}
.guide{margin:12px 0 4px;padding:12px 16px;background:var(--card);border:1px solid var(--line);border-radius:6px}
.guide h3{font-size:13.5px;margin:0 0 8px}.guide dl{display:grid;grid-template-columns:auto 1fr;gap:6px 12px;margin:0}
.guide dt{white-space:nowrap;font-size:13px}.guide dd{margin:0;font-size:12.5px;overflow-wrap:anywhere}
.guide .more{margin-top:10px;font-size:12px;color:var(--ink2)}
@media (max-width:560px){.brand .tool{text-align:left}.guide dl{grid-template-columns:1fr}.guide dt{margin-top:6px}}
"""


def _logo():
    from .branding import logo_file
    p = logo_file()
    if not p:
        return ''
    mime = {'svg': 'image/svg+xml', 'png': 'image/png', 'jpg': 'image/jpeg', 'jpeg': 'image/jpeg',
            'webp': 'image/webp'}[p.suffix[1:].lower()]
    with open(p, 'rb') as f:
        return f"<img src='data:{mime};base64,{base64.b64encode(f.read()).decode()}' alt='{_e(INSTITUTION)}'>"


def _brand():
    people = ' · '.join(f'{_e(role)}: <b>{_e(name)}</b>' for role, name in PEOPLE)
    return (f"<header class='brand'>{_logo()}<div class='org'><div class='inst'>{_e(INSTITUTION)}</div>"
            f"<div class='dept'>{_e(DEPARTMENT)} · {_e(LAB)}</div></div>"
            f"<div class='tool'><b>DMDSeek</b>v{__version__} · DMD exon DEL/DUP · PacBio HiFi</div></header>"
            f"<div class='credit'>{people}<br>{_e(AFFILIATION)}</div>")


GUIDE = """<section class='guide'><h3>How to read this report</h3><dl>
<dt><span class='tier tA'>A</span> confirmed</dt><dd>Exon(s) lost or gained by depth <b>and</b> an exact breakpoint
junction (split reads) that matches the depth change. Highest confidence.</dd>
<dt><span class='tier tB'>B</span> review</dt><dd>Exonic depth change without a matching junction (or with a complex
junction). Likely real when large; check IGV and confirm by MLPA / qPCR.</dd>
<dt><span class='tier tC'>C</span> intronic</dt><dd>Breakpoint-supported deletion / duplication that does not include an
exon. Exons are intact; usually not the cause of disease.</dd>
<dt><span class='tier tD'>D</span> weak</dt><dd>Small or noisy depth signal. Not reported here (events TSV only).</dd>
<dt>Exon review flag</dt><dd>A single exon with unusual depth but no call - e.g. a small heterozygous deletion in a
woman, or a sequence variant under an MLPA probe. Check before reporting.</dd>
<dt>Ratio</dt><dd>Sample depth / reference depth (1.0 = normal). Male: deletion ~0, duplication ~2.
Female: heterozygous deletion ~0.5, duplication ~1.5.</dd>
<dt>Reading frame</dt><dd>Exon-length rule: out-of-frame usually Duchenne, in-frame usually Becker (exceptions exist).
Events touching exon 1 or 79 are not classified.</dd>
</dl><div class='more'>Research use only. Confirm reportable findings by an orthogonal method (MLPA, qPCR) and
review the reads in IGV.</div></section>"""


def html_report(path, res, gene_png, exon_png):
    ab = [e for e in res['events'] if e['tier'] in ('A', 'B')]
    if ab:
        cls = 'pos' if ab[0]['type'] == 'DEL' else 'dup'
        lines = []
        for e in ab:
            zyg = ('hemizygous' if res['sex'] == 'M' else
                   ('heterozygous' if e['hmm_CN'] in (1, 3) else f"CN {e['est_CN']:.1f}"))
            lines.append(f"<div class='big'>{_e(e['type'])} {_e(e['exons_long'])} — {_e(zyg)}, "
                         f"{_e(e['frame'].replace('_', ' ').lower())} "
                         f"<span class='tier t{e['tier']}'>{e['tier']}</span></div>"
                         f"<div class='sub'>chrX:{e['start'] + 1:,}-{e['end']:,} ({e['size_bp']:,} bp, depth) · "
                         f"breakpoint: {_e(e['bp_text'])} · ratio {e['median_ratio']:.2f}</div>")
        result = f"<div class='result {cls}'>{''.join(lines)}</div>"
    else:
        result = ("<div class='result'><div class='big'>No exonic DEL/DUP detected in DMD</div>"
                  "<div class='sub'>Depth and breakpoint evidence did not support an exonic copy-number change. "
                  "Small single-exon heterozygous events and sequence variants are outside the scope of this test.</div></div>")
    xf = [e for e in res['exon_table'] if e['flag']]
    if xf:
        result = result[:-6] + (f"<div class='sub flag'>Exon-level review: exon "
                                f"{', '.join(str(e['num']) for e in xf)} "
                                f"{'has' if len(xf) == 1 else 'have'} unusual depth without a called event "
                                f"(see below) — check IGV / MLPA.</div></div>")
    flags = res['flags']
    qc_rows = [
        ('Input', res['bam']), ('Data type', f"{'WGS (whole genome)' if res['mode'] == 'wgs' else 'targeted capture'} "
              f"({res.get('mode_source', 'auto-detected')}; reads in DMD / mapped = {res['mode_fraction']:.4f})"),
        ('Sex used', f"{res['sex']} ({res['sex_source']}; X dosage {res['x_dosage']:.2f}, chrY fraction "
                     f"{res['y_fraction']:.5f})"),
        ('Reference panel', f"{res['panel']} (profile correlation {res['panel_corr']:.3f})"
         if res['mode'] == 'targeted' else res['panel']),
        ('Depth scale', f"{res['scale']:.1f}x ({'flank median' if res['mode'] == 'wgs' else 'on-target level'})"),
        ('Noise', f"robust SD of ratio {res['sample_sd']:.3f}; k = {res['noise_k']:.2f} x reference spread"),
        ('Usable bins', f"{res['n_valid']:,} / {res['n_bins']:,}"),
        ('Breakpoint caller', res['sawfish_status']),
        ('QC flags', ', '.join(flags) if flags else 'none'),
    ]
    qc = ''.join(f"<tr><td>{_e(k)}</td><td>{_e(v)}</td></tr>" for k, v in qc_rows)
    ev_rows = ''.join(
        f"<tr><td><span class='tier t{e['tier']}'>{e['tier']}</span></td><td>{_e(e['type'])}</td>"
        f"<td>{_e(e['exons_long'])}</td><td>{_e(e['frame'])}</td><td>chrX:{e['start'] + 1:,}-{e['end']:,}</td>"
        f"<td>{e['size_bp']:,}</td><td>{e['median_ratio']:.2f}</td><td>{e['est_CN']:.2f}</td>"
        f"<td>{_e(e['bp_text'])}</td><td>{_e(e['notes'])}</td></tr>"
        for e in res['events'] if e['tier'] in ('A', 'B', 'C'))
    events_html = (f"<div class='scroll'><table><tr><th>Tier</th><th>Type</th><th>Exons</th><th>Frame</th>"
                   f"<th>Depth interval (GRCh38, 1-based)</th><th>Size bp</th><th>Ratio</th><th>Est. CN</th>"
                   f"<th>Breakpoint (junction)</th><th>Notes</th></tr>{ev_rows}</table></div>"
                   if ev_rows else "<p class='sub'>No A/B/C events.</p>")
    xflags = ('<ul>' + ''.join(f"<li class='flag'>exon {e['num']}: ratio {e['ratio']:.2f} (z {e['z_sample']:.1f} vs this sample's exon noise) — "
                               f"{_e(e['flag'])}</li>" for e in xf) + '</ul>') if xf else \
        "<p class='sub'>None.</p>"
    ex_rows = ''.join(f"<tr><td>{e['num']}</td><td>chrX:{e['start'] + 1:,}-{e['end']:,}</td><td>{e['end'] - e['start']}</td>"
                      f"<td>{e['depth']:.1f}</td><td>{e['ratio']:.2f}</td><td>{e['z']:.1f}</td><td>{_e(e['state'])}</td></tr>"
                      for e in res['exon_table'])
    sv_only = res.get('sv_only', [])
    svo = ('<ul>' + ''.join(f"<li>{_e(s['type'])} {_e(s['exons'])} chrX:{s['start']:,}-{s['end']:,} "
                            f"(depth ratio {s['depth_ratio']:.2f}) — junction without depth change</li>" for s in sv_only)
           + '</ul>') if sv_only else "<p class='sub'>None.</p>"
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>DMD del/dup — {_e(res['sample'])}</title>
<style>{CSS}</style></head><body><main>
{_brand()}
<h1>DMD exon deletion / duplication — {_e(res['sample'])}</h1>
<div class="sub">{'WGS' if res['mode'] == 'wgs' else 'Targeted'} · reference {_e(res['reference_label'])} ({_e(res['reference_created'])}) ·
PacBio HiFi long-read · GRCh38 · MANE NM_004006.3 (79 exons)</div>
{result}
{GUIDE}
<h2>Exon-level view</h2><img src="{_img(exon_png)}" alt="exon-level depth ratio">
<h2>Gene view</h2><img src="{_img(gene_png)}" alt="DMD depth ratio">
<h2>Events</h2>{events_html}
<h2>Exon-level review flags (not calls)</h2>{xflags}
<h2>Breakpoints without a depth change (exonic)</h2>{svo}
<h2>Input and QC</h2><table class="kv">{qc}</table>
<details><summary>All 79 exons</summary><div class="scroll"><table><tr><th>Exon</th><th>GRCh38</th><th>bp</th>
<th>Depth</th><th>Ratio</th><th>z</th><th>State</th></tr>{ex_rows}</table></div></details>
<footer>DMDSeek {__version__} — developed at the {_e(LAB)}, {_e(DEPARTMENT)}, {_e(INSTITUTION)}.
{' · '.join(f'{_e(r)}: {_e(n)}' for r, n in PEOPLE)}. Research use only.</footer>
</main></body></html>"""
    with open(path, 'w', encoding='utf8') as f:
        f.write(doc)
