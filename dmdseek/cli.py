"""dmdseek command line: call | build-reference | validate | info"""
import argparse
import csv
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

from . import __version__, ui
from .branding import bundle_dir
from .model import Params


def default_reference():
    env = os.environ.get('DMDSEEK_REFERENCE')
    if env:
        return Path(env)
    return bundle_dir() / 'reference' / 'dmd_reference.npz'


def find_sawfish(user=None):
    for c in [user, os.environ.get('DMDSEEK_SAWFISH'), str(bundle_dir() / 'bin' / 'sawfish'), shutil.which('sawfish')]:
        if c and Path(c).is_file() and os.access(c, os.X_OK):
            return str(c)
    return None


def die(msg):
    print(f"{ui.red('ERROR')} {msg}  {ui.kao('err')}", file=sys.stderr)
    sys.exit(2)


def cmd_call(a):
    from .bamio import InputError
    from .pipeline import call_sample, sample_name
    from .reference import Reference
    mode = 'wgs' if a.wgs else 'targeted' if a.target else a.mode
    jobs = [(None, b) for b in a.bam]
    if a.list:
        from .tables import read_sample_sheet
        try:
            rows_in = read_sample_sheet(a.list)
        except (OSError, ValueError) as exc:
            die(str(exc))
        jobs += [(r['sample_id'], r['bam']) for r in rows_in]
    if not jobs:
        a.parser.print_help()
        die('give BAM file(s) or --list')
    if len(jobs) > 1 and a.sample:
        die('--sample can only be used with a single BAM')
    if not Path(a.reference).is_file():
        die(f'no reference panel at {a.reference}\n'
            f'      build one from your own cohort:  dmdseek build-reference --target-bams list.tsv -o my_ref.npz\n'
            f'      then:  dmdseek call ... --reference my_ref.npz   (or copy it to {bundle_dir() / "reference" / "dmd_reference.npz"})')
    try:
        ref = Reference(a.reference)
    except (OSError, ValueError, KeyError) as exc:
        die(f'reference panel {a.reference} not usable ({exc})')
    exe = find_sawfish(a.sawfish) if a.ref else None
    p = Params()
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    print(ui.banner())
    print(ui.item('box', 'reference', f"{ref.meta.get('label')}  {ui.dim('panels: ' + ', '.join(ref.meta['panels']))}", 2))
    if exe:
        from .sawfish import sawfish_version
        saw = sawfish_version(exe) or exe
    else:
        saw = ui.yellow('not used') + ui.dim(' (give --ref GRCh38.fa for breakpoints)' if not a.ref else ' (not found)')
    print(ui.item('fish', 'breakpoints', saw, 2))
    print(ui.item('dna', 'data type', {'auto': 'detected per sample (WGS or targeted)', 'wgs': 'WGS (set by --wgs)',
                                       'targeted': 'targeted capture (set by --target)'}[mode], 2))
    print(ui.item('folder', 'output', str(out), 2))
    print(f"\n  {ui.kao('start')} {len(jobs)} sample{'s' if len(jobs) > 1 else ''} to go\n")

    rows, table, failed, t_all = [], [], 0, time.time()
    for k, (given, bam) in enumerate(jobs, 1):
        t0 = time.time()
        name = a.sample or given or sample_name(bam)
        print(f"  {ui.bold(f'[{k}/{len(jobs)}]')} {ui.bold(name)}")
        try:
            sdir = out / name if (len(jobs) > 1 or a.list) else out
            res = call_sample(bam, ref, p, sdir, sample=name, mode=mode, sex=a.sex, fasta=a.ref,
                              sawfish_exe=exe, threads=a.threads, cram_ref=a.ref,
                              forced_mode=mode != 'auto', log=print)
            good = res['final_call'].startswith(('NO_EXONIC', 'DEL', 'DUP')) and not any(
                f.startswith(('SAWFISH_FAILED', 'NO_DEPTH', 'NOISY', 'CHECK_X', 'SEX_UND', 'NO_ONTARGET'))
                for f in res['flags'])
            review = [e['num'] for e in res['exon_table'] if e['flag']]
            print(ui.item('ok' if good and not review else 'warn', 'result', ui.call_color(res['final_call'])))
            if review:
                print(ui.item('flag', 'exon review', ui.yellow(f"exon {', '.join(map(str, review))}: unusual depth, "
                                                               f"not called - check IGV / MLPA") + '  ' + ui.kao('think')))
            if res['flags']:
                print(ui.item('warn', 'QC flags', ui.yellow(', '.join(res['flags'])) + '  ' + ui.kao('warn')))
            report = sdir / f'{name}.dmdseek.report.html'
            print(ui.item('doc', 'report', str(report)))
            print(ui.item('clock', 'time', ui.elapsed(t0)) + '\n')
            rows.append([name, res['mode'], res['sex'], res['final_call'], ';'.join(res['flags']) or 'PASS',
                         str(sdir / f'{name}.dmdseek.events.tsv'), str(report), bam])
            table.append([name, res['mode'], res['sex'], res['final_call'] + (
                f"  (review: ex{','.join(map(str, review))})" if review else '')])
        except (InputError, FileNotFoundError, ValueError) as exc:
            failed += 1
            print(ui.item('err', 'ERROR', ui.red(str(exc)) + '  ' + ui.kao('err')) + '\n', file=sys.stderr)
            rows.append([name, '.', '.', f'ERROR: {exc}', 'ERROR', '.', '.', bam])
            table.append([name, '.', '.', 'ERROR'])
        except Exception:
            failed += 1
            traceback.print_exc()
            rows.append([name, '.', '.', 'ERROR (see log)', 'ERROR', '.', '.', bam])
            table.append([name, '.', '.', 'ERROR'])
    with open(out / 'dmdseek_summary.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['sample', 'mode', 'sex', 'final_call', 'flags', 'events_tsv', 'report_html', 'bam'])
        w.writerows(rows)
    print(ui.rule('summary'))
    print(ui.table(table, ['sample', 'mode', 'sex', 'result']))
    print(ui.rule())
    face = ui.kao('happy') if not failed else ui.kao('err')
    print(f"  {face}  {len(jobs) - failed}/{len(jobs)} done in {ui.elapsed(t_all)} {icon_arrow()} "
          f"{out / 'dmdseek_summary.tsv'}")
    print(ui.dim('  Tiers: A = depth + breakpoint agree · B = review · C = intronic · see the report for details\n'))
    sys.exit(1 if failed else 0)


def icon_arrow():
    return '→' if ui.FANCY else '->'


def cmd_build(a):
    from .reference import build_reference
    if not (a.wgs_pon or a.target_pon or a.wgs_bams or a.target_bams):
        die('give --wgs-bams / --target-bams (or --wgs-pon / --target-pon)')
    if (a.wgs_pon and a.wgs_bams) or (a.target_pon and a.target_bams):
        die('use either BAM lists or PoN files per mode, not both')
    if not a.exon_bed:
        die('--exon-bed is required (no bundled data/DMD_MANE_79exons.bed)')
    build_reference(a)


def cmd_validate(a):
    from .validate import validate
    validate(a.truth, a.summary, a.out)


def cmd_info(a):
    from .branding import logo_file
    from .reference import Reference
    print(ui.banner())
    print(ui.item('box', 'bundle', str(bundle_dir()), 2))
    exe = find_sawfish()
    if exe:
        from .sawfish import sawfish_version
        print(ui.item('fish', 'sawfish', f'{exe} ({sawfish_version(exe)})', 2))
    else:
        print(ui.item('warn', 'sawfish', 'not found (breakpoint step will be skipped)', 2))
    lf = logo_file()
    print(ui.item('sparkle', 'logo', str(lf) if lf else 'none (text header only)', 2))
    try:
        ref = Reference(a.reference)
        print(ui.item('chart', 'reference', a.reference, 2))
        for k, v in ref.meta.items():
            if k != 'panels':
                print(f'      {ui.dim(k + ":")} {v}')
        for k, v in ref.meta['panels'].items():
            print(f'      {ui.dim("panel " + k + ":")} ' + ', '.join(f'{x}={y}' for x, y in v.items() if x != 'normalisation'))
        print(f'      {ui.dim("recurrent intronic events:")} {len(ref.recurrent)}')
    except Exception as exc:
        print(ui.item('err', 'reference', f'not usable ({exc})', 2))


# ------------------------------------------------------------------ argument parsing
class Parser(argparse.ArgumentParser):
    def format_help(self):
        return ui.banner() + '\n' + super().format_help()

    def error(self, message):
        self.print_usage(sys.stderr)
        die(f'{message}   (see: {self.prog} -help)')


def _help(p):
    p.add_argument('-h', '-help', '--help', action='help', help='show this help and exit')


EXAMPLES = """\
quick start:
  dmdseek call sample.bam --ref GRCh38.fa -o out/              WGS or targeted is detected from the reads
  dmdseek call sample.bam --ref GRCh38.fa -o out/ -target      say it is a targeted (capture) BAM
  dmdseek call sample.bam --ref GRCh38.fa -o out/ -wgs         say it is a whole-genome BAM
  dmdseek call --list samples.csv --ref GRCh38.fa -o out/ -t 8 many samples (CSV/TSV/XLSX: sample_id, bam)
  dmdseek info                                                 version, bundled sawfish, reference panel
  dmdseek <command> -help                                      all options of one command

output per sample: <sample>.dmdseek.report.html (open in a browser), events / exons TSV, summary JSON.
Tiers: A = depth and breakpoint agree · B = needs review · C = intronic · D = weak (TSV only).
Research use only - confirm findings by MLPA / qPCR and IGV.

public release: first install sawfish (conda install -c bioconda sawfish) and build a reference panel
from your own cohort (dmdseek build-reference; WGS >=8 males and >=8 females, targeted >=8 per capture
batch, >=20 recommended) - see README "Before you start".
"""

CALL_EXAMPLES = """\
examples:
  dmdseek call S1.bam --ref GRCh38.fa -o out/
  dmdseek call S1.bam S2.bam --ref GRCh38.fa -o out/ -target -t 8
  dmdseek call --list samples.xlsx --ref GRCh38.fa -o out/

sample sheet (--list, build-reference): CSV, TSV or Excel (.xlsx, first sheet) with columns
sample_id, bam (+ optional sex M/F, batch). Paths: absolute, relative to the sheet, or Windows
paths (E:\\data\\x.bam -> /mnt/e/data/x.bam on WSL). Templates: templates/ in the bundle.

data type: by default each BAM is classified from its reads (targeted if >=1% of mapped reads fall in
DMD or the BAM has <1M mapped reads); -wgs / -target force it. The chosen type is printed per sample
and in the report.
"""


def main(argv=None):
    ap = Parser(prog='dmdseek', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter,
                description='DMD exon deletion/duplication calling from PacBio HiFi BAMs (WGS or targeted), GRCh38.',
                epilog=EXAMPLES)
    _help(ap)
    ap.add_argument('-v', '-version', '--version', action='version', version=f'DMDSeek {__version__}')
    sub = ap.add_subparsers(dest='cmd', metavar='<command>', title='commands')

    c = sub.add_parser('call', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter,
                       help='call DMD exon DEL/DUP in one or more BAMs  (start here)',
                       description='Call DMD exon deletions / duplications and write an HTML report per sample.',
                       epilog=CALL_EXAMPLES)
    g = c.add_argument_group('input')
    g.add_argument('bam', nargs='*', help='indexed GRCh38 BAM(s) (or CRAM with --ref)')
    g.add_argument('--list', '-list', metavar='TSV', help='sample sheet (CSV, TSV or XLSX) with columns sample_id, bam')
    g.add_argument('--ref', '-ref', metavar='FASTA', help='GRCh38 FASTA (+ .fai): enables breakpoints (recommended)')
    g.add_argument('-o', '--out', '-out', required=True, metavar='DIR', help='output folder')
    t = c.add_argument_group('data type (default: detected from the reads)')
    tx = t.add_mutually_exclusive_group()
    tx.add_argument('-wgs', '--wgs', action='store_true', help='whole-genome BAM')
    tx.add_argument('-target', '--target', '--targeted', action='store_true', help='targeted capture BAM (e.g. Twist DMD panel)')
    tx.add_argument('--mode', choices=['auto', 'wgs', 'targeted'], default='auto', help=argparse.SUPPRESS)
    o = c.add_argument_group('options')
    o.add_argument('--sample', '-sample', help='sample name (single BAM only); default from the file / folder name')
    o.add_argument('--sex', '-sex', choices=['M', 'F'], help='override the sex inferred from the reads')
    o.add_argument('--reference', '-reference', metavar='NPZ', default=str(default_reference()),
                   help='reference panel (default: the bundled one)')
    o.add_argument('--sawfish', metavar='EXE', help='sawfish executable (default: bundled, then PATH)')
    o.add_argument('-t', '--threads', '-threads', type=int, default=4, metavar='N', help='threads (default 4)')
    _help(o)
    c.set_defaults(func=cmd_call, parser=c)

    b = sub.add_parser('build-reference', add_help=False, help='build an anonymised reference panel from a cohort')
    bed = bundle_dir() / 'data' / 'DMD_MANE_79exons.bed'
    b.add_argument('--exon-bed', default=str(bed) if bed.is_file() else None,
                   help='79-exon BED of MANE NM_004006.3 on GRCh38 (exon_01..exon_79); default: data/ of the bundle')
    b.add_argument('--wgs-bams', metavar='SHEET', help='sample sheet (CSV/TSV/XLSX): sample_id, bam [, sex] - WGS cohort')
    b.add_argument('--target-bams', metavar='SHEET', help='sample sheet (CSV/TSV/XLSX): sample_id, bam [, sex, batch] - '
                   'targeted cohort (BAM incl. off-target reads)')
    b.add_argument('--target-window', default='31097676-33339609', help='capture interval for --target-bams (0-based start-end)')
    b.add_argument('-t', '--threads', type=int, default=4, help='BAMs read in parallel')
    b.add_argument('--wgs-pon', help='alternative: CJH_DMD_PoN.npz (dmd_pon.py)')
    b.add_argument('--wgs-pon-exons', help='CJH_DMD_PoN_exons.npz')
    b.add_argument('--wgs-events', nargs='*', help='events TSV of the WGS cohort (carrier bins masked)')
    b.add_argument('--target-pon', help='DMD_target_PoN.npz (dmd_target_prep.py)')
    b.add_argument('--target-pon-exons', help='DMD_target_PoN_exons.npz')
    b.add_argument('--target-qc', help='DMD_target_input_QC.tsv (off-target reads for X dosage)')
    b.add_argument('--target-events', nargs='*', help='events TSV of the targeted cohort')
    b.add_argument('--exclude-samples', nargs='*', help='text files with sample IDs to leave out')
    b.add_argument('--label', default='dmdseek reference', help='free-text label stored in the reference (no identifiers)')
    b.add_argument('-o', '--out', required=True)
    _help(b)
    b.set_defaults(func=cmd_build)

    v = sub.add_parser('validate', add_help=False, help='compare calls with a truth table (e.g. MLPA)')
    v.add_argument('--truth', required=True, help='CSV/TSV/XLSX: sample_id, type (DEL/DUP/NONE/EXCLUDE), exons (e.g. 45~47)')
    v.add_argument('--summary', required=True, help='dmdseek_summary.tsv from "dmdseek call"')
    v.add_argument('-o', '--out', required=True)
    _help(v)
    v.set_defaults(func=cmd_validate)

    i = sub.add_parser('info', add_help=False, help='version, bundled sawfish, logo and reference panel')
    i.add_argument('--reference', default=str(default_reference()))
    _help(i)
    i.set_defaults(func=cmd_info)

    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        ap.print_help()
        sys.exit(0)
    subs = {'call': c, 'build-reference': b, 'validate': v}
    if len(argv) == 1 and argv[0] in subs:           # bare sub-command -> its help
        subs[argv[0]].print_help()
        sys.exit(0)
    a = ap.parse_args(argv)
    if not getattr(a, 'func', None):
        ap.print_help()
        sys.exit(0)
    try:
        a.func(a)
    except (ValueError, OSError, KeyError) as exc:      # input problems: message, not a traceback
        die(str(exc))


if __name__ == '__main__':
    main()
