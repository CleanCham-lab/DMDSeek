"""Breakpoint evidence: Sawfish on a regional chrX BAM (pysam extraction, SA tags restricted)."""
import gzip
import re
import shutil
import subprocess
from pathlib import Path

import pysam

from .constants import CHROM, DMD_END, DMD_START, SV_WINDOW

MATE_RE = re.compile(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]')


def sawfish_version(exe):
    try:
        out = subprocess.run([exe, '--version'], capture_output=True, text=True, timeout=30)
        return (out.stdout or out.stderr).strip().split('\n')[0]
    except (OSError, subprocess.SubprocessError):
        return None


def regional_bam(bam_path, out_bam, ref=None):
    """chrX:29-35.5 Mb reads into a BAM whose header has only chrX; SA segments outside the window
    (other chromosomes, ALT contigs, distal chrX) are removed so Sawfish accepts the file."""
    lo, hi = SV_WINDOW
    src = pysam.AlignmentFile(str(bam_path), 'rc' if str(bam_path).endswith('.cram') else 'rb',
                              reference_filename=ref)
    hdr = src.header.to_dict()
    xlen = dict(zip(src.references, src.lengths))[CHROM]
    hdr['SQ'] = [{'SN': CHROM, 'LN': xlen}]
    hdr.setdefault('PG', []).append({'ID': 'dmdseek-regional', 'PN': 'dmdseek'})
    n = changed = 0
    with pysam.AlignmentFile(str(out_bam), 'wb', header=hdr) as out:
        for read in src.fetch(CHROM, lo, hi):
            d = read.to_dict()
            if d.get('next_ref_name') not in ('=', '*', CHROM):
                d['next_ref_name'], d['next_ref_pos'] = '*', '0'
            tags = []
            for t in d.get('tags', []):
                if not t.startswith('SA:Z:'):
                    tags.append(t)
                    continue
                keep = []
                for seg in t[5:].split(';'):
                    f = seg.split(',')
                    if len(f) >= 6 and f[0] == CHROM and lo <= int(f[1]) <= hi:
                        keep.append(seg)
                new = 'SA:Z:' + ';'.join(keep) + ';' if keep else None
                if new != t:
                    changed += 1
                if new:
                    tags.append(new)
            d['tags'] = tags
            out.write(pysam.AlignedSegment.from_dict(d, out.header))
            n += 1
    src.close()
    pysam.index(str(out_bam))
    return n, changed


def run(bam_path, workdir, ref_fasta, exe, threads=4, cram_ref=None, log=print):
    workdir = Path(workdir)
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    small = workdir / 'regional.chrX_29M_35.5M.bam'
    n, changed = regional_bam(bam_path, small, cram_ref)
    from . import ui
    log(ui.item('link', 'breakpoints', f'Sawfish on {n:,} alignments (chrX:29-35.5 Mb) ...'))
    if n == 0:
        return []
    disc, jc = workdir / 'discover', workdir / 'joint_call'
    with open(workdir / 'sawfish.log', 'w') as lf:
        # --cov-regex '.': the regional BAM has only chrX (Sawfish's default depth regex wants autosomes)
        for cmd in ([exe, 'discover', '--threads', str(threads), '--ref', str(ref_fasta), '--bam', str(small),
                     '--cov-regex', '.', '--disable-cnv', '--output-dir', str(disc)],
                    [exe, 'joint-call', '--threads', str(threads), '--sample', str(disc), '--output-dir', str(jc)]):
            lf.write('$ ' + ' '.join(cmd) + '\n')
            lf.flush()
            r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
            if r.returncode != 0:
                raise RuntimeError(f'sawfish failed ({r.returncode}); see {workdir / "sawfish.log"}')
    vcf = jc / 'genotyped.sv.vcf.gz'
    if not vcf.is_file():
        raise RuntimeError(f'sawfish produced no {vcf}')
    return parse_vcf(vcf)


def parse_vcf(vcf):
    """Records overlapping DMD, or BND with a mate inside DMD; coordinates 1-based POS -> as given."""
    out = []
    op = gzip.open if str(vcf).endswith('.gz') else open
    with op(vcf, 'rt', encoding='utf8', errors='replace') as f:
        for line in f:
            if line.startswith('#'):
                continue
            x = line.rstrip('\n').split('\t')
            if len(x) < 8 or x[0] != CHROM:
                continue
            pos = int(x[1])
            info = dict(kv.split('=', 1) if '=' in kv else (kv, '1') for kv in x[7].split(';'))
            typ = info.get('SVTYPE', '.')
            try:
                end = int(info.get('END', pos))
            except ValueError:
                end = pos
            ref, alt = x[3], x[4]
            if typ != 'BND' and end == pos and len(ref) > 1:
                end = pos + len(ref) - 1
            mates = [(m.group(1), int(m.group(2))) for m in MATE_RE.finditer(alt)]
            in_dmd = min(pos, end) <= DMD_END and max(pos, end) >= DMD_START
            mate_in = any(c == CHROM and DMD_START <= q <= DMD_END for c, q in mates)
            if not (in_dmd or mate_in):
                continue
            gt = '.'
            if len(x) >= 10 and 'GT' in x[8].split(':'):
                gt = x[9].split(':')[x[8].split(':').index('GT')]
            ins_len = len(alt) - 1 if typ == 'INS' and not alt.startswith('<') else 0
            out.append(dict(caller='Sawfish', type=typ, start=pos, end=max(end, pos + 1), filt=x[6],
                            qual=x[5], id=x[2], ins_len=ins_len, mates=mates, gt=gt, alt=alt,
                            svclaim=info.get('SVCLAIM', '.')))
    return out
