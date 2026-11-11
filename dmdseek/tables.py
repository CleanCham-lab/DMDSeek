"""Sample sheets and truth tables: CSV, TSV/TXT or Excel (.xlsx, first sheet).

Column names are matched without regard to case, spaces or hyphens, and common alternatives are
accepted (e.g. "Sample ID", "BAM path", "Gender"). BAM paths may be absolute, relative to the sheet's
folder, or Windows paths (E:\\data\\x.bam), which are converted to /mnt/e/data/x.bam under WSL.
"""
import csv
import io
import os
import re
from pathlib import Path

ALIASES = {
    'sample_id': ('sample_id', 'sample', 'sampleid', 'id', 'name', 'sample_name'),
    'bam': ('bam', 'bam_path', 'bampath', 'path', 'file', 'bam_file', 'cram', 'mapped_bam'),
    'sex': ('sex', 'gender', '성별'),
    'batch': ('batch', 'capture_batch', 'run', 'batch_id'),
    'type': ('type', 'mlpa', 'call', 'truth_type'),
    'exons': ('exons', 'mlpa_exons', 'exon', 'truth_exons'),
}
_KEYMAP = {v: k for k, vs in ALIASES.items() for v in vs}


def _norm_key(k):
    k = re.sub(r'[\s\-]+', '_', str(k or '').strip().lower())
    return _KEYMAP.get(k, _KEYMAP.get(k.replace('_', ''), k))


def _cell(v):
    if v is None:
        return ''
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v).strip()


def _read_text(path):
    raw = Path(path).read_bytes()
    for enc in ('utf-8-sig', 'cp949', 'latin-1'):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode('utf-8', 'replace')


def read_table(path):
    """-> list of dicts with normalised column names (header = first non-empty row)."""
    p = Path(path)
    ext = p.suffix.lower()
    if ext in ('.xlsx', '.xlsm'):
        try:
            import openpyxl
        except ImportError as exc:                       # pragma: no cover
            raise ValueError(f'{path}: reading Excel needs openpyxl ({exc}); save the sheet as CSV') from exc
        wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
        grid = [[_cell(v) for v in row] for row in wb.worksheets[0].iter_rows(values_only=True)]
        wb.close()
    elif ext == '.xls':
        raise ValueError(f'{path}: old .xls format - save it as .xlsx or .csv')
    else:
        text = _read_text(p)
        lines = [l for l in text.splitlines() if l.strip() and not l.lstrip().startswith('#')]
        if ext == '.csv':
            delim = ','
        elif ext in ('.tsv', '.txt'):
            delim = '\t'
        else:
            delim = '\t' if '\t' in (lines[0] if lines else '') else ','
        grid = [[c.strip() for c in row] for row in csv.reader(io.StringIO('\n'.join(lines)), delimiter=delim)]
    grid = [r for r in grid if any(c for c in r)]
    if not grid:
        return []
    raw = [re.sub(r'[\s\-]+', '_', str(h or '').strip().lower()) for h in grid[0]]
    header = [_norm_key(h) for h in grid[0]]
    # several columns mapping to one name (e.g. "sample_id" and "name"): keep the exact one, else the first
    keep = {}
    for i, h in enumerate(header):
        if h and (h not in keep or (raw[i] == h and raw[keep[h]] != h)):
            keep[h] = i
    rows = []
    for r in grid[1:]:
        r = r + [''] * (len(header) - len(r))
        rows.append({h: r[i] for h, i in keep.items()})
    return rows


def _resolve_path(p, base):
    p = os.path.expanduser(p.strip().strip('"').strip("'"))
    m = re.match(r'^([A-Za-z]):[\\/](.*)$', p)
    if m and not Path(p).exists():                       # Windows path on WSL
        p = f"/mnt/{m.group(1).lower()}/{m.group(2).replace(chr(92), '/')}"
    if not os.path.isabs(p):
        p = os.path.abspath(os.path.join(base, p))       # no symlink resolution (index may sit next to the link)
    return p


def _sex(v):
    v = (v or '').strip().lower()
    if v in ('m', 'male', 'man', 'xy', '남', '남자', '남성'):
        return 'M'
    if v in ('f', 'female', 'woman', 'xx', '여', '여자', '여성'):
        return 'F'
    return ''


def read_sample_sheet(path):
    """Sample sheet with columns sample_id, bam [, sex, batch] -> list of dicts (paths resolved)."""
    rows = read_table(path)
    if not rows or 'sample_id' not in rows[0] or 'bam' not in rows[0]:
        cols = ', '.join(rows[0]) if rows else 'none'
        raise ValueError(f'{path}: needs columns sample_id and bam (optional: sex, batch); found: {cols}')
    base = os.path.dirname(os.path.abspath(path))
    out, seen = [], set()
    for i, r in enumerate(rows, 2):
        sid, bam = r.get('sample_id', ''), r.get('bam', '')
        if not sid and not bam:
            continue
        if not sid or not bam:
            raise ValueError(f'{path}: row {i} has an empty sample_id or bam')
        if sid in seen:
            raise ValueError(f'{path}: sample_id {sid} appears twice')
        seen.add(sid)
        r = dict(r, sample_id=sid, bam=_resolve_path(bam, base), sex=_sex(r.get('sex')),
                 batch=(r.get('batch') or '').strip())
        out.append(r)
    return out
