"""Terminal look: banner, colours, icons. Colours only on a terminal (not in log files);
NO_COLOR=1 turns colours off, DMDSEEK_PLAIN=1 also drops icons and kaomoji."""
import os
import sys
import time

from . import __version__
from .branding import DEPARTMENT, INSTITUTION, LAB

PLAIN = bool(os.environ.get('DMDSEEK_PLAIN'))
_enc = (getattr(sys.stdout, 'encoding', None) or '').lower()
FANCY = not PLAIN and 'utf' in _enc
COLOR = sys.stdout.isatty() and not os.environ.get('NO_COLOR') and not PLAIN

LOGO = r"""
 ___  __  __ ___  ___          _
|   \|  \/  |   \/ __| ___ ___| |__
| |) | |\/| | |) \__ \/ -_) -_) / /
|___/|_|  |_|___/|___/\___\___|_\_\
"""


def _c(code, s):
    return f'\033[{code}m{s}\033[0m' if COLOR else str(s)


def bold(s): return _c('1', s)
def dim(s): return _c('2', s)
def red(s): return _c('31;1', s)
def green(s): return _c('32;1', s)
def yellow(s): return _c('33;1', s)
def blue(s): return _c('34;1', s)
def magenta(s): return _c('35;1', s)
def cyan(s): return _c('36', s)


ICON = dict(dna='🧬', person='👤', link='🔗', chart='📊', ok='✅', warn='🔶', err='❌', doc='📄', box='📦',
            fish='🐟', folder='📂', search='🔍', sparkle='✨', clock='⏳', flag='🚩')
KAO = dict(start='(•̀ᴗ•́)و', happy='\\(^o^)/', warn='(^_^;)', err='(>_<)', think='(・_・ )')
TIER_COLOR = {'A': green, 'B': yellow, 'C': magenta, 'D': dim}


def icon(name):
    return ICON.get(name, '') if FANCY else ''


def kao(name):
    return KAO.get(name, '') if FANCY else {'happy': ':)', 'warn': ':/', 'err': ':(', 'start': '', 'think': ''}[name]


def banner(subtitle='DMD exon deletion / duplication caller for PacBio HiFi'):
    out = [cyan(LOGO.rstrip('\n'))] if FANCY else []
    out.append(f"  {bold('/// DMDSeek')} {icon('search')} {dim('v' + __version__)}   {subtitle}")
    out.append(dim(f'      {LAB} · {DEPARTMENT} · {INSTITUTION}'))
    return '\n'.join(out) + '\n'


def item(name, label, value, indent=5):
    ic = icon(name)
    return f"{' ' * indent}{ic + ' ' if ic else ''}{dim(f'{label:<12}')} {value}"


def call_color(text):
    t = str(text)
    if t.startswith(('DEL', 'DUP')):
        tier = t[-2] if t.endswith(']') else ''
        return TIER_COLOR.get(tier, bold)(t)
    if t.startswith('NO_EXONIC'):
        return green(t)
    if t.startswith(('NO_DEPTH', 'ERROR')):
        return red(t)
    return t


def elapsed(t0):
    s = int(time.time() - t0)
    return f'{s // 60}m{s % 60:02d}s' if s >= 60 else f'{s}s'


def rule(title=''):
    w = 72
    if not title:
        return dim('─' * w if FANCY else '-' * w)
    side = (w - len(title) - 2) // 2
    ch = '─' if FANCY else '-'
    return dim(ch * side) + ' ' + bold(title) + ' ' + dim(ch * side)


def table(rows, header):
    """Plain-width table (colour codes added after padding)."""
    widths = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
    widths[-1] = 0
    out = ['  ' + '  '.join(dim(f'{h:<{w}}') for h, w in zip(header, widths))]
    for r in rows:
        cells = [f'{str(v):<{w}}' for v, w in zip(r, widths)]
        cells[-1] = call_color(r[-1])
        out.append('  ' + '  '.join(cells))
    return '\n'.join(out)
