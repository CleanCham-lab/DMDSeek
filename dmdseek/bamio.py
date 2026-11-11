"""BAM access with pysam: depth per bin / exon (mosdepth-equivalent), index stats, mode and sex."""
import numpy as np
import pysam

from .constants import (CHROM, CHRX_LEN_GRCH38, DMD_END, DMD_START, EXCLUDE_FLAGS, MIN_MAPQ,
                        TARGETED_MIN_FRACTION)

AUTOSOMES = {f'chr{i}' for i in range(1, 23)} | {str(i) for i in range(1, 23)}


class InputError(Exception):
    pass


def open_bam(path, ref=None):
    try:
        bam = pysam.AlignmentFile(str(path), 'rc' if str(path).endswith('.cram') else 'rb',
                                  reference_filename=ref)
    except (ValueError, OSError) as exc:
        raise InputError(f'cannot open {path}: {exc}')
    if not bam.has_index():
        raise InputError(f'{path} has no index (.bai/.csi). Run: samtools index {path}')
    lengths = dict(zip(bam.references, bam.lengths))
    if CHROM not in lengths:
        raise InputError(f'{path}: no {CHROM} contig (UCSC "chr" names on GRCh38 required)')
    if lengths[CHROM] != CHRX_LEN_GRCH38:
        raise InputError(f'{path}: {CHROM} length {lengths[CHROM]} != GRCh38 {CHRX_LEN_GRCH38}')
    return bam


def index_stats(bam):
    """Mapped reads per contig from the index (samtools idxstats equivalent)."""
    return {s.contig: s.mapped for s in bam.get_index_statistics()}


def coverage(bam, start, end, chrom=CHROM):
    """Per-base depth over [start, end): MAPQ >= 20, flags 1796 excluded, only aligned blocks
    (M/=/X; deletions and N not counted), supplementary counted - as mosdepth default."""
    n = end - start
    diff = np.zeros(n + 1, dtype=np.int64)
    a_list, b_list = [], []
    nreads = 0
    for read in bam.fetch(chrom, start, end):
        if read.flag & EXCLUDE_FLAGS or read.mapping_quality < MIN_MAPQ:
            continue
        nreads += 1
        for bs, be in read.get_blocks():
            if be <= start or bs >= end:
                continue
            a_list.append(max(bs, start) - start)
            b_list.append(min(be, end) - start)
    if a_list:
        np.add.at(diff, np.asarray(a_list), 1)
        np.add.at(diff, np.asarray(b_list), -1)
    return np.cumsum(diff[:-1]).astype(np.float64), nreads


def interval_means(cov, cov_start, starts, ends):
    """Mean depth of each [start, end) interval (all inside the coverage window)."""
    cs = np.concatenate([[0.0], np.cumsum(cov)])
    a = np.asarray(starts) - cov_start
    b = np.asarray(ends) - cov_start
    return (cs[b] - cs[a]) / np.maximum(b - a, 1)


def detect_mode(bam, stats):
    """targeted if >=1% of mapped reads fall in DMD, or if the BAM is too small to be a WGS
    (<1 M mapped reads, i.e. <~5x HiFi) - the latter catches a male whole-gene deletion in capture data."""
    total = sum(stats.values())
    in_dmd = bam.count(CHROM, DMD_START, DMD_END)
    frac = in_dmd / max(total, 1)
    targeted = frac >= TARGETED_MIN_FRACTION or total < 1_000_000
    return ('targeted' if targeted else 'wgs'), frac, total, in_dmd


def sex_from_stats(bam, stats, mode):
    """Returns (call, details). WGS: X dosage (chrX vs autosome reads per Mb) and chrY fraction.
    Targeted: chrY fraction of off-target reads (X dosage is added later from the capture scale)."""
    lengths = dict(zip(bam.references, bam.lengths))
    y = stats.get('chrY', 0) + stats.get('Y', 0)
    if mode == 'wgs':
        auto_reads = sum(v for c, v in stats.items() if c in AUTOSOMES)
        auto_len = sum(lengths[c] for c in stats if c in AUTOSOMES)
        x_dos = (stats.get(CHROM, 0) / lengths[CHROM]) / max(auto_reads / max(auto_len, 1), 1e-12)
        y_frac = y / max(sum(stats.values()), 1)
        by_x = 'U' if auto_reads == 0 else ('M' if x_dos < 0.7 else ('F' if x_dos > 0.8 else 'U'))
        by_y = 'U' if auto_reads == 0 else ('M' if y_frac >= 1e-3 else ('F' if y_frac <= 3e-4 else 'U'))
        call = by_x if by_x != 'U' else by_y
        flag = '' if (by_x == by_y or 'U' in (by_x, by_y)) else f'SEX_SIGNALS_DISAGREE(X={by_x},Y={by_y})'
        return call, dict(x_dosage=x_dos, y_fraction=y_frac, by_x=by_x, by_y=by_y, flag=flag)
    off = sum(stats.values()) - stats.get(CHROM, 0) - y
    y_frac = y / max(off, 1)
    by_y = 'M' if y_frac >= 5e-4 else ('F' if y_frac <= 2e-4 else 'U')
    return by_y, dict(x_dosage=float('nan'), y_fraction=y_frac, by_x='.', by_y=by_y, flag='',
                      offtarget_reads=off)
