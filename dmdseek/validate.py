"""Compare dmdseek calls with a truth table (MLPA / orthogonal results)."""
import csv
import re
from collections import Counter
from pathlib import Path


def parse_exons(s):
    out = set()
    for part in re.split(r'[,;]', (s or '').strip()):
        part = part.strip().lower().replace('ex', '').replace('exon', '')
        m = re.match(r'^(\d+)\s*[~\-]\s*(\d+)$', part)
        if m:
            a, b = sorted((int(m.group(1)), int(m.group(2))))
            out |= set(range(a, b + 1))
        elif part.isdigit():
            out.add(int(part))
    return out


def label(nums):
    nums = sorted(nums)
    if not nums:
        return '.'
    runs, s, q = [], nums[0], nums[0]
    for n in nums[1:]:
        if n == q + 1:
            q = n
            continue
        runs.append((s, q))
        s = q = n
    runs.append((s, q))
    return ','.join(f'{a}' if a == b else f'{a}~{b}' for a, b in runs)


def validate(truth_path, summary_path, out_path):
    calls = {}
    with open(summary_path, encoding='utf8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            calls[r['sample']] = r
    events, review = {}, {}
    for s, r in calls.items():
        ev_file = Path(r['events_tsv'])
        if ev_file.is_file():
            with open(ev_file, encoding='utf8') as f:
                events[s] = [e for e in csv.DictReader(f, delimiter='\t') if e['tier'] in ('A', 'B')]
            js = Path(str(ev_file).replace('.events.tsv', '.summary.json'))
            if js.is_file():
                import json
                review[s] = set(json.load(open(js)).get('exon_review_flags') or [])
    rows, stat = [], Counter()
    from .tables import read_table
    for t in read_table(truth_path):
        s = t['sample_id']
        typ = (t.get('type') or 'none').strip().upper()
        if typ in ('EXCLUDE', 'NA', 'SKIP'):
            stat['EXCLUDED'] += 1
            continue
        typ = typ if typ in ('DEL', 'DUP') else 'NONE'
        tset = parse_exons(t.get('exons') or '')
        if s not in calls:
            status, got = 'NOT_RUN', None
        elif calls[s].get('flags') == 'ERROR':
            status, got = 'ERROR', None
        else:
            evs = events.get(s, [])
            hit = [e for e in evs if e['type'] == typ and
                   (parse_exons(e['exons']) | parse_exons(e['partial_exons'])) & tset]
            if typ == 'NONE':
                status, got = ('FP' if evs else 'TN'), (evs[0] if evs else None)
            elif hit:
                got = max(hit, key=lambda e: len(parse_exons(e['exons']) & tset))
                cset = parse_exons(got['exons']) | parse_exons(got['partial_exons'])
                status = 'EXACT' if cset == tset else 'BOUNDARY_DIFF'
            else:
                got = evs[0] if evs else None
                status = 'OTHER_EVENT' if evs else (
                    'MISSED(exon_review_flag)' if tset & review.get(s, set()) else 'MISSED')
        stat[status] += 1
        cset = (parse_exons(got['exons']) | parse_exons(got['partial_exons'])) if got else set()
        rows.append([s, typ, label(tset), got['type'] if got else '.', label(cset), got['tier'] if got else '.',
                     got['junction_match'] if got else '.', status, calls.get(s, {}).get('flags', '.')])
    with open(out_path, 'w', newline='', encoding='utf8') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['sample_id', 'truth_type', 'truth_exons', 'call_type', 'call_exons', 'tier', 'junction',
                    'status', 'qc_flags'])
        w.writerows(sorted(rows, key=lambda r: (r[7] in ('EXACT', 'TN'), r[7], r[0])))
    pos = sum(stat[k] for k in ('EXACT', 'BOUNDARY_DIFF', 'MISSED', 'MISSED(exon_review_flag)', 'OTHER_EVENT'))
    if stat['ERROR'] or stat['NOT_RUN']:
        print(f'[VALIDATE] {stat["ERROR"]} sample(s) failed and {stat["NOT_RUN"]} not run - not counted')
    neg = stat['TN'] + stat['FP']
    print(f'[VALIDATE] {dict(stat)}')
    if pos:
        print(f'[VALIDATE] positives detected {stat["EXACT"] + stat["BOUNDARY_DIFF"]}/{pos}; '
              f'exact exon range {stat["EXACT"]}/{pos}; missed but flagged for exon review '
              f'{stat["MISSED(exon_review_flag)"]}')
    if neg:
        print(f'[VALIDATE] negatives without A/B call {stat["TN"]}/{neg}')
    print(f'[VALIDATE] -> {out_path}')
