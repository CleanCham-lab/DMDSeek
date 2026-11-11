#!/usr/bin/env bash
set -Eeuo pipefail
# =====================================================================
# K-fold cross-validation of the targeted reference.
# Families are kept in the same fold (stratified by batch). For each fold the reference is built
# WITHOUT that fold's samples, then those samples are called with it -> honest estimate for new samples.
#
#   bash crossval_target.sh <dmdseek executable> <tgt_bams.tsv> <truth.tsv> <out_dir> [K=5]
#
# env (optional): EXON_BED (default: the bundled 79-exon BED), FASTA (breakpoints; strongly recommended), THREADS (4), TARGET_EVENTS (carrier masks),
#                 TARGET_PON + TARGET_PON_EXONS + TARGET_QC  -> fast reference build from an existing PoN
#                 (otherwise the training BAMs are re-read for every fold)
# =====================================================================
EXE=$1; LIST=$2; TRUTH=$3; OUT=$4; K=${5:-5}
THREADS=${THREADS:-4}
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)

python3 - "$LIST" "$OUT" "$K" <<'PY'
import csv, random, sys, collections
lst, out, K = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = list(csv.DictReader(open(lst), delimiter='\t'))
fam = collections.defaultdict(list)
for r in rows:
    fam[(r.get('batch') or '', r.get('family') or r['sample_id'])].append(r)
fold = {}
for b in sorted({k[0] for k in fam}):
    keys = sorted(k for k in fam if k[0] == b)
    random.Random(1).shuffle(keys)
    for i, k in enumerate(keys):
        fold[k] = i % K
for k in range(K):
    test = [r for key, rs in fam.items() if fold[key] == k for r in rs]
    train = [r for key, rs in fam.items() if fold[key] != k for r in rs]
    for name, rs in (('test', test), ('train', train)):
        with open(f'{out}/fold{k + 1}.{name}.tsv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(rs)
    with open(f'{out}/fold{k + 1}.test.ids', 'w') as f:
        f.write(''.join(r['sample_id'] + '\n' for r in test))
    print(f'[cv] fold {k + 1}: test {len(test)} samples, train {len(train)}')
PY

: > "$OUT/cv_summary.tsv"
for k in $(seq 1 "$K"); do
    echo "[cv] ===== fold $k/$K"
    if [[ -n "${TARGET_PON:-}" ]]; then
        SRC=(--target-pon "$TARGET_PON" --target-pon-exons "$TARGET_PON_EXONS" --target-qc "$TARGET_QC"
             --exclude-samples "$OUT/fold$k.test.ids")
    else
        SRC=(--target-bams "$OUT/fold$k.train.tsv")
    fi
    "$EXE" build-reference ${EXON_BED:+--exon-bed $EXON_BED} "${SRC[@]}" ${TARGET_EVENTS:+--target-events $TARGET_EVENTS} \
        --label "cv-fold$k" -t "$THREADS" -o "$OUT/fold$k.reference.npz" 2>&1 | grep -E '^\[(TARGETED|REFERENCE|WARN)' || true
    "$EXE" call --list "$OUT/fold$k.test.tsv" --reference "$OUT/fold$k.reference.npz" \
        ${FASTA:+--ref $FASTA} -t "$THREADS" -o "$OUT/fold$k" > "$OUT/fold$k.log" 2>&1 || echo "[cv] fold $k: some samples failed (see fold$k.log)"
    if [[ ! -s "$OUT/cv_summary.tsv" ]]; then head -1 "$OUT/fold$k/dmdseek_summary.tsv" > "$OUT/cv_summary.tsv"; fi
    tail -n +2 "$OUT/fold$k/dmdseek_summary.tsv" >> "$OUT/cv_summary.tsv"
done
"$EXE" validate --truth "$TRUTH" --summary "$OUT/cv_summary.tsv" -o "$OUT/cv_validation.tsv"
