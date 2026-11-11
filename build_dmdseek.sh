#!/usr/bin/env bash
set -Eeuo pipefail
# =====================================================================
# Build the self-contained dmdseek folder for Linux x86_64 (WSL included):
#
#   dmdseek/
#     dmdseek                      executable (Python + numpy + pysam + matplotlib inside)
#     _internal/                  runtime libraries
#     bin/sawfish                 breakpoint caller (copied from your installation)
#     reference/dmd_reference.npz reference panel (aggregated statistics only)
#     data/DMD_MANE_79exons.bed   MANE NM_004006.3 exons on GRCh38 (copied from --exon-bed)
#     README.md
#
#   1) bash build_dmdseek.sh --exon-bed DMD_MANE_79exons.bed [--logo logo.png|svg] [--sawfish /path/sawfish] [--python python3]
#   2) build/dist/dmdseek/dmdseek build-reference ... -o build/dist/dmdseek/reference/dmd_reference.npz
#   3) bash build_dmdseek.sh --package [--public]                                 # smoke test + tar.gz
#      (--public: without sawfish binary, logo and reference panel -> dmdseek-<ver>-linux-x86_64-public.tar.gz)
#   (or give --reference /path/dmd_reference.npz in step 1 to do everything at once)
#
# Build on the OLDEST Linux you want to support: the executable runs on that glibc and newer
# (built on Ubuntu 22.04 -> runs on 22.04/24.04, RHEL 9, ...).
# Behind an SSL-inspecting proxy:  PIP_OPTS="--trusted-host pypi.org --trusted-host files.pythonhosted.org"
# =====================================================================
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PY=python3
REF=""
SAWFISH=""
EXON_BED=""
LOGO=""
PACKAGE_ONLY=0
PUBLIC=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --package) PACKAGE_ONLY=1; shift ;;
        --public) PUBLIC=1; shift ;;
        --reference) REF="$2"; shift 2 ;;
        --sawfish) SAWFISH="$2"; shift 2 ;;
        --exon-bed) EXON_BED="$2"; shift 2 ;;
        --logo) LOGO="$2"; shift 2 ;;
        --python) PY="$2"; shift 2 ;;
        *) echo "unknown option $1" >&2; exit 2 ;;
    esac
done
SAWFISH="${SAWFISH:-$(command -v sawfish || true)}"
VERSION=$(sed -n "s/^__version__ = '\(.*\)'/\1/p" "$HERE/dmdseek/__init__.py")
B="$HERE/build"
OUT="$B/dist/dmdseek"
mkdir -p "$B"

package() {
    [[ -s "$OUT/reference/dmd_reference.npz" ]] || {
        echo "[NEXT] no reference in the bundle yet. Build it with the bundled executable, e.g.:"
        echo "       $OUT/dmdseek build-reference --exon-bed DMD_MANE_79exons.bed --wgs-bams wgs.tsv --target-bams tgt.tsv \\"
        echo "           --wgs-events ... --target-events ... -o $OUT/reference/dmd_reference.npz"
        echo "       then: bash $0 --package"
        exit 0; }
    echo "[smoke test]"
    "$OUT/dmdseek" --version
    "$OUT/dmdseek" info | sed 's/^/      /'
    if [[ $PUBLIC == 1 ]]; then
        # public release: no sawfish binary (own licence; users install it from bioconda), no logo,
        # no reference panel (users build one with build-reference, or add a released panel later)
        TAR="$HERE/dmdseek-${VERSION}-linux-x86_64-public.tar.gz"
        tar -czf "$TAR" -C "$B/dist" --exclude='dmdseek/bin/sawfish' --exclude='dmdseek/bin/sawfish-license' \
            --exclude='dmdseek/branding' --exclude='dmdseek/reference/*.npz' dmdseek
    else
        TAR="$HERE/dmdseek-${VERSION}-linux-x86_64.tar.gz"
        tar -czf "$TAR" -C "$B/dist" dmdseek
    fi
    echo "[DONE] $TAR ($(du -h "$TAR" | cut -f1)) - unpack anywhere and run ./dmdseek/dmdseek call ..."
}
if [[ $PACKAGE_ONLY == 1 ]]; then package; exit 0; fi

echo "[1/5] build environment (python: $($PY --version 2>&1))"
if [[ ! -x "$B/venv/bin/python" ]]; then
    "$PY" -m venv "$B/venv"
fi
"$B/venv/bin/python" -m pip install ${PIP_OPTS:-} -q --upgrade pip
"$B/venv/bin/python" -m pip install ${PIP_OPTS:-} -q "numpy>=1.24" "matplotlib>=3.7" "pysam>=0.22" "openpyxl>=3.1" "pyinstaller>=6.0"

echo "[2/5] PyInstaller"
KEEP="$B/previous_reference.npz"
if [[ -z "$REF" && -s "$OUT/reference/dmd_reference.npz" ]]; then
    cp "$OUT/reference/dmd_reference.npz" "$KEEP"       # rebuild keeps the reference panel already built
    REF="$KEEP"; echo "      keeping the reference panel of the previous build"
fi
rm -rf "$B/dist" "$B/work"
"$B/venv/bin/pyinstaller" --noconfirm --clean --onedir --name dmdseek \
    --distpath "$B/dist" --workpath "$B/work" --specpath "$B" \
    --collect-all pysam --hidden-import openpyxl --exclude-module tkinter --exclude-module IPython \
    --paths "$HERE" "$HERE/dmdseek_main.py" > "$B/pyinstaller.log" 2>&1 \
    || { echo "[ERROR] PyInstaller failed - see $B/pyinstaller.log" >&2; exit 1; }

echo "[3/5] sawfish"
mkdir -p "$OUT/bin" "$OUT/reference"
if [[ -n "$SAWFISH" && -x "$SAWFISH" ]]; then
    if command -v ldd >/dev/null && ldd "$SAWFISH" 2>/dev/null | grep -q "not found"; then
        echo "[ERROR] $SAWFISH needs libraries that are not on this system:" >&2
        ldd "$SAWFISH" | grep "not found" >&2; exit 1
    fi
    cp -L "$SAWFISH" "$OUT/bin/sawfish"
    chmod 755 "$OUT/bin/sawfish"
    # licence text: conda env share/, then the conda package cache (bioconda keeps info/licenses)
    pre="$(cd "$(dirname "$(readlink -f "$SAWFISH")")/.." && pwd)"
    for lic in "$pre/share/licenses/sawfish" $(ls -d "$pre"/../../pkgs/sawfish-*/info/licenses 2>/dev/null | sort -r) \
               $(ls -d "$pre"/pkgs/sawfish-*/info/licenses 2>/dev/null | sort -r); do
        if [[ -d "$lic" ]] && [[ -n "$(ls -A "$lic")" ]]; then cp -r "$lic" "$OUT/bin/sawfish-license"; break; fi
    done
    meta=$(ls "$pre"/conda-meta/sawfish-*.json 2>/dev/null | head -1 || true)
    [[ -n "$meta" ]] && python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print('      sawfish licence (conda):', d.get('license','?'))" "$meta" || true
    if [[ -d "$OUT/bin/sawfish-license" ]]; then echo "      licence text: bin/sawfish-license/"
    else echo "      [WARN] sawfish licence text not found - add it to bin/sawfish-license/ before distributing"; fi
    echo "      bundled: $("$OUT/bin/sawfish" --version 2>&1 | head -1)"
else
    echo "      [WARN] sawfish not found - the bundle will run depth-only unless sawfish is on PATH"
fi

echo "[4/5] reference + docs"
[[ -n "$REF" ]] && cp "$REF" "$OUT/reference/dmd_reference.npz"
for f in README.md README_KO.md LICENSE CITATION.cff; do [[ -f "$HERE/$f" ]] && cp "$HERE/$f" "$OUT/"; done
rm -rf "$OUT/templates" && cp -r "$HERE/templates" "$OUT/templates"
# licences of the bundled Python packages (numpy, matplotlib, pysam/htslib, ...)
rm -rf "$OUT/licenses"; mkdir -p "$OUT/licenses"
for d in "$B"/venv/lib/python*/site-packages/*.dist-info; do
    pkg=$(basename "$d" .dist-info)
    case "$pkg" in pip-*|setuptools-*|pyinstaller-*|pyinstaller_hooks_contrib-*|altgraph-*|wheel-*) continue ;; esac
    files=$(find "$d" \( -iname 'LICEN[CS]E*' -o -iname 'COPYING*' -o -iname 'NOTICE*' \) -type f 2>/dev/null || true)
    [[ -z "$files" ]] && continue
    mkdir -p "$OUT/licenses/$pkg"
    for f in $files; do cp "$f" "$OUT/licenses/$pkg/"; done
done
echo "      third-party licences: $(ls "$OUT/licenses" | wc -l) packages -> licenses/"
if [[ -n "$LOGO" ]]; then
    ext="${LOGO##*.}"; ext="${ext,,}"
    case "$ext" in png|svg|jpg|jpeg|webp) ;; *) echo "[ERROR] --logo must be png/svg/jpg/webp" >&2; exit 1 ;; esac
    rm -rf "$OUT/branding"; mkdir -p "$OUT/branding" && cp "$LOGO" "$OUT/branding/logo.$ext"
    echo "      logo: $LOGO -> branding/logo.$ext (embedded in every report)"
fi
if [[ -n "$EXON_BED" ]]; then
    [[ $(grep -cv -e '^#' -e '^[[:space:]]*$' "$EXON_BED") == 79 ]] || { echo "[ERROR] $EXON_BED: expected 79 exon lines" >&2; exit 1; }
    mkdir -p "$OUT/data" && cp "$EXON_BED" "$OUT/data/DMD_MANE_79exons.bed"
    echo "      exon BED: $EXON_BED -> data/DMD_MANE_79exons.bed"
else
    echo "      [WARN] no --exon-bed: build-reference will need --exon-bed"
fi
echo "$VERSION" > "$OUT/VERSION"

echo "[5/5] package"
package
