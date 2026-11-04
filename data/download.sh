#!/bin/bash
# Download (parts of) the simulation outputs from the public Nextcloud share and verify their sha256 checksums.
#
#   bash data/download.sh TARGET [TARGET ...]          download into $DATA_DIR (default: data/sim)
#   bash data/download.sh -n TARGET ...                dry run: list the files and the total size only
#
# TARGET is one of the names below, or a glob pattern over the paths in data/MANIFEST.tsv
# (quote it), e.g. 'J_sweep/bq/6.?/freq.npy' or 'D_sweep/no_bq_5.5/320/*'.
#
#   fig1 fig3 ... figS1    what one figure needs (fig1, figS3 and Tables S2/S4 need nothing: they use data/derived/)
#     fig2                 J_sweep/{bq/6.0,no_bq/5.5}/data.h5                                     2.8 GB
#     fig3                 .npy of J_sweep/bq/6.0, J_sweep/no_bq/5.5                              4 MB
#     fig4                 .npy of J_sweep/bq/6.0, 6.1, J_sweep/no_bq/5.5, 5.6                    6 MB
#     fig5                 .npy of J_sweep/bq/12.0, J_sweep/no_bq/12.0                            3 MB
#     figS1                freq.npy + fft_avg.npy of all 122 diameter-sweep runs                  2 MB
#   figures                all of the above                                                       2.8 GB
#   ref                    the two reference data.h5 (= fig2; also launch_front.py, frame_check.py)
#   npy                    every .npy file of every run (wave_extract.py, j_vs_freq.py)           0.8 GB
#   h5-J, h5-D             every data.h5 of the current / diameter sweep (power_balance_J.py / power_balance.py)
#                                                                                         381 GB / 223 GB
#   all                    everything in the manifest                                             604 GB
#
# Share layout (same as $DATA_DIR):  J_sweep/{bq,no_bq}/<J>/   D_sweep/{bq,no_bq_5.5}/<D>/
# Files already present with the right checksum are skipped, so an interrupted download can simply be restarted.
# Requires bash >= 4, curl, sha256sum.
set -euo pipefail

SHARE_URL="${SHARE_URL:-https://genkai-nextcloud.hpc.kyushu-u.ac.jp/s/iHs5kSBoLckLR6a}"
HOST="${SHARE_URL%%/s/*}"; TOKEN="${SHARE_URL##*/s/}"; TOKEN="${TOKEN%%/*}"
WEBDAV="$HOST/public.php/webdav"            # public-share WebDAV: user = share token, empty password

HERE="$(cd "$(dirname "$0")" && pwd)"
DEST="${DATA_DIR:-$HERE/sim}"
TSV="$HERE/MANIFEST.tsv"                    # path  bytes  sha256

NPY_MAPS="final heat_map freq fft"          # files read by load_maps() in figures/common.py
declare -A SETS=(
  [fig1]="" [figS3]="" [tables]=""
  [fig2]="J_sweep/bq/6.0/data.h5 J_sweep/no_bq/5.5/data.h5"
  [ref]="J_sweep/bq/6.0/data.h5 J_sweep/no_bq/5.5/data.h5"
  [figS1]="D_sweep/*/freq.npy D_sweep/*/fft_avg.npy"
  [npy]="*.npy"
  [h5-J]="J_sweep/*/data.h5"
  [h5-D]="D_sweep/*/data.h5"
  [all]="*"
)
maps() { local r f out=""; for r in "$@"; do for f in $NPY_MAPS; do out+=" J_sweep/$r/$f.npy"; done; done; echo "$out"; }
SETS[fig3]="$(maps bq/6.0 no_bq/5.5)"
SETS[fig4]="$(maps bq/6.0 bq/6.1 no_bq/5.5 no_bq/5.6)"
SETS[fig5]="$(maps bq/12.0 no_bq/12.0)"
SETS[figures]="${SETS[fig2]} ${SETS[fig3]} ${SETS[fig4]} ${SETS[fig5]} ${SETS[figS1]}"

usage() { sed -n '2,25p' "$0"; exit 1; }
DRY=0
[ "${1:-}" = "-n" ] && { DRY=1; shift; }
[ $# -ge 1 ] || usage

# expand targets into glob patterns
pats=()
for t in "$@"; do
    if [[ -v SETS[$t] ]]; then
        [ -z "${SETS[$t]}" ] && echo "$t: needs no download (uses data/derived/)"
        read -ra p <<< "${SETS[$t]}"; pats+=("${p[@]}")
    elif [[ "$t" == *[/*?]* ]]; then pats+=("$t")
    else echo "unknown target '$t'" >&2; usage; fi
done
[ ${#pats[@]} -gt 0 ] || exit 0

# select manifest entries matching any pattern
paths=(); sizes=(); sums=(); total=0
while IFS=$'\t' read -r p b s; do
    [ "$p" = path ] && continue
    for pat in "${pats[@]}"; do
        # shellcheck disable=SC2053
        if [[ "$p" == $pat ]]; then paths+=("$p"); sizes+=("$b"); sums+=("$s"); total=$((total + b)); break; fi
    done
done < "$TSV"
[ ${#paths[@]} -gt 0 ] || { echo "no file in data/MANIFEST.tsv matches: ${pats[*]}" >&2; exit 1; }
human() { awk -v b="$1" 'BEGIN{split("B kB MB GB TB",u); i=1; while (b>=1000 && i<5) {b/=1000; i++}; printf "%.1f %s", b, u[i]}'; }
echo "${#paths[@]} files, $(human $total) -> $DEST"
if [ $DRY = 1 ]; then printf '%s\n' "${paths[@]}"; exit 0; fi

ok=0; skip=0; bad=0
for i in "${!paths[@]}"; do
    p="${paths[$i]}"; out="$DEST/$p"; want="${sums[$i]}"
    if [ -f "$out" ] && [ "$(stat -c %s "$out")" = "${sizes[$i]}" ] && \
       [ "$(sha256sum "$out" | cut -d' ' -f1)" = "$want" ]; then skip=$((skip + 1)); continue; fi
    # keep a shorter file as a partial download (resumed with -C -); a complete but wrong file is removed
    [ -f "$out" ] && [ "$(stat -c %s "$out")" -ge "${sizes[$i]}" ] && rm -f "$out"
    mkdir -p "$(dirname "$out")"
    printf '[%d/%d] %s (%s)\n' $((i + 1)) ${#paths[@]} "$p" "$(human "${sizes[$i]}")"
    if ! curl -fsSL --retry 5 --retry-delay 10 -C - -u "$TOKEN:" -o "$out" "$WEBDAV/$p"; then
        # resume may fail if the server ignores ranges: retry from scratch once
        rm -f "$out"; curl -fsSL --retry 5 --retry-delay 10 -u "$TOKEN:" -o "$out" "$WEBDAV/$p"
    fi
    if [ "$(sha256sum "$out" | cut -d' ' -f1)" = "$want" ]; then ok=$((ok + 1))
    else echo "  CHECKSUM MISMATCH: $p (file removed)" >&2; rm -f "$out"; bad=$((bad + 1)); fi
done
echo "downloaded $ok, already present $skip, failed $bad  (data in $DEST)"
[ $bad = 0 ]
