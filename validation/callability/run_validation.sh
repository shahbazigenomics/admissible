#!/usr/bin/env bash
# Run admissible's callability check and the independent bedtools baseline for every
# row of a kits table, and compare them.
#
#   run_validation.sh KITS.tsv OUTDIR
#
# NOT RUN by the author: this environment had neither bedtools nor mosdepth, and no
# capture-kit data. Nothing in 0.2.0 has been validated by this script yet.
# Exit codes: 0 every sample agrees and >= 2 kits were covered; 1 a disagreement;
# 2 INCOMPLETE (fewer than 2 kits, or a row could not be run); 3 a tool is missing.
set -u
kits=${1:-}; outdir=${2:-}
if [ -z "$kits" ] || [ -z "$outdir" ] || [ ! -r "$kits" ]; then
  echo "usage: run_validation.sh KITS.tsv OUTDIR" >&2; exit 2
fi
here=$(cd "$(dirname "$0")" && pwd)
for tool in admissible python3 bedtools; do
  command -v "$tool" >/dev/null 2>&1 || { echo "run_validation.sh: $tool not found" >&2; exit 3; }
done
mkdir -p "$outdir" || exit 2
{ echo "admissible: $(admissible --version 2>&1)"
  echo "bedtools:   $(bedtools --version 2>&1)"
  echo "date:       $(date -u +%FT%TZ)"; } > "$outdir/versions.txt"

n_rows=0; n_ok=0; n_diff=0; n_fail=0; kit_list=""
while IFS=$'\t' read -r kit sample cov target excl; do
  case "$kit" in ''|\#*) continue ;; esac
  n_rows=$((n_rows + 1))
  d="$outdir/$kit/$sample"; mkdir -p "$d"
  printf 'FAM\t%s\t0\t0\t0\t2\n' "$sample" > "$d/one.ped"

  if ! "$here/baseline.sh" "$cov" "$target" "$d/baseline" ${excl:+"$excl"}; then
    echo "$kit/$sample: baseline failed" >&2; n_fail=$((n_fail + 1)); continue
  fi
  args=(callability --ped "$d/one.ped" --coverage "$sample=$cov" --target "$target"
        --by-region --json "$d/admissible.json")
  flag=()
  if [ -n "${excl:-}" ]; then args+=(--exclude-bed "$excl"); flag=(--with-exclusion); fi
  # admissible exits non-zero for FAIL/UNKNOWN verdicts, which is not an error here
  admissible "${args[@]}" > "$d/admissible.txt" 2>&1
  python3 "$here/compare.py" "$d/baseline" "$d/admissible.json" ${flag[@]+"${flag[@]}"} \
      --out "$d/result.tsv" > "$d/compare.txt"
  rc=$?
  cat "$d/compare.txt" | sed "s|^|$kit/$sample: |"
  case $rc in 0) n_ok=$((n_ok + 1)) ;; 1) n_diff=$((n_diff + 1)) ;; *) n_fail=$((n_fail + 1)) ;; esac
  case " $kit_list " in *" $kit "*) ;; *) kit_list="$kit_list $kit" ;; esac
done < "$kits"

# shellcheck disable=SC2086  # kit_list is deliberately split into one kit per line
n_kits=$(printf '%s\n' $kit_list | grep -c .)
echo "rows: $n_rows  identical: $n_ok  different: $n_diff  could-not-run: $n_fail  kits:$kit_list"
[ "$n_diff" -eq 0 ] || exit 1
if [ "$n_kits" -lt 2 ] || [ "$n_fail" -gt 0 ] || [ "$n_ok" -eq 0 ]; then
  echo "INCOMPLETE: needs at least two kits and every row to run" >&2; exit 2
fi
echo "all compared numbers identical across $n_kits kits"
exit 0
