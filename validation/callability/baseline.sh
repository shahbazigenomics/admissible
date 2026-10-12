#!/usr/bin/env bash
# Independent baseline for ONE sample: callable bases in the target, computed with
# bedtools and awk only. It shares no code with admissible, which is the point.
#
#   baseline.sh COVERAGE TARGET OUTDIR [EXCLUDE]
#
# Writes into OUTDIR:
#   baseline.tsv   metric <TAB> value        (target_bp, callable_in_target_bp,
#                                             callable_in_target_minus_exclude_bp)
#   regions.tsv    region <TAB> bp <TAB> callable_bp <TAB> callable_minus_exclude_bp
#
# NOT RUN by the author: this environment had neither bedtools nor mosdepth.
# Exit codes: 0 done, 3 a required tool is missing, 4 an input is unusable.
set -u
export LC_ALL=C   # sort and join must agree on collation

cov=${1:-}; target=${2:-}; out=${3:-}; excl=${4:-}
if [ -z "$cov" ] || [ -z "$target" ] || [ -z "$out" ]; then
  echo "usage: baseline.sh COVERAGE TARGET OUTDIR [EXCLUDE]" >&2; exit 4
fi
for tool in bedtools awk sort; do
  command -v "$tool" >/dev/null 2>&1 || { echo "baseline.sh: $tool not found" >&2; exit 3; }
done
for f in "$cov" "$target"; do
  [ -r "$f" ] || { echo "baseline.sh: cannot read $f" >&2; exit 4; }
done
[ -z "$excl" ] || [ -r "$excl" ] || { echo "baseline.sh: cannot read $excl" >&2; exit 4; }
mkdir -p "$out" || exit 4

# Contig names are normalised the way admissible does it (strip "chr", M -> MT), so
# the two sides are compared on the same footing. This is the one assumption the
# baseline shares with admissible; it is stated in the README.
NORM='function norm(c){sub(/^[Cc][Hh][Rr]/,"",c); u=toupper(c);
  if(u=="M"||u=="MT")return "MT"; if(u=="X"||u=="Y")return u; return c}'

reader() { case "$1" in *.gz|*.bgz) gzip -dc "$1" ;; *) cat "$1" ;; esac; }

# callable = rows labelled 10:inf by mosdepth --quantize 0:10:
reader "$cov" | awk -F'\t' -v OFS='\t' "$NORM"'
  /^#|^track|^browser/ || NF<4 {next} $4=="10:inf" && $3>$2 {print norm($1),$2,$3}' \
  | sort -k1,1 -k2,2n | bedtools merge -i - > "$out/callable.bed" || exit 4

# the target with names (one line per BED line) and its merged union
reader "$target" | awk -F'\t' -v OFS='\t' "$NORM"'
  /^#|^track|^browser/ || NF<3 {next} $2>=0 && $3>$2 {
    n=(NF>=4 && $4!="")?$4:(norm($1)":"$2"-"$3); print norm($1),$2,$3,n}' \
  | sort -k1,1 -k2,2n > "$out/target.named.bed" || exit 4
cut -f1-3 "$out/target.named.bed" | bedtools merge -i - > "$out/target.merged.bed" || exit 4

if [ -n "$excl" ]; then
  reader "$excl" | awk -F'\t' -v OFS='\t' "$NORM"'
    /^#|^track|^browser/ || NF<3 {next} $2>=0 && $3>$2 {print norm($1),$2,$3}' \
    | sort -k1,1 -k2,2n | bedtools merge -i - > "$out/exclude.bed" || exit 4
else
  : > "$out/exclude.bed"
fi

bp() { awk '{s+=$3-$2} END{printf "%d", s+0}' "$1"; }

bedtools intersect -a "$out/target.merged.bed" -b "$out/callable.bed" \
  > "$out/callable_in_target.bed" || exit 4
bedtools subtract -a "$out/callable_in_target.bed" -b "$out/exclude.bed" \
  > "$out/callable_in_target_minus_exclude.bed" || exit 4

{
  printf 'target_bp\t%s\n' "$(bp "$out/target.merged.bed")"
  printf 'callable_in_target_bp\t%s\n' "$(bp "$out/callable_in_target.bed")"
  printf 'callable_in_target_minus_exclude_bp\t%s\n' \
    "$(bp "$out/callable_in_target_minus_exclude.bed")"
} > "$out/baseline.tsv"

# Per region. A region is every line sharing a name, as in admissible; the two agree
# only if lines of one name do not overlap each other (compare.py reports it if not).
cut -f1-4 "$out/target.named.bed" > "$out/regions.in.bed"
bedtools intersect -a "$out/regions.in.bed" -b "$out/callable.bed" -wao \
  | awk -F'\t' '{c[$4]+=$NF} END{for(k in c) printf "%s\t%d\n", k, c[k]}' | sort > "$out/regions.callable.tsv"
bedtools subtract -a "$out/callable.bed" -b "$out/exclude.bed" > "$out/callable.x.bed" || exit 4
bedtools intersect -a "$out/regions.in.bed" -b "$out/callable.x.bed" -wao \
  | awk -F'\t' '{c[$4]+=$NF} END{for(k in c) printf "%s\t%d\n", k, c[k]}' \
  | sort > "$out/regions.callable_x.tsv"
awk -F'\t' '{b[$4]+=$3-$2} END{for(k in b) printf "%s\t%d\n", k, b[k]}' \
  "$out/regions.in.bed" | sort > "$out/regions.bp.tsv"
join -t "$(printf '\t')" "$out/regions.bp.tsv" "$out/regions.callable.tsv" \
  | join -t "$(printf '\t')" - "$out/regions.callable_x.tsv" > "$out/regions.tsv" || exit 4
exit 0
