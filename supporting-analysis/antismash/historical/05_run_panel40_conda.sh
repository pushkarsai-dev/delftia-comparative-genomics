#!/bin/bash
# ============================================================================
#  antiSMASH 8.0.4 on the frozen 40-genome Delftia panel  (conda route)
#
#  Run 04_setup_antismash_conda.sh first.
#
#  Usage:   bash "05_run_panel40_conda.sh"
#
#  Safe to stop (Ctrl-C) and re-run: finished genomes are skipped.
# ============================================================================

set -u
ENVNAME="antismash8"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IN="$HERE/panel40_fasta"
OUT="$HERE/results"
LOGS="$HERE/logs"

CPUS="${CPUS:-$(( $(sysctl -n hw.ncpu 2>/dev/null || echo 4) - 1 ))}"
[ "$CPUS" -lt 1 ] && CPUS=1

# Identical to the Docker route -- see README.md for the rationale.
OPTS=(
  --genefinding-tool prodigal
  --hmmdetection-strictness relaxed
  --cc-mibig
  --cb-knownclusters
  --cb-subclusters
  --asf
  --cpus "$CPUS"
)

say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
fail() { printf '\n\033[31mERROR: %s\033[0m\n\n' "$*" >&2; exit 1; }

# ------------------------------------------------------------------ conda
# Must verify conda actually runs, not merely that conda.sh exists.
CONDA_SH=""
for c in "/opt/miniforge3" "$HOME/miniforge3" "/opt/homebrew/Caskroom/miniforge/base" \
         "$HOME/anaconda3" "$HOME/opt/anaconda3" "/opt/anaconda3" \
         "$HOME/miniconda3" "$HOME/opt/miniconda3" "/opt/miniconda3"; do
  [ -f "$c/etc/profile.d/conda.sh" ] || continue
  if [ -x "$c/bin/conda" ] && "$c/bin/conda" --version >/dev/null 2>&1; then
    CONDA_SH="$c/etc/profile.d/conda.sh"; break
  fi
done
[ -n "$CONDA_SH" ] || fail "No working conda found. Run 07_install_miniforge.sh, then 04_setup_antismash_conda.sh."
# shellcheck disable=SC1090
source "$CONDA_SH"

# Activate by full prefix rather than by name. The Anaconda uninstall took
# ~/.conda with it, so conda's environments.txt registry may be missing and
# name lookup can fail even though the environment is perfectly fine.
CONDA_BASE="$(cd "$(dirname "$(dirname "$CONDA_SH")")/.." && pwd)"
ENVPREFIX="$CONDA_BASE/envs/$ENVNAME"

if [ -d "$ENVPREFIX" ]; then
  conda activate "$ENVPREFIX" || fail "Could not activate $ENVPREFIX"
else
  conda activate "$ENVNAME" || fail \
"Environment '$ENVNAME' not found (looked for $ENVPREFIX).
 Run 04_setup_antismash_conda.sh first."
fi
command -v antismash >/dev/null 2>&1 || fail "antismash not on PATH inside '$ENVNAME'."

[ -d "$IN" ] || fail "Input folder not found: $IN"
mkdir -p "$OUT" "$LOGS"

# macOS ships bash 3.2, which has no `mapfile`.
FILES=()
while IFS= read -r _line; do
  [ -n "$_line" ] && FILES+=("$_line")
done < <(ls "$IN"/*.fna 2>/dev/null | sort)
N=${#FILES[@]}
[ "$N" -gt 0 ] || fail "No .fna files in $IN"

say "antiSMASH $(antismash --version 2>&1 | head -1) — ${N} genomes, ${CPUS} cores each"
echo "    input:   $IN"
echo "    output:  $OUT"
echo "    options: ${OPTS[*]}"
echo
echo "    Expect roughly 5-15 min per genome (longer under Rosetta)."
echo "    You can Ctrl-C and re-run to resume."

# --------------------------------------------------------- smoke test first
SMOKE="$IN/BAA554_hybrid.fna"
if [ -f "$SMOKE" ] && [ ! -f "$OUT/BAA554_hybrid/index.html" ]; then
  say "Smoke test on BAA554_hybrid.fna before running the rest"
  antismash "$SMOKE" --output-dir "$OUT/BAA554_hybrid" \
      --output-basename BAA554_hybrid "${OPTS[@]}" 2>&1 \
      | tee "$LOGS/BAA554_hybrid.log"
  if [ ! -f "$OUT/BAA554_hybrid/index.html" ]; then
    fail "Smoke test failed. Read $LOGS/BAA554_hybrid.log before going further.
 If it reports an unrecognised option, check that option against the installed version."
  fi
  say "Smoke test passed. Continuing with the remaining genomes."
fi

# ------------------------------------------------------------------- main
i=0; done_n=0; skip_n=0; fail_n=0; FAILED=()
START=$(date +%s)

for f in "${FILES[@]}"; do
  i=$((i+1))
  base="$(basename "$f" .fna)"

  if [ -f "$OUT/$base/index.html" ]; then
    printf '[%2d/%d] %-34s already done, skipping\n' "$i" "$N" "$base"
    skip_n=$((skip_n+1)); continue
  fi

  printf '[%2d/%d] %-34s ' "$i" "$N" "$base"
  t0=$(date +%s)
  antismash "$f" --output-dir "$OUT/$base" --output-basename "$base" \
      "${OPTS[@]}" > "$LOGS/$base.log" 2>&1
  t1=$(date +%s)

  if [ -f "$OUT/$base/index.html" ]; then
    printf 'ok (%dm%02ds)\n' $(( (t1-t0)/60 )) $(( (t1-t0)%60 ))
    done_n=$((done_n+1))
  else
    printf 'FAILED  -- see logs/%s.log\n' "$base"
    fail_n=$((fail_n+1)); FAILED+=("$base")
  fi
done

END=$(date +%s)

# ---------------------------------------------------------------- summary
say "Building region-count summary"
SUM="$HERE/antismash_region_counts.tsv"
{
  printf 'genome\tregions\tstatus\n'
  for f in "${FILES[@]}"; do
    base="$(basename "$f" .fna)"
    if [ -f "$OUT/$base/index.html" ]; then
      n=$(ls "$OUT/$base"/*.region*.gbk 2>/dev/null | wc -l | tr -d ' ')
      printf '%s\t%s\tok\n' "$base" "$n"
    else
      printf '%s\t\tFAILED\n' "$base"
    fi
  done
} > "$SUM"

say "Finished"
cat <<EOF

  completed this run : $done_n
  already done       : $skip_n
  failed             : $fail_n
  wall clock         : $(( (END-START)/3600 ))h $(( ((END-START)%3600)/60 ))m

  Region counts      : antismash_region_counts.tsv
  Per-genome reports : results/<genome>/index.html   (open in a browser)
  Logs               : logs/

EOF

if [ "$fail_n" -gt 0 ] && [ ${#FAILED[@]} -gt 0 ]; then
  echo "  Failed genomes:"
  printf '    %s\n' "${FAILED[@]}"
  echo
fi

column -t -s $'\t' "$SUM" 2>/dev/null || cat "$SUM"
