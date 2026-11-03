#!/usr/bin/env bash
# ============================================================================
#  Re-run antiSMASH 8 on the two panel genomes that failed on 2 Aug 2026.
#
#  WHY THEY FAILED
#    antiSMASH calls Prodigal 2.6.3 in metagenomic mode. Prodigal crashed hard
#    (no error message, truncated log) on one low-GC contig in each genome:
#      GCA_029835695.1  JAOEAI010000069.1   5,326 bp, 51.7% GC
#      GCA_037198095.1  JBBBVL010000006.1  11,762 bp, 51.4% GC
#    Both are clean ACGT with no ambiguity codes, against a genome average of
#    ~66.5% GC. The input is fine; this is the Prodigal 2.6.3 meta-mode crash.
#
#  THE FIX
#    Skip antiSMASH's gene finding and supply CDS calls made with pyrodigal
#    3.7.1, a reimplementation of Prodigal 2.6.3 without the crash. It called
#    5,965 and 5,865 CDS on these genomes without error, using per-genome
#    training (meta=False) -- the same settings as the panel-wide analyses in
#    methods sections 48 to 50.
#
#  Flags below are identical to 05_run_panel40_conda.sh except that
#  --genefinding-tool prodigal is replaced by --genefinding-gff3, so the
#  results remain comparable with the other 38 genomes.
#
#  RUN ON THE MAC, where the antismash8 conda env already exists.
#    bash run_rerun_two.sh
# ============================================================================

set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

CPUS="${CPUS:-$(( $(sysctl -n hw.ncpu 2>/dev/null || echo 4) - 1 ))}"
[ "$CPUS" -lt 1 ] && CPUS=1

for C in "$HOME/miniforge3/etc/profile.d/conda.sh" \
         "/opt/miniforge3/etc/profile.d/conda.sh" \
         "$HOME/miniconda3/etc/profile.d/conda.sh"; do
  [ -f "$C" ] && { source "$C"; break; }
done
conda activate antismash8 || { echo "could not activate antismash8"; exit 1; }
antismash --version

for G in GCA_029835695.1 GCA_037198095.1; do
  echo "=== $G"
  antismash \
    --genefinding-gff3 "${G}.gff" \
    --hmmdetection-strictness relaxed \
    --cc-mibig \
    --cb-knownclusters \
    --cb-subclusters \
    --asf \
    --cpus "$CPUS" \
    --output-dir "out_${G}" \
    "${G}.fna" 2>&1 | tee "${G}.rerun.log"
done

echo
echo "=== regions found ==="
for G in GCA_029835695.1 GCA_037198095.1; do
  printf "%s\t%s regions\n" "$G" "$(ls out_${G}/*region*.gbk 2>/dev/null | wc -l | tr -d ' ')"
done

echo
echo "=== delftibactin among KnownClusterBlast hits? ==="
for G in GCA_029835695.1 GCA_037198095.1; do
  if grep -rqls delftibactin "out_${G}" 2>/dev/null; then
    echo "$G  YES"
    grep -rhs -m1 -A3 "Significant hits" "out_${G}"/knownclusterblast/*.txt 2>/dev/null | head -4
  else
    echo "$G  no"
  fi
done
