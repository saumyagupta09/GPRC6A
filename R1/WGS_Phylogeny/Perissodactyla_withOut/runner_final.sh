#!/usr/bin/env bash
set -euo pipefail

# Start wall-clock timer
SECONDS=0

# (Optional) name a log file; comment out if you just want stdout
LOG_FILE="${LOG_FILE:-run_metrics.log}"

# Helper to format seconds -> HH:MM:SS
fmt_hms() {
  local s=$1
  printf '%02d:%02d:%02d' $((s/3600)) $(((s%3600)/60)) $((s%60))
}

# Helper to read bash `times` output and return total CPU seconds (shell + children)
# `times` prints two lines: "user system" for the shell, then for all children,
# with each value like "0m1.23s". We parse and sum them.
cpu_seconds_now() {
  ( times ) 2>&1 | awk '
    function t(x, m, s) {
      sub(/s$/, "", x)
      split(x, a, "m")
      m=a[1]+0
      s=a[2]+0
      return m*60 + s
    }
    NR==1 { su=$1; ss=$2 }
    NR==2 { cu=$1; cs=$2 }
    END   { print t(su)+t(ss)+t(cu)+t(cs) }'
}

# Capture CPU seconds at (near) start, so we can subtract if needed
CPU_START=$(cpu_seconds_now)



echo "Running actual code now."

clade="Perissodactyla"
python ../orthofinder_synteny_supermatrix_final.py --saturation --species-list "$clade".txt --outdir out_supermatrix_"$clade" --threads 8 --allow-missing --min-species 0.9 --min-aa-len 20 --min-col-occupancy 0.6 --max-gap-frac 0.5 --min-codon-len 60


cd out_supermatrix_"$clade"

iqtree2 -s supermatrix_codon.fna -p partitions.txt -m TESTMERGEONLY -rclusterf 5 --merge-model 1 --merge-rate 1 -T AUTO --prefix premerge

iqtree2 -s supermatrix_codon.fna -p premerge.best_scheme.nex -m MFP -T AUTO --prefix concat_fast  -B 1000 -bnni

ls per_gene_alignments_codon/*.aln.fna | parallel -j 70 'iqtree2 -s {} -m MFP -T 1 -B 1000 --prefix {.}'

cat per_gene_alignments_codon/*.treefile > gene_trees.best.tre

iqtree2 -t concat_fast.treefile -s supermatrix_codon.fna --gcf gene_trees.best.tre --scf 1000 --prefix cf_results --threads 30








# Wall-clock elapsed (in seconds)
ELAPSED_WALL_SEC=$SECONDS

# Total CPU seconds consumed by THIS script and its children
CPU_END=$(cpu_seconds_now)
CPU_SEC=$(( CPU_END - CPU_START ))

# Compute hours = total CPU core-seconds / 3600
# (e.g., a 4-core job for 30 minutes ˜ 2.0 compute-hours)
# Use bc if present for fractional hours; otherwise integer division fallback.
if command -v bc >/dev/null 2>&1; then
  COMPUTE_HOURS=$(echo "scale=4; $CPU_SEC / 3600" | bc)
else
  COMPUTE_HOURS="$(( CPU_SEC / 3600 ))"
fi

# Pretty print
WALL_HMS=$(fmt_hms "$ELAPSED_WALL_SEC")

{
  echo "---------------------------------------------"
  echo "Host                : $(hostname)"
  echo "Ended at            : $(date -Is)"
  echo "Wall-clock runtime  : $WALL_HMS (${ELAPSED_WALL_SEC}s)"
  echo "CPU time (seconds)  : ${CPU_SEC}"
  echo "Total compute hours : ${COMPUTE_HOURS}"
  echo "---------------------------------------------"
} | tee -a "${LOG_FILE:-/dev/null}"
