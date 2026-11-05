#!/usr/bin/env bash
# Reproduce the synthetic-scene benchmarks used in the report.
#   experiments/run_benchmarks.sh <pointcloud_tool> <output_dir>
set -euo pipefail
TOOL=${1:-build/pointcloud/pointcloud_tool}
OUT=${2:-experiments/results/$(date +%F)}
mkdir -p "$OUT"
{
  echo "date: $(date -Iseconds)"
  echo "commit: $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs || echo unknown)"
  echo "threads: ${OMP_NUM_THREADS:-$(nproc)}"
  echo "compiler: $(${CXX:-g++} --version | head -1)"
} > "$OUT/machine.txt"

# Main sweep: N x theta x kappa, Gaussian KDE h = 5 cm, leaf 8, plus mean shift and WLOP.
"$TOOL" bench --sizes 10000,20000,50000,100000 --thetas 0,0.3,0.5,0.7,1.0 --kappas 0,1,0.5 \
        --h 0.05 --leaf 8 --repeats 2 --csv "$OUT/bench_main.csv" | tee "$OUT/bench_main_log.txt"

# Leaf-capacity sensitivity at N = 50k.
for L in 4 16 32; do
  "$TOOL" bench --sizes 50000 --thetas 0,0.5,0.7 --kappas 0,1 --h 0.05 --leaf "$L" --repeats 2 \
          --csv "$OUT/bench_leaf$L.csv" | tee "$OUT/bench_leaf${L}_log.txt"
done

# End-to-end demo on the full synthetic frame (timings quoted in the report).
D="$OUT/demo"; mkdir -p "$D"
"$TOOL" synthetic -o "$D/synthetic_labels.ply" --pgm "$D/synthetic_depth.pgm" --labels
"$TOOL" depth2cloud "$D/synthetic_depth.pgm" --scale 1000 -o "$D/frame.ply"
{
  "$TOOL" downsample "$D/frame.ply" --voxel 0.01 -o "$D/ds_1cm.ply"
  "$TOOL" sor "$D/frame.ply" --k 20 --std 2.0 -o "$D/sor.ply"
  "$TOOL" normals "$D/sor.ply" --k 20 --viewpoint 0,0,0 --curvature -o "$D/normals.ply"
  "$TOOL" density "$D/sor.ply" --h 0.02 --theta 0.5 -o "$D/density_full.ply"
  "$TOOL" density "$D/ds_1cm.ply" --h 0.03 --theta 0.5 --compare -o "$D/density_ds.ply"
  "$TOOL" density "$D/ds_1cm.ply" --h 0.03 --theta 0.5 --kappa 1 --compare -o "$D/density_ds_k1.ply"
  "$TOOL" meanshift "$D/ds_1cm.ply" --h 0.08 --theta 0.5 --min-size 100 -o "$D/meanshift_ds.ply"
  "$TOOL" resample "$D/ds_1cm.ply" --count 4000 --h 0.08 --iters 15 --theta 0.5 -o "$D/wlop_ds.ply"
} | tee "$OUT/demo_log.txt"
echo "results in $OUT"
