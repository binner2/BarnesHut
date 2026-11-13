#!/usr/bin/env bash
# Real-data measurements on a TUM RGB-D sequence (docs/GERCEK_VERI_PLANI.md, sections 0, 1.1, 3).
#   experiments/run_tum.sh <pointcloud_tool> <sequence_dir> <output_dir> [n_frames=20] [stride=2] [repeats=5]
# <sequence_dir> is the extracted rgbd_dataset_freiburgN_<name>/ folder (needs depth/*.png).
# Produces, in <output_dir>:
#   machine.txt                                 CPU, threads, compiler, commit
#   tum_<seq>_frames.csv                        benchframe rows (schema of experiments/README.md)
#   tum_<seq>_info.csv, tum_<seq>_info_stride<S>.csv   per-frame cloud statistics (full res / benchmark res)
#   profiles/<frame>.csv                        NN-distance vs depth profile per full-resolution frame
#   tum_<seq>_log.txt                           console log
# Resumable: frames already present in the frames CSV are skipped.
set -euo pipefail
TOOL=${1:-build/pointcloud/pointcloud_tool}
SEQ=${2:?sequence dir}
OUT=${3:-experiments/results/$(date +%F)}
NF=${4:-20}
STRIDE=${5:-2}
REPEATS=${6:-5}
NAME=$(basename "$SEQ" | sed 's/^rgbd_dataset_//')
CLOUDS=$(dirname "$SEQ")/${NAME}_clouds
mkdir -p "$OUT/profiles" "$CLOUDS/s1" "$CLOUDS/s$STRIDE"
{
  echo "date: $(date -Iseconds)"
  echo "commit: $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs || echo unknown)"
  echo "cores: $(nproc)"
  echo "threads: ${OMP_NUM_THREADS:-$(nproc)}"
  echo "compiler: $(${CXX:-g++} --version | head -1)"
  echo "flags: -O3 -march=native (CMake Release)"
  echo "sequence: $NAME"
  echo "frames: $NF equally spaced of $(ls "$SEQ/depth" | wc -l)"
  echo "benchmark stride: $STRIDE  repeats: $REPEATS"
  echo "intrinsics: fx=fy=525 cx=319.5 cy=239.5 scale=5000 (ROS default, TUM file_formats page)"
} > "$OUT/machine.txt"
LOG="$OUT/tum_${NAME}_log.txt"
FR="$OUT/tum_${NAME}_frames.csv"

# 1. All depth frames -> PLY (full resolution and benchmark stride) with the documented converter.
if [ "$(ls "$CLOUDS/s1" | wc -l)" -lt "$(ls "$SEQ/depth" | wc -l)" ]; then
  tools/tum_depth_to_ply.sh "$SEQ" "$CLOUDS/s1" "$TOOL" 1 | tee -a "$LOG"
fi
if [ "$(ls "$CLOUDS/s$STRIDE" | wc -l)" -lt "$(ls "$SEQ/depth" | wc -l)" ]; then
  tools/tum_depth_to_ply.sh "$SEQ" "$CLOUDS/s$STRIDE" "$TOOL" "$STRIDE" | tee -a "$LOG"
fi

# 2. Equally spaced frame selection (first and last included).
mapfile -t ALL < <(ls "$SEQ/depth" | sed 's/\.png$//' | sort)
TOTAL=${#ALL[@]}
SEL=()
for ((i = 0; i < NF; i++)); do
  idx=$(( (i * (TOTAL - 1) + (NF - 1) / 2) / (NF - 1) ))
  SEL+=("${ALL[$idx]}")
done
printf '%s\n' "${SEL[@]}" > "$OUT/tum_${NAME}_selected_frames.txt"

# 3. Per frame: statistics (full res + stride), depth profile, benchmark.
for f in "${SEL[@]}"; do
  if [ -f "$FR" ] && grep -q "frame=$f;" "$FR"; then
    echo "skip $f (already measured)" | tee -a "$LOG"; continue
  fi
  echo "### frame $f  $(date -Iseconds)" | tee -a "$LOG"
  "$TOOL" info "$CLOUDS/s1/$f.ply" --tag "$f" --csv "$OUT/tum_${NAME}_info.csv" \
          --depth-profile "$OUT/profiles/$f.csv" --bin 0.25 --fx 525 --stride 1 | tee -a "$LOG"
  "$TOOL" info "$CLOUDS/s$STRIDE/$f.ply" --tag "$f" --csv "$OUT/tum_${NAME}_info_stride$STRIDE.csv" | tee -a "$LOG"
  "$TOOL" benchframe "$CLOUDS/s$STRIDE/$f.ply" --tag "$f" --csv "$FR" --repeats "$REPEATS" \
          --hs 0.02,0.05 --thetas 0.5,0.7 --kappas 0,1 --cutoff 4 \
          --ms-h 0.08 --ms-thetas 0,0.5,0.7,1.0 --ms-iters 30 \
          --wlop-h 0.08 --wlop-thetas 0,0.5,0.7,1.0 --wlop-iters 10 --wlop-frac 0.1 | tee -a "$LOG"
done
echo "done $(date -Iseconds)" | tee -a "$LOG"
