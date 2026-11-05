#!/usr/bin/env bash
# Convert every depth PNG of a TUM RGB-D sequence to a PLY cloud.
#   tools/tum_depth_to_ply.sh <sequence_dir> <out_dir> [pointcloud_tool] [stride]
# TUM depth PNGs are 16-bit with scale 5000 (raw/5000 = metres); ROS default
# intrinsics fx=fy=525.0, cx=319.5, cy=239.5 (file_formats page, accessed 2026-11-04).
set -euo pipefail
SEQ=$1; OUT=$2; TOOL=${3:-build/pointcloud/pointcloud_tool}; STRIDE=${4:-1}
mkdir -p "$OUT"
for f in "$SEQ"/depth/*.png; do
  b=$(basename "$f" .png)
  "$TOOL" depth2cloud "$f" --fx 525 --fy 525 --cx 319.5 --cy 239.5 --scale 5000 --stride "$STRIDE" -o "$OUT/$b.ply" > /dev/null
done
echo "$(ls "$OUT" | wc -l) clouds in $OUT"
