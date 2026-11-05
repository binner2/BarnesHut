#!/usr/bin/env python3
"""Capture depth frames from an Intel RealSense camera as 16-bit PNGs plus an
intrinsics file readable by pointcloud_tool (--intrinsics).

    python3 tools/realsense_capture.py out_dir --frames 30 --width 640 --height 480

Requires pyrealsense2, numpy, Pillow. NOT tested in the authoring session (no
camera). Depth scale is read from the device and written to the intrinsics file
as raw units per metre, so --scale on the command line is not needed."""
import argparse, os, time
import numpy as np
from PIL import Image
import pyrealsense2 as rs

ap = argparse.ArgumentParser()
ap.add_argument("out"); ap.add_argument("--frames", type=int, default=30)
ap.add_argument("--width", type=int, default=640); ap.add_argument("--height", type=int, default=480)
ap.add_argument("--fps", type=int, default=30); ap.add_argument("--warmup", type=int, default=30)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
pipe = rs.pipeline(); cfg = rs.config()
cfg.enable_stream(rs.stream.depth, a.width, a.height, rs.format.z16, a.fps)
profile = pipe.start(cfg)
sensor = profile.get_device().first_depth_sensor()
depth_scale_m = sensor.get_depth_scale()                 # metres per raw unit
intr = profile.get_stream(rs.stream.depth).as_video_stream_profile().get_intrinsics()
with open(os.path.join(a.out, "intrinsics.txt"), "w") as f:
    f.write("# fx fy cx cy width height depth_scale(raw units per metre)\n")
    f.write(f"{intr.fx} {intr.fy} {intr.ppx} {intr.ppy} {intr.width} {intr.height} {1.0/depth_scale_m}\n")
try:
    for _ in range(a.warmup): pipe.wait_for_frames()
    for i in range(a.frames):
        fr = pipe.wait_for_frames().get_depth_frame()
        img = np.asanyarray(fr.get_data()).astype(np.uint16)
        Image.fromarray(img).save(os.path.join(a.out, f"{time.time():.6f}.png"))
finally:
    pipe.stop()
print(f"wrote {a.frames} frames and intrinsics.txt to {a.out}")
