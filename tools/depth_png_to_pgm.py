#!/usr/bin/env python3
"""Convert a 16-bit depth PNG (Kinect/RealSense/TUM RGB-D style) to a 16-bit PGM
that pointcloud_tool reads without libpng, or directly to an XYZ/PLY cloud.

    python3 tools/depth_png_to_pgm.py depth.png depth.pgm
    python3 tools/depth_png_to_pgm.py depth.png cloud.ply --fx 525 --fy 525 --cx 319.5 --cy 239.5 --scale 5000

Only numpy and Pillow are required.
"""
import argparse
import struct
import sys

import numpy as np
from PIL import Image


def load_depth(path):
    img = Image.open(path)
    arr = np.array(img)
    if arr.ndim == 3:
        arr = arr[..., 0]
    if arr.dtype != np.uint16:
        arr = arr.astype(np.uint16)
    return arr


def write_pgm(arr, path):
    h, w = arr.shape
    with open(path, "wb") as f:
        f.write(f"P5\n{w} {h}\n65535\n".encode())
        f.write(arr.astype(">u2").tobytes())


def back_project(arr, fx, fy, cx, cy, scale, zmin, zmax, stride):
    v, u = np.mgrid[0:arr.shape[0]:stride, 0:arr.shape[1]:stride]
    z = arr[::stride, ::stride].astype(np.float64) / scale
    ok = (arr[::stride, ::stride] > 0) & (z >= zmin) & (z <= zmax)
    x = (u - cx) * z / fx
    y = (v - cy) * z / fy
    return np.stack([x[ok], y[ok], z[ok]], axis=1)


def write_ply(pts, path):
    with open(path, "wb") as f:
        f.write(("ply\nformat binary_little_endian 1.0\n"
                 f"element vertex {len(pts)}\n"
                 "property float x\nproperty float y\nproperty float z\nend_header\n").encode())
        f.write(pts.astype("<f4").tobytes())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("output", help=".pgm, .ply or .xyz")
    ap.add_argument("--fx", type=float, default=525.0)
    ap.add_argument("--fy", type=float, default=525.0)
    ap.add_argument("--cx", type=float, default=319.5)
    ap.add_argument("--cy", type=float, default=239.5)
    ap.add_argument("--scale", type=float, default=1000.0, help="raw units per metre (TUM PNG: 5000)")
    ap.add_argument("--min", type=float, default=0.05)
    ap.add_argument("--max", type=float, default=10.0)
    ap.add_argument("--stride", type=int, default=1)
    a = ap.parse_args()

    arr = load_depth(a.input)
    out = a.output.lower()
    if out.endswith(".pgm"):
        write_pgm(arr, a.output)
    else:
        pts = back_project(arr, a.fx, a.fy, a.cx, a.cy, a.scale, a.min, a.max, a.stride)
        if out.endswith(".ply"):
            write_ply(pts, a.output)
        else:
            np.savetxt(a.output, pts, fmt="%.6f")
        print(f"{len(pts)} points", file=sys.stderr)
    print(f"wrote {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
