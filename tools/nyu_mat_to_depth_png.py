#!/usr/bin/env python3
"""Export NYU Depth v2 labeled .mat (MATLAB v7.3 / HDF5) to 16-bit depth PNGs
(millimetres, scale 1000) and 16-bit label PNGs.

    python3 tools/nyu_mat_to_depth_png.py nyu_depth_v2_labeled.mat data/nyu [--limit N]

The dataset page states depth values are in metres. Requires numpy, h5py, Pillow.
NOT tested in the authoring session (dataset not downloaded); verify array
orientation on the first frame (HDF5 stores MATLAB arrays transposed)."""
import argparse, os
import numpy as np
import h5py
from PIL import Image

ap = argparse.ArgumentParser()
ap.add_argument("mat"); ap.add_argument("out"); ap.add_argument("--limit", type=int, default=0)
a = ap.parse_args()
os.makedirs(os.path.join(a.out, "depth"), exist_ok=True)
os.makedirs(os.path.join(a.out, "label"), exist_ok=True)
with h5py.File(a.mat, "r") as f:
    depths, labels = f["depths"], f["labels"]      # stored as N x W x H
    n = depths.shape[0] if a.limit == 0 else min(a.limit, depths.shape[0])
    for i in range(n):
        d = np.asarray(depths[i]).T                 # -> H x W, metres
        l = np.asarray(labels[i]).T
        Image.fromarray(np.clip(np.round(d * 1000.0), 0, 65535).astype(np.uint16)).save(os.path.join(a.out, "depth", f"{i:04d}.png"))
        Image.fromarray(l.astype(np.uint16)).save(os.path.join(a.out, "label", f"{i:04d}.png"))
print(f"wrote {n} frames to {a.out} (depth scale 1000)")
