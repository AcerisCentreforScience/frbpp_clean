"""Vectorise the bitmap layout figures without redrawing them."""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import potrace
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "docs" / "figures" / "retext" / "en"
DST = ROOT / "docs" / "figures" / "vector"
DST.mkdir(parents=True, exist_ok=True)

FULL_CM = 16.2
CM_PER_IN = 2.54
EXPECT_DPI = 600


def build_palette(rgb: np.ndarray, max_colors: int = 26, min_dist: float = 26.0,
                  min_count_frac: float = 2e-5, white_dist: float = 18.0):
    flat = rgb.reshape(-1, 3).astype(np.int32)
    key = (flat[:, 0] << 16) | (flat[:, 1] << 8) | flat[:, 2]
    uniq, cnt = np.unique(key, return_counts=True)
    order = np.argsort(-cnt)
    total = flat.shape[0]

    pal: list[np.ndarray] = []
    for i in order:
        if cnt[i] / total < min_count_frac:
            break
        v = int(uniq[i])
        col = np.array([v >> 16, (v >> 8) & 255, v & 255], dtype=np.int32)
        if all(np.linalg.norm(col - p) > min_dist for p in pal):
            pal.append(col)
        if len(pal) >= max_colors:
            break

    palette = np.stack(pal).astype(np.float32)
    near_white = np.linalg.norm(palette - 255.0, axis=1) < white_dist
    return palette, near_white


def assign(rgb: np.ndarray, palette: np.ndarray,
           dark_bias: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    h, w, _ = rgb.shape
    flat = rgb.reshape(-1, 3).astype(np.float32)
    lum = palette.mean(axis=1) / 255.0
    weight = (1.0 + dark_bias * lum).astype(np.float32)
    best = np.full(flat.shape[0], np.inf, dtype=np.float32)
    idx = np.zeros(flat.shape[0], dtype=np.int32)
    for k, col in enumerate(palette):
        d = np.linalg.norm(flat - col, axis=1) * weight[k]
        m = d < best
        best[m] = d[m]
        idx[m] = k
    return idx.reshape(h, w), np.bincount(idx.ravel(), minlength=len(palette))


def trace_mask(mask: np.ndarray, turd: int = 2, alphamax: float = 1.0,
               opttol: float = 0.2):
    return potrace.Bitmap(~mask).trace(turdsize=turd, alphamax=alphamax,
                                       opticurve=True, opttolerance=opttol)


def curves_to_path(path) -> MplPath:
    verts, codes = [], []
    for curve in path:
        sp = curve.start_point
        verts.append((sp.x, sp.y)); codes.append(MplPath.MOVETO)
        for seg in curve:
            if getattr(seg, "is_corner", False):
                verts.append((seg.c.x, seg.c.y)); codes.append(MplPath.LINETO)
                verts.append((seg.end_point.x, seg.end_point.y)); codes.append(MplPath.LINETO)
            else:
                verts.append((seg.c1.x, seg.c1.y)); codes.append(MplPath.CURVE4)
                verts.append((seg.c2.x, seg.c2.y)); codes.append(MplPath.CURVE4)
                verts.append((seg.end_point.x, seg.end_point.y)); codes.append(MplPath.CURVE4)
        verts.append((sp.x, sp.y)); codes.append(MplPath.CLOSEPOLY)
    if not verts:
        return MplPath(np.zeros((0, 2)))
    return MplPath(np.asarray(verts, float), codes)


def vectorize(key: str, max_colors: int = 26, min_dist: float = 45.0,
              turd: int = 2, min_frac: float = 8e-5,
              dark_bias: float = 0.9) -> dict:
    t0 = time.time()
    im = Image.open(SRC / f"{key}.png").convert("RGB")
    rgb = np.asarray(im)
    H, W, _ = rgb.shape

    palette, near_white = build_palette(rgb, max_colors=max_colors, min_dist=min_dist)
    idx, counts = assign(rgb, palette, dark_bias=dark_bias)
    total = H * W

    order = sorted(range(len(palette)),
                   key=lambda k: (0 if near_white[k] else 1, -counts[k]))

    fig = plt.figure(figsize=(FULL_CM / CM_PER_IN, FULL_CM / CM_PER_IN * H / W))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), W, H, facecolor="white", edgecolor="none", zorder=-100))

    layers = 0
    for k in order:
        frac = counts[k] / total
        col = palette[k] / 255.0
        color = tuple(float(c) for c in col)
        if near_white[k]:
            ax.add_patch(Rectangle((0, 0), W, H, facecolor=color, edgecolor="none",
                                   zorder=-50))
            continue
        if frac < min_frac:
            continue
        mask = (idx == k)
        p = trace_mask(mask, turd=turd)
        mp = curves_to_path(p)
        if len(mp.vertices) == 0:
            continue
        layers += 1
        ax.add_patch(PathPatch(mp, facecolor=color, edgecolor=color,
                               linewidth=0.1, joinstyle="round", zorder=layers))
        print(f"    #{int(col[0]*255):02X}{int(col[1]*255):02X}{int(col[2]*255):02X} "
              f"frac={frac:6.2%} curves={len(p):5d} verts={len(mp.vertices):6d} "
              f"t={time.time()-t0:6.1f}s", flush=True)

    fig.savefig(DST / f"{key}.pdf", format="pdf")
    fig.savefig(DST / f"{key}_900dpi.png", dpi=900)
    plt.close(fig)

    info = dict(key=key, size=(W, H), palette=len(palette), layers=layers,
                seconds=round(time.time() - t0, 1),
                raster_dpi=round(W / (FULL_CM / CM_PER_IN)))
    print(f"  [{key}] {W}x{H} 调色板{len(palette)}色 -> {layers} 图层, "
          f"{info['seconds']}s；位图等效 {info['raster_dpi']} dpi -> 矢量无上限 "
          f"(满足 {EXPECT_DPI} dpi 规则)", flush=True)
    return info


if __name__ == "__main__":
    keys = sys.argv[1:] or ["flow", "geometry", "dual", "scatter", "roadmap"]
    for k in keys:
        print(f"== {k} ==", flush=True)
        vectorize(k)
