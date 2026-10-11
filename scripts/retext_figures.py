"""Pixel-level helper used by build_retext to erase and redraw text boxes."""
from __future__ import annotations
import json
import os
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
REF_DIR = Path(r"C:\Users\Admin\Downloads")
OUT = ROOT / "docs" / "figures" / "retext"
OUT.mkdir(parents=True, exist_ok=True)

_PLACEHOLDER: list = []


def flush_rotated(im: Image.Image):
    for tmp, (cx, cy) in _PLACEHOLDER:
        w, h = tmp.size
        im.paste(tmp, (int(cx - w / 2), int(cy - h / 2)), tmp)
    _PLACEHOLDER.clear()
    return im

FONT_DIRS = [Path(r"C:\Windows\Fonts")]
CJK_BOLD = "msyhbd.ttc"
CJK_REG = "msyh.ttc"
LAT_BOLD = "timesbd.ttf"
LAT_REG = "times.ttf"


def _find(name: str) -> str | None:
    for d in FONT_DIRS:
        p = d / name
        if p.exists():
            return str(p)
    return None


def font(size: int, cjk: bool = False, bold: bool = False):
    cand = []
    if cjk:
        cand += [CJK_BOLD if bold else CJK_REG]
    else:
        cand += [LAT_BOLD if bold else LAT_REG, "arialbd.ttf" if bold else "arial.ttf"]
    cand += [CJK_BOLD if bold else CJK_REG]
    for c in cand:
        f = _find(c)
        if f:
            try:
                return ImageFont.truetype(f, size)
            except Exception:
                continue
    return ImageFont.load_default()


REFS = {
    "composite": "一幅适用于Natur-图片-CHIME-FRB同源关联学术图.jpg",
    "flow": "一幅适用于Natur-图片-CHIME-FRB学术关联流程图.jpg",
    "geometry": "一幅适用于Natur-图片-CHIME望远镜几何周期示意图.jpg",
    "scatter": "一幅适用于Natur-图片-FRB配对学术散点图.jpg",
    "roadmap": "一幅适用于Natur-图片-前瞻引力波学术路线图.jpg",
    "dual": "一幅适用于Natur-图片-双面板学术对比图.jpg",
}


def sample_bg(im: Image.Image, box, pad=6) -> tuple[int, int, int]:
    a = np.asarray(im).astype(np.int16)
    x0, y0, x1, y1 = box
    H, W, _ = a.shape
    xs0, ys0 = max(0, x0 - pad), max(0, y0 - pad)
    xs1, ys1 = min(W, x1 + pad), min(H, y1 + pad)
    ring = np.concatenate([
        a[ys0:ys1, xs0:x0].reshape(-1, 3),
        a[ys0:ys1, x1:xs1].reshape(-1, 3),
        a[ys0:y0, xs0:xs1].reshape(-1, 3),
        a[y1:ys1, xs0:xs1].reshape(-1, 3),
    ]) if (x1 > x0 and y1 > y0) else a.reshape(-1, 3)
    if ring.size == 0:
        ring = a.reshape(-1, 3)
    med = np.median(ring, axis=0).astype(int)
    return tuple(int(v) for v in med)


def erase(im: Image.Image, box, pad=2, filled: bool = True, bg=None):
    x0, y0, x1, y1 = [int(v) for v in box]
    x0, y0 = max(0, x0 - pad), max(0, y0 - pad)
    x1, y1 = min(im.width, x1 + pad), min(im.height, y1 + pad)
    if x1 <= x0 or y1 <= y0:
        return
    color = bg if bg is not None else sample_bg(im, (x0, y0, x1, y1))
    ImageDraw.Draw(im).rectangle([x0, y0, x1 - 1, y1 - 1], fill=color)


def put(d: ImageDraw.ImageDraw, box, text, fnt: ImageFont.FreeTypeFont,
        color=(0, 0, 0), align="center", valign="center", anchor=None,
        rotate=0):
    x0, y0, x1, y1 = box
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    if rotate:
        tmp = Image.new("RGBA", (max(2, x1 - x0) * 3, max(2, y1 - y0) * 3),
                        (0, 0, 0, 0))
        td = ImageDraw.Draw(tmp)
        W2, H2 = tmp.size
        kw = {"anchor": "mm"} if anchor is None else {"anchor": anchor}
        td.multiline_text((W2 / 2, H2 / 2), text, font=fnt, fill=color,
                          align="center", spacing=4, **kw)
        tmp = tmp.rotate(rotate, expand=True, resample=Image.BICUBIC)
        _PLACEHOLDER.append((tmp, (int(cx), int(cy))))
        return
    kw = {}
    if anchor:
        kw["anchor"] = anchor
    if "\n" in text:
        lines = text.split("\n")
        hs = []
        for ln in lines:
            bb = d.textbbox((0, 0), ln, font=fnt)
            hs.append(bb[3] - bb[1])
        lh = max(hs) + 3
        total = lh * len(lines)
        ty = cy - total / 2
        for ln in lines:
            bb = d.textbbox((0, 0), ln, font=fnt)
            w = bb[2] - bb[0]
            tx = cx - w / 2 if align == "center" else (x0 if align == "left" else x1 - w)
            d.text((tx, ty), ln, font=fnt, fill=color)
            ty += lh
        return
    bb = d.textbbox((0, 0), text, font=fnt)
    w, h = bb[2] - bb[0], bb[3] - bb[1]
    if anchor is None:
        tx = cx - w / 2 if align == "center" else (x0 if align == "left" else x1 - w)
        ty = cy - h / 2 if valign == "center" else (y0 if valign == "top" else y1 - h)
        d.text((tx - bb[0], ty - bb[1]), text, font=fnt, fill=color)
    else:
        d.text((cx, cy), text, font=fnt, fill=color, anchor=anchor)


def erase_pill(im: Image.Image, box, pad=4):
    x0, y0, x1, y1 = [int(v) for v in box]
    x0, y0 = max(0, x0 - pad), max(0, y0 - pad)
    x1, y1 = min(im.width, x1 + pad), min(im.height, y1 + pad)
    a = np.asarray(im).astype(np.int16)
    corners = np.concatenate([
        a[y0:y0 + 6, x0:x0 + 6].reshape(-1, 3),
        a[y0:y0 + 6, x1 - 6:x1].reshape(-1, 3),
        a[y1 - 6:y1, x0:x0 + 6].reshape(-1, 3),
        a[y1 - 6:y1, x1 - 6:x1].reshape(-1, 3),
    ])
    color = tuple(int(v) for v in np.median(corners, axis=0))
    ImageDraw.Draw(im).rectangle([x0, y0, x1 - 1, y1 - 1], fill=color)
