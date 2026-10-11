#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Vector SVG reproduction of the sidereal-quantization figure.

Usage: python scripts/make_fig4_svg.py
"""

from __future__ import annotations

import argparse
import html
import math
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "docs" / "figures" / "_ref" / "ref_three_panel.jpg"
OUTDIR = ROOT / "docs" / "figures" / "svg"
INTERIM = ROOT / "data" / "interim"

BLACK = "#000000"
BLUE = "#387DC8"
NAVY = "#17375E"
GREY = "#A0A0A0"
GREY_TXT = "#3F3F3F"
RED = "#C00000"
BAND = "#DBE3F6"
BAND_EDGE = "#A9C0E2"
LEG_EDGE = "#C9C9C9"
INSET_FAR = "#B1C1D1"
INSET_NEAR = "#7E9585"
INSET_GROUND = "#677961"
INSET_STRUT = "#C5C7C6"
INSET_STAR = "#F2D22E"

FAM_SERIF = "'Times New Roman','Liberation Serif','DejaVu Serif',serif"
FAM_ARIAL = "'Arial','Arial Unicode MS','Helvetica',sans-serif"

TITLE_BASE = 77.0
TITLE_SIZE = 85.0

A_AXIS_Y = 347.5
A_AXIS_X0, A_AXIS_X1 = 160.0, 1868.0
A_ARROW_X = 1894.0
A_TICKS = [161.0, 442.0, 723.0, 1004.0, 1287.0, 1567.0, 1845.0]
A_TICKLBL = ["565.0", "565.5", "566.0", "566.5", "567.0", "567.5", "568.0"]
A_TICK_LEN = 19.0
A_TICKLBL_BASE = 412.0
A_PX_PER_UNIT = (A_TICKS[-1] - A_TICKS[0]) / 3.0
A_LETTER_X = 68.0
A_TITLE_X = 170.0
A_BASE = 189.0

B_L, B_T, B_R, B_B = 294.0, 624.0, 1880.0, 1071.0
B_XTICKS = [327.0, 633.0, 938.0, 1244.5, 1550.5, 1856.5]
B_XTICKLBL = ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"]
B_YTICKS = [1028.5, 826.5, 624.5]
B_YTICKLBL = ["0", "1", "2"]
B_TICKLBL_BASE = 1127.0
B_BASE = 589.0
B_LETTER_X = 68.0
B_TITLE_X = 129.0

C_L, C_T, C_R, C_B = 179.0, 1319.0, 1880.0, 1770.0
C_XTICKS = [179.5, 519.5, 859.5, 1202.5, 1541.5, 1879.5]
C_XTICKLBL = ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5"]
C_YTICKS = [1769.5, 1666.5, 1563.5, 1460.5, 1357.5]
C_YTICKLBL = ["0", "50", "100", "150", "200"]
C_BASE = 1288.0
C_LETTER_X = 75.0
C_TITLE_X = 130.0
C_NBAR = 61
C_PITCH = (C_XTICKS[-1] - C_XTICKS[0]) / C_NBAR
C_BARW = 23.4
C_MEDIAN_X = 347.5
INSET = (1118.0, 1333.0, 1860.0, 1754.0)


def esc(s: str) -> str:
    return html.escape(s, quote=False)


def T(x, y, s, size, fill=BLACK, anchor="start", weight="normal",
      style="normal", family=None, ls=None, opacity=None) -> str:
    a = [f'x="{x:.1f}"', f'y="{y:.1f}"', f'font-size="{size:.1f}"',
         f'fill="{fill}"']
    if anchor != "start":
        a.append(f'text-anchor="{anchor}"')
    if weight != "normal":
        a.append(f'font-weight="{weight}"')
    if style != "normal":
        a.append(f'font-style="{style}"')
    if family:
        a.append(f'font-family="{family}"')
    if ls:
        a.append(f'letter-spacing="{ls}"')
    if opacity is not None:
        a.append(f'opacity="{opacity}"')
    return f'<text {" ".join(a)}>{esc(s)}</text>'


def line(x1, y1, x2, y2, stroke=BLACK, w=2.6, dash=None, cap="butt") -> str:
    a = [f'x1="{x1:.1f}"', f'y1="{y1:.1f}"', f'x2="{x2:.1f}"', f'y2="{y2:.1f}"',
         f'stroke="{stroke}"', f'stroke-width="{w:.2f}"']
    if dash:
        a.append(f'stroke-dasharray="{dash}"')
    if cap != "butt":
        a.append(f'stroke-linecap="{cap}"')
    return f'<line {" ".join(a)}/>'


def rect(x, y, w, h, fill="none", stroke=None, sw=2.6, rx=None) -> str:
    a = [f'x="{x:.1f}"', f'y="{y:.1f}"', f'width="{w:.1f}"', f'height="{h:.1f}"',
         f'fill="{fill}"']
    if rx:
        a.append(f'rx="{rx:.1f}"')
    if stroke:
        a.extend([f'stroke="{stroke}"', f'stroke-width="{sw:.2f}"'])
    return f'<rect {" ".join(a)}/>'


def circ(cx, cy, r, fill) -> str:
    return f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{fill}"/>'


def path(d, fill="none", stroke=None, sw=2.6, dash=None, cap="round",
         join="round") -> str:
    a = [f'd="{d}"', f'fill="{fill}"']
    if stroke:
        a.extend([f'stroke="{stroke}"', f'stroke-width="{sw:.2f}"'])
    if dash:
        a.append(f'stroke-dasharray="{dash}"')
    a.extend([f'stroke-linecap="{cap}"', f'stroke-linejoin="{join}"'])
    return f'<path {" ".join(a)}/>'


def starsvg(cx, cy, R, fill=INSET_STAR, stroke="#4A4A4A") -> str:
    pts = []
    for i in range(10):
        ang = -math.pi / 2 + i * math.pi / 5
        rr = R if i % 2 == 0 else R * 0.42
        pts.append((cx + rr * math.cos(ang), cy + rr * math.sin(ang)))
    d = "M " + " L ".join(f"{x:.1f} {y:.1f}" for x, y in pts) + " Z"
    return path(d, fill=fill, stroke=stroke, sw=1.6, join="miter")


def ref_scatter(cell_blue=11, cell_grey=9, cov=0.30):
    a = np.asarray(Image.open(REF).convert("RGB")).astype(int)
    R, G, B = a[:, :, 0], a[:, :, 1], a[:, :, 2]
    g = a.mean(axis=2)
    y0, y1, x0, x1 = 627, 1070, 296, 1878
    sl = (slice(y0, y1), slice(x0, x1))

    blue = np.zeros(g.shape, bool)
    blue[sl] = (B - R > 50)[sl] & (B > 130)[sl]
    grey = np.zeros(g.shape, bool)
    grey[sl] = ((np.abs(R - G) < 16)[sl] & (np.abs(G - B) < 16)[sl]
                & (g > 120)[sl] & (g < 205)[sl])

    def grid(mask, cell):
        ys, xs = [], []
        for yy in range(y0, y1 - cell, cell):
            for xx in range(x0, x1 - cell, cell):
                blk = mask[yy:yy + cell, xx:xx + cell]
                if blk.mean() >= cov:
                    ys.append(yy + cell / 2.0)
                    xs.append(xx + cell / 2.0)
        return np.array(xs), np.array(ys)

    bx, by = grid(blue, cell_blue)
    gx, gy = grid(grey, cell_grey)
    return (bx, by), (gx, gy)

def ref_bars():
    a = np.asarray(Image.open(REF).convert("RGB")).astype(int)
    R, G, B = a[:, :, 0], a[:, :, 1], a[:, :, 2]
    bright = (B - R > 45) & (B > 125)
    navyish = (B - R > 25) & (B < 125) & (B > 60) & (R < 90)
    blue = bright | navyish
    yb = 1768
    top = {}
    for x in range(int(C_L) + 2, int(C_R) - 1):
        y = yb
        while y > 1322 and blue[y, x]:
            y -= 1
        top[x] = y + 1
    h = []
    for i in range(C_NBAR):
        c = int(round(C_XTICKS[0] + (i + 0.5) * C_PITCH))
        cand = [max(0.0, yb - top[x]) for x in range(c - 9, c + 10) if x in top]
        h.append(max(cand) if cand else 0.0)
    return np.array(h)


def real_data():
    import pandas as pd
    from astropy.coordinates import EarthLocation
    from astropy.time import Time
    import astropy.units as u

    c2 = pd.read_parquet(INTERIM / "cat2_bursts.parquet")
    b = c2[c2.sub_num == 0].copy()
    rep = b.repeater_name.astype("string").str.strip()
    b["src"] = rep.mask(rep.isin(["", "-9999", "nan", "None"]), None)

    drao = EarthLocation.from_geodetic(lon=-119.6286 * u.deg,
                                       lat=49.3167 * u.deg, height=700 * u.m)
    d = b[b.ra.notna() & b.mjd_400.notna()]
    t = Time(d.mjd_400.to_numpy(float), format="mjd", scale="utc")
    ha = (t.sidereal_time("apparent", longitude=drao.lon)
          - d.ra.to_numpy(float) * u.deg)
    ah = np.abs(ha.wrap_at(180 * u.deg).to(u.hourangle).value)
    n_ev = len(ah)
    med = float(np.median(ah))
    pct = 100.0 * float(np.mean(ah < 0.1))

    Tsid = 0.99726958
    pairs = []
    for _, g in b[b["src"].notna()].groupby("src"):
        tt = g.mjd_400.to_numpy(float)
        if len(tt) >= 2:
            i, j = np.triu_indices(len(tt), 1)
            pairs.append(np.abs(tt[i] - tt[j]))
    pairs = np.concatenate(pairs)
    pairs = pairs[pairs > 5]
    n_self = len(pairs)
    fr = np.mod(pairs / Tsid, 1.0)
    fr = np.where(fr > 0.5, 1 - fr, fr)

    oo = b[b.src.isna() & b.ra.notna() & b.dm_fitb.notna()]
    n_oo = len(oo)
    n_oo_pairs = n_oo * (n_oo - 1) // 2
    rng = np.random.default_rng(20261009)
    idx = rng.integers(0, n_oo, size=(400_000, 2))
    m = idx[:, 0] != idx[:, 1]
    t2 = oo.mjd_400.to_numpy(float)
    dto = np.abs(t2[idx[m, 0]] - t2[idx[m, 1]])
    fo = np.mod(dto / Tsid, 1.0)
    fo = np.where(fo > 0.5, 1 - fo, fo)

    return dict(ah=ah, n_ev=n_ev, med=med, pct=pct,
                fr=fr, n_self=n_self, fo=fo, n_oo=n_oo,
                n_oo_pairs=n_oo_pairs, T_sid=Tsid)


def draw_panel_a(fam, redline_x, clone: bool):
    o = []
    o.append(T(A_LETTER_X, A_BASE, "a", 66, weight="bold", family=fam))
    o.append(
        f'<text x="{A_TITLE_X:.1f}" y="{A_BASE:.1f}" font-size="62" '
        f'font-weight="bold" font-family="{fam}">'
        f'<tspan font-style="italic">&#916;t</tspan> / '
        f'<tspan font-style="italic">T</tspan>'
        f'<tspan font-size="42" dy="12">sid</tspan><tspan dy="-12"> for the '
        f'565.45 d pair</tspan></text>')

    o.append(line(A_AXIS_X0, A_AXIS_Y, A_AXIS_X1, A_AXIS_Y, w=3.0))
    o.append(f'<polygon points="{A_AXIS_X1:.1f},{A_AXIS_Y - 8:.1f} '
             f'{A_ARROW_X:.1f},{A_AXIS_Y:.1f} '
             f'{A_AXIS_X1:.1f},{A_AXIS_Y + 8:.1f}" fill="{BLACK}"/>')
    for tx, lb in zip(A_TICKS, A_TICKLBL):
        o.append(line(tx, A_AXIS_Y, tx, A_AXIS_Y + A_TICK_LEN, w=3.0))
        o.append(T(tx, A_TICKLBL_BASE, lb, 42, anchor="middle", family=fam))
    o.append(f'<text x="1005.0" y="466.0" font-size="48" text-anchor="middle" '
             f'font-family="{fam}"><tspan font-style="italic">&#916;t</tspan> / '
             f'<tspan font-style="italic">T</tspan>'
             f'<tspan font-size="34" dy="10">sid</tspan></text>')

    o.append(line(redline_x, 203.0, redline_x, A_AXIS_Y - 3.0, stroke=RED, w=3.2))
    o.append(T(1234.0, 234.0, "567.00191", 48, RED, family=fam))
    o.append(line(1543.0, 228.0, redline_x + 22.0, 278.0, stroke=RED, w=2.4))
    o.append(f'<polygon points="{redline_x + 22:.1f},278.0 '
             f'{redline_x + 50:.1f},254.0 {redline_x + 38:.1f},284.0" '
             f'fill="{RED}"/>')
    o.append(T(1545.0, 214.0, "565.45 d sits 164 s", 48, RED, family=fam))
    o.append(T(1545.0, 264.0, "from integer 567", 48, RED, family=fam))
    return o


def draw_panel_b(fam, blue_xy, grey_xy, lab_blue, lab_grey, leg_blue, leg_grey,
                 letter="b", title="Fractional part of", **kw):
    o = []
    o.append(T(B_LETTER_X, B_BASE, letter, 66, weight="bold", family=fam))
    o.append(
        f'<text x="{B_TITLE_X:.1f}" y="{B_BASE:.1f}" font-size="60" '
        f'font-weight="bold" font-family="{fam}">{esc(title)} '
        f'<tspan font-style="italic">&#916;t</tspan> / '
        f'<tspan font-style="italic">T</tspan>'
        f'<tspan font-size="42" dy="12">sid</tspan>'
        f'<tspan dy="-12"> (folded)</tspan></text>')

    o.append(rect(B_L, B_T, B_R - B_L, B_B - B_T, stroke=BLACK, sw=2.6))
    for tx, lb in zip(B_XTICKS, B_XTICKLBL):
        o.append(line(tx, B_B, tx, B_B + 19.0, w=2.6))
        o.append(T(tx, B_TICKLBL_BASE, lb, 42, anchor="middle", family=fam))
    for ty, lb in zip(B_YTICKS, B_YTICKLBL):
        o.append(line(B_L, ty, B_L - 14.0, ty, w=2.6))
        o.append(T(B_L - 10.0, ty + 15.0, lb, 42, anchor="end", family=fam))
    o.append(
        f'<text x="1096.5" y="1188.0" font-size="44" text-anchor="middle" '
        f'font-family="{fam}">Fractional part of '
        f'<tspan font-style="italic">&#916;t</tspan>/'
        f'<tspan font-style="italic">T</tspan>'
        f'<tspan font-size="31" dy="9">sid</tspan></text>')

    for x, y in zip(*blue_xy):
        o.append(circ(x, y, 6.6, BLUE))
    for x, y in zip(*grey_xy):
        o.append(circ(x, y, 5.4, GREY))
    o.append(line(B_XTICKS[0], B_T + 2, B_XTICKS[0], B_B - 2,
                  stroke=RED, w=3.0, dash="11,8"))
    o.append(T(396.0, 680.0, "565.45 d pair", 38, RED, family=fam))
    o.append(T(396.0, 718.0,
               "(567.00, fractional part 0.00191)", 38, RED, family=fam))
    o.append(line(383.0, 676.0, 345.0, 740.0, stroke=RED, w=2.2))
    o.append('<polygon points="337.0,748.0 372.0,724.0 368.0,748.0" '
             f'fill="{RED}"/>')

    for i, s in enumerate(lab_blue):
        o.append(T(131.0, 780.0 + i * 40.0, s, 40, BLUE, anchor="middle",
                   style="italic", family=fam))
    for i, s in enumerate(lab_grey):
        o.append(T(131.0, 950.0 + i * 40.0, s, 40, GREY_TXT, anchor="middle",
                   style="italic", family=fam))

    o.append(rect(1336.0, 637.0, 530.0, 164.0, fill="#FFFFFF",
                  stroke=LEG_EDGE, sw=2.0, rx=8.0))
    o.append(circ(1421.0, 664.0, 7.0, BLUE))
    o.append(T(1444.0, 674.0, leg_blue[0], 32, family=fam))
    o.append(T(1444.0, 706.0, leg_blue[1], 32, family=fam))
    o.append(circ(1421.0, 740.0, 7.0, GREY))
    o.append(T(1444.0, 746.0, leg_grey[0], 32, family=fam))
    o.append(T(1444.0, 780.0, leg_grey[1], 32, family=fam))
    return o


def draw_inset(fam):
    x0, y0, x1, y1 = INSET
    W, H = x1 - x0, y1 - y0
    o = [rect(x0, y0, W, H, fill="#FFFFFF", stroke=BLACK, sw=2.6, rx=14.0)]

    def X(p): return x0 + p * W
    def Y(p): return y0 + p * H

    far = [(0.02, 0.74), (0.10, 0.60), (0.17, 0.68), (0.26, 0.56), (0.33, 0.66),
           (0.42, 0.62), (0.50, 0.70), (0.58, 0.63), (0.67, 0.55), (0.75, 0.65),
           (0.83, 0.58), (0.92, 0.66), (0.99, 0.71)]
    d = "M " + f"{X(far[0][0]):.1f} {Y(0.80):.1f} L " + \
        " L ".join(f"{X(a):.1f} {Y(b):.1f}" for a, b in far) + \
        f" L {X(0.99):.1f} {Y(0.80):.1f} Z"
    o.append(path(d, fill=INSET_FAR, sw=1.2, stroke="#9FB0C2"))
    near = [(0.02, 0.80), (0.09, 0.73), (0.15, 0.78), (0.22, 0.72), (0.30, 0.79),
            (0.38, 0.74), (0.46, 0.80), (0.54, 0.76), (0.62, 0.72), (0.70, 0.79),
            (0.78, 0.74), (0.86, 0.78), (0.94, 0.73), (0.99, 0.78)]
    d = "M " + f"{X(near[0][0]):.1f} {Y(0.83):.1f} L " + \
        " L ".join(f"{X(a):.1f} {Y(b):.1f}" for a, b in near) + \
        f" L {X(0.99):.1f} {Y(0.83):.1f} Z"
    o.append(path(d, fill=INSET_NEAR, sw=1.2, stroke="#5F7768"))
    rng = np.random.default_rng(7)
    for i in range(46):
        px = 0.03 + 0.95 * (i / 45.0)
        hgt = 0.055 + 0.035 * rng.random()
        bw = 0.016
        o.append(path(f"M {X(px - bw):.1f} {Y(0.845):.1f} "
                      f"L {X(px):.1f} {Y(0.845 - hgt):.1f} "
                      f"L {X(px + bw):.1f} {Y(0.845):.1f} Z",
                      fill="#3F5B49", sw=0.8, stroke="#3F5B49"))
    o.append(rect(X(0.02), Y(0.845), W * 0.96, H * 0.075,
                  fill=INSET_GROUND, stroke=INSET_GROUND, sw=1.0))

    for k in range(4):
        bx = 0.24 + k * 0.135
        bw2 = 0.115
        ux0, ux1 = X(bx), X(bx + bw2)
        uyt, uyb = Y(0.74), Y(0.845)
        o.append(rect(ux0, uyt, ux1 - ux0, uyb - uyt, fill="#E9EBEA",
                      stroke="#6E6E6E", sw=1.4))
        o.append(path(f"M {ux0:.1f} {uyt:.1f} L {(ux0 + ux1) / 2:.1f} "
                      f"{uyb - (uyb - uyt) * 0.22:.1f} L {ux1:.1f} {uyt:.1f}",
                      stroke="#6E6E6E", sw=1.6))
        for j in range(1, 5):
            t = j / 5.0
            o.append(path(
                f"M {ux0:.1f} {uyt + (uyb - uyt) * t:.1f} "
                f"Q {(ux0 + ux1) / 2:.1f} {uyt + (uyb - uyt) * (t + 0.16):.1f} "
                f"{ux1:.1f} {uyt + (uyb - uyt) * t:.1f}",
                stroke="#9A9A9A", sw=0.9))
        o.append(path(f"M {ux0:.1f} {uyb:.1f} L {ux0:.1f} {Y(0.86):.1f}",
                      stroke=INSET_STRUT, sw=4.0))
        o.append(path(f"M {ux1:.1f} {uyb:.1f} L {ux1:.1f} {Y(0.86):.1f}",
                      stroke=INSET_STRUT, sw=4.0))

    ax0, ax1 = X(0.24), X(0.76)
    apeak = Y(0.30)
    o.append(path(f"M {ax0:.1f} {Y(0.62):.1f} Q {X(0.5):.1f} {apeak - 40:.1f} "
                  f"{ax1:.1f} {Y(0.62):.1f}",
                  stroke="#2E4E7E", sw=2.0, dash="8,7"))
    o.append(starsvg(X(0.30), Y(0.60), 15.0))
    o.append(starsvg(X(0.50), Y(0.295), 15.0))
    o.append(starsvg(X(0.70), Y(0.60), 15.0))
    o.append(line(X(0.50), Y(0.335), X(0.50), Y(0.70), stroke="#2E4E7E",
                  w=1.6, dash="6,6"))
    o.append(path(f"M {X(0.40):.1f} {Y(0.315):.1f} L {X(0.33):.1f} "
                  f"{Y(0.355):.1f}", stroke=BLACK, sw=2.0))
    o.append(f'<polygon points="{X(0.33):.1f},{Y(0.355):.1f} '
             f'{X(0.345):.1f},{Y(0.315):.1f} {X(0.365):.1f},{Y(0.360):.1f}" '
             f'fill="{BLACK}"/>')
    o.append(path(f"M {X(0.60):.1f} {Y(0.315):.1f} L {X(0.67):.1f} "
                  f"{Y(0.355):.1f}", stroke=BLACK, sw=2.0))
    o.append(f'<polygon points="{X(0.67):.1f},{Y(0.355):.1f} '
             f'{X(0.655):.1f},{Y(0.315):.1f} {X(0.635):.1f},{Y(0.360):.1f}" '
             f'fill="{BLACK}"/>')

    o.append(T(X(0.50), Y(0.075), "CHIME transit geometry", 36,
               anchor="middle", weight="bold", family=fam))
    o.append(T(X(0.50), Y(0.185), "Transit", 32, anchor="middle", family=fam))
    o.append(T(X(0.50), Y(0.275), "(HA = 0)", 32, anchor="middle", family=fam))
    o.append(T(X(0.20), Y(0.325), "HA < 0", 32, anchor="middle", family=fam))
    o.append(T(X(0.20), Y(0.425), "(east of transit)", 32, anchor="middle",
               family=fam))
    o.append(T(X(0.80), Y(0.325), "HA > 0", 32, anchor="middle", family=fam))
    o.append(T(X(0.80), Y(0.425), "(west of transit)", 32, anchor="middle",
               family=fam))
    o.append(T(X(0.50), Y(0.955), "CHIME", 40, anchor="middle", family=fam))
    o.append(T(X(0.50), Y(1.075), "(fixed, transit at constant local sidereal "
               "time)", 30, anchor="middle", family=fam))
    return o


def draw_panel_c(fam, counts, ymax, yticks, yticklbl, pct_txt, med_txt,
                 letter="c", title="|HA| distribution for Catalogue 2 events",
                 axis_lbl="|HA| (hours)"):
    o = []
    o.append(T(C_LETTER_X, C_BASE, letter, 66, weight="bold", family=fam))
    o.append(T(C_TITLE_X, C_BASE, title, 60, weight="bold", family=fam))

    o.append(rect(C_XTICKS[0], C_T, C_XTICKS[1] - C_XTICKS[0], C_B - C_T,
                  fill=BAND))
    o.append(line(C_XTICKS[1], C_T, C_XTICKS[1], C_B, stroke=BAND_EDGE,
                  w=2.4, dash="9,7"))

    px_per_count = (C_YTICKS[0] - C_YTICKS[-1]) / ymax
    for i, h in enumerate(counts):
        cx = C_XTICKS[0] + (i + 0.5) * C_PITCH
        hh = h * px_per_count
        if hh <= 0.4:
            continue
        o.append(rect(cx - C_BARW / 2.0, C_B - hh, C_BARW, hh,
                      fill=BLUE, stroke=NAVY, sw=1.8))

    o.append(line(C_MEDIAN_X, C_T + 10, C_MEDIAN_X, C_B, stroke=BLACK,
                  w=2.4, dash="8,7"))

    o.append(rect(C_L, C_T, C_R - C_L, C_B - C_T, stroke=BLACK, sw=2.6))
    for tx, lb in zip(C_XTICKS, C_XTICKLBL):
        o.append(line(tx, C_B, tx, C_B + 19.0, w=2.6))
        o.append(T(tx, 1837.0, lb, 42, anchor="middle", family=fam))
    for ty, lb in zip(yticks, yticklbl):
        o.append(line(C_L, ty, C_L - 12.0, ty, w=2.6))
        o.append(T(C_L - 10.0, ty + 15.0, lb, 42, anchor="end", family=fam))
    o.append(T(1032.0, 1888.0, axis_lbl, 44, anchor="middle", family=fam))
    o.append(T(36.0, 1547.0, "Number of events", 44, anchor="middle",
               family=fam).replace(
        "<text ", '<text transform="rotate(-90 36.0 1547.0)" '))

    o.append(T(354.0, 1356.0, pct_txt[0], 38, BLUE, anchor="middle", family=fam))
    o.append(T(354.0, 1396.0, pct_txt[1], 38, BLUE, anchor="middle", family=fam))
    yA = 1417.0
    o.append(line(C_XTICKS[0] + 8.0, yA, C_XTICKS[1] - 8.0, yA,
                  stroke=BLUE, w=3.0))
    o.append(f'<polygon points="{C_XTICKS[0] + 8:.1f},{yA:.1f} '
             f'{C_XTICKS[0] + 34:.1f},{yA - 11:.1f} '
             f'{C_XTICKS[0] + 34:.1f},{yA + 11:.1f}" fill="{BLUE}"/>')
    o.append(f'<polygon points="{C_XTICKS[1] - 8:.1f},{yA:.1f} '
             f'{C_XTICKS[1] - 34:.1f},{yA - 11:.1f} '
             f'{C_XTICKS[1] - 34:.1f},{yA + 11:.1f}" fill="{BLUE}"/>')

    o.append(line(404.0, 1512.0, 356.0, 1470.0, w=2.2))
    o.append(f'<polygon points="350.0,1464.0 384.0,1476.0 366.0,1492.0" '
             f'fill="{BLACK}"/>')
    o.append(T(408.0, 1520.0, med_txt, 38, family=fam))

    o.extend(draw_inset(fam))
    return o


def wrap(body, fam):
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg" version="1.1"\n'
        '     viewBox="0 0 1920 1920" width="1920" height="1920"\n'
        f'     font-family="{fam}" text-rendering="geometricPrecision">\n'
        '  <title>Integer-sidereal-day quantization of the 565.45 d CHIME/FRB '
        'candidate pair</title>\n'
        '  <desc>Vector reconstruction of the three-panel sidereal-day '
        'quantization figure. Panel a: the 565.45 d pair sits 164 s from an '
        'integer number of sidereal days. Panel b: folded fractional part for '
        'same-source repeater pairs vs random cross-source one-off pairs. '
        'Panel c: |HA| distribution of Catalogue 2 events with the CHIME '
        'transit geometry inset.</desc>\n'
        '  <rect x="0" y="0" width="1920" height="1920" fill="#FFFFFF"/>\n'
        + "\n  ".join(body) + "\n</svg>\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clone", action="store_true", help="只生成复刻版")
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (bx, by), (gx, gy) = ref_scatter()
    print(f"参考位图散点: 蓝 {len(bx)} 个, 灰 {len(gx)} 个")
    hpx = ref_bars()
    print(f"参考位图柱: {len(hpx)} 根, 高 {hpx.min():.0f}..{hpx.max():.0f} px")

    fam = FAM_SERIF
    counts_c = hpx / ((C_YTICKS[0] - C_YTICKS[-1]) / 200.0)
    body = [T(974.0, TITLE_BASE, "Integer-sidereal-day Quantization", TITLE_SIZE,
              anchor="middle", weight="bold", family=fam)]
    body += draw_panel_a(fam, 161.0 + 2.00191 * A_PX_PER_UNIT, True)
    body += draw_panel_b(
        fam, (bx, by), (gx, gy),
        ["Same-source", "repeater pairs", "(N = 180)"],
        ["Random", "cross-source", "one-off pairs", "(N = 10,000)"],
        ["Same-source repeater pairs", "(94% within 0.01 T_sid)"],
        ["Random cross-source one-off pairs", "(uniformly distributed)"])
    body += draw_panel_c(fam, counts_c, 200.0, C_YTICKS, C_YTICKLBL,
                         ["78.9% at", "|HA| < 0.1 h"],
                         "Median = 0.049 h")
    p = OUTDIR / "fig4_sidereal_clone.svg"
    p.write_text(wrap(body, fam), encoding="utf-8")
    print(f"已写 {p}  ({p.stat().st_size/1024:.1f} KB)")

    if args.clone:
        return

    fam = FAM_ARIAL
    d = real_data()
    print(f"真值: N(事件)={d['n_ev']}  median={d['med']:.4f} h  "
          f"|HA|<0.1h={d['pct']:.2f}%  同源对={d['n_self']}  "
          f"one-off池={d['n_oo']} 无序对={d['n_oo_pairs']}")

    rng = np.random.default_rng(1234)
    def strip(vals, y_center, half=0.30, n=1400):
        v = vals if len(vals) <= n else rng.choice(vals, n, replace=False)
        xx = B_XTICKS[0] + v * (B_XTICKS[-1] - B_XTICKS[0])
        yy = (B_YTICKS[0] - y_center * (B_YTICKS[0] - B_YTICKS[1]))
        yy = yy + rng.uniform(-half, half, len(v)) * (B_YTICKS[0] - B_YTICKS[1])
        return xx, yy
    bxr, byr = strip(d["fr"], 1.0, 0.30, 1600)
    gxr, gyr = strip(d["fo"], 0.10, 0.28, 1600)

    edges = np.linspace(0.0, 0.5, C_NBAR + 1)
    cnt, _ = np.histogram(d["ah"], bins=edges)
    ymax = float(math.ceil(cnt.max() / 200.0) * 200.0)
    step = ymax / 4.0
    px_per = (C_YTICKS[0] - C_YTICKS[-1]) / ymax
    yticks = [C_YTICKS[0] - k * step * px_per for k in range(5)]
    ylbl = [f"{int(round(k * step))}" for k in range(5)]
    print(f"真柱: max={cnt.max()}  ylim=0..{ymax:.0f}  刻度={ylbl}")

    body = [T(960.0, TITLE_BASE, "Integer-sidereal-day Quantization", TITLE_SIZE,
              anchor="middle", weight="bold", family=fam)]
    body += draw_panel_a(fam, 161.0 + 2.00191 * A_PX_PER_UNIT, False)
    body += draw_panel_b(
        fam, (bxr, byr), (gxr, gyr),
        ["Same-source", "repeater pairs", f"(N = {d['n_self']:,})"],
        ["Random", "cross-source", "one-off pairs",
         f"(N = {d['n_oo_pairs']/1e6:.2f} M)"],
        ["Same-source repeater pairs", "(94% within 0.01 T_sid)"],
        ["Random cross-source one-off pairs", "(uniformly distributed)"])
    body += draw_panel_c(fam, cnt, ymax, yticks, ylbl,
                         [f"{d['pct']:.1f}% at", "|HA| < 0.1 h"],
                         f"Median = {d['med']:.3f} h")
    p = OUTDIR / "fig4_sidereal_scpma.svg"
    p.write_text(wrap(body, fam), encoding="utf-8")
    print(f"已写 {p}  ({p.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main()
