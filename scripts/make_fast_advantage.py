#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Two-panel figure quantifying what a pointed FAST follow-up adds.

Usage: python scripts/make_fast_advantage.py [--dpi 900] [--outdir DIR]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch
from fig_gate import report_collisions

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "figures"

PRINT_W = 15.0 / 2.54
MIN_FS = 7.0


def FS(v: float) -> float:
    return max(MIN_FS, float(v))


def fit(ax, x, y, s, size, max_w=0.30, xy=None, **kw):
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    axw = ax.get_window_extent(renderer=r).width
    axh = ax.get_window_extent(renderer=r).height
    kw.setdefault("clip_on", False)
    if xy is not None:
        t = ax.annotate(s, xy=xy, xytext=(x, y), fontsize=size, **kw)
    else:
        t = ax.text(x, y, s, fontsize=size, **kw)
    lim = max_w * axw
    for _ in range(28):
        bb = t.get_window_extent(renderer=r)
        if bb.width <= lim and bb.height <= axh:
            break
        size = max(MIN_FS * 0.85, size - 0.3)
        t.set_fontsize(size)
    return t


plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.linewidth"] = 0.8

NAVY = "#1F3A5F"
LBLUE = "#5B8DBF"
TEAL = "#1B9E9E"
ORANGE = "#D9822B"
RED = "#C62828"
GREY = "#8894A0"
UGREY = "#C7D0D9"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dpi", type=float, default=900.0)
    ap.add_argument("--outdir", default=str(OUT))
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(PRINT_W, 3.15 * PRINT_W / 7.4),
                                   gridspec_kw=dict(width_ratios=[1.15, 1.0], wspace=0.36))

    F = np.logspace(-1, 3.2, 400)
    alpha = 1.5
    N0 = 1e4
    N = N0 * (F / 1.0) ** (-alpha)
    ax1.plot(F, N, color=NAVY, lw=1.4, zorder=5)
    ax1.fill_between(F, N * 0.5, N * 2.0, color=LBLUE, alpha=0.18, lw=0, zorder=2)

    F0_chime, F0_fast = 1.58, 0.028
    ax1.plot([F0_chime, F0_chime], [1e-1, 3e4], color=ORANGE, lw=1.3, ls="--", zorder=6)
    ax1.plot([F0_fast, F0_fast], [1e-1, 3e4], color=TEAL, lw=1.3, ls="--", zorder=6)
    WBOX = dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85)
    fit(ax1, F0_fast * 1.6, 2.0e2, f"FAST {F0_fast*1e3:.0f} mJy ms",
        FS(6.4), max_w=0.26, color=TEAL, ha="left", va="center",
        fontweight="bold", zorder=9, bbox=WBOX)
    fit(ax1, F0_chime * 1.15, 3.5e1, f"CHIME {F0_chime:.1f} Jy ms",
        FS(6.4), max_w=0.26, color=ORANGE, ha="left", va="center",
        fontweight="bold", zorder=9, bbox=WBOX)

    ax1.scatter([3.4, 2.7], [N0 * (3.4 / 1.0) ** -alpha, N0 * (2.7 / 1.0) ** -alpha],
                s=30, marker="*", color=RED, zorder=10,
                edgecolors="white", linewidths=0.4)
    fit(ax1, 0.97, 0.96,
        "this pair (measured 3.4 / 2.7 Jy ms):\n"
        "edge of the CHIME threshold, deep inside FAST", FS(6.4), max_w=0.55,
        xy=(3.6, N0 * (3.6) ** -alpha), textcoords=ax1.transAxes,
        color=RED, ha="right", va="top", zorder=11, linespacing=1.35, bbox=WBOX,
        arrowprops=dict(arrowstyle="-|>", color=RED, lw=0.7,
                        connectionstyle="arc3,rad=0.15"))

    ax1.annotate("", xy=(F0_chime, 0.25), xytext=(F0_fast, 0.25),
                 arrowprops=dict(arrowstyle="<|-|>", color=GREY, lw=0.8), zorder=7)
    fit(ax1, 0.30, 0.105, "threshold leverage $\\approx \\times 50$",
        FS(6.4), max_w=0.40, transform=ax1.transAxes,
        color=GREY, ha="left", va="bottom", fontweight="bold", bbox=WBOX)

    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlim(0.05, 1.5e3); ax1.set_ylim(1e-1, 3e4)
    ax1.set_xlabel("fluence density $F$ (Jy ms)", fontsize=FS(7.2))
    ax1.set_ylabel("cumulative event rate (arb. units)", fontsize=FS(7.2))
    ax1.tick_params(labelsize=FS(6.4))
    ax1.grid(True, which="major", ls=":", lw=0.4, color=UGREY, alpha=0.7)
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)
    ax1.set_title("Detection-threshold leverage: rate $\\propto F^{-3/2}$",
                  fontsize=FS(7.4), color=NAVY, weight="bold", pad=5)

    devs = ["CHIME", "VLA", "MeerKAT", "FAST"]
    sefd = {"CHIME": 100.0, "VLA": 24.0, "MeerKAT": 32.0, "FAST": 2.0}
    t_fast_ref = 0.4
    t = [t_fast_ref * (sefd[d] / sefd["FAST"]) ** 2 for d in devs]
    cols = [ORANGE, GREY, GREY, TEAL]
    y = np.arange(len(devs))[::-1]
    ax2.barh(y, t, height=0.52, color=cols, zorder=5)
    for yy, tt, dd in zip(y, t, devs):
        lab = f"{tt:g} h" if tt >= 1 else f"{tt*60:.0f} min"
        if dd == "FAST":
            ax2.text(tt * 1.5, yy, lab, va="center", ha="left",
                     fontsize=FS(6.8), fontweight="bold", color="#0E6B6B", zorder=12)
        else:
            ax2.text(tt * 1.15, yy, lab, va="center", ha="left",
                     fontsize=FS(6.8), fontweight="bold", color=NAVY, zorder=9)
    ax2.set_yticks(y); ax2.set_yticklabels(devs, fontsize=FS(7.0))
    ax2.set_xscale("log")
    ax2.set_xlim(0.15, 4.0e3); ax2.set_ylim(-0.65, 3.6)
    ax2.set_xticks([0.25, 1, 10, 100, 1000])
    ax2.set_xticklabels(["0.25", "1", "10", "100", "1000"])
    ax2.set_xlabel("integration time for S/N = 20, 1 ms pulse (h, log)",
                   fontsize=FS(7.2))
    ax2.tick_params(labelsize=FS(6.4))
    ax2.grid(True, axis="x", ls=":", lw=0.4, color=UGREY, alpha=0.7)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    ax2.set_title("Time-resource efficiency ($t \\propto$ SEFD$^2$)",
                  fontsize=FS(7.4), color=NAVY, weight="bold", pad=5)
    fit(ax2, 0.98, 0.99, "FAST verdict $\\approx 24$ min; others: $10^1$-$10^3$ h",
        FS(6.2), max_w=0.46, transform=ax2.transAxes, color="#0E6B6B",
        ha="right", va="top", fontweight="bold", zorder=9)
    ax2.add_patch(FancyBboxPatch((0.15, -0.45), 3.0, 0.9,
                                 boxstyle="round,pad=0.02,rounding_size=0.06",
                                 fc="#EAF6F5", ec=TEAL, lw=0.6, zorder=0))

    report_collisions(fig, "fig8_fast_advantage")
    fig.savefig(out / "fig8_fast_advantage.png", bbox_inches="tight", dpi=a.dpi)
    plt.close(fig)
    print(f"saved -> {out / 'fig8_fast_advantage.png'}  dpi={a.dpi:g}  "
          f"设计宽 {PRINT_W*2.54:.1f} cm  字号地板 {MIN_FS:g} pt")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    raise SystemExit(main())
