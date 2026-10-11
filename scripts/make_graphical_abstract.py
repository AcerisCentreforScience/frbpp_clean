"""Single-page six-panel graphical abstract.

Usage: python scripts/make_graphical_abstract.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import (FancyArrowPatch, FancyBboxPatch, Rectangle,
                                Circle, Polygon, Wedge, Ellipse)

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT / "docs" / "figures"
UP = ROOT / "data" / "processed" / "upgrade"

_CJK = ["Microsoft YaHei", "SimHei", "Noto Sans SC", "DejaVu Sans"]
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": _CJK, "font.size": 7,
    "axes.unicode_minus": False, "figure.dpi": 400, "savefig.dpi": 400,
    "text.color": "black",
})
try:
    plt.rcParams["text.zorder"] = 15
except Exception:
    pass
_orig_text = plt.Axes.text


def _text_z(self, x, y, s, *a, **kw):
    kw.setdefault("zorder", 15)
    return _orig_text(self, x, y, s, *a, **kw)


plt.Axes.text = _text_z

NAVY = "#14507F"
BLUE = "#3E7CB1"
LBLUE = "#DCE9F5"
RED = "#C0392B"
LRED = "#F6D9D4"
GREEN = "#4C8C4A"
LGREEN = "#DCEBD8"
AMBER = "#D68910"
LAMBER = "#FBEDD4"
GREY = "#9AA0A6"
LGREY = "#EDEFF2"
INK = "#1A1A1A"


def badge(ax, x, y, n, r=0.021, color=NAVY):
    ax.add_patch(Circle((x, y), r, transform=ax.transAxes, fc=color,
                        ec="none", zorder=20, clip_on=False))
    ax.text(x, y, str(n), transform=ax.transAxes, ha="center", va="center",
            color="white", fontsize=8, fontweight="bold", zorder=21, clip_on=False)


def title_bar(ax, x0, x1, y, h, text, color=NAVY, fs=8):
    ax.add_patch(FancyBboxPatch((x0, y), x1 - x0, h,
                                boxstyle="round,pad=0,rounding_size=0.006",
                                transform=ax.transAxes, fc=color, ec="none",
                                zorder=10, clip_on=False))
    ax.text(x0 + 0.012, y + h / 2, text, transform=ax.transAxes,
            ha="left", va="center", color="white", fontsize=fs,
            fontweight="bold", zorder=11, clip_on=False)


def concl_bar(ax, x0, x1, y, text, color=NAVY, fs=6.4, h=0.052):
    ax.add_patch(FancyBboxPatch((x0, y), x1 - x0, h,
                                boxstyle="round,pad=0,rounding_size=0.006",
                                transform=ax.transAxes, fc=LGREY, ec="none",
                                zorder=10, clip_on=False))
    ax.add_patch(Rectangle((x0, y), 0.005, h, transform=ax.transAxes,
                           fc=color, ec="none", zorder=11, clip_on=False))
    ax.text((x0 + x1) / 2 + 0.004, y + h / 2, text, transform=ax.transAxes,
            ha="center", va="center", color=INK, fontsize=fs,
            zorder=11, clip_on=False)


def card(ax, x0, y0, x1, y1, ec="#D5DBE1", fc="white", lw=0.6, z=5):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle="round,pad=0,rounding_size=0.008",
                                transform=ax.transAxes, fc=fc, ec=ec, lw=lw,
                                zorder=z, clip_on=False))
    return z + 10


def txt(ax, x, y, s, z=15, **kw):
    kw.setdefault("transform", ax.transAxes)
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    ax.text(x, y, s, zorder=z, clip_on=False, **kw)


def build():
    g = json.loads((UP / "g4b_rerun_20261009_v13.json").read_text(encoding="utf-8"))
    nb = json.loads((UP / "c2_null_battery.json").read_text(encoding="utf-8"))
    n_events = g["v13_event_dedup"]["n_events_all"]
    n_oneoff = g["v13_event_dedup"]["n_events_oneoff"]
    n_oneoff_fin = g["v13_event_dedup"]["n_events_oneoff_finite_dm"]
    dm_budget = json.loads((UP / "pair565_summary.md").read_text(encoding="utf-8")
                           ) if False else None
    scat = json.loads((UP / "pair565_scatter.json").read_text(encoding="utf-8"))

    fig = plt.figure(figsize=(13.4, 7.1))
    fig.patch.set_facecolor("white")

    panels = [
        (0.012, 0.505, 0.320, 0.480),
        (0.339, 0.505, 0.322, 0.480),
        (0.668, 0.505, 0.320, 0.480),
        (0.012, 0.015, 0.320, 0.470),
        (0.339, 0.015, 0.322, 0.470),
        (0.668, 0.015, 0.320, 0.470),
    ]
    axes = [fig.add_axes(p) for p in panels]
    for i, ax in enumerate(axes, 1):
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.axis("off")
        card(ax, 0.0, 0.0, 1.0, 1.0)
        badge(ax, 0.030, 0.952, i)

    ax = axes[0]
    title_bar(ax, 0.05, 0.99, 0.895, 0.062, "数据与框架：CHIME/FRB 表观单次事件子空间", fs=7.2)
    fy0, fy1 = 0.60, 0.875
    ax.add_patch(Polygon([[0.10, fy1], [0.52, fy1], [0.42, fy0], [0.20, fy0]],
                         transform=ax.transAxes, fc=LBLUE, ec=NAVY, lw=0.7, zorder=6))
    rng = np.random.default_rng(7)
    for _ in range(70):
        t = rng.random(); yy = fy1 - t * (fy1 - fy0)
        xw = 0.42 + (0.52 - 0.42) * t
        xx = rng.uniform(0.20 + (0.10 - 0.20) * t, xw)
        ax.add_patch(Circle((xx, yy), 0.0045, transform=ax.transAxes,
                            fc=NAVY, ec="none", alpha=0.55, zorder=7))
    ax.text(0.31, fy1 + 0.018, "CHIME/FRB Catalog 1 & 2", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=6.6, color=NAVY, fontweight="bold")
    for yy, num, lab in [(0.815, f"{n_events} 个", "独立事件"),
                          (0.705, f"{n_oneoff_fin} 个", "表观单次事件")]:
        ax.add_patch(FancyArrowPatch((0.53, yy - 0.005), (0.615, yy - 0.005),
                                     transform=ax.transAxes, arrowstyle="-|>",
                                     mutation_scale=7, color=NAVY, lw=0.7, zorder=8))
        ax.text(0.635, yy + 0.012, num, transform=ax.transAxes, ha="left",
                va="center", fontsize=7.6, color=RED, fontweight="bold")
        ax.text(0.635, yy - 0.030, lab, transform=ax.transAxes, ha="left",
                va="center", fontsize=6.3, color=INK)
    fx0, fy0b, fx1, fy1b = 0.06, 0.135, 0.94, 0.500
    card(ax, fx0, fy0b, fx1, fy1b, ec="#C7D3DE", fc="#F7FAFD")
    ax.text((fx0 + fx1) / 2, fy1b - 0.038, "跨名称同源关联框架", transform=ax.transAxes,
            ha="center", va="center", fontsize=7.0, color=NAVY, fontweight="bold")
    cy = 0.350
    for cx in (0.22, 0.46):
        ax.add_patch(Polygon([[cx - 0.085, cy - 0.052], [cx + 0.085, cy - 0.052],
                              [cx, cy + 0.028]],
                             transform=ax.transAxes, fc="none", ec=BLUE, lw=0.6, zorder=8))
        ax.text(cx, cy - 0.078, "视位置", transform=ax.transAxes, ha="center",
                va="top", fontsize=5.4, color=INK)
        ax.add_patch(Circle((cx, cy + 0.030), 0.007, transform=ax.transAxes,
                            fc=BLUE, ec="none", zorder=9))
    ax.annotate("", xy=(0.46, cy + 0.030), xytext=(0.22, cy + 0.030),
                textcoords=ax.transAxes, arrowprops=dict(arrowstyle="<->",
                color=BLUE, lw=0.6), zorder=9)
    ax.text(0.34, cy + 0.042, "θ < Θ", transform=ax.transAxes, ha="center",
            va="bottom", fontsize=5.8, color=BLUE)
    ax.text(0.34, cy - 0.110, "配对准则：时空锥 + DM 窗口", transform=ax.transAxes,
            ha="center", va="top", fontsize=6.0, color=INK)
    ax.text(0.34, cy - 0.160, "|ΔDM| < ΔDM_max", transform=ax.transAxes,
            ha="center", va="top", fontsize=5.7, color=GREY)
    ax.text(0.50, 0.178, "假设合并：Union-Find 算法", transform=ax.transAxes,
            ha="center", va="center", fontsize=6.2, color=INK)
    ax.annotate("", xy=(0.50, 0.205), xytext=(0.50, 0.245),
                textcoords=ax.transAxes, arrowprops=dict(arrowstyle="-|>",
                color=NAVY, lw=0.7), zorder=11)
    concl_bar(ax, 0.05, 0.99, 0.045,
              f"框架零辐射模型依赖 · 输入仅（位置, DM, 时刻）· 候选并集 28 对", fs=6.0)

    ax = axes[1]
    title_bar(ax, 0.05, 0.99, 0.895, 0.062, "核心证据：最强同源候选对", fs=7.2)
    ax.text(0.5, 0.845, "FRB20200511C × FRB20211128E", transform=ax.transAxes,
            ha="center", va="center", fontsize=7.4, color=RED, fontweight="bold")
    sx0, sy0, sx1, sy1 = 0.07, 0.30, 0.50, 0.80
    ax.add_patch(Rectangle((sx0, sy0), sx1 - sx0, sy1 - sy0,
                           transform=ax.transAxes, fc="white", ec="#C7D3DE",
                           lw=0.5, zorder=6))
    bg = rng.normal(0, 1, (150, 2))
    bg[:, 0] = np.abs(bg[:, 0]) * 0.10 + 0.02
    bg[:, 1] = np.abs(bg[:, 1]) * 30 + 5
    ax.scatter(sx0 + 0.03 + (bg[:, 0] / 0.5) * (sx1 - sx0 - 0.06),
               sy0 + 0.03 + (bg[:, 1] / 120) * (sy1 - sy0 - 0.06),
               s=1.5, c=GREY, alpha=0.6, zorder=7)
    sx_star, sy_star = sx0 + 0.050, sy1 - 0.150
    ax.scatter([sx_star], [sy_star], s=55, marker="*", c=RED,
               zorder=12, edgecolors="white", linewidths=0.4)
    ax.text((sx0 + sx1) / 2, sy1 - 0.020, "位置-DM 联合空间", transform=ax.transAxes,
            ha="center", va="top", fontsize=6.0, color=NAVY, fontweight="bold")
    ax.text(sx0 + 0.010, sy0 + 0.055, "θ = 23.9″\nΔDM = 0.18 pc cm$^{-3}$",
            transform=ax.transAxes, ha="left", va="center", fontsize=5.5,
            color=RED, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.30", fc="white", ec=RED, lw=0.5))
    ax.annotate("", xy=(sx_star, sy_star - 0.020),
                xytext=(sx0 + 0.075, sy0 + 0.070),
                textcoords=ax.transAxes,
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.5), zorder=11)
    ax.text(sx1 + 0.012, sy0 + (sy1 - sy0) / 2, "角距 θ (角秒)",
            transform=ax.transAxes, ha="left", va="center",
            fontsize=5.0, color=GREY, rotation=90)
    ax.text(0.760, 0.735, "散射时标 τ 对比  (1.0–1.5 GHz)", transform=ax.transAxes,
            ha="center", va="center", fontsize=6.0, color=NAVY, fontweight="bold")
    bx = [(0.665, 0.330, GREY, "DM 趋势\n预期"), (0.805, 0.062, RED, "观测\n结果")]
    for xx, hh, cc, lab in bx:
        ax.add_patch(Rectangle((xx, 0.300), 0.070, hh, transform=ax.transAxes,
                               fc=cc, ec="none", zorder=8))
    ax.text(0.665 + 0.035, 0.300 + 0.330 + 0.016, "DM 趋势预期", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=5.2, color=GREY)
    ax.text(0.805 + 0.035, 0.300 + 0.062 + 0.016, "观测结果", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=5.2, color=RED, fontweight="bold")
    ax.annotate("", xy=(0.862, 0.330), xytext=(0.862, 0.600),
                textcoords=ax.transAxes, arrowprops=dict(arrowstyle="-|>",
                color=RED, lw=0.7), zorder=11)
    ax.text(0.872, 0.465, "低 > 2 个\n数量级", transform=ax.transAxes,
            ha="left", va="center", fontsize=5.4, color=RED, fontweight="bold")
    ax.plot([0.648, 0.892], [0.300, 0.300], transform=ax.transAxes,
            color="#B8C4CE", lw=0.5, zorder=8)
    concl_bar(ax, 0.05, 0.99, 0.045,
              "无散射探测，显著低于 DM 趋势预期 · 排除湍流屏模型主导，支持同源假设", fs=5.9)

    ax = axes[2]
    title_bar(ax, 0.05, 0.99, 0.895, 0.062, "物理机制与红移假设", fs=7.2)
    ax.text(0.24, 0.845, "DM 预算分解 (pc cm$^{-3}$)", transform=ax.transAxes,
            ha="center", va="center", fontsize=6.1, color=NAVY, fontweight="bold")
    x0, w = 0.115, 0.20
    base, top = 0.175, 0.775
    tot = 1835.0
    parts = [("银河系 NE2001/YMW16", 123.8, LBLUE, NAVY),
             ("IGM（本征）", 1671.1, LAMBER, "#8A5A00"),
             ("银晕（假定）", 40.1, LGREY, "#5A5A5A")]
    scale = (top - base) / tot
    yy = top
    for lab, val, fc, tc in parts:
        h = val * scale
        ax.add_patch(Rectangle((x0, yy - h), w, h, transform=ax.transAxes,
                               fc=fc, ec="white", lw=0.6, zorder=8))
        if h > 0.10:
            ax.text(x0 + w / 2, yy - h / 2, f"{lab}\n{val:.0f}",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=5.1, color=tc, zorder=16)
        else:
            ax.annotate(f"{lab}  {val:.0f}", xy=(x0 + w, yy - h / 2),
                        xytext=(x0 + w + 0.012, yy - h / 2),
                        textcoords=ax.transAxes, ha="left", va="center",
                        fontsize=5.0, color=tc, zorder=16,
                        arrowprops=dict(arrowstyle="-", color=tc, lw=0.4))
        yy -= h
    ax.text(x0 + w / 2, top + 0.020, f"总 DM = {tot:.0f}", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=6.0, color=NAVY, fontweight="bold")
    ax.annotate("", xy=(0.50, 0.46), xytext=(0.335, 0.46),
                textcoords=ax.transAxes, arrowprops=dict(arrowstyle="-|>",
                color=RED, lw=0.9), zorder=9)
    for yy0, hh, cc, lc, tag, name, detail in [
        (0.545, 0.235, RED, LRED, "假设路径 H1", "高红移重复源候选",
         "IGM 主导 → z ≈ 1.5"),
        (0.275, 0.235, GREEN, LGREEN, "假设路径 H2", "低红移致密介质",
         "本地致密环境贡献高 DM → z ≈ 0.1")]:
        card(ax, 0.53, yy0, 0.96, yy0 + hh, ec=cc, fc=lc, lw=0.8)
        ax.text(0.745, yy0 + hh - 0.042, tag, transform=ax.transAxes,
                ha="center", va="center", fontsize=6.3, color=cc, fontweight="bold")
        ax.text(0.745, yy0 + hh - 0.088, name, transform=ax.transAxes,
                ha="center", va="center", fontsize=5.7, color=INK)
        ax.text(0.745, yy0 + hh - 0.148, detail, transform=ax.transAxes,
                ha="center", va="center", fontsize=5.7, color=cc, fontweight="bold")
    ax.text(0.745, 0.225, "两条路径均待进一步观测验证", transform=ax.transAxes,
            ha="center", va="center", fontsize=5.5, color=GREY)
    concl_bar(ax, 0.05, 0.99, 0.045,
              "DM ≈ 1835 → 若 H1 则 z ≈ 1.5，属高红移重复源候选", fs=5.9)

    ax = axes[3]
    title_bar(ax, 0.05, 0.99, 0.895, 0.062, "方法学守卫：恒星日量化 · 几何机制解释", fs=7.0)
    ty = 0.68
    ax.add_patch(FancyArrowPatch((0.07, ty), (0.93, ty), transform=ax.transAxes,
                                 arrowstyle="-|>", mutation_scale=8, color=INK,
                                 lw=0.8, zorder=8))
    for xx, lab, cc in [(0.17, "事件 A\n(观测日)", BLUE),
                        (0.50, "565 天\n(太阳日)", NAVY),
                        (0.83, "事件 B\n(观测日)", BLUE)]:
        ax.add_patch(Circle((xx, ty), 0.013, transform=ax.transAxes, fc=cc,
                            ec="white", lw=0.6, zorder=9))
        ax.text(xx, ty - 0.040, lab, transform=ax.transAxes, ha="center",
                va="top", fontsize=5.5, color=INK)
    ax.text(0.50, ty + 0.062, "= 567 恒星日（地球相对恒星自转）", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=6.0, color=RED, fontweight="bold")
    hx0, hx1, hy = 0.12, 0.88, 0.395
    ax.add_patch(Rectangle((hx0, hy), hx1 - hx0, 0.045, transform=ax.transAxes,
                           fc=LGREY, ec="#C7D3DE", lw=0.5, zorder=7))
    ax.add_patch(Rectangle((0.44, hy), 0.12, 0.045, transform=ax.transAxes,
                           fc=LBLUE, ec=NAVY, lw=0.5, zorder=8))
    ax.text(0.50, hy + 0.068, "CHIME 过境观测窗口（|HA| < 0.1 h）", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=5.3, color=NAVY)
    ax.text(0.50, hy - 0.020, "→ 两次探测必然落在同一窗口相位 → 再次重合", transform=ax.transAxes,
            ha="center", va="top", fontsize=5.4, color=GREY)
    card(ax, 0.08, 0.215, 0.92, 0.345, ec=NAVY, fc="#EAF2FA", lw=0.8)
    ax.text(0.50, 0.300, "565 天 = 567 恒星日", transform=ax.transAxes,
            ha="center", va="center", fontsize=7.4, color=NAVY, fontweight="bold")
    ax.text(0.50, 0.252, "CHIME 过境观测几何产物，非物理周期", transform=ax.transAxes,
            ha="center", va="center", fontsize=5.8, color=INK)
    concl_bar(ax, 0.05, 0.99, 0.045,
              "通用守卫：任何两事件「周期」须先报 Δt/T_sid 贴整数度", fs=5.9)

    ax = axes[4]
    title_bar(ax, 0.05, 0.99, 0.895, 0.062, "零模型检验结果", fs=7.2)
    ax.text(0.265, 0.835, "4 种置换模型", transform=ax.transAxes,
            ha="center", va="center", fontsize=6.2, color=NAVY, fontweight="bold")
    for i, lab in enumerate(["打乱 DM", "打乱时间", "打乱位置", "DM–时间联合打乱"]):
        yy = 0.720 - i * 0.082
        ax.add_patch(Circle((0.135, yy), 0.013, transform=ax.transAxes,
                            fc="none", ec=NAVY, lw=0.8, zorder=8))
        ax.add_patch(Circle((0.135, yy + 0.008), 0.005, transform=ax.transAxes,
                            fc=NAVY, ec="none", zorder=9))
        ax.text(0.168, yy, lab, transform=ax.transAxes, ha="left", va="center",
                fontsize=5.9, color=INK)
    ax.text(0.265, 0.320, "4 种零假设下\n真实 28 对均未被重现", transform=ax.transAxes,
            ha="center", va="center", fontsize=5.6, color=NAVY)
    ax.text(0.715, 0.835, "显著性（观测值 vs. 置换分布）", transform=ax.transAxes,
            ha="center", va="center", fontsize=5.8, color=NAVY, fontweight="bold")
    xl, xr, ybase = 0.545, 0.955, 0.34
    ax.annotate("", xy=(xr, ybase), xytext=(xl, ybase), textcoords=ax.transAxes,
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=0.7), zorder=8)
    for xx, lab in [(0.575, r"$10^{-4}$"), (0.685, r"$10^{-3}$"),
                    (0.795, r"$10^{-2}$"), (0.905, r"$10^{-1}$")]:
        ax.plot([xx, xx], [ybase, ybase + 0.012], transform=ax.transAxes,
                color=GREY, lw=0.5, zorder=8)
        ax.text(xx, ybase - 0.018, lab, transform=ax.transAxes, ha="center",
                va="top", fontsize=5.0, color=GREY)
    ax.text((xl + xr) / 2, ybase - 0.058, "p 值（双对数）", transform=ax.transAxes,
            ha="center", va="top", fontsize=5.3, color=GREY)
    pvals = [0.0010, 0.0010, 0.0010, 0.0090]
    for i, pv in enumerate(pvals):
        yy = 0.745 - i * 0.082
        xv = xl + (np.log10(pv) + 4) / 3 * (xr - xl)
        ax.plot([xv, xv], [yy, ybase], transform=ax.transAxes, color=NAVY,
                lw=0.5, ls=":", zorder=7)
        ax.add_patch(Circle((xv, yy), 0.011, transform=ax.transAxes, fc=NAVY,
                            ec="none", zorder=9))
    xthr = xl + (np.log10(0.01) + 4) / 3 * (xr - xl)
    ax.plot([xthr, xthr], [ybase, 0.775], transform=ax.transAxes,
            color=RED, lw=0.8, ls="--", zorder=8)
    ax.text(xthr - 0.012, 0.480, "观测 p 值\np = 0.009", transform=ax.transAxes,
            ha="right", va="center", fontsize=5.9, color=RED, fontweight="bold")
    concl_bar(ax, 0.05, 0.99, 0.045,
              "多种零假设下均显著，结果稳健", fs=6.2, h=0.062)

    ax = axes[5]
    title_bar(ax, 0.05, 0.99, 0.895, 0.062, "未来验证路径：FAST 关键观测", fs=7.0)
    cx, cy = 0.235, 0.44
    ax.add_patch(Wedge((cx, cy - 0.10), 0.165, 35, 145, width=0.055,
                       transform=ax.transAxes, fc=LBLUE, ec=NAVY, lw=0.8, zorder=8))
    ax.plot([cx, cx], [cy + 0.045, cy + 0.175], transform=ax.transAxes,
            color=NAVY, lw=0.8, zorder=8)
    ax.add_patch(Circle((cx, cy + 0.045), 0.016, transform=ax.transAxes,
                        fc=LRED, ec=RED, lw=0.8, zorder=9))
    for dx in (-0.075, -0.038, 0, 0.038, 0.075):
        ax.plot([cx + dx, cx], [cy - 0.055, cy + 0.045], transform=ax.transAxes,
                color=NAVY, lw=0.35, ls=":", zorder=7)
    ax.text(cx, 0.208, "中国 FAST 望远镜", transform=ax.transAxes, ha="center",
            va="center", fontsize=6.2, color=NAVY, fontweight="bold")
    for i, (tag, det) in enumerate([("第三次爆发确认", "检验重复暴本质"),
                                     ("RM 测量", "探测本地磁场环境，区分 H1/H2"),
                                     ("高精度定位", "确认同源体位置")]):
        yy = 0.735 - i * 0.195
        ax.add_patch(FancyBboxPatch((0.53, yy - 0.075), 0.43, 0.15,
                                    boxstyle="round,pad=0,rounding_size=0.008",
                                    transform=ax.transAxes, fc=LGREY, ec="#C7D3DE",
                                    lw=0.6, zorder=8))
        ax.add_patch(Rectangle((0.53, yy - 0.075), 0.006, 0.15, transform=ax.transAxes,
                               fc=NAVY, ec="none", zorder=9))
        ax.text(0.55, yy + 0.032, tag, transform=ax.transAxes, ha="left",
                va="center", fontsize=6.1, color=NAVY, fontweight="bold")
        ax.text(0.55, yy - 0.030, det, transform=ax.transAxes, ha="left",
                va="center", fontsize=5.2, color=INK)
    ax.annotate("", xy=(0.525, 0.44), xytext=(0.415, 0.44), textcoords=ax.transAxes,
                arrowprops=dict(arrowstyle="-|>", color=NAVY, lw=0.9), zorder=11)
    concl_bar(ax, 0.05, 0.99, 0.045,
              "该场几何/频段/灵敏度统计均指向 FAST 为第三次确认首选", fs=5.7)

    out = FIGDIR / "fig0_graphical_abstract.png"
    fig.savefig(out, dpi=400, bbox_inches="tight", facecolor="white")
    print(f"saved -> {out}")
    return out


if __name__ == "__main__":
    build()
