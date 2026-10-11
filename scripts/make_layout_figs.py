#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Native vector redraw of the dual and geometry layout figures.

Usage: python scripts/make_layout_figs.py [--outdir DIR]
Output: docs/figures/vector2/{dual,geometry}.pdf
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch
from fig_gate import report_collisions as _shared_gate


def report_collisions(fig, name):
    """Shared layout gate, made fatal: an overprinted figure is not a draft figure."""
    n = _shared_gate(fig, name)
    if n:
        raise SystemExit(f"{name}: {n} 处叠字/越出画布/压数据，图不出文件")
    return n

ROOT = Path(__file__).resolve().parents[1]
CAT2 = ROOT / "data/interim/cat2_bursts.parquet"
SCAT = ROOT / "data/processed/upgrade/pair565_scatter.json"
OUT_DEFAULT = ROOT / "docs/figures/vector2"

CM = 2.54
TEXT_W_IN = 15.0 / CM
H_IN = 8.2 / CM

C_BLUE = "#2b5b98"
C_NAVY = "#1f3a5f"
C_LBLUE = "#5b8dbf"
C_GREY = "#8a8a8a"
C_RED = "#c8282e"
C_GOLD = "#c8862a"
C_SLATE = "#5f7082"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7.9,
    "axes.labelsize": 7.5,
    "axes.titlesize": 8.0,
    "xtick.labelsize": 7.3,
    "ytick.labelsize": 7.3,
    "legend.fontsize": 7.3,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.6,
    "ytick.major.size": 2.6,
    "savefig.bbox": None,
})


def load_scat():
    j = json.load(open(SCAT, encoding="utf-8"))
    c2 = pd.read_parquet(CAT2)
    b = c2[c2.sub_num == 0].copy()
    meas = b[(b.scat_time > 0) & b.scat_time_err.notna()]
    exp = {"burst_level": len(b), "tau_measured": len(meas)}
    got = {"burst_level": j["population"]["burst_level"],
           "tau_measured": j["population"]["tau_measured"]}
    assert exp == got, f"样本选择与产物不一致：本脚本 {exp} vs 产物 {got}"
    sc = j["scaling"]
    assert abs(sc["b"] - 0.5886244847276179) < 1e-9
    return j, b, meas


def panel_scattering(ax, j, meas):
    sc = j["scaling"]
    a, bb, rstd = sc["a"], sc["b"], sc["resid_std_dex"]
    lo, hi = sc["dm_range_of_fit"]
    x = np.logspace(np.log10(lo), np.log10(hi), 200)
    y = 10.0 ** (a + bb * np.log10(x)) * 1e3

    ax.scatter(meas.dm_fitb.astype(float), meas.scat_time.astype(float) * 1e3,
               s=5, color=C_GREY, alpha=0.42, linewidths=0,
               label=f"measured $\\tau$ (n = {len(meas)})")
    ax.fill_between(x, y * 10 ** (-rstd), y * 10 ** rstd, color=C_LBLUE,
                    alpha=0.20, linewidth=0, label=f"$\\pm1\\sigma$ ({rstd:.2f} dex)")
    ax.plot(x, y, color="black", lw=1.0, label="catalogue fit")

    dm_pair = 1835.0
    floor_ms = sc["tau_floor_used_ms"]
    pred_ms = sc["tau_predicted_at_1835_ms"]
    ax.plot([dm_pair], [pred_ms], marker="o", ms=4.5, color=C_GOLD, zorder=6,
            label=f"fit predicts {pred_ms:.1f} ms")
    ax.annotate(f"$R \\approx {sc['ratio_R_predicted_over_floor']:.0f}\\times$ above floor",
                xy=(dm_pair, pred_ms), xytext=(0.03, 0.60),
                textcoords=ax.transAxes, fontsize=7.55, color=C_GOLD,
                ha="left", va="top",
                arrowprops=dict(arrowstyle="-|>", lw=0.7, color=C_GOLD,
                                connectionstyle="arc3,rad=-0.16"))
    ax.annotate("", xy=(dm_pair, floor_ms * 0.55), xytext=(dm_pair, floor_ms * 3.2),
                arrowprops=dict(arrowstyle="-|>", lw=1.1, color=C_RED))
    ax.plot([dm_pair], [floor_ms], marker="*", ms=11, color=C_RED, zorder=7,
            label=f"this pair: $\\tau < {floor_ms:.3f}$ ms")
    ax.annotate("no measurable scattering:\nturbulence-screen H2 disfavoured",
                xy=(dm_pair, floor_ms * 1.1), xytext=(0.03, 0.38),
                textcoords=ax.transAxes,
                fontsize=7.55, color=C_RED, ha="left", va="top", linespacing=1.4,
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none", alpha=0.85),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=C_RED,
                                connectionstyle="arc3,rad=0.16"))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1.0, 5.0e3); ax.set_ylim(1.2e-2, 4.0e4)
    ax.set_xlabel("DM (pc cm$^{-3}$)")
    ax.set_ylabel("$\\tau$ (ms)")
    ax.set_title("Scattering time vs. DM", color=C_NAVY, pad=4)
    leg = ax.legend(loc="upper left", frameon=False, handletextpad=0.4,
                    borderaxespad=0.3, labelspacing=0.32)
    leg.set_zorder(10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(direction="out")


def panel_budget(ax, j):
    dm_mw, dm_halo, dm_igm, dm_obs = 123.8, 40.0, 1671.1, 1834.95
    total = dm_mw + dm_halo + dm_igm
    assert abs(total - dm_obs) < 0.2, f"预算分解与总 DM 不符：{total} vs {dm_obs}"
    segs = [("Galactic", dm_mw, C_BLUE), ("Halo", dm_halo, C_SLATE),
            ("IGM", dm_igm, C_GOLD)]
    y0 = 0.0
    for name, val, col in segs:
        ax.bar([0.0], [val], bottom=y0, width=0.42, color=col, linewidth=0)
        yl = y0 + val / 2
        ax.annotate(f"{name} {val:.0f} ({val / dm_obs * 100:.1f}%)",
                    xy=(0.215, yl), xytext=(0.28, yl), ha="left", va="center",
                    fontsize=7.45, color=C_NAVY,
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=C_SLATE))
        y0 += val
    ax.text(0.0, dm_obs * 1.045, f"total DM {dm_obs:.2f}", ha="center", va="bottom",
            fontsize=7.55, color=C_NAVY)

    ax.text(0.86, dm_obs * 0.93, "H1: $z \\approx 1.5$\nhigh-redshift repeater",
            ha="center", va="center", fontsize=7.55, color=C_NAVY, linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.42", fc="#e8eef7", ec=C_LBLUE, lw=0.7))
    ax.text(0.86, dm_obs * 0.66, "H2: low-redshift\ndense local medium",
            ha="center", va="center", fontsize=7.55, color=C_NAVY, linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.42", fc="#eef0f2", ec=C_SLATE, lw=0.7))
    ax.add_patch(FancyArrowPatch((0.22, dm_obs * 0.78), (0.44, dm_obs * 0.90),
                                 arrowstyle="-|>", mutation_scale=7,
                                 lw=0.8, color=C_LBLUE, linestyle="--"))
    ax.add_patch(FancyArrowPatch((0.22, dm_obs * 0.66), (0.44, dm_obs * 0.66),
                                 arrowstyle="-|>", mutation_scale=7,
                                 lw=0.8, color=C_SLATE, linestyle="--"))
    ax.set_xlim(-1.05, 1.62); ax.set_ylim(0, dm_obs * 1.16)
    ax.set_xticks([])
    ax.set_ylabel("DM (pc cm$^{-3}$)")
    ax.set_title("DM Budget Decomposition", color=C_NAVY, pad=4)
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(direction="out")


def draw_dual(j, meas, outdir: Path) -> Path:
    fig = plt.figure(figsize=(TEXT_W_IN, H_IN))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 0.82],
                          left=0.075, right=0.985, top=0.945, bottom=0.115, wspace=0.30)
    panel_scattering(fig.add_subplot(gs[0, 0]), j, meas)
    panel_budget(fig.add_subplot(gs[0, 1]), j)
    report_collisions(fig, "dual")
    p = outdir / "dual.pdf"
    fig.savefig(p, format="pdf")
    plt.close(fig)
    return p


def panel_time(ax):
    dt_days, n_sid = 565.4537, 567.00191
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.02, 0.525), 0.96, 0.46,
                                boxstyle="round,pad=0.012,rounding_size=0.03",
                                fc="#eef4fb", ec=C_LBLUE, lw=0.8))
    ax.text(0.50, 0.925, "Time relation", ha="center", color=C_NAVY,
            fontsize=7.5, weight="bold")
    ax.text(0.50, 0.815, f"$\\Delta t = {dt_days:.4f}$", ha="center",
            fontsize=7.9, color=C_NAVY)
    ax.text(0.50, 0.745, "solar days", ha="center", fontsize=7.45, color=C_SLATE)
    ax.text(0.50, 0.645, f"$= {n_sid:.5f}$", ha="center",
            fontsize=7.9, color=C_NAVY, weight="bold")
    ax.text(0.50, 0.605, "sidereal days", ha="center", fontsize=7.45, color=C_SLATE)
    ax.text(0.50, 0.556, "(164 s from an exact integer)", ha="center",
            fontsize=7.3, color=C_SLATE)
    ax.add_patch(FancyBboxPatch((0.02, 0.02), 0.96, 0.47,
                                boxstyle="round,pad=0.012,rounding_size=0.03",
                                fc="#f6f7f8", ec=C_SLATE, lw=0.7))
    ax.text(0.50, 0.465, "Geometric lock", ha="center", color=C_NAVY,
            fontsize=7.5, weight="bold")
    ax.text(0.50, 0.245,
            "A drift-scan catalogue\nrecords every event within\n"
            "$|\\mathrm{HA}| < 0.5$ h of its\nown transit, so two\n"
            "detections differ by an\ninteger number of sidereal\n"
            "days regardless of\nthe source.",
            ha="center", va="center", fontsize=7.45, color=C_SLATE, linespacing=1.45)


def panel_geometry(ax):
    ax.set_xlim(-1.70, 1.70); ax.set_ylim(-1.42, 2.95)
    ax.set_aspect("equal"); ax.axis("off")
    ax.text(0, 2.80, "CHIME Fixed Zenith Telescope", ha="center",
            color=C_NAVY, fontsize=7.8, weight="bold")
    ax.text(0, 2.53, "geometry with source declination $\\delta$", ha="center",
            fontsize=7.7, color=C_NAVY)

    earth = Circle((0, 0), 0.62, fc="#dceaf6", ec=C_BLUE, lw=0.9)
    ax.add_patch(earth)
    ax.add_patch(Ellipse((0, 0), 1.24, 0.34, fc="none", ec=C_SLATE, lw=0.7, ls="--"))
    ax.annotate("", xy=(0, 1.16), xytext=(0, -1.16),
                arrowprops=dict(arrowstyle="-", lw=0.8, color=C_NAVY, ls=(0, (4, 3))))
    ax.text(0.0, 1.22, "N", ha="center", fontsize=7.9, color=C_NAVY, weight="bold")
    ax.text(0.0, -1.34, "S", ha="center", fontsize=7.9, color=C_NAVY, weight="bold")
    ax.add_patch(FancyArrowPatch((-0.22, 1.06), (0.22, 1.06), connectionstyle="arc3,rad=-0.5",
                                 arrowstyle="-|>", mutation_scale=7, lw=0.9, color=C_NAVY))
    ax.add_patch(FancyArrowPatch((-1.34, -0.14), (-0.92, -0.28), connectionstyle="arc3,rad=0.3",
                                 arrowstyle="-|>", mutation_scale=8, lw=0.9, color=C_SLATE))
    ax.add_patch(FancyArrowPatch((0.92, 0.28), (1.34, 0.14), connectionstyle="arc3,rad=0.3",
                                 arrowstyle="-|>", mutation_scale=8, lw=0.9, color=C_SLATE))
    ax.text(-1.66, -0.86, "Earth's rotation\n(1 sidereal day)", ha="left",
            fontsize=7.45, color=C_SLATE, linespacing=1.4)

    obs = np.array([-0.30, 0.545])
    zen = obs + np.array([0.0, 0.60])
    ang = np.radians(30.0)
    src_dir = np.array([np.sin(ang), np.cos(ang)])
    src = zen + src_dir * 0.86
    ax.plot(*obs, marker="o", ms=3.2, color=C_NAVY, zorder=6)
    ax.annotate("", xy=zen, xytext=obs,
                arrowprops=dict(arrowstyle="-|>", lw=0.9, color=C_NAVY))
    ax.text(zen[0] - 0.06, zen[1] + 0.06, "Zenith", ha="right", fontsize=7.45, color=C_SLATE)
    ax.annotate("", xy=src, xytext=zen,
                arrowprops=dict(arrowstyle="-|>", lw=1.0, color=C_BLUE))
    ax.plot(*src, marker="*", ms=10, color=C_BLUE, zorder=6)
    ax.text(src[0] + 0.10, src[1] + 0.02, "Source\n(declination $\\delta$)",
            ha="left", va="center", fontsize=7.45, color=C_SLATE, linespacing=1.35)
    th = np.linspace(0.0, ang, 24)
    ax.plot(zen[0] + 0.34 * np.sin(th), zen[1] + 0.34 * np.cos(th), color=C_LBLUE, lw=0.9)
    ax.text(zen[0] + 0.15, zen[1] + 0.36, "$\\delta$", fontsize=7.5, color=C_LBLUE)
    ax.annotate("fixed zenith beam\n(transit scan)", xy=(-0.05, 0.30),
                xytext=(1.66, 0.72), fontsize=7.45, color=C_SLATE, ha="right",
                arrowprops=dict(arrowstyle="-|>", lw=0.6, color=C_SLATE,
                                connectionstyle="arc3,rad=-0.2"))


def panel_discriminants(ax):
    rows = [
        ("Phase stability", "drifts", "stable"),
        ("Multi-source", "source-dep.", "universal"),
        ("Period 565 vs 567", "not integer", "integer"),
        ("Rotation link", "none", "direct"),
        ("Predictive power", "limited", "high"),
        ("Instrument robustness", "sensitive", "robust"),
        ("Sidereal/solar link", "no", "yes"),
    ]
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.0, 0.985, "Key Discriminants", color=C_NAVY, fontsize=7.8,
            weight="bold", va="top")
    xs = (0.0, 0.455, 0.735)
    heads = ("Feature", "Period", "Geometric")
    top, rowh = 0.905, 0.106
    ax.add_patch(plt.Rectangle((0.0, top - rowh), 1.0, rowh, fc="#eef4fb",
                               ec="none", linewidth=0))
    for x, h in zip(xs, heads):
        ax.text(x + 0.014, top - rowh / 2, h, fontsize=7.45, color=C_NAVY,
                weight="bold", va="center")
    for i, (feat, phys, geo) in enumerate(rows):
        y = top - rowh * (i + 1)
        if i % 2 == 0:
            ax.add_patch(plt.Rectangle((0.0, y - rowh), 1.0, rowh, fc="#f8fafc",
                                       ec="none", linewidth=0))
        ax.text(xs[0] + 0.014, y - rowh / 2, feat, fontsize=7.35, va="center",
                color=C_NAVY)
        ax.text(xs[1] + 0.014, y - rowh / 2, phys, fontsize=7.35, va="center",
                color=C_SLATE)
        ax.text(xs[2] + 0.014, y - rowh / 2, "\u2713 " + geo, fontsize=7.35,
                va="center", color=C_BLUE)
    ax.plot([0, 1, 1, 0, 0], [top, top, top - rowh * 8, top - rowh * 8, top],
            color=C_LBLUE, lw=0.7)
    for x in xs[1:]:
        ax.plot([x, x], [top, top - rowh * 8], color=C_LBLUE, lw=0.5)


def draw_geometry(outdir: Path) -> Path:
    fig = plt.figure(figsize=(TEXT_W_IN, H_IN))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.85, 1.10, 1.45],
                          left=0.018, right=0.988, top=0.985, bottom=0.015, wspace=0.10)
    panel_time(fig.add_subplot(gs[0, 0]))
    panel_geometry(fig.add_subplot(gs[0, 1]))
    panel_discriminants(fig.add_subplot(gs[0, 2]))
    report_collisions(fig, "geometry")
    p = outdir / "geometry.pdf"
    fig.savefig(p, format="pdf")
    plt.close(fig)
    return p


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(OUT_DEFAULT))
    a = ap.parse_args()
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    j, b, meas = load_scat()
    p1 = draw_dual(j, meas, outdir)
    p2 = draw_geometry(outdir)
    for p in (p1, p2):
        print(f"  {p.name:14s} {p.stat().st_size/1024:6.0f} KB  画布 "
              f"{TEXT_W_IN*CM:.1f} x {H_IN*CM:.1f} cm（设计宽；RAA 版心 14.5 cm ⇒ 印 "
              f"{TEXT_W_IN*CM*14.5/15.0:.1f} x {H_IN*CM*14.5/15.0:.1f} cm，缩放 0.967）")
    print(f"\nτ–DM 点云用真实数据：n = {len(meas)} 可测暴（占暴级 {j['population']['burst_level']} 的 "
          f"{len(meas)/j['population']['burst_level']*100:.1f}%），拟合 a/b/残差回指 pair565_scatter.json")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    raise SystemExit(main())
