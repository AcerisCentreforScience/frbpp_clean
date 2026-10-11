"""Native vector redraw of the scatter, flow and roadmap layout figures.

Usage: python scripts/make_layout_figs2.py
Output: docs/figures/vector2/{scatter,flow,roadmap}.pdf
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon
from fig_gate import report_collisions

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from make_layout_figs import (CM, H_IN, TEXT_W_IN, C_BLUE, C_GOLD, C_LBLUE,  # noqa: E402
                              C_NAVY, C_RED, C_SLATE)

ROOT = Path(__file__).resolve().parents[1]
ISO = ROOT / "data/processed/upgrade/c2_isolation_figure_v1.json"
OUTDIR = ROOT / "docs/figures/vector2"


def load_iso():
    return json.loads(ISO.read_text(encoding="utf-8"))


def panel_nn(ax, j):
    ec = np.array([[d["theta"], d["frac"]] for d in j["nn_ecdf"]], float)
    ax.plot(ec[:, 0], ec[:, 1] * 100, color=C_BLUE, lw=1.1,
            label="nearest other event")
    nn = j["nearest_neighbour"]
    pair, second = nn["member_1nn_arcsec"], nn["member_2nn_arcsec"]
    ax.axvline(pair, color=C_RED, lw=0.8, ls="--")
    ax.axvline(second, color=C_GOLD, lw=0.8, ls=":")
    f_pair = np.interp(pair, ec[:, 0], ec[:, 1]) * 100
    ax.scatter([pair], [f_pair], s=26, color=C_RED, zorder=6)
    ax.annotate(f"this pair\n{pair:.1f}″", xy=(pair, f_pair),
                xytext=(7.0, 34), fontsize=7.45, color=C_RED, ha="right",
                arrowprops=dict(arrowstyle="-|>", lw=0.7, color=C_RED))
    ax.annotate(f"next-nearest for each member\n{second:.0f}″ "
                f"({nn['ratio_2nn_over_pair']:.0f}$\\times$ the pair)",
                xy=(second, 99), xytext=(second * 0.10, 66), fontsize=7.45,
                color=C_GOLD, ha="center",
                arrowprops=dict(arrowstyle="-|>", lw=0.7, color=C_GOLD))
    ax.text(0.98, 0.06, f"{nn['frac_events_nn_gt_pair']*100:.1f}% of one-off events have a\n"
                        f"nearest neighbour farther than\nthis pair's separation",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.45, color=C_SLATE,
            linespacing=1.4,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none", alpha=0.82))
    ax.set_xscale("log"); ax.set_xlim(1, 3e5); ax.set_ylim(0, 100)
    ax.set_xlabel("nearest-neighbour separation (arcsec)")
    ax.set_ylabel("cumulative fraction (%)")
    ax.set_title("Position isolation", color=C_NAVY, pad=4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_count(ax, j):
    cu = np.array([[d["theta"], d["n_all"], d["n_cross"], d["n_cross_dm"]]
                   for d in j["cross_day_curve"]], float)
    ax.plot(cu[:, 0], np.maximum(cu[:, 1], 0.5), color=C_SLATE, lw=0.9, ls="--",
            label="all pairs")
    ax.plot(cu[:, 0], np.maximum(cu[:, 2], 0.5), color=C_BLUE, lw=1.1,
            label="cross-day ($\\Delta t > 1$ d)")
    ax.plot(cu[:, 0], np.maximum(cu[:, 3], 0.5), color=C_RED, lw=1.1,
            label="cross-day and $|\\Delta\\mathrm{DM}| < 0.5$")
    box = j["isolation_box"]
    ax.scatter([30], [box["n_pairs_in_box_all"]], s=26, color=C_SLATE, zorder=6)
    ax.scatter([30], [box["n_pairs_in_box_cross_day"]], s=26, color=C_RED, zorder=6)
    ax.annotate(f"{box['n_pairs_in_box_all']} pairs", xy=(30, box["n_pairs_in_box_all"]),
                xytext=(56, 120), fontsize=7.45, color=C_SLATE,
                arrowprops=dict(arrowstyle="-", lw=0.6, color=C_SLATE))
    ax.annotate(f"{box['n_pairs_in_box_cross_day']} pairs\n(this pair and its\n"
                "113.7 d companion)", xy=(30, box["n_pairs_in_box_cross_day"]),
                xytext=(60, 0.72), fontsize=7.45, color=C_RED, linespacing=1.4,
                ha="left", va="bottom",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none", alpha=0.82),
                arrowprops=dict(arrowstyle="-|>", lw=0.7, color=C_RED))
    ax.text(0.02, 0.94, f"$p \\leq 0.009$\n(4 permutation models,\nBonferroni)",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.45, color=C_NAVY,
            linespacing=1.4,
            bbox=dict(boxstyle="round,pad=0.34", fc="#f4f6f8", ec=C_SLATE, lw=0.6))
    ax.text(240, 700, "all pairs", fontsize=7.35, color=C_SLATE)
    ax.text(1750, 110, "cross-day\n($\\Delta t > 1$ d)", fontsize=7.35, color=C_BLUE,
            ha="right", va="center", linespacing=1.35)
    ax.text(1900, 6.5, "cross-day and\n$|\\Delta\\mathrm{DM}| < 0.5$", fontsize=7.35,
            color=C_RED, ha="right", va="bottom", linespacing=1.35)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1, 2e3); ax.set_ylim(0.4, 5e3)
    ax.set_xlabel("pair-count threshold $\\theta < \\theta_{\\rm cut}$ (arcsec)")
    ax.set_ylabel("number of pairs")
    ax.set_title("Counting experiment", color=C_NAVY, pad=4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def draw_scatter(outdir: Path) -> Path:
    j = load_iso()
    fig = plt.figure(figsize=(TEXT_W_IN, H_IN))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0],
                          left=0.072, right=0.985, top=0.945, bottom=0.125, wspace=0.30)
    panel_nn(fig.add_subplot(gs[0, 0]), j)
    panel_count(fig.add_subplot(gs[0, 1]), j)
    report_collisions(fig, "scatter")
    report_fit("scatter")
    p = outdir / "scatter.pdf"
    fig.savefig(p, format="pdf"); plt.close(fig)
    return p


def _funnel(ax, cx, y_top, y_bot, w_top, w_bot, text, y_num, subs, y_sub, fc, tc):
    ax.add_patch(Polygon([(cx - w_top / 2, y_top), (cx + w_top / 2, y_top),
                          (cx + w_bot / 2, y_bot), (cx - w_bot / 2, y_bot)],
                         closed=True, fc=fc, ec="none"))
    put_fit(ax, cx, y_num, text, size=10, max_w=0.8 * min(w_top, w_bot) / 100.0,
            ha="center", va="center", color=tc, weight="bold")
    for k, s in enumerate(subs):
        put_fit(ax, cx, y_sub - 4.3 * k, s, size=6.2,
                max_w=1.30 * min(w_top, w_bot) / 100.0,
                ha="center", va="center", color=tc)


MIN_PLOT_PT = 7.3
FIT_OVERFLOW: list[tuple[str, float, float]] = []


def put_fit(ax, x, y, s, size=MIN_PLOT_PT, max_w=0.20, **kw):
    r = ax.figure.canvas.get_renderer()
    axw = ax.get_window_extent(renderer=r).width
    size = max(float(size), MIN_PLOT_PT)
    t = ax.text(x, y, s, fontsize=size, **kw)
    for _ in range(30):
        if t.get_window_extent(renderer=r).width <= max_w * axw:
            break
        size -= 0.25
        t.set_fontsize(size)
    if size < MIN_PLOT_PT:
        size = MIN_PLOT_PT
        t.set_fontsize(size)
        FIT_OVERFLOW.append((str(s)[:34], round(max_w * axw, 1),
                             round(t.get_window_extent(renderer=r).width, 1)))
    return t


def report_fit(name: str) -> None:
    if FIT_OVERFLOW:
        print(f"!! {name}: {len(FIT_OVERFLOW)} 处文字在 {MIN_PLOT_PT} pt 地板下放不进"
              f"（已保持地板字号，需缩短文字或放宽框）")
        for s, have, need in FIT_OVERFLOW:
            print(f"     {s!r:36s} 可用 {have} 需要 {need}")
        FIT_OVERFLOW.clear()
    else:
        print(f"{name}: 全部文字 >= {MIN_PLOT_PT} pt 地板")


def draw_flow(outdir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(TEXT_W_IN, H_IN))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    fig.canvas.draw()

    put_fit(ax, 13, 97.6, "5608 rows", size=6.4, max_w=0.22,
            ha="center", va="center", color=C_SLATE)
    put_fit(ax, 13, 90.4, "CHIME/FRB\nCatalogs 1 and 2", size=7.6, max_w=0.28,
            ha="center", va="center", color=C_NAVY, weight="bold", linespacing=1.12)
    _funnel(ax, 12.5, 85.0, 64.0, 23, 17, "5014", 79.5, ["independent events"], 73.0,
            "#3f6798", "white")
    put_fit(ax, 12.5, 59.5, "594 duplicate rows", size=6.2, max_w=0.26,
            ha="center", color=C_SLATE)
    _funnel(ax, 12.5, 55.0, 25.0, 17, 23, "3665", 49.0,
            ["one-off subspace", "(finite DM, pos.)"], 43.0,
            "#c9d7e6", C_NAVY)

    ax.add_patch(FancyBboxPatch((26, 20), 44, 74, boxstyle="round,pad=0.6,rounding_size=1.2",
                                fc="white", ec=C_NAVY, lw=0.9))
    put_fit(ax, 48, 89.5, "Cross-Name Association Framework", size=8.2, max_w=0.445,
            ha="center", color=C_NAVY, weight="bold")
    steps = [
        (33.0, "Step 1", "Cone + DM\npairing"),
        (47.0, "Step 2", "Union-find\nmerging"),
        (61.0, "Step 3", "Null-model\ntest"),
    ]
    for cx, head, body in steps:
        ax.add_patch(FancyBboxPatch((cx - 6, 76), 12, 6,
                                    boxstyle="round,pad=0.3,rounding_size=1.0",
                                    fc="#9db8d6", ec="none"))
        put_fit(ax, cx, 79, head, size=6.8, max_w=0.128, ha="center", va="center",
                color="white", weight="bold")
        put_fit(ax, cx, 70, body, size=6.2, max_w=0.128, ha="center", va="center",
                color=C_NAVY, linespacing=1.5)
        if head == "Step 1":
            ax.add_patch(Polygon([(cx - 4.6, 58), (cx + 4.6, 58), (cx, 46)],
                                 closed=True, fc="#dbe7f4", ec=C_LBLUE, lw=0.7))
            for dx, dy in [(-1.6, 54.0), (1.0, 52.0), (0.1, 55.4)]:
                ax.add_patch(Circle((cx + dx, dy), 0.5, fc=C_BLUE, ec="none"))
            ax.add_patch(FancyBboxPatch((cx - 5.0, 38), 10, 3.4,
                                        boxstyle="round,pad=0.2,rounding_size=0.4",
                                        fc="white", ec=C_SLATE, lw=0.6))
            ax.add_patch(Polygon([(cx - 1.3, 38.4), (cx + 1.3, 38.4), (cx + 1.3, 41.0),
                                  (cx - 1.3, 41.0)], closed=True, fc=C_BLUE, ec="none"))
            put_fit(ax, cx, 34.5, "DM", size=6.4, max_w=0.128, ha="center", color=C_NAVY)
            put_fit(ax, cx, 30.5, "$\\Delta t$, $\\Delta$DM", size=6.2, max_w=0.128,
                    ha="center", color=C_SLATE)
        elif head == "Step 2":
            for dx, dy, c in [(-3.4, 55.0, "#bcd0e6"), (-1.0, 50.5, C_BLUE),
                              (2.0, 54.2, "#8fb0d4"), (4.6, 47.0, C_NAVY),
                              (2.6, 42.6, "#8fb0d4"), (6.0, 43.4, C_BLUE)]:
                ax.add_patch(Circle((cx + dx, dy), 1.25, fc=c, ec="none", zorder=5))
            for a, b in [((-3.4, 55.0), (-1.0, 50.5)), ((-1.0, 50.5), (2.0, 54.2)),
                         ((2.0, 54.2), (4.6, 47.0)), ((4.6, 47.0), (2.6, 42.6)),
                         ((2.6, 42.6), (6.0, 43.4))]:
                ax.plot([cx + a[0], cx + b[0]], [a[1], b[1]], color=C_SLATE, lw=0.6, zorder=4)
            for ox, oy, w, h in [(-1.8, 52.4, 6.6, 9.6), (4.0, 45.4, 6.6, 9.6)]:
                ax.add_patch(matplotlib.patches.Ellipse((cx + ox, oy), w, h, angle=18,
                                                        fc="none", ec=C_SLATE,
                                                        lw=0.6, ls="--"))
        else:
            rng = np.random.default_rng(1)
            h = np.abs(rng.normal(1.5, 1.0, 11))
            for i, v in enumerate(h):
                ax.add_patch(Polygon([(cx - 5.4 + i * 1.0, 40), (cx - 4.55 + i * 1.0, 40),
                                      (cx - 4.55 + i * 1.0, 40 + v * 2.1),
                                      (cx - 5.4 + i * 1.0, 40 + v * 2.1)],
                                     closed=True, fc="#bcd0e6", ec="none"))
            ax.plot([cx + 5.0, cx + 5.0], [40, 47], color=C_NAVY, lw=1.3)
            put_fit(ax, cx, 35.5, "statistic", size=6.2, max_w=0.128,
                    ha="center", color=C_SLATE)
            put_fit(ax, cx, 52.0, "$p \\leq 0.009$", size=6.8, max_w=0.128,
                    ha="center", color=C_RED, weight="bold")
    for x0 in (39.5, 53.5):
        ax.add_patch(FancyArrowPatch((x0, 62), (x0 + 2.0, 62), arrowstyle="-|>",
                                     mutation_scale=8, lw=1.0, color=C_SLATE))
    ax.add_patch(FancyArrowPatch((21.5, 55), (25.5, 55), arrowstyle="-|>",
                                 mutation_scale=10, lw=1.2, color=C_SLATE))
    ax.add_patch(FancyArrowPatch((71.5, 55), (76.5, 55), arrowstyle="-|>",
                                 mutation_scale=10, lw=1.2, color=C_SLATE))

    ax.add_patch(FancyBboxPatch((77.5, 36), 21, 38, boxstyle="round,pad=0.6,rounding_size=1.4",
                                fc="#dbe7f4", ec=C_LBLUE, lw=0.9))
    put_fit(ax, 88, 72.0, "HIGHEST-PRIORITY", size=6.4, max_w=0.24, va="top",
            ha="center", color=C_SLATE, weight="bold")
    put_fit(ax, 88, 66.2, "SURVIVING\nCANDIDATE", size=6.4, max_w=0.21, va="top", linespacing=1.15,
            ha="center", color=C_SLATE, weight="bold")
    ax.add_patch(FancyBboxPatch((78.6, 36.4), 18.8, 18.0, boxstyle="round,pad=0.4,rounding_size=0.9",
                                fc="#3f6798", ec="none"))
    put_fit(ax, 88, 49.6, "FRB20200511C", size=7.0, max_w=0.20, ha="center",
            color="white", weight="bold")
    ax.text(88, 45.6, "$\\times$", ha="center", va="center", fontsize=9, color="white")
    put_fit(ax, 88, 40.0, "FRB20211128E", size=7.0, max_w=0.20, ha="center",
            color="white", weight="bold")
    report_collisions(fig, "flow")
    report_fit("flow")
    p = outdir / "flow.pdf"
    fig.savefig(p, format="pdf"); plt.close(fig)
    return p



def draw_roadmap(outdir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(TEXT_W_IN, H_IN))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    cols = [
        (17, C_BLUE, "#e8eef7", "FAST Telescope\nThird burst",
         "Repeatability\nConfirmation", "dish"),
        (50, "#0f9b8e", "#e3f4f1", "Polarization RM\nMagnetic environment",
         "Magnetized Medium\nDiagnosis", "rm"),
        (83, "#6a4fb0", "#efeaf8", "VLBI / Optical\nPrecise position\nH1 vs H2", "Redshift\nAnchor", "array"),
    ]
    for cx, ec, fc, head, out, icon in cols:
        if icon == "dish":
            ax.add_patch(matplotlib.patches.Ellipse((cx, 78), 20, 7.5, angle=0,
                                                    fc="#dbe7f4", ec=ec, lw=0.9))
            for dx in (-9.5, -5, 0, 5, 9.5):
                ax.plot([cx + dx, cx + dx * 0.55], [79.5, 84.5], color=ec, lw=0.6)
            ax.plot([cx - 9.5, cx - 7.5], [76.5, 71], color=ec, lw=0.8)
            ax.plot([cx + 9.5, cx + 7.5], [76.5, 71], color=ec, lw=0.8)
            ax.add_patch(Circle((cx, 76.4), 0.85, fc=ec, ec="none"))
        elif icon == "rm":
            for r, a0 in [(3.0, 0), (6.0, 40), (9.0, 80)]:
                th = np.linspace(np.radians(a0), np.radians(a0 + 300), 40)
                ax.plot(cx + r * np.cos(th) * 1.5, 78 + r * np.sin(th) * 0.42,
                        color="#8a6fc0", lw=0.6)
            for ang in range(0, 360, 45):
                dx, dy = 7.6 * np.cos(np.radians(ang)), 3.4 * np.sin(np.radians(ang))
                ax.add_patch(FancyArrowPatch((cx + dx * 0.45, 78 + dy * 0.45),
                                             (cx + dx, 78 + dy), arrowstyle="-|>",
                                             mutation_scale=5, lw=0.6, color=ec))
            ax.text(cx + 13.5, 82, "RM", fontsize=7.45, color=C_NAVY)
            ax.add_patch(FancyBboxPatch((cx + 11.5, 70), 4.2, 13,
                                        boxstyle="round,pad=0.1,rounding_size=0.3",
                                        fc="none", ec="none"))
            for i, c in enumerate(["#3b6fb5", "#6f8fd0", "#a58fd0", "#d08fb5", "#e05a5a"]):
                ax.add_patch(Polygon([(cx + 12, 81 - i * 2.4), (cx + 14.6, 81 - i * 2.4),
                                      (cx + 14.6, 78.6 - i * 2.4), (cx + 12, 78.6 - i * 2.4)],
                                     closed=True, fc=c, ec="none"))
        else:
            for dx, dy in [(-8, 6), (8, 6), (-8, -4), (8, -4)]:
                ax.add_patch(matplotlib.patches.Ellipse((cx + dx, 78 + dy), 5.6, 2.1,
                                                        angle=-18, fc="#dbe7f4",
                                                        ec=ec, lw=0.7))
            ax.add_patch(Circle((cx, 78), 1.9, fc="#8a6fc0", ec="none"))
            for dx, dy in [(-8, 6), (8, 6), (-8, -4), (8, -4)]:
                ax.plot([cx, cx + dx * 0.8], [78, 78 + dy * 0.8], color=ec, lw=0.5, ls="--")
            for k, (lab, c) in enumerate([("H1", C_GOLD), ("H2", "#3fb8d8")]):
                ax.add_patch(Circle((cx + 14.5, 81.5 - k * 4.2), 1.0, fc=c, ec="none"))
                ax.text(cx + 16.3, 81.5 - k * 4.2, lab, va="center", fontsize=7.35,
                        color=C_NAVY)
        ax.add_patch(FancyBboxPatch((cx - 15.5, 53), 31, 13,
                                    boxstyle="round,pad=0.8,rounding_size=6",
                                    fc=ec, ec="none"))
        put_fit(ax, cx, 59.5, head, size=6.4, max_w=0.28, ha="center", va="center",
                color="white", weight="bold", linespacing=1.5)
        n = int((cx - 17) // 33) + 1
        ax.add_patch(Circle((cx - 11, 46), 2.6, fc=ec, ec="none"))
        ax.text(cx - 11, 46, str(n), ha="center", va="center", fontsize=7.6,
                color="white", weight="bold")
        ax.text(cx, 34, "Key Scientific Output", ha="center", fontsize=7.35, color=C_SLATE)
        ax.text(cx, 26.5, out, ha="center", va="center", fontsize=7.6, color=C_NAVY,
                weight="bold", linespacing=1.5)
    for x0 in (23.5, 56.5):
        ax.add_patch(FancyArrowPatch((x0, 46), (x0 + 26.5, 46), arrowstyle="-|>",
                                     mutation_scale=11, lw=1.4, color=C_SLATE))
    ax.text(50, 12.5, "Current uncertainty window: $z \\in [0.3, 2.5]$ (H1); "
                      "H2 lies at $z < 0.1$", ha="center", fontsize=7.55, color=C_SLATE)
    report_collisions(fig, "roadmap")
    report_fit("roadmap")
    p = outdir / "roadmap.pdf"
    fig.savefig(p, format="pdf"); plt.close(fig)
    return p


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for fn in (draw_scatter, draw_flow, draw_roadmap):
        p = fn(OUTDIR)
        print(f"  {p.name:14s} {p.stat().st_size/1024:6.0f} KB")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    raise SystemExit(main())
