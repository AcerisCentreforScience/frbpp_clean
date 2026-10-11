"""Mechanism figures fig1-fig6 for the evidence chain.

Usage: python scripts/make_figs.py --only all --dpi 900 --outdir docs/figures/hires
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
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
from matplotlib.text import Text
from matplotlib.lines import Line2D
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from frbpp.upgrade.crossmatch import (build_event_table, geometry_pairs,  # noqa: E402
                                      LOGIC_VERSION, same_event_mask)
from fig_gate import report_collisions as _shared_gate  # noqa: E402

def report_collisions(fig, name):
    """Shared layout gate, made fatal: an overprinted raster is not a draft figure."""
    n = _shared_gate(fig, name)
    if n:
        raise SystemExit(f"{name}: {n} 处叠字/越出画布/压数据，位图不出图")
    return n


FIGDIR = ROOT / "docs" / "figures"

MIN_FS = 7.30
PRINT_W = 15.0 / 2.54


def SZ(w, h):
    return (PRINT_W, h * PRINT_W / w)


FIT_OVERFLOW: list[tuple[str, float, float]] = []

def FS(v):
    return float(max(MIN_FS, float(v)))


def fit(ax, x, y, s, size, max_w=0.30, xy=None, **kw):
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    axw = ax.get_window_extent(renderer=r).width
    kw.setdefault("clip_on", False)
    t = (ax.annotate(s, xy=xy, xytext=(x, y), fontsize=size, **kw) if xy is not None
         else ax.text(x, y, s, fontsize=size, **kw))
    floor = float(MIN_FS)

    def text_w():
        return Text.get_window_extent(t, renderer=r).width

    for _ in range(28):
        if text_w() <= max_w * axw:
            break
        if floor and size - 0.3 < floor:
            t.set_fontsize(floor)
            FIT_OVERFLOW.append((str(s)[:34], round(max_w * axw, 1), round(text_w(), 1)))
            return t
        size -= 0.3
        t.set_fontsize(size)
    return t


def report_fit(name: str) -> None:
    if not FIT_OVERFLOW:
        print(f"{name}: 全部文字 >= {MIN_FS:g} pt 地板")
        return
    print(f"!! {name}: {len(FIT_OVERFLOW)} 处文字在 {MIN_FS:g} pt 地板下放不进（已保持地板字号，"
          f"需缩短文字或放宽框）")
    for s, have, need in FIT_OVERFLOW:
        print(f"     {s!r:36s} 可用 {have} 需要 {need}")
    bad = len(FIT_OVERFLOW)
    FIT_OVERFLOW.clear()
    raise SystemExit(f"{name}: {bad} 处文字放不进 {MIN_FS:g} pt 地板，位图不出图")


UP = ROOT / "data" / "processed" / "upgrade"
INTERIM = ROOT / "data" / "interim"

C_GREY, C_BLUE, C_RED, C_ORANGE = "#B3B3B3", "#2166AC", "#B2182B", "#E08214"
C_GREEN, C_TEAL, C_PURPLE = "#1B7837", "#01665E", "#762A83"
C_INK = "#1A1A1A"
_CJK = ["DejaVu Sans", "Noto Sans SC", "Microsoft YaHei", "SimHei"]
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": _CJK, "font.size": 7,
    "axes.unicode_minus": False,
    "axes.labelsize": 7.5, "axes.titlesize": 7.5, "axes.linewidth": 0.6,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 400, "savefig.dpi": 400, "savefig.bbox": "tight",
})


def panel_label(ax, s, dx=-0.13, dy=1.05):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=FS(9), fontweight="bold",
            ha="left", va="top")


def sep_arcsec(ra1, dec1, ra2, dec2):
    c1 = SkyCoord(np.asarray(ra1) * u.deg, np.asarray(dec1) * u.deg)
    c2 = SkyCoord(np.asarray(ra2) * u.deg, np.asarray(dec2) * u.deg)
    return c1.separation(c2).arcsec


def load(v13=True):
    c1 = pd.read_parquet(INTERIM / "cat1_bursts.parquet")
    c2 = pd.read_parquet(INTERIM / "cat2_bursts.parquet")
    return c1, c2, build_event_table(c1, c2, event_dedup=v13)


def fig1():
    inv = json.loads((UP / "g4b_rerun_20261009_v13.json").read_text(encoding="utf-8"))
    fig, ax = plt.subplots(figsize=SZ(7.4, 7.4))
    ax.set_axis_off()
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)

    steps = [
        ("1  catalog merge",
         "Cat1 600 rows + Cat2 5008 rows",
         "5608 rows / 4242 one-off",
         C_GREY, "same event counted twice at row level"),
        ("2  event-layer dedup  (v1.3)",
         "drop Cat1 rows with same name & |dt| < 0.05 d",
         "5014 events / 3742 one-off",
         C_TEAL, "event pool is not an independent-event set"),
        ("3  pairing-layer dedup  (iron rule)",
         "same_event_mask drops same-event pairs",
         "44 -> 44 candidates (unchanged)",
         C_BLUE, "~6e4 self-pairs contaminated by dups"),
        ("4  9-cell (theta, dDM) scan",
         "theta in {60,180,600}\" x dDM in {1,5,20}",
         "9/9 cells: obs > mean + 3 sigma",
         C_BLUE, "single-cell tuning (theta selection)"),
        ("5  candidate union",
         "de-dup by name, keep tightest theta",
         "28 pairs (13 cross-day + 15 same-day)",
         C_ORANGE, "pure chance coincidence"),
        ("6  four permutation nulls",
         "RA-strat / DM / RA+DM / full, 1000 each",
         "0/1000 reach 28 pairs (p_g = 9.0e-3)",
         C_PURPLE, "joint position-and-DM chance"),
        ("7  connected components",
         "28 pairs -> 21 events -> 8 components",
         "[7,2,2,2,2,2,2,2] -> 7 hypotheses",
         C_GREEN, "28 pairs = 28 independent finds (inflated)"),
        ("8  primary-candidate isolation MC",
         "error convolution, 1000 realizations",
         "min theta(|dDM|<0.5) median 1.2e4\"",
         C_RED, "23.9\" is unsurprising given ~700\" errors"),
    ]
    n = len(steps)
    top, bot = 92.0, 4.0
    cell = (top - bot) / n
    bh = cell * 0.80
    x_l, w_l = 2.0, 50.0
    x_r, w_r = 54.5, 43.5
    for k, (title, op, out, col, kills) in enumerate(steps):
        y = top - (k + 0.5) * cell
        yc = y + bh / 2
        ax.add_patch(FancyBboxPatch((x_l, y), w_l, bh,
                                    boxstyle="round,pad=0.15,rounding_size=1.1",
                                    facecolor="white", edgecolor=col, lw=1.05, zorder=2))
        ax.text(x_l + 1.6, yc + bh * 0.30, title, fontsize=FS(7.1), fontweight="bold",
                color=col, va="center")
        ax.text(x_l + 1.6, yc + bh * 0.03, op, fontsize=FS(5.9), color="0.28", va="center")
        ax.text(x_l + 1.6, yc - bh * 0.28, "-> " + out, fontsize=FS(5.9), color=C_INK,
                va="center", fontweight="bold")
        ax.add_patch(FancyBboxPatch((x_r, y), w_r, bh,
                                    boxstyle="round,pad=0.15,rounding_size=1.1",
                                    facecolor=col, alpha=0.08, edgecolor="none", zorder=1))
        ax.text(x_r + 1.6, yc, "x  " + kills, fontsize=FS(5.9), color="0.22", va="center",
                linespacing=1.35)
        if k < n - 1:
            ax.add_patch(FancyArrowPatch((x_l + w_l / 2, y - 0.15),
                                         (x_l + w_l / 2, y - cell + bh + 0.15),
                                         arrowstyle="-|>", mutation_scale=7,
                                         color="0.5", lw=0.8, zorder=3))
    fit(ax, 2.0, 99.0, "Observation-to-null evidence chain   (logic v%s, 2026-10-09)"
        % LOGIC_VERSION, FS(8.2), max_w=0.96, fontweight="bold", color=C_INK, va="top")
    fit(ax, 2.0, 1.0, "Key invariant: step 2 shrinks only the denominator (5608->5014); the criterion\n"
        "set is frozen at 28 pairs, which is what every null below is measured against.",
        FS(5.7), max_w=0.96, color="0.32", va="bottom")
    report_collisions(fig, 'fig1_evidence_chain')
    fig.savefig(FIGDIR / "fig1_evidence_chain.png")
    plt.close(fig)


def fig2():
    audit = json.loads((UP / "event_table_audit.json").read_text(encoding="utf-8"))
    nb = json.loads((UP / "c2_null_battery.json").read_text(encoding="utf-8"))
    rl, dd = audit["row_level"], audit["dedup"]
    ind = nb["independence"]
    n_cat1 = rl["cat1_rows"]
    n_cat2 = rl["cat2_rows_ok"]
    n_rows = rl["all"]
    n_cat1_only = dd["cat1_rows_only_in_cat1"]
    n_dup = dd["cat1_rows_matched_in_cat2"]
    n_all = dd["all"]
    n_rep = dd["repeater"]
    n_oo = dd["one_off"]
    n_dom = dd["one_off_finite_dm"]
    n_pairs = ind["n_pairs"]; n_ev = ind["n_unique_events"]; n_comp = ind["n_components"]
    sizes = ind["component_sizes"]
    assert n_cat2 + n_cat1_only == n_all, (n_cat2, n_cat1_only, n_all)
    assert n_rep + n_oo == n_all, (n_rep, n_oo, n_all)
    assert n_rows - n_dup == n_all, (n_rows, n_dup, n_all)
    assert n_dom <= n_oo
    assert (n_pairs, n_ev, n_comp) == (28, 21, 8), (n_pairs, n_ev, n_comp)
    assert sum(sizes) == n_ev, (sizes, n_ev)
    n_screen = n_dom * (n_dom - 1) // 2

    fig = plt.figure(figsize=SZ(7.4, 3.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.14)

    def cascade(ax, steps, head):
        ax.set_axis_off(); ax.set_xlim(0, 100); ax.set_ylim(0, 100)
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        axw = ax.get_window_extent(renderer=r).width
        fit(ax, 1, 98.5, head, FS(7.6), max_w=0.97 * axw / axw,
            fontweight="bold", color=C_INK, va="top")
        y = 88.0
        h = 20.0
        for k, (lab, num, note) in enumerate(steps):
            last = k == len(steps) - 1
            ax.add_patch(FancyBboxPatch((1, y - h), 98, h - 3.0,
                                        boxstyle="round,pad=0.5,rounding_size=2.0",
                                        fc="#F2F6FA" if not last else "#E7F1E9",
                                        ec=C_BLUE if not last else C_GREEN, lw=0.9))
            bw = 0.98 * axw
            fit(ax, 4, y - 4.0, lab, FS(6.6), max_w=0.74 * bw / axw, va="top",
                color="0.12")
            fit(ax, 96, y - 4.0, f"{num:,}", FS(9.0), max_w=0.30 * bw / axw, ha="right",
                va="top", fontweight="bold", color=C_BLUE if not last else C_GREEN)
            if note:
                fit(ax, 4, y - 11.0, note, FS(5.9), max_w=0.96 * bw / axw, va="top",
                    color="0.35", linespacing=1.35)
            if not last:
                ax.add_patch(FancyArrowPatch((50, y - h + 0.6), (50, y - h - 2.4),
                                             arrowstyle="-|>", mutation_scale=7,
                                             color="0.45", lw=0.9))
            y -= h + 3.0
        bad = []
        for t in ax.texts:
            tb = t.get_window_extent(renderer=r)
            if tb.x1 > ax.get_window_extent(renderer=r).x1 + 1:
                bad.append(t.get_text()[:28])
        if bad:
            print(f"   !! 文字越出面板右缘: {bad}")
        return len(bad)

    ax = fig.add_subplot(gs[0, 0])
    cascade(ax, [
        ("catalogue rows concatenated", n_rows,
         f"{n_cat1} Cat1 + {n_cat2:,} localized Cat2 rows"),
        ("distinct events", n_all,
         f"{n_dup} Cat1 rows also in Cat2 dropped"),
        ("apparent one-off events", n_oo,
         f"{n_rep:,} repeater-linked set aside"),
        ("declared search domain", n_dom,
         "finite dispersion measure, RA/Dec"),
    ], "From rows to the search domain")
    panel_label(ax, "a", dx=-0.10)

    ax = fig.add_subplot(gs[0, 1])
    cascade(ax, [
        ("name pairs screened", n_screen,
         f"all pairs in the {n_dom:,}-event domain"),
        ("candidate pairs", n_pairs,
         "passing the $\\theta$ and $\\Delta$DM vetoes"),
        ("unique events involved", n_ev,
         "the pairs share bursts"),
        ("independent hypotheses", n_comp,
         "sizes " + ", ".join(map(str, sizes))),
    ], "What the search returns")
    panel_label(ax, "b", dx=-0.10)

    report_collisions(fig, 'fig2_event_dedup_mechanism')
    fig.savefig(FIGDIR / "fig2_event_dedup_mechanism.png")
    plt.close(fig)

def fig3():
    c1, c2, ev = load(v13=True)
    gi, gj = geometry_pairs(ev, 600.0)
    keep = ~same_event_mask(ev, gi, gj)
    gi, gj = gi[keep], gj[keep]
    oo = np.array([(not ev.src_rep[i]) and (not ev.src_rep[j]) for i, j in zip(gi, gj)])
    ia, ib = gi[oo], gj[oo]
    th = sep_arcsec(ev.ra[ia], ev.dec[ia], ev.ra[ib], ev.dec[ib])
    ddm = np.abs(ev.dm[ia] - ev.dm[ib])
    ok = np.isfinite(ddm) & (th > 0)
    cand = pd.read_csv(UP / "c2_candidates_v1.csv")
    iso = json.loads((UP / "pair565_isolation.json").read_text(encoding="utf-8"))
    samples = np.array(iso["samples"]["min_theta05_arcsec"], float)

    fig = plt.figure(figsize=SZ(7.6, 3.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 0.9], wspace=0.32)

    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(th[ok], np.maximum(ddm[ok], 8e-4), s=6, c=C_GREY, alpha=0.6,
               edgecolors="none", rasterized=True, zorder=1)
    for _, r in cand.iterrows():
        x, y = r["theta_arcsec"], max(abs(r["ddm"]), 8e-4)
        if r["group"] == "cross_day":
            ax.scatter(x, y, s=26, facecolors="none", edgecolors=C_BLUE, lw=0.9, zorder=3)
        else:
            ax.scatter(x, y, s=12, marker="s", facecolors="none",
                       edgecolors="#9AA7B4", lw=0.7, zorder=2)
    ax.axvspan(8, 23.86, color=C_RED, alpha=0.05, lw=0)
    ax.axhline(0.5, color=C_ORANGE, lw=0.8, ls="--")
    ax.axvline(600, color="0.4", lw=0.7, ls=":")
    px, py = 23.86, 0.1769
    ax.scatter(px, py, s=110, marker="*", c=C_RED, edgecolors="k", lw=0.4, zorder=6)
    fit(ax, 0.03, 0.26,
        "FRB20200511C\nx FRB20211128E\n23.9'', 0.18 = 0.12$\\sigma$ (floor)",
        FS(5.8), max_w=0.92, xy=(px, py), textcoords=ax.transAxes,
        color=C_RED, ha="left", va="top", linespacing=1.35,
        arrowprops=dict(arrowstyle="-|>", color=C_RED, lw=0.7, shrinkA=0, shrinkB=3))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(8, 900); ax.set_ylim(5e-4, 3e3)
    ax.set_xlabel("Angular separation $\\theta$ (arcsec)")
    ax.set_ylabel("|$\\Delta$DM| (pc cm$^{-3}$)")
    ax.grid(True, which="major", ls=":", lw=0.3, alpha=0.4)
    ax.legend(handles=[
        Line2D([], [], ls="", marker="o", ms=3, mfc=C_GREY, mec="none", label="one-off pairs (n=44)"),
        Line2D([], [], ls="", marker="o", ms=4, mfc="none", mec=C_BLUE, mew=0.9, label="cross-day"),
        Line2D([], [], ls="", marker="s", ms=3.2, mfc="none", mec="#9AA7B4", mew=0.7, label="same-day"),
        Line2D([], [], ls="", marker="*", ms=7, mfc=C_RED, mec="k", mew=0.4, label="primary"),
    ], loc="upper left", frameon=False, handletextpad=0.3, borderaxespad=0.2)
    panel_label(ax, "a", dx=-0.16)

    ax = fig.add_subplot(gs[0, 1])
    lo, hi = 30.0, 20000.0
    bins = np.logspace(np.log10(lo), np.log10(hi), 34)
    ax.hist(samples, bins=bins, color=C_BLUE, alpha=0.85, edgecolor="white", lw=0.25)
    ax.axvline(20.1, color=C_RED, lw=1.0)
    ax.axvline(23.9, color=C_RED, lw=1.0, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("min $\\theta$ (|$\\Delta$DM|<0.5), arcsec")
    ax.set_ylabel("MC realizations")
    ax.set_ylim(0, 240)
    fit(ax, 0.03, 0.62, "obs 20.1''\n113.7 d pair", FS(5.8), max_w=0.62,
        xy=(20.1, 45), textcoords=ax.transAxes, color=C_RED, ha="left", va="center",
        linespacing=1.35, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85),
        arrowprops=dict(arrowstyle="-|>", color=C_RED, lw=0.7))
    fit(ax, 0.97, 0.97, "error-convolved null\n1,000 draws\nmedian 1.2$\\times$10$^4$''\n"
        "0 reach the box", FS(5.6), max_w=0.88, transform=ax.transAxes,
        ha="right", va="top", color="0.15", linespacing=1.4)
    ax.tick_params(axis="y", labelsize=FS(6.0))
    panel_label(ax, "b", dx=-0.22)

    ax = fig.add_subplot(gs[0, 2])
    ax.set_axis_off(); ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    fit(ax, 0, 99, "Error chain", FS(7.0), max_w=0.94,
        fontweight="bold", color=C_INK, va="top")
    items = [
        ("localization", "$\\sigma$ ~ 705''/745''\n(ra / dec median)", C_GREY),
        ("independence", "$\\sigma_{comb}$ ~ 1189''\n-> 23.9'' trivial", C_GREY),
        ("same beam", "systematic errors\nare correlated", C_GREEN),
        ("precision", "beats 1189''\nMC conservative", C_GREEN),
    ]
    y = 94
    for k, (t, b, col) in enumerate(items):
        ax.add_patch(Rectangle((1, y - 21), 98, 20,
                               facecolor=col, alpha=0.10, edgecolor=col, lw=0.9))
        fit(ax, 5, y - 2.5, t, FS(5.9), max_w=0.94, fontweight="bold",
            color="0.15", va="top")
        fit(ax, 5, y - 9.5, b, FS(5.4), max_w=0.94, va="top", color="0.25",
            linespacing=1.4)
        if k < 3:
            ax.add_patch(FancyArrowPatch((88, y - 21.4), (88, y - 23.6),
                                         arrowstyle="-|>", mutation_scale=6,
                                         color="0.45", lw=0.8))
        y -= 23.0
    panel_label(ax, "c", dx=-0.05)
    report_collisions(fig, 'fig3_joint_isolation')
    fig.savefig(FIGDIR / "fig3_joint_isolation.png")
    plt.close(fig)


def fig4():
    c2 = pd.read_parquet(INTERIM / "cat2_bursts.parquet")
    b = c2[c2.sub_num == 0].copy()
    T_sid = 0.99726958
    a = b[b.tns_name == "FRB20200511C"].iloc[0]
    e = b[b.tns_name == "FRB20211128E"].iloc[0]
    dt = float(e.mjd_400 - a.mjd_400)
    n_sid = dt / T_sid

    fig = plt.figure(figsize=SZ(7.4, 3.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.15, 1.0], wspace=0.42)

    ax = fig.add_subplot(gs[0, 0])
    k = np.arange(565, 570)
    ax.vlines(k, 0, 1, color="0.85", lw=6)
    ax.vlines(n_sid, 0, 1, color=C_RED, lw=1.6, zorder=5)
    ax.set_xlim(565, 569)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("$\\Delta t \\,/\\, T_{sid}$")
    ax.annotate("567.00191\n(164 s from integer)",
                xy=(n_sid, 0.55), xytext=(565.05, 0.90), fontsize=FS(5.8), color=C_RED,
                arrowprops=dict(arrowstyle="-", color=C_RED, lw=0.7))
    fit(ax, 0.97, 0.04, "P(chance within\n0.002 of an integer)\n~ 0.4%",
        FS(5.6), max_w=0.92, transform=ax.transAxes, color="0.2", va="bottom",
        ha="right")
    panel_label(ax, "a", dx=-0.20)

    ax = fig.add_subplot(gs[0, 1])
    rng = np.random.default_rng(20261009)
    rep = b.repeater_name.astype("string").str.strip()
    b["src"] = rep.mask(rep.isin(["", "-9999", "nan", "None"]), None)
    oo = b[b.src.isna() & b.ra.notna() & b.dm_fitb.notna()]
    n = len(oo); t = oo.mjd_400.to_numpy(float)
    idx = rng.integers(0, n, size=(500_000, 2))
    m = idx[:, 0] != idx[:, 1]
    dto = np.abs(t[idx[m, 0]] - t[idx[m, 1]]); dto = dto[dto > 5]
    fo = np.mod(dto / T_sid, 1.0); fo = np.where(fo > 0.5, 1 - fo, fo)
    p = []
    for s, g in b[b["src"].notna()].groupby("src"):
        tt = g.mjd_400.to_numpy(float)
        if len(tt) >= 2:
            i, j = np.triu_indices(len(tt), 1)
            p.append(np.abs(tt[i] - tt[j]))
    p = np.concatenate(p); p = p[p > 5]
    fr = np.mod(p / T_sid, 1.0); fr = np.where(fr > 0.5, 1 - fr, fr)
    bins = np.linspace(0, 0.5, 61)
    ax.hist(fo, bins=bins, density=True, color=C_GREY, alpha=0.8, edgecolor="white",
            lw=0.2, label="cross-source (random)")
    ax.hist(fr, bins=bins, density=True, color=C_BLUE, alpha=0.72, edgecolor="white",
            lw=0.2, label="same-source (repeaters)")
    pf = abs(np.mod(dt / T_sid + 0.5, 1.0) - 0.5)
    ax.axvline(pf, color=C_RED, lw=1.2, zorder=6, label="this pair = 567.00 $T_{sid}$")
    ax.set_xlim(0, 0.09); ax.set_xlabel("fractional part of $\\Delta t/T_{sid}$ (folded)")
    ax.set_ylabel("probability density")
    ax.legend(loc="upper right", frameon=False, handletextpad=0.3, fontsize=FS(5.8))
    panel_label(ax, "b", dx=-0.16)

    DRAO = EarthLocation.from_geodetic(lon=-119.6286 * u.deg, lat=49.3167 * u.deg,
                                       height=700 * u.m)
    d = b[b.ra.notna() & b.mjd_400.notna()]
    tt = Time(d.mjd_400.to_numpy(float), format="mjd", scale="utc")
    ha = (tt.sidereal_time("apparent", longitude=DRAO.lon) - d.ra.to_numpy(float) * u.deg)
    ah = np.abs(ha.wrap_at(180 * u.deg).to(u.hourangle).value)
    ax = fig.add_subplot(gs[0, 2])
    ax.hist(ah, bins=np.linspace(0, 1.2, 45), color=C_ORANGE, alpha=0.9,
            edgecolor="white", lw=0.25)
    ax.axvspan(0, 0.1, color=C_BLUE, alpha=0.10, lw=0)
    ax.axvline(np.median(ah), color="k", lw=0.8, ls="--")
    ax.set_xlabel("|hour angle| (h)"); ax.set_ylabel("Cat2 events")
    fit(ax, 0.97, 0.28, "median 0.049 h\n78.9% < 0.1 h\nCHIME window",
        FS(5.8), max_w=0.92, transform=ax.transAxes, ha="right", va="top", color="0.15")
    axi = ax.inset_axes([0.50, 0.38, 0.40, 0.42])
    theta = np.linspace(-1.2, 1.2, 100)
    axi.plot(theta, np.cos(theta * np.pi / 2.4), color=C_BLUE, lw=1.0)
    axi.axvspan(-0.1, 0.1, color=C_RED, alpha=0.20, lw=0)
    axi.set_xticks([]); axi.set_yticks([])
    for s in ("top", "right"):
        axi.spines[s].set_visible(False)
    axi.set_xlabel("HA", fontsize=FS(5.0), labelpad=0.5)
    axi.set_title("transit", fontsize=FS(5.0), pad=1)
    panel_label(ax, "c", dx=-0.20)
    report_collisions(fig, 'fig4_sidereal_quantization')
    fig.savefig(FIGDIR / "fig4_sidereal_quantization.png")
    plt.close(fig)


def fig5_null_battery():
    j = json.loads((UP / "c2_null_battery.json").read_text(encoding="utf-8"))
    models = [("ra", "RA permutation\n(Dec-stratified)"),
              ("dm", "DM permutation\n(distribution-preserving)"),
              ("ra_dm", "RA + DM permutation"),
              ("full", "Fully randomized")]
    fig, axes = plt.subplots(1, 4, figsize=SZ(7.6, 2.7))
    for ax, (key, title) in zip(axes.ravel(), models):
        p = j["perm"][key]
        ref = int(p["real_cand"])
        hist = np.array(p["hist"], float)
        edges = np.array(p["hist_edges"], float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        ax.bar(centers, hist, width=0.82, color=C_BLUE, alpha=0.85, edgecolor="white",
               lw=0.3)
        ax.axvline(ref, color=C_RED, lw=1.1)
        ax.set_yscale("log")
        ax.set_ylim(0.6, 3000)
        ax.set_xlim(-1, ref + 1.5)
        ax.set_xlabel("candidate pairs / realization", fontsize=FS(6.0))
        if key == "ra":
            ax.set_ylabel("realizations")
        ax.set_title(title, fontsize=FS(6.2), pad=3)
        ax.text(0.04, 0.97, f"0/1000 $\\geq$ {ref}",
                transform=ax.transAxes, fontsize=FS(6.2), ha="left", va="top",
                color="0.15")
        panel_label(ax, dict(zip([k for k, _ in models], "abcd"))[key], dx=-0.30, dy=1.16)
    fig.text(0.5, 0.005, "single p <= 0.001 (B=1000 floor)  |  Bonferroni x9 global p <= 0.009  |  all four nulls fail alike",
             ha="center", fontsize=FS(5.8), color="0.30")
    fig.tight_layout(w_pad=1.1, rect=(0, 0.045, 1, 0.86))
    report_collisions(fig, 'fig5_null_battery')
    fig.savefig(FIGDIR / "fig5_null_battery.png")
    plt.close(fig)


def fig6():
    c2 = pd.read_parquet(INTERIM / "cat2_bursts.parquet")
    b = c2[c2.sub_num == 0]
    meas = b[(b.scat_time > 0) & b.scat_time_err.notna()].copy()
    sc = json.loads((UP / "pair565_scatter.json").read_text(encoding="utf-8"))
    cell = sc["control_cell"]["tau_ms_measured"]
    acoef, bcoef = -3.548, 0.589
    fig, ax = plt.subplots(figsize=SZ(4.5, 3.5))
    ax.scatter(meas.dm_fitb, meas.scat_time * 1e3, s=3.5, c=C_GREY, alpha=0.5,
               edgecolors="none", rasterized=True, label="measured $\\tau$ (n=1,811)")
    xs = np.logspace(2, np.log10(4000), 100)
    ax.plot(xs, 10 ** (acoef + bcoef * np.log10(xs)) * 1e3, color=C_BLUE, lw=1.0,
            label="catalog scaling\nlog $\\tau$ = $-$3.55 + 0.59 log DM")
    ax.add_patch(Rectangle((1500, 0.015), 500, 3000, facecolor=C_ORANGE, alpha=0.13,
                           edgecolor="none", zorder=0))
    ax.text(1750, 1.05e3, "matched cell\nSNR 8-15\nmedian 31 ms", fontsize=FS(5.5),
            ha="center", color="0.2", linespacing=1.35, va="top")
    for val, ls in ((cell["p25"], ":"), (cell["median"], "-"), (cell["p75"], ":")):
        ax.plot([1500, 2000], [val, val], color=C_ORANGE, lw=0.8, ls=ls, zorder=4)
    ax.annotate("", xy=(1835, 5.2), xytext=(1835, 3000),
                arrowprops=dict(arrowstyle="->", color=C_RED, lw=1.2))
    ax.scatter([1835], [5.2], s=42, marker="v", c=C_RED, edgecolors="k", lw=0.4, zorder=6)
    ax.scatter([1835], [cell["median"]], s=30, facecolors="none", edgecolors=C_RED,
               lw=0.9, zorder=6)
    ax.annotate("FRB20200511C x FRB20211128E\nno measurable $\\tau$ (<~5.2 ms)\nvs expected 24 ms / cell 31 ms",
                xy=(1835, 4.6), xytext=(112, 3.0e-1), fontsize=FS(5.4), color=C_RED,
                arrowprops=dict(arrowstyle="-", color=C_RED, lw=0.6), linespacing=1.45,
                va="top")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(100, 4200); ax.set_ylim(0.015, 3000)
    ax.set_xlabel("DM (pc cm$^{-3}$)"); ax.set_ylabel("scattering time $\\tau$ (ms)")
    ax.grid(True, which="major", ls=":", lw=0.3, alpha=0.4)
    ax.legend(loc="lower left", frameon=False, handletextpad=0.4, borderaxespad=0.2,
              fontsize=FS(5.6))
    ax.text(0.035, 0.97, ">~2 orders below both the trend\nand the matched control cell\n-> DM not from a turbulent screen\n-> independent H1/H2 discriminant",
            transform=ax.transAxes, fontsize=FS(5.4), ha="left", va="top", color="0.15",
            linespacing=1.45)
    fig.tight_layout()
    report_collisions(fig, 'fig6_scattering_discriminant')
    fig.savefig(FIGDIR / "fig6_scattering_discriminant.png")
    plt.close(fig)


ALL = {"fig1": fig1, "fig2": fig2, "fig3": fig3, "fig4": fig4,
       "fig5": fig5_null_battery, "fig6": fig6}
SUBMISSION = ["fig2", "fig3", "fig4", "fig5", "fig6"]
NAMES = {"fig1_evidence_chain.png": "fig1", "fig2_event_dedup_mechanism.png": "fig2",
         "fig3_joint_isolation.png": "fig3", "fig4_sidereal_quantization.png": "fig4",
         "fig5_null_battery.png": "fig5", "fig6_scattering_discriminant.png": "fig6"}

if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="重出机制图")
    ap.add_argument("--dpi", type=float, default=400.0, help="输出 dpi（hires/ 用 900）")
    ap.add_argument("--outdir", default=str(FIGDIR), help="输出目录（默认 docs/figures）")
    ap.add_argument("--only", default="all",
                    help=f"all（={','.join(SUBMISSION)}）或逗号分隔，如 fig3,fig5；"
                         f"未入稿的 fig1 需显式点名")
    ap.add_argument("--min-fs", type=float, default=7.30,
                    help="字号地板 pt（默认 7.3，145mm 版心下印 7.06）；0 表示保持旧行为")
    a = ap.parse_args()

    globals()["MIN_FS"] = a.min_fs
    matplotlib.rcParams.update({"font.size": FS(7.0), "axes.labelsize": FS(7.5),
                                "axes.titlesize": FS(7.5), "xtick.labelsize": FS(6.5),
                                "ytick.labelsize": FS(6.5), "legend.fontsize": FS(6.0)})
    globals()["FIGDIR"] = Path(a.outdir)
    matplotlib.rcParams["figure.dpi"] = a.dpi
    matplotlib.rcParams["savefig.dpi"] = a.dpi
    FIGDIR.mkdir(parents=True, exist_ok=True)

    want = list(SUBMISSION) if a.only == "all" else [x.strip()
                                                     for x in a.only.split(",") if x.strip()]
    unknown = [w for w in want if w not in ALL]
    if unknown:
        raise SystemExit(f"未知图名 {unknown}；默认清单 {SUBMISSION}，可显式点名的全部 {list(ALL)}")
    for w in want:
        ALL[w]()
        report_fit(w)
        fname = next(k for k, v in NAMES.items() if v == w)
        print(f"{w} ok  ->  {FIGDIR / fname}  dpi={a.dpi:g}  min-fs={a.min_fs:g}")
