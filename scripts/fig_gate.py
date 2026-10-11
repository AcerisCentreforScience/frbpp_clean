"""Shared layout gate: auto-fit of text blocks and overprint detection.
Used by the figure generators.
"""
from __future__ import annotations

import numpy as np
from matplotlib.text import Text

MIN_FS_DEFAULT = 7.0


def fit(ax, x, y, s, size, max_w=0.30, xy=None, min_fs=MIN_FS_DEFAULT, **kw):
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    axw = ax.get_window_extent(renderer=r).width
    kw.setdefault("clip_on", False)
    if xy is not None:
        t = ax.annotate(s, xy=xy, xytext=(x, y), fontsize=size, **kw)
    else:
        t = ax.text(x, y, s, fontsize=size, **kw)
    for _ in range(30):
        if t.get_window_extent(renderer=r).width <= max_w * axw:
            break
        size = max(min_fs * 0.85, size - 0.25)
        t.set_fontsize(size)
    return t


def all_axes(axes):
    """Every axes in a figure, including insets.

    matplotlib puts axes created by Axes.inset_axes() in the parent's child_axes and
    not in fig.axes, so walking fig.axes alone leaves inset titles and labels unchecked.
    """
    for ax in axes:
        yield ax
        yield from all_axes(ax.child_axes)


def _items(fig, renderer):
    out = []
    for ax in all_axes(fig.axes):
        lg = ax.get_legend()
        tagged = ([("text", t) for t in ax.texts]
                  + ([("tick", t) for t in list(ax.get_xticklabels())
                                     + list(ax.get_yticklabels())] if ax.axison else [])
                  + ([("axis", t) for t in (ax.xaxis.label, ax.yaxis.label, ax.title)]
                     if ax.axison else [])
                  + [("legend", t) for t in (list(lg.get_texts()) if lg else [])])
        for kind, t in tagged:
            s = t.get_text().strip() if hasattr(t, "get_text") else ""
            if not s:
                continue
            try:
                out.append((kind, s.replace("\n", " ")[:36],
                            t.get_window_extent(renderer=renderer)))
            except Exception:
                pass
    return out


def report_text_over_data(fig, name):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    hits = []
    for ax in all_axes(fig.axes):
        pts = []
        for ln in ax.get_lines():
            if not ln.get_linestyle().startswith("None") and ln.get_marker() == "None":
                xy = ln.get_xydata()
                if len(xy) > 1:
                    pts.append(ax.transData.transform(xy))
        for pc in ax.collections:
            off = pc.get_offsets()
            if len(off):
                pts.append(ax.transData.transform(off))
        if not pts:
            continue
        allpts = np.vstack(pts)
        for t in ax.texts:
            s = t.get_text().strip()
            if len(s) <= 4:
                continue
            bp = t.get_bbox_patch()
            if bp is not None:
                fc = bp.get_facecolor()
                if (fc[3] >= 0.7 and not (abs(fc[0] - 1) < 1e-6
                                          and abs(fc[1] - 1) < 1e-6
                                          and abs(fc[2] - 1) < 1e-6 and fc[3] == 0)):
                    continue
            tb = Text.get_window_extent(t, renderer=r)
            inside = ((allpts[:, 0] > tb.x0 + 1) & (allpts[:, 1] < tb.x1 - 1)
                      & (allpts[:, 1] > tb.y0 + 1) & (allpts[:, 1] < tb.y1 - 1))
            if inside.sum() > 8:
                hits.append((s.replace("\n", " ")[:34], int(inside.sum())))
    for s, n in hits:
        print(f"   !! {name} 文字被数据穿过: {s!r}（{n} 个数据点落在字框内）")
    return len(hits)


def report_collisions(fig, name, also_check_axes_edge=True):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    items = _items(fig, r)
    hits = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            ka, sa, a = items[i]
            kb, sb, b = items[j]
            if ka in ("tick", "axis") and kb in ("tick", "axis"):
                continue
            if len(sa.strip()) <= 1 or len(sb.strip()) <= 1:
                continue
            ox = max(0.0, min(a.x1, b.x1) - max(a.x0, b.x0))
            oy = max(0.0, min(a.y1, b.y1) - max(a.y0, b.y0))
            if ox * oy > 0.15 * min(a.width * a.height, b.width * b.height):
                hits.append((sa, sb))
    for x, y in hits:
        print(f"   !! {name} 叠字: {x!r} <-> {y!r}")

    if also_check_axes_edge:
        fb = fig.get_window_extent(renderer=r)
        seen = set()
        for kind, s, bb in items:
            if kind != "text" or len(s) <= 2 or s in seen:
                continue
            if bb.x1 > fb.x1 + 1 or bb.x0 < fb.x0 - 1 or bb.y1 > fb.y1 + 1 or bb.y0 < fb.y0 - 1:
                seen.add(s)
                print(f"   !! {name} 文字越出画布: {s!r}")
                hits.append((s, "<canvas edge>"))
    return len(hits) + report_text_over_data(fig, name)
