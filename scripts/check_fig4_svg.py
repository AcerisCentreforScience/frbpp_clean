#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Structural and overflow checks on the SVG produced by make_fig4_svg.py.

Usage: python scripts/check_fig4_svg.py
"""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SVG = ROOT / "docs" / "figures" / "svg"
NS = "{http://www.w3.org/2000/svg}"

W = {"default": 0.50, "weight": 0.55}
NARROW = set("iltfj.,;:'|!()[]")
WIDE = set("mwMW@")


def text_w(s: str, size: float, bold: bool) -> float:
    tot = 0.0
    for ch in s:
        if ch in NARROW:
            f = 0.30
        elif ch in WIDE:
            f = 0.85
        elif ch.isdigit():
            f = 0.50
        elif ch == " ":
            f = 0.26
        else:
            f = 0.50
        tot += f * size
    return tot * (1.04 if bold else 1.0)


def load(p):
    root = ET.parse(p).getroot()
    els = []
    for e in root.iter():
        if e.tag.startswith(NS):
            els.append(e)
    return root, els


def by_tag(els, t):
    return [e for e in els if e.tag == NS + t]


def f(e, k, d=0.0):
    v = e.get(k)
    return float(v) if v is not None else d


def main():
    problems: list[str] = []
    for name in ("fig4_sidereal_clone.svg", "fig4_sidereal_scpma.svg"):
        p = SVG / name
        print("=" * 74)
        print(name)
        if not p.exists():
            problems.append(f"{name} 不存在")
            continue
        try:
            root, els = load(p)
        except ET.ParseError as e:
            problems.append(f"{name} XML 不合法: {e}")
            print("  XML 解析失败:", e)
            continue
        vb = root.get("viewBox")
        print(f"  viewBox           : {vb}")
        if vb != "0 0 1920 1920":
            problems.append(f"{name} viewBox 不是 0 0 1920 1920")

        rects = by_tag(els, "rect")
        lines = by_tag(els, "line")
        circs = by_tag(els, "circle")
        texts = by_tag(els, "text")
        polys = by_tag(els, "polygon")
        paths = by_tag(els, "path")
        print(f"  rect={len(rects)} line={len(lines)} circle={len(circs)} "
              f"text={len(texts)} polygon={len(polys)} path={len(paths)}")

        frames = [(f(r, "x"), f(r, "y"), f(r, "width"), f(r, "height"))
                  for r in rects if r.get("fill") == "none"]
        want = {"B": (294.0, 624.0, 1586.0, 447.0),
                "C": (179.0, 1319.0, 1701.0, 451.0)}
        for k, w in want.items():
            hit = [fr for fr in frames
                   if abs(fr[0] - w[0]) < 1.5 and abs(fr[1] - w[1]) < 1.5]
            if not hit:
                problems.append(f"{name}: 找不到面板 {k} 的外框 {w}")
            else:
                h = hit[0]
                ok = (abs(h[2] - w[2]) < 2 and abs(h[3] - w[3]) < 2)
                print(f"  面板 {k} 框        : {h}  {'OK' if ok else '尺寸不符'}")

        ins = [r for r in rects if abs(f(r, "x") - 1118.0) < 1.5
               and abs(f(r, "y") - 1333.0) < 1.5]
        print(f"  内嵌框            : {'找到' if ins else '缺失'}"
              f"{'' if ins else ''}")
        if not ins:
            problems.append(f"{name}: 内嵌插画框缺失")

        bars = [r for r in rects
                if re.match(r"^#387DC8$", (r.get("fill") or ""))
                and abs(f(r, "width") - 23.4) < 0.6]
        print(f"  柱                : {len(bars)} 根")
        if len(bars) < 55:
            problems.append(f"{name}: 柱只有 {len(bars)} 根，应为 61")

        bl = [c for c in circs if c.get("fill") == "#387DC8"]
        gr = [c for c in circs if c.get("fill") == "#A0A0A0"]
        print(f"  散点              : 蓝 {len(bl)}  灰 {len(gr)}")
        bad = [c for c in bl + gr
               if not (294 <= f(c, "cx") <= 1880 and 624 <= f(c, "cy") <= 1071)]
        if bad:
            problems.append(f"{name}: {len(bad)} 个散点落在面板 b 框外")

        over = []
        for t in texts:
            if t.get("transform"):
                continue
            size = f(t, "font-size", 20.0)
            s = "".join(t.itertext())
            if not s.strip():
                continue
            bold = (t.get("font-weight") == "bold")
            wpx = text_w(s, size, bold)
            x = f(t, "x")
            anc = t.get("text-anchor", "start")
            x0 = x if anc == "start" else (x - wpx / 2 if anc == "middle"
                                           else x - wpx)
            x1 = x0 + wpx
            y = f(t, "y")
            if x0 < -150 or x1 > 2070 or y > 1990:
                over.append((s[:38], round(x0, 1), round(x1, 1), round(y, 1)))
        if over:
            print("  文字越界（估算）:")
            for s, a, b, yy in over:
                print(f"      x {a:7.1f}..{b:7.1f}  y={yy:6.1f}  «{s}»")
        else:
            print("  文字越界（估算）  : 无")

        left = [t for t in texts
                if t.get("text-anchor") == "middle" and abs(f(t, "x") - 131.0) < 1.0]
        if len(left) != 7:
            problems.append(f"{name}: 面板 b 左侧标签应为 7 行，实为 {len(left)}")
        for t in left:
            s = "".join(t.itertext())
            size = f(t, "font-size", 40.0)
            half = text_w(s, size, False) / 2
            if 131.0 - half < 2:
                problems.append(f"{name}: 左侧标签 «{s}» 超出画布左缘")
            if 131.0 + half > 292:
                problems.append(f"{name}: 左侧标签 «{s}» 压到面板 b 边框")
        print(f"  面板 b 左侧标签    : {len(left)} 行，居中 x=131，"
              f"最宽半宽 {max((text_w(''.join(t.itertext()), f(t,'font-size',40.0), False)/2) for t in left):.1f} px")

        print()

    print("=" * 74)
    if problems:
        print(f"发现 {len(problems)} 处问题：")
        for x in problems:
            print("  -", x)
        return 1
    print("结构自检全部通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
