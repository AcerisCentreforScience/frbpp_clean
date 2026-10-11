#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Rewrite in-text citations as natbib citep/citet/citealt by grammatical role."""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEX = ROOT / "docs" / "submission" / "raa_manuscript.tex"

EXACT = {
    "NE2001": "citealt",
    "Jiang et al.\\ 2020": "citealt",
    "Outriggers 2025": "citealt",
    "Chronoscope 2026": "citealt",
    "sample in": "citet",
    "^\\cite{pandhi2024}": "citet",
}
LOOKBEHIND = 60


def convert(src: str) -> tuple[str, dict]:
    tally = {"citep": 0, "citet": 0, "citealt": 0}
    out = []
    pos = 0
    for m in re.finditer(r"\\cite\{", src):
        ctx = src[max(0, m.start() - LOOKBEHIND):m.start()]
        mode = "citep"
        for feat, mm in EXACT.items():
            if feat.startswith("^"):
                line_head = ctx.rfind("\n") + 1
                if ctx[line_head:].strip() == "":
                    mode = mm
                    break
            elif feat in ctx:
                mode = mm
                break
        tally[mode] += 1
        out.append(src[pos:m.start()])
        out.append("\\" + mode + "{")
        pos = m.end()
    out.append(src[pos:])
    return "".join(out), tally


def main() -> int:
    src = TEX.read_text(encoding="utf-8")
    before = len(re.findall(r"\\cite\{", src))
    new, tally = convert(src)
    after = len(re.findall(r"\\cite\{", new))
    print(f"\\cite 站点 {before} 处 -> 剩余 {after} 处")
    print(f"分派：citep {tally['citep']} / citet {tally['citet']} / citealt {tally['citealt']}"
          f"   合计 {sum(tally.values())}")
    total = sum(tally.values())
    if total != before or after != 0:
        print("!! 分派数与站点数不一致，或有站点未被替换，拒绝写入")
        return 1
    if "--write" in sys.argv:
        bak = TEX.with_suffix(".tex.bak_before_raa_cite")
        if not bak.exists():
            shutil.copy2(TEX, bak)
            print(f"已备份 -> {bak.name}")
        TEX.write_text(new, encoding="utf-8")
        print("已写入")
    else:
        print("（预览模式，未写入；加 --write 生效）")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    raise SystemExit(main())
