#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Point the six submission audits at any manuscript.

Usage: python scripts/audit_submission.py [manuscript_stem]
"""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUB = ROOT / "docs" / "submission"
sys.path.insert(0, str(ROOT / "scripts"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import raa_readiness as rr  # noqa: E402


def main() -> int:
    stem = sys.argv[1] if len(sys.argv) > 1 else "raa_manuscript"
    tex = SUB / f"{stem}.tex"
    pdf = SUB / "_build" / f"{stem}.pdf"
    if not tex.exists():
        raise SystemExit(f"找不到稿件 {tex}")
    if not pdf.exists():
        raise SystemExit(f"找不到编译产物 {pdf} —— 先跑 xelatex 两遍，再审计")
    rr.TEX = tex
    rr.BUILD_PDF = pdf
    print(f"== 审计目标：{tex.name} / {pdf.name} ==")
    return rr.main()


if __name__ == "__main__":
    raise SystemExit(main())
