#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Extract text from a compiled PDF for manual reading.

Usage: python scripts/dump_pdf_text.py <pdf>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pypdfium2 as pdfium
except ImportError:  # pragma: no
    sys.exit("需要 pypdfium2：pip install -i https://pypi.org/simple pypdfium2")


def page_texts(pdf_path: Path) -> list[str]:
    doc = pdfium.PdfDocument(str(pdf_path))
    return [doc[i].get_textpage().get_text_range() for i in range(len(doc))]


def main() -> int:
    ap = argparse.ArgumentParser(description="抽取 PDF 文字用于人工复核")
    ap.add_argument("--pdf", required=True, type=Path)
    ap.add_argument("--front", action="store_true",
                    help="只打印首页前置部分（标题/作者/摘要/关键词/PACS）")
    ap.add_argument("--front-chars", type=int, default=1600,
                    help="--front 打印的字符数（默认 1600）")
    ap.add_argument("--find", action="append", default=[],
                    help="打印含该子串的页码（可重复）")
    ap.add_argument("--max-pages", type=int, default=0,
                    help=">0 时只打印前 N 页的全文")
    args = ap.parse_args()

    if not args.pdf.exists():
        sys.exit(f"找不到 {args.pdf}")
    pages = page_texts(args.pdf)
    print(f"[{args.pdf.name}] 共 {len(pages)} 页\n")

    if args.front:
        print("=" * 70)
        print("首页前置部分（请逐字核对标题、作者、摘要、关键词、PACS）")
        print("=" * 70)
        print(pages[0][: args.front_chars])
        print()
        return 0

    if args.find:
        for needle in args.find:
            hits = [i + 1 for i, t in enumerate(pages) if needle in t]
            print(f"  {needle!r} -> 页码 {hits if hits else '未找到'}")
        print()
        return 0

    upto = args.max_pages or len(pages)
    for i, t in enumerate(pages[:upto]):
        print("=" * 70)
        print(f"=== p/{len(pages)}  (第 {i + 1} 页) ===")
        print("=" * 70)
        print(t)
    return 0


if __name__ == "__main__":
    sys.exit(main())
