"""Audit of figure and table float references."""
from __future__ import annotations

import re
import sys
from pathlib import Path

DEFAULT_TEX = Path(__file__).resolve().parents[1] / "docs" / "submission" / "raa_manuscript.tex"


def strip_comments(text: str) -> str:
    out = []
    for line in text.splitlines():
        idx = None
        for i, ch in enumerate(line):
            if ch == "%" and (i == 0 or line[i - 1] != "\\"):
                idx = i
                break
        out.append(line if idx is None else line[:idx])
    return "\n".join(out)


def line_of(text: str, pos: int) -> int:
    return text[:pos].count("\n") + 1


def audit(text: str, prefix: str) -> tuple[list[tuple], int]:
    labels = [(m.group(1), line_of(text, m.start()))
              for m in re.finditer(r"\\label\{(" + re.escape(prefix) + r"[^}]*)\}", text)]
    rows, issues = [], 0
    prev_cite = 0
    for lab, env_line in labels:
        m = re.search(r"\\ref\{" + re.escape(lab) + r"\}", text)
        if m is None:
            rows.append((lab, env_line, None, "UNCITED"))
            issues += 1
            continue
        cite_line = line_of(text, m.start())
        notes = []
        if cite_line > env_line:
            notes.append("cite-after-env")
        if cite_line < prev_cite:
            notes.append("ORDER")
        if notes:
            issues += 1
        rows.append((lab, env_line, cite_line, " ".join(notes)))
        prev_cite = cite_line
    return rows, issues


def main(argv: list[str]) -> int:
    tex = Path(argv[1]) if len(argv) > 1 else DEFAULT_TEX
    if not tex.exists():
        print(f"NOT FOUND: {tex}")
        return 1
    text = strip_comments(tex.read_text(encoding="utf-8"))

    total = 0
    for kind, prefix in (("FIGURES", "fig:"), ("TABLES", "tab:")):
        rows, issues = audit(text, prefix)
        print(f"===== {kind} ({len(rows)}) =====")
        for lab, env_line, cite_line, status in rows:
            print("  %-20s env L%-6d cite L%-8s %s"
                  % (lab, env_line, cite_line if cite_line else "-", status))
        print(f"  issues: {issues}\n")
        total += issues

    print("=" * 60)
    print("FLOAT AUDIT PASSED" if total == 0 else f"FLOAT AUDIT: {total} issue(s)")
    return 0 if total == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
