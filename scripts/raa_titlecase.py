#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert subsection headings and table captions to word-by-word title case."""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEX = ROOT / "docs" / "submission" / "raa_manuscript.tex"

SMALL = {"a", "an", "and", "as", "at", "but", "by", "for", "from", "in", "into",
         "of", "on", "or", "the", "to", "with", "vs", "via", "per", "nor", "yet"}


def cap_word(w: str, capital: bool) -> str:
    if not w or not re.search(r"[A-Za-z]", w):
        return w
    segs = w.split("-")
    out = []
    for k, s in enumerate(segs):
        m = re.match(r"([^A-Za-z]*)([A-Za-z]*)(.*)", s)
        pre, tok, rest = m.group(1), m.group(2), m.group(3)
        if not tok:
            out.append(s)
            continue
        if tok.isupper() and len(tok) >= 2:
            out.append(s)
            continue
        if tok.lower() in SMALL and not capital:
            out.append(pre + tok.lower() + rest)
        else:
            out.append(pre + tok[0].upper() + tok[1:] + rest)
    return "-".join(out)


def titlecase(text: str) -> str:
    out, i, n = [], 0, len(text)
    buf = ""
    need_cap = {"start": True, "after": False}

    def flush() -> None:
        nonlocal buf
        if not buf:
            return
        parts = re.split(r"(\s+)", buf)
        for j, p in enumerate(parts):
            if p == "" or p.isspace():
                continue
            capital = need_cap["start"] or need_cap["after"]
            parts[j] = cap_word(p, capital)
            need_cap["start"] = False
            need_cap["after"] = bool(re.search(r"[:;.!]\s*$|[—–-]{2,}\s*$", p))
        out.append("".join(parts))
        buf = ""
        if re.search(r"[:;.]\s*$", "".join(out)):
            need_cap["after"] = True

    while i < n:
        c = text[i]
        if c == "$":
            flush()
            j = text.find("$", i + 1)
            j = n - 1 if j == -1 else j
            out.append(text[i:j + 1])
            i = j + 1
            continue
        if c == "\\":
            flush()
            m = re.match(r"\\[A-Za-z]+|\\.", text[i:])
            tok = m.group(0)
            out.append(tok)
            i += len(tok)
            if i < n and text[i] == "{":
                depth, j = 0, i
                while j < n:
                    if text[j] == "{":
                        depth += 1
                    elif text[j] == "}":
                        depth -= 1
                        if depth == 0:
                            break
                    j += 1
                out.append(text[i:j + 1])
                i = j + 1
            continue
        if c == "{":
            flush()
            depth, j = 0, i
            while j < n:
                if text[j] == "{":
                    depth += 1
                elif text[j] == "}":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            out.append(text[i:j + 1])
            i = j + 1
            continue
        if c == "%":
            flush()
            j = text.find("\n", i)
            j = n if j == -1 else j
            out.append(text[i:j])
            i = j
            continue
        buf += c
        i += 1
    flush()
    return "".join(out)


def spans_of(text: str, macro: str) -> list[tuple[int, int, str]]:
    res = []
    for m in re.finditer(r"\\" + macro + r"\s*\{", text):
        i = m.end() - 1
        depth, j = 0, i
        while j < len(text):
            if text[j] in "{}":
                depth += 1 if text[j] == "{" else -1
                if depth == 0:
                    break
            j += 1
        res.append((i + 1, j, text[i + 1:j]))
    return res


def main() -> int:
    src = TEX.read_text(encoding="utf-8")
    body_start = src.find("\\begin{document}")
    body = src[body_start:]

    targets = [(s, e, c) for s, e, c in spans_of(body, "subsection")]
    targets += [(s, e, c) for s, e, c in spans_of(body, "subsubsection")]
    caps = []
    for s, e, c in spans_of(body, "caption"):
        t = body.rfind("\\begin{table", 0, s)
        f = body.rfind("\\begin{figure", 0, s)
        if t > f:
            caps.append((s, e, c))
    targets += caps

    edits = []
    problems: list[str] = []
    print(f"待改 {len(targets)} 条（subsection+subsubsection {len(targets)-len(caps)} / 表题 {len(caps)}）\n"
          + "=" * 100)
    for s, e, content in sorted(targets):
        new = titlecase(content)
        if new != content:
            a = len(re.findall(r"[A-Za-z]", content))
            b = len(re.findall(r"[A-Za-z]", new))
            if b < a:
                problems.append(f"字母数从 {a} 掉到 {b}：{content[:50]!r}")
                continue
            edits.append((body_start + s, body_start + e, content, new))
            print(f"-  {content}")
            print(f"+  {new}")
    print("=" * 100)
    if problems:
        print(f"!! 自检失败 {len(problems)} 条，拒绝写入：")
        for p in problems:
            print(f"   {p}")
        return 1
    if "--write" in sys.argv and edits:
        bak = TEX.with_suffix(".tex.bak_before_raa_titlecase")
        if not bak.exists():
            shutil.copy2(TEX, bak)
            print(f"已备份 -> {bak.name}")
        chars = list(src)
        for s, e, old, new in sorted(edits, reverse=True):
            assert src[s:e] == old, f"位置漂移：{src[s:e][:40]!r} != {old[:40]!r}"
            src = src[:s] + new + src[e:]
        TEX.write_text(src, encoding="utf-8")
        print(f"已写入 {len(edits)} 条")
    elif not edits:
        print("无需修改")
    else:
        print("（预览模式，未写入；加 --write 生效）")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    raise SystemExit(main())
