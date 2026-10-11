"""Consistency check for the RAA LaTeX manuscript.

Usage: python scripts/check_tex.py [path/to/manuscript.tex]
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

TEX = Path(__file__).resolve().parents[1] / "docs/submission/raa_manuscript.tex"
if len(sys.argv) > 1:
    TEX = Path(sys.argv[1])
FIGDIRS = [
    Path(__file__).resolve().parents[1] / "docs/submission",
    Path(__file__).resolve().parents[1] / "docs/figures",
    Path(__file__).resolve().parents[1] / "docs/figures/retext/en",
]


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


def main() -> int:
    raw = TEX.read_text(encoding="utf-8")
    body = strip_comments(raw)
    problems: list[str] = []

    cites: set[str] = set()
    for m in re.finditer(r"\\cite\{([^}]*)\}", body):
        for key in m.group(1).split(","):
            key = key.strip()
            if key:
                cites.add(key)
    bibs = set(re.findall(r"\\bibitem\{([^}]*)\}", body))

    missing_bib = sorted(cites - bibs)
    unused_bib = sorted(bibs - cites)
    if missing_bib:
        problems.append(f"[1] \\cite keys with no \\bibitem: {missing_bib}")
    if unused_bib:
        problems.append(f"[1] \\bibitem never cited (ok if intentional): {unused_bib}")

    refs: set[str] = set()
    for m in re.finditer(r"\\(?:ref|eqref)\{([^}]*)\}", body):
        refs.add(m.group(1).strip())
    labels = set(re.findall(r"\\label\{([^}]*)\}", body))
    missing_lbl = sorted(refs - labels)
    unused_lbl = sorted(labels - refs)
    if missing_lbl:
        problems.append(f"[2] \\ref keys with no \\label: {missing_lbl}")
    if unused_lbl:
        print(f"INFO [2] {len(unused_lbl)} labels not referenced (normal, for numbering only)")

    ph = re.findall(r"Placeholder[^\n]*", raw, flags=re.IGNORECASE)
    if ph:
        problems.append(f"[3] leftover placeholders: {ph}")

    tm = re.search(r"\\title\{(.*?)\n\}", raw, flags=re.DOTALL)
    if not tm:
        problems.append("[4] could not locate \\title{}")
    else:
        t = re.sub(r"[%\\{}]", " ", tm.group(1))
        words = [w for w in re.split(r"\s+", t) if w]
        print(f"TITLE words: {len(words)}  (internal limit 20; RAA only asks 'concise yet informative')")
        print(f"TITLE has colon: {':' in tm.group(1)}  has question: {'?' in tm.group(1)}")
        print(f"TITLE -> {' '.join(words)}")
        if len(words) > 20:
            problems.append(f"[4] title exceeds 20 words ({len(words)})")
        if ":" in tm.group(1) or "?" in tm.group(1):
            problems.append("[4] title contains a colon or question mark")

    am = None
    km = re.search(r"\\abstract\s*\{", raw)
    if km:
        i = km.end() - 1
        depth, j = 0, i
        while j < len(raw):
            if raw[j] == "{":
                depth += 1
            elif raw[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if depth == 0:
            am = raw[i + 1:j]
    if am is None:
        m2 = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", raw, flags=re.DOTALL)
        am = m2.group(1) if m2 else None
    if am is None:
        problems.append("[5] could not locate the abstract")
    else:
        a = re.sub(r"\\keywords\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}", " ", am)
        a = strip_comments(a)
        a = re.sub(r"\$[^$]*\$", " MATH ", a)
        a_clean = re.sub(r"\\[a-zA-Z]+\*?", " ", a)
        a_clean = re.sub(r"[{}$\\]", " ", a_clean)
        words = [w for w in re.split(r"\s+", a_clean) if re.search(r"[A-Za-z0-9]", w)]
        print(f"ABSTRACT words: {len(words)}  (source lower bound; RAA limit 200)")
        if len(words) > 200:
            problems.append(f"[5] abstract lower bound already exceeds 200 words ({len(words)})")

    figs = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", body)
    for f in figs:
        p = Path(f)
        found = None
        if p.is_absolute():
            if p.exists():
                found = p
        else:
            cand = TEX.parent / p
            if cand.exists():
                found = cand
        if found is None:
            for d in FIGDIRS:
                cand = d / p.name
                if cand.exists():
                    found = cand
                    print(f"    (fallback by filename: {f} -> {cand})")
                    break
        status = "OK" if found else "MISSING"
        print(f"FIG {f}: {status}" + (f" -> {found}" if found else ""))
        if not found:
            problems.append(f"[6] figure not found: {f}")

    cjk = re.findall(r"[\u4e00-\u9fff]+", body)
    if cjk:
        problems.append(f"[7] CJK characters present in body: {cjk}")

    begins = re.findall(r"\\begin\{([^}]*)\}", body)
    ends = re.findall(r"\\end\{([^}]*)\}", body)
    from collections import Counter
    cb, ce = Counter(begins), Counter(ends)
    for env in set(cb) | set(ce):
        if cb[env] != ce[env]:
            problems.append(f"[8] environment unbalanced: {env} begin={cb[env]} end={ce[env]}")
    print(f"ENVIRONMENTS: {dict(cb)}")

    lines = raw.splitlines()
    glued: list[str] = []
    for i in range(len(lines) - 1):
        a = lines[i].rstrip()
        if not re.search(r"[A-Za-z0-9\u4e00-\u9fff\)\}\$\.\,\;]%$", a):
            continue
        nxt = lines[i + 1].strip()
        if not nxt or nxt.startswith("%"):
            continue
        if re.match(r"^[A-Za-z0-9\u4e00-\u9fff]", nxt):
            glued.append(f"L{i+1}: ...{a.strip()[-45:]} ++ {nxt[:45]}")
    if glued:
        problems.append(f"[9] trailing '%' will glue words together: {glued}")
    else:
        print("GLUE-CHECK: no trailing-'%' word gluing")

    print("\n" + "=" * 60)
    if problems:
        print(f"PROBLEMS ({len(problems)}):")
        for p in problems:
            print("  -", p)
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
