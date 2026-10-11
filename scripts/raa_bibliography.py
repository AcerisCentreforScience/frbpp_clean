#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert the embedded bibliography from numeric to author-year style."""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEX = ROOT / "docs" / "submission" / "raa_manuscript.tex"

JOURNAL_MAP = {
    "Astrophys J Suppl Ser": r"\apjs",
    "Astrophys J Lett": r"\apjl",
    "Astron Astrophys Rev": r"\aapr",
    "Astron Astrophys": r"\aap",
    "Astrophys J": r"\apj",
    "Mon Not R Astron Soc": r"\mnras",
    "Res Astron Astrophys": "RAA",
    "Nature": r"\nat",
    "Science": "Science",
    "Universe": "Universe",
}
JOURNAL_KEYS = sorted(JOURNAL_MAP, key=len, reverse=True)

COLLABS = ("CHIME/FRB Collaboration", "Planck Collaboration")

PARTICLES = {"van", "de", "der", "den", "di", "da", "del", "la", "le", "dos", "du", "st"}


def format_initials(chunk: str) -> str:
    """'Amiri M' -> 'Amiri M.'; 'Macquart J-P' -> 'Macquart J.-P.'"""
    for coll in COLLABS:
        if chunk.startswith(coll):
            return chunk
    parts = chunk.split()
    if len(parts) < 2:
        return chunk
    surname = " ".join(parts[:2]) if parts[0].lower() in PARTICLES and len(parts) > 2 else parts[0]
    rest = parts[2:] if surname != parts[0] else parts[1:]
    ini = []
    for tok in rest:
        ini.append("-".join(p + "." for p in tok.split("-") if p)
                   if "-" in tok else tok + ".")
    return (surname + " " + " ".join(ini)).strip() if ini else surname


def author_block(authors: str) -> str:
    out = []
    for chunk in authors.split(", "):
        c = chunk.strip()
        if c.lower() in ("et al", "et al."):
            out.append("et al.")
        else:
            out.append(format_initials(c))
    s = ", ".join(out)
    return s.replace(", et al.", " et al.")


def surname_of(authors: str) -> str:
    for coll in COLLABS:
        if authors.startswith(coll):
            return coll
    parts = authors.split(",")[0].split()
    if parts and parts[0].lower() in PARTICLES and len(parts) > 2:
        return " ".join(parts[:2])
    return parts[0]


def n_authors(authors: str) -> int:
    toks = [t.strip() for t in authors.split(",")]
    return len([t for t in toks if t.lower() not in ("et al", "et al.")])


def parse_entry(raw: str) -> dict:
    raw = " ".join(raw.split())
    m = re.search(r"(?:,\s*)?arXiv:([A-Za-z0-9\.\-/]+),\s*(\d{4})\s*$", raw)
    if m:
        head = raw[: m.start()].rstrip(", ")
        authors, _, title = head.partition(". ")
        return {"kind": "preprint", "authors": authors.strip(), "title": title.strip(),
                "arxiv": m.group(1), "year": m.group(2), "raw": raw}
    for j in JOURNAL_KEYS:
        pat = re.compile(re.escape(j) + r",\s*(\d{4}),\s*([\w\-\.]+):\s*([\w\-\.\, ]+?)\s*$")
        m = pat.search(raw)
        if m:
            head = raw[: m.start()].rstrip()
            authors, _, title = head.partition(". ")
            return {"kind": "journal", "authors": authors.strip(), "title": title.strip(),
                    "journal": j, "year": m.group(1), "vol": m.group(2),
                    "pages": m.group(3).strip(), "raw": raw}
    raise ValueError(f"无法解析条目：{raw[:90]}")


def make_label(e: dict, first_surname: str) -> str:
    et = " et al." if n_authors(e["authors"]) >= 3 and not first_surname.endswith("Collaboration") \
        else ""
    if first_surname in COLLABS:
        et = ""
    if n_authors(e["authors"]) == 2 and not first_surname in COLLABS:
        second = e["authors"].split(", ")[1]
        sec_sur = surname_of(second) if "," not in second else second.split(",")[0]
        return f"{first_surname} \\& {sec_sur}({e['year']})"
    return f"{first_surname}{et}({e['year']})"


def build(entries: list[dict]) -> list[dict]:
    recs = []
    for e in entries:
        sur = surname_of(e["authors"])
        recs.append({**e, "surname": sur, "label": make_label(e, sur)})
    recs.sort(key=lambda r: (r["surname"].lower(), r["year"], r["title"].lower()))
    seen: dict[str, list[dict]] = {}
    for r in recs:
        seen.setdefault(re.sub(r"[()]", "", r["label"]).lower(), []).append(r)
    for base, group in seen.items():
        if len(group) > 1:
            for i, r in enumerate(group):
                suffix = "abcdefgh"[i]
                r["label"] = re.sub(r"\((\d{4})\)$", rf"(\1{suffix})", r["label"])
    return recs


def render(r: dict) -> str:
    a = author_block(r["authors"])
    if r["kind"] == "preprint":
        return f"{a}, {r['year']}, arXiv:{r['arxiv']}"
    jr = JOURNAL_MAP[r["journal"]]
    return f"{a}, {r['year']}, {jr}, {r['vol']}, {r['pages']}"


def numbers_preserved(r: dict, new_body: str) -> list[str]:
    want = [r["year"]]
    if r["kind"] == "preprint":
        want.append(r["arxiv"])
    else:
        want += [r["vol"], r["pages"]]
    return sorted(x for x in want if x and x not in new_body)


def main() -> int:
    src = TEX.read_text(encoding="utf-8")
    mb = re.search(r"\\begin\{thebibliography\}\{[^}]*\}(.*?)\\end\{thebibliography\}", src, re.S)
    if not mb:
        print("找不到 thebibliography")
        return 2
    bib = mb.group(1)
    pat = re.compile(r"\\bibitem(\[[^\]]*\])?\{([^}]*)\}")
    hits = list(pat.finditer(bib))
    entries = []
    for i, m in enumerate(hits):
        end = hits[i + 1].start() if i + 1 < len(hits) else len(bib)
        raw = bib[m.end():end].strip()
        entries.append({"key": m.group(2), "raw": " ".join(raw.split())})
        entries[-1].update(parse_entry(" ".join(raw.split())))
    recs = build(entries)

    print(f"条目 {len(recs)} 条；键序 = 字母序结果\n" + "=" * 100)
    problems: list[str] = []
    lines: list[str] = []
    for r in recs:
        body = render(r)
        missing = numbers_preserved(r, body)
        if missing:
            problems.append(f"{r['key']}: 新条目里丢了 {missing}")
        for chunk in r["authors"].split(", "):
            sur = chunk.split()[0]
            if sur.lower() not in body.lower() and sur.lower() not in ("et", "van"):
                problems.append(f"{r['key']}: 作者 {sur!r} 在新条目里丢了")
        lines.append(f"\\bibitem[{r['label']}]{{{r['key']}}} {body}")
        print(f"  {r['key']:18s} [{r['label']}]")
        print(f"      新 {body}")
        print(f"      丢 {r['title'][:76]}")
    labels = [re.sub(r"\s+", " ", re.sub(r"[()]", "", l)) for l in lines]
    dup = {l for l in labels if labels.count(l) > 1}
    if dup:
        problems.append(f"标签碰撞：{sorted(dup)}")

    print("=" * 100)
    if problems:
        print(f"断言失败 {len(problems)} 项，拒绝写入：")
        for p in problems:
            print(f"   !! {p}")
        return 1
    print("断言：数字与作者全部保留、标签无碰撞 -> 通过")

    new_block = ("\\begin{thebibliography}{99}\n\n"
                 + "\n\n".join(lines) + "\n\n\\end{thebibliography}")
    if "--write" in sys.argv:
        bak = TEX.with_suffix(".tex.bak_before_raa_bib")
        if not bak.exists():
            shutil.copy2(TEX, bak)
            print(f"已备份 -> {bak.name}")
        head = src[: mb.start()]
        tail = src[mb.end():]
        TEX.write_text(head + new_block + tail, encoding="utf-8")
        print(f"已写入 {TEX.name}（{len(lines)} 条）")
    else:
        print("（预览模式，未写入；加 --write 生效）")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    raise SystemExit(main())
