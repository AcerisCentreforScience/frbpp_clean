"""Bibliography audit: bibitem-to-citation cross-check, ordering, field extraction.

Usage: python scripts/audit_refs.py [path/to/manuscript.tex]
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
CITE_RE = re.compile(r"\\cite[tpalt]*\{([^}]*)\}")
BIB_RE = re.compile(r"^\s*\\bibitem(?:\[([^\]]*)\])?\{([^}]+)\}")

TEX_DEFAULT = ROOT / "docs/submission/raa_manuscript.tex"
tex_path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else TEX_DEFAULT
tex = tex_path.read_text(encoding="utf-8")
lines = tex.splitlines()

bib, labels = [], {}
for i, ln in enumerate(lines, 1):
    m = BIB_RE.match(ln)
    if m:
        note, key = m.group(1), m.group(2)
        bib.append(key)
        labels[key] = note or ""

cited = []
for m in CITE_RE.finditer(tex):
    for k in m.group(1).split(","):
        k = k.strip()
        if k and k not in cited:
            cited.append(k)

print(f"bibitem 条数 = {len(bib)}    正文引用键数 = {len(cited)}")
print("重复 bibitem 键:", sorted({k for k in bib if bib.count(k) > 1}) or "无")
print("定义了但从未被引用:", [k for k in bib if k not in cited] or "无")
print("引用但无对应 bibitem:", [k for k in cited if k not in bib] or "无")


def sortkey(note: str) -> str:
    head = re.split(r"\(", note)[0]
    head = re.sub(r"\\&", "and", head)
    head = re.sub(r"\\[a-zA-Z]+", "", head)
    return " ".join(head.split()).lower()


order = [sortkey(labels[k]) for k in bib if labels[k]]
inv = [(bib[i - 1], order[i - 1], bib[i], order[i])
       for i in range(1, len(order)) if order[i] < order[i - 1]]
print(f"\n字母序：{'已排好' if not inv else '倒序 ' + str(len(inv)) + ' 处'}")
for a, ka, b, kb in inv[:10]:
    print(f"   {a!r} ({ka}) 排在 {b!r} ({kb}) 之前")
missing_label = [k for k in bib if not labels[k]]
if missing_label:
    print("   无标签的条目（作者-年份制不允许）:", missing_label)

print("\n===== 条目字段（供逐条外部核对）=====")
JOURNAL_MACROS = {
    "apj": "ApJ", "apjl": "ApJL", "apjs": "ApJS", "aap": "A&A", "aaps": "A&AS",
    "aapr": "A&ARv", "mnras": "MNRAS", "nat": "Nature", "pasp": "PASP",
    "ssr": "Space Sci. Rev.", "raa": "RAA", "sci": "Science", "aj": "AJ",
    "araa": "ARA&A", "prd": "Phys. Rev. D", "prl": "Phys. Rev. Lett.",
}
for i, k in enumerate(bib, 1):
    m = next((j for j, ln in enumerate(lines, 1) if BIB_RE.match(ln) and BIB_RE.match(ln).group(2) == k), 0)
    txt = lines[m - 1] if m else ""
    txt = re.sub(r"^\s*\\bibitem(\[[^\]]*\])?\{[^}]+\}\s*", "", txt)
    txt = re.sub(r"\\(" + "|".join(JOURNAL_MACROS) + r")\b",
                 lambda x: JOURNAL_MACROS[x.group(1)], txt)
    txt = re.sub(r"[{}\\]", "", txt).replace("  ", " ").strip()
    print(f"[{i:>2}] ({k}) L{m}: {txt[:150]}")
if "--dump-only" in sys.argv:
    pass
