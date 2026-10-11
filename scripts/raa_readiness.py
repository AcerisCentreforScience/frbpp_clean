#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Submission-compliance auditor for the journal's requirements."""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEX = ROOT / "docs" / "submission" / "raa_manuscript.tex"

ABSTRACT_MAX_WORDS = 200
KEYWORDS_MAX = 6
BUILD_PDF = ROOT / "docs" / "submission" / "_build" / "raa_manuscript.pdf"


def rendered_abstract_words(pdf=BUILD_PDF):
    import glob
    import os
    path = None
    if os.path.exists(str(pdf)):
        path = str(pdf)
    else:
        cand = glob.glob(str(pdf.parent / "*.pdf"))
        path = cand[0] if len(cand) == 1 else None
    if path is None:
        return None, "编译产物不存在"
    try:
        import pymupdf
    except ImportError:
        return None, "本机无 pymupdf"
    doc = pymupdf.open(path)
    t = doc[0].get_text()
    doc.close()
    i = t.find("Abstract")
    if i < 0:
        return None, "首页找不到 Abstract"
    j = t.find("Key words", i)
    if j < 0:
        j = t.find("1. INTRODUCTION", i)
    ab = t[i + len("Abstract"):j if j > 0 else len(t)]
    raw = [w for w in re.split(r"\s+", ab) if re.search(r"[A-Za-z0-9]", w)]
    words: list[str] = []
    for w in raw:
        if words and words[-1].endswith("-") and re.match(r"[a-z]", w):
            words[-1] = words[-1][:-1] + w
        else:
            words.append(w)
    return len(words), f"（断行复原前 {len(raw)}）"

SMALL_WORDS = {
    "a", "an", "and", "as", "at", "but", "by", "for", "from", "in", "into", "of",
    "on", "or", "the", "to", "with", "vs", "via", "per", "nor", "yet",
}


def strip_comments(body: str) -> str:
    out = []
    for line in body.splitlines():
        i, cut = 0, len(line)
        while i < len(line):
            if line[i] == "%" and (i == 0 or line[i - 1] != "\\"):
                cut = i
                break
            i += 1
        out.append(line[:cut])
    return "\n".join(out)


def brace_args(text: str, macro: str) -> list[str]:
    return [c for _, c in brace_arg_spans(text, macro)]


def brace_arg_spans(text: str, macro: str) -> list[tuple[int, str]]:
    out: list[tuple[int, str]] = []
    for m in re.finditer(r"\\" + macro + r"\s*\{", text):
        i = m.end() - 1
        depth, j = 0, i
        while j < len(text):
            c = text[j]
            if j > i and text[j - 1] == "\\":
                j += 1
                continue
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    out.append((i + 1, text[i + 1:j]))
                    break
            j += 1
    return out


def latex_words(text: str) -> list[str]:
    t = re.sub(r"\\cite[tpalt\*]*\{[^}]*\}", "CIT", text)
    t = re.sub(r"\\ref\{[^}]*\}|\\eqref\{[^}]*\}", "REF", t)
    t = re.sub(r"\\texttt\{([^}]*)\}", r"\1", t)
    t = re.sub(r"\\[a-zA-Z]+\*?", " ", t)
    t = re.sub(r"\$[^$]*\$", "MATH", t)
    t = re.sub(r"[{}\\&%$]", " ", t)
    return [w for w in t.split() if w.strip("-,.;:")]


def get_body(src: str) -> tuple[str, str, str]:
    m = re.search(r"\\begin\{document\}", src)
    pre, rest = src[: m.start()], src[m.start():]
    mb = re.search(r"\\begin\{thebibliography\}\{[^}]*\}(.*?)\\end\{thebibliography\}", rest, re.S)
    if not mb:
        return pre, rest, ""
    return pre, rest[: mb.start()] + rest[mb.end():], mb.group(1)


def bib_entries(bib: str) -> list[dict]:
    entries = []
    pat = re.compile(r"\\bibitem(\[[^\]]*\])?\{([^}]*)\}")
    for m in re.finditer(pat, bib):
        note = m.group(1)[1:-1] if m.group(1) else None
        nxt = bib.find("\\bibitem", m.end())
        raw = bib[m.end(): nxt if nxt != -1 else len(bib)]
        text = " ".join(x.strip() for x in raw.strip().splitlines())
        entries.append({"key": m.group(2), "note": note, "text": text,
                        "pos": m.start()})
    return entries


def note_sortkey(note: str) -> str:
    head = re.split(r"[(，,]", note)[0]
    head = re.sub(r"\\&", "and", head)
    head = re.sub(r"\\[a-zA-Z]+", "", head)
    return " ".join(head.split()).lower()


def note_display(note: str) -> str:
    s = re.sub(r"\\&", "&", note)
    s = re.sub(r"\\o", "ø", s)
    s = re.sub(r"\\[a-zA-Z]+\s?", "", s)
    s = s.replace("(", " ").replace(")", "")
    return " ".join(s.split()).lower()


def main() -> int:
    src = TEX.read_text(encoding="utf-8")
    hard: list[str] = []
    soft: list[str] = []

    pre, body, bib = get_body(strip_comments(src))

    print(f"稿件：{TEX.relative_to(ROOT)}   （{len(src.splitlines())} 行）")
    print("=" * 96)

    dc = re.search(r"\\documentclass(?:\[[^\]]*\])?\{([^}]*)\}", src)
    cls = dc.group(1) if dc else "?"
    cls_ok = cls == "raa"
    print(f"[1] documentclass = {cls}  ->  {'OK' if cls_ok else '需为 raa'}")
    if not cls_ok:
        hard.append(f"documentclass 是 {cls}，RAA 要求 raa.cls")
    has_natbib = re.search(r"\\usepackage(?:\[[^\]]*\])?\{natbib\}", src)
    bibpunct = re.search(r"\\bibpunct\{([^}]*)\}\{([^}]*)\}\{([^}]*)\}\{([^}]*)\}", src)
    print(f"[2] natbib={bool(has_natbib)}  bibpunct={bibpunct.groups() if bibpunct else '缺失'}")
    if not has_natbib:
        hard.append("未载入 natbib（作者-年份制必需）")
    if bibpunct and bibpunct.group(4) != "a":
        hard.append(f"\\bibpunct 第 4 参数是 {bibpunct.group(4)!r}，作者-年份制需 'a'")
    if re.search(r"\\usepackage(?:\[[^\]]*\])?\{cite\}", src):
        hard.append("cite 宏包与 natbib 冲突，必须删除")

    bare = re.findall(r"(?<![\w\\])~?\[\d{1,2}\](?:--\[\d{1,2}\])?(?!\{)", body)
    bare = [b for b in bare if not b.strip().endswith("]")]
    n_cite = len(re.findall(r"\\cite\{", body))
    n_p = len(re.findall(r"\\citep\{", body))
    n_t = len(re.findall(r"\\citet\{", body))
    n_alt = len(re.findall(r"\\citealt\{", body))
    print(f"[3] 裸编号 {len(bare)} 处  |  \\cite {n_cite} / \\citep {n_p} / \\citet {n_t} / \\citealt {n_alt}")
    if bare:
        hard.append(f"正文残留 {len(bare)} 处裸编号引用，例如 {bare[:3]}")

    ab = None
    for _, content in brace_arg_spans(body, "abstract"):
        ab = content
        break
    if ab is None:
        m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", body, re.S)
        ab = m.group(1) if m else None
    n_src = None
    if ab is not None:
        ab = re.sub(r"\\keywords\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}", " ", ab)
        n_src = len(latex_words(ab))
    n_rendered, note = rendered_abstract_words()
    if ab is None:
        print("[4] 摘要 未找到")
        hard.append("找不到摘要环境")
    else:
        low = n_src
        verdict = n_rendered if n_rendered is not None else low
        basis = "渲染真值" if n_rendered is not None else f"源文件下界（{note}）"
        print(f"[4] 摘要词数：{verdict}（判据={basis}；源文件下界 {low}，上限 {ABSTRACT_MAX_WORDS}）")
        if n_rendered is not None and n_src is not None and n_rendered - n_src > 10:
            soft.append(f"摘要渲染值 {n_rendered} 比源文件下界 {n_src} 高 {n_rendered - n_src} 词，"
                        f"下界已不可信，以渲染值为准")
        if n_rendered is None:
            soft.append("摘要词数只用到源文件下界：先编译再看这一项")
        if verdict > ABSTRACT_MAX_WORDS:
            hard.append(f"摘要 {verdict} 词 > {ABSTRACT_MAX_WORDS}（{basis}）")

    kw = re.search(r"\\keywords\{([^}]*)\}", body)
    if kw:
        parts = [x.strip() for x in re.split(r"\s---\s", kw.group(1)) if x.strip()]
        words: list[str] = []
        for p in parts:
            words += [w.strip() for w in re.split(r",(?![^()]*\))", p) if w.strip()]
        print(f"[5] 关键词 {len(words)} 个（上限 {KEYWORDS_MAX}）：{words}")
        if len(words) > KEYWORDS_MAX:
            hard.append(f"关键词 {len(words)} 个 > {KEYWORDS_MAX}")
        joined = " ".join(words)
        if "---" not in joined and len(words) > 1 and ":" not in joined:
            soft.append("关键词未用 RAA/ApJ 的 'genus: species --- genus: species' 体式")
    else:
        print("[5] 关键词 未找到")
        hard.append("找不到 \\keywords")

    ents = bib_entries(bib)
    print(f"[6] 文献条目 {len(ents)} 条")
    no_note = [e["key"] for e in ents if not e["note"]]
    if no_note:
        hard.append(f"{len(no_note)} 条 \\bibitem 缺 [Author et al.(YYYY)] 标签：{no_note[:6]}")
    keyed = [e for e in ents if e["note"]]
    bare_amp = [e["key"] for e in keyed if re.search(r"(?<!\\)&", e["note"])]
    if bare_amp:
        hard.append(f"标签里有未转义的 &（应为 \\&），会在表格里当成对齐符：{bare_amp}")
    keys = [note_display(e["note"]) for e in keyed]
    bad_format = [e["key"] for e in keyed if not re.search(r"\(\d{4}[a-z]?\)", e["note"])]
    if bad_format:
        hard.append(f"{len(bad_format)} 条标签里没有 (YYYY)：{bad_format[:6]}")
    dup = {k for k in keys if keys.count(k) > 1}
    if dup:
        hard.append(f"作者-年份标签碰撞（引用会歧义）：{sorted(dup)[:6]}")
    order = [note_sortkey(e["note"]) for e in keyed]
    unsorted = [(keyed[i]["key"], order[i - 1], order[i])
                for i in range(1, len(order)) if order[i] < order[i - 1]]
    print(f"    字母序：{'已排好' if not unsorted else f'{len(unsorted)} 处倒序'}")
    if unsorted:
        hard.append(f"文献表非字母序，{len(unsorted)} 处倒序，首 3 处：{unsorted[:3]}")
    with_title = []
    for e in ents:
        m = re.search(r"(?<=[a-z\)\.])\s+(19|20)\d{2},", e["text"])
        if m and len(e["text"][: m.start()]) > 60 and re.search(r"[A-Z][a-z]{3,}", e["text"][: m.start()]):
            head = e["text"][: m.start()]
            if re.search(r"\b(?:the|and|of|for|with)\b", head) and "Collaboration" not in head:
                with_title.append(e["key"])
    if with_title:
        hard.append(f"{len(with_title)} 条期刊条目仍带论文标题（RAA 样式不含标题）：{with_title[:5]} ...")

    tbl_labels = re.findall(r"\\label\{(tab:[^}]*)\}", body)
    first_at = {}
    for lab in tbl_labels:
        m = re.search(r"Table~?\\ref\{" + re.escape(lab) + r"\}", body)
        if m:
            first_at[lab] = m.start()
    never = [l for l in tbl_labels if l not in first_at]
    if never:
        hard.append(f"{len(never)} 张表正文未按号提及：{never}")
    else:
        nums = {l: i + 1 for i, l in enumerate(tbl_labels)}
        seq = sorted(first_at, key=first_at.get)
        inv = [(seq[i - 1], nums[seq[i - 1]], seq[i], nums[seq[i]])
               for i in range(1, len(seq)) if nums[seq[i]] < nums[seq[i - 1]]]
        print(f"[7] 表格 {len(tbl_labels)} 张全部被提及；首次提及顺序 {'单调' if not inv else f'{len(inv)} 处倒序'}")
        if inv:
            hard.append(f"表格首次提及非数值升序：{inv[:3]}")

    def cap_problems(text: str) -> list[str]:
        t = re.sub(r"\\(?:label|ref|eqref|cite[tpalt\*]*)\{[^}]*\}", " ", text)
        for _ in range(3):
            t = re.sub(r"\\(?:emph|textit|textbf|texttt|textrm|text|mbox)\{[^{}]*\}", " ", t)
        bad = []
        plain = re.sub(r"\\[a-zA-Z]+|[{}$]", " ", t)
        for w in plain.split():
            core = re.sub(r"[^A-Za-z]", "", w)
            if len(core) < 2 or core.lower() in SMALL_WORDS:
                continue
            if core.isupper():
                continue
            if core[0].islower():
                bad.append(core)
        return bad

    subs = brace_args(body, "subsection") + brace_args(body, "subsubsection")
    caps = brace_arg_spans(body, "caption")
    tab_caps, fig_caps = [], []
    for pos, content in caps:
        t = body.rfind("\\begin{table", 0, pos)
        f = body.rfind("\\begin{figure", 0, pos)
        (tab_caps if t > f else fig_caps).append(content)
    sub_bad = [(t, cap_problems(t)) for t in subs if cap_problems(t)]
    cap_bad = [(t, cap_problems(t)) for t in tab_caps if cap_problems(t)]
    print(f"[8] 已检查 subsection {len(subs)} 条、表题 {len(tab_caps)} 条、"
          f"图题 {len(fig_caps)} 条（图题按 RAA 示例为句首大写，不判逐词大写）；"
          f"小写实词命中 {len(sub_bad)} / {len(cap_bad)} 条")
    for t, p in sub_bad:
        soft.append(f"subsection 需逐词大写：{t[:60]} -> {p[:5]}")
    for t, p in cap_bad:
        soft.append(f"caption 需逐词大写：{t[:60]} -> {p[:5]}")

    need = {"\\volnopage": "编辑保留项 \\volnopage",
            "\\setcounter{page}": "编辑保留项 \\setcounter{page}",
            "\\authorrunning": "\\authorrunning",
            "\\titlerunning": "\\titlerunning",
            "\\maketitle": "\\maketitle",
            "\\begin{acknowledgements}": "acknowledgements 环境",
            "\\label{lastpage}": "\\label{lastpage}"}
    miss = [v for k, v in need.items() if k not in src]
    print(f"[9] raa 构件缺失 {len(miss)} 项")
    for v in miss:
        soft.append(f"缺 {v}（raa.cls 模板要求）")

    print("=" * 96)
    print(f"硬违规 {len(hard)} 项")
    for h in hard:
        print(f"   !! {h}")
    print(f"软提醒 {len(soft)} 项")
    for s in soft:
        print(f"   ~  {s}")
    if not hard and not soft:
        print("RAA 合规：全部检查通过")
    return 1 if hard else 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    raise SystemExit(main())
