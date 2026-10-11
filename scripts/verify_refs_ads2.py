"""Verify each bibitem against ADS records using title-phrase matching."""
import json
import os
import re
import time
import pathlib

import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
os.environ.setdefault("NO_PROXY", "*")
TOKEN = json.loads((ROOT / "config/ads.local.json").read_text(encoding="utf-8"))["api_token"]
H = {"Authorization": f"Bearer {TOKEN}"}
URL = "https://api.adsabs.harvard.edu/v1/search/query"
FL = "bibcode,year,v,page,pub,doctypes,title"


def entries_from(tex_path):
    lines = pathlib.Path(tex_path).read_text(encoding="utf-8").splitlines()
    out, i = [], 0
    while i < len(lines):
        m = re.match(r"\s*\\bibitem\{([^}]+)\}", lines[i])
        if not m:
            i += 1
            continue
        key = m.group(1)
        buf = [re.sub(r"^\s*\\bibitem\{[^}]+\}\s*", "", lines[i])]
        j = i + 1
        while j < len(lines) and not re.match(r"\s*\\bibitem\{", lines[j]) \
                and not re.match(r"\s*\\end\{thebibliography\}", lines[j]):
            buf.append(lines[j].strip())
            j += 1
            if len(buf) > 7:
                break
        raw = " ".join(b for b in buf if b)
        raw = re.sub(r"\\cite[tp]?\{[^}]*\}", "", raw)
        raw = re.sub(r"\\(textbf|emph|url|doi|href)\{([^}]*)\}?", r"\2", raw)
        raw = re.sub(r"[{}]", "", raw).replace("\\&", "&").replace("  ", " ").strip()
        out.append((key, raw))
        i = j
    return out


def split_fields(raw):
    parts = re.split(r"(?<=[a-z0-9\)\.])\.\s+(?=[A-Z])", raw, maxsplit=2)
    authors = parts[0] if parts else ""
    title = parts[1] if len(parts) > 1 else ""
    tail = parts[2] if len(parts) > 2 else ""
    yr = re.search(r"\b(19[5-9]\d|20[0-4]\d)\b", tail) or re.search(r"\b(19[5-9]\d|20[0-4]\d)\b", raw)
    vp = re.search(r"(\d{4})\s*,\s*([A-Za-z\.\- ]*?)\s*(\d+)\s*:\s*([0-9A-Za-z\u2013\-]+)", tail)
    journal = (vp.group(2).strip() if vp else (re.search(r",\s*(\d{4})\s*,", tail) and "" or ""))
    jm = re.search(r"\.\s+(.+?),\s*(?:19|20)\d\d", tail)
    return dict(authors=authors.strip(), title=title.strip(), tail=tail.strip(),
                year=(yr.group(1) if yr else None),
                journal=(jm.group(1).strip() if jm else ""),
                vol=(vp.group(3) if vp else None), page=(vp.group(4) if vp else None),
                arxiv=(re.search(r"arXiv:([0-9a-zA-Z\.\-/]+)", raw) or [None, None])[1])


def ads(q, rows=4):
    r = requests.get(URL, headers=H, params={"q": q, "fl": FL, "rows": rows}, timeout=45)
    if r.status_code != 200:
        return None, f"HTTP{r.status_code}:{r.text[:120]}"
    try:
        j = r.json()
    except Exception:
        return None, f"NON-JSON:{r.text[:120]}"
    if "error" in j:
        return None, f"ADS-ERR:{j['error'].get('message','')[:110]}"
    resp = j.get("response") or {}
    docs = resp.get("docs", resp.get("documents"))
    if docs is None:
        return None, f"NODOCS:{str(j)[:120]}"
    return docs, None


ABB = {"ApJS": "Astrophys J Suppl Ser", "ApJ": "Astrophys J", "Nature": "Nature",
       "A&A": "Astron Astrophys", "ARA&A": "Annu Rev Astron Astrophys",
       "Universe": "Universe", "RAA": "Res Astron Astrophys", "A&ARv": "Astron Astrophys Rev"}

items = entries_from(ROOT / "docs/submission/raa_manuscript.tex")
print(f"条目 {len(items)} 条\n")
bad, ok, nf = [], 0, 0
for key, raw in items:
    f = split_fields(raw)
    probe = f["title"] or raw
    probe = re.sub(r"[^A-Za-z0-9 \-:/\(\)]", " ", probe)
    words = probe.split()
    phrase = " ".join(words[:9]) if len(words) > 4 else " ".join(words)
    docs, err = ads(f'abs:"{phrase}"')
    if docs is None:
        bad.append((key, "QUERY-FAIL", err))
        print(f"[{key:16s}] QUERY-FAIL {err}")
        continue
    cand = None
    for d in docs:
        t = " ".join(d.get("title") or [])
        if key == "panstarrs" or "arxiv" in f:
            pass
        cand = cand or d
        if str(d.get("year") or "") == (f["year"] or ""):
            cand = d
            break
        if (d.get("v") and str(d.get("v")) == (f["vol"] or "")):
            cand = d
    if cand is None:
        nf += 1
        bad.append((key, "NOT-FOUND", f'probe="{phrase[:50]}" | 稿内 {f["journal"][:16]} '
                                      f'{f["vol"]}:{f["page"]} {f["year"]}'))
        print(f"[{key:16s}] NOT-FOUND   probe='{phrase[:46]}'")
        continue
    ay, av = str(cand.get("year") or ""), str(cand.get("v") or "")
    ap = cand.get("page")
    ap = (ap[0] if isinstance(ap, list) and ap else (ap or ""))
    notes = []
    if f["year"] and ay and f["year"] != ay:
        notes.append(f"年 稿{f['year']}/ADS{ay}")
    if f["vol"] and av and f["vol"] != av:
        notes.append(f"卷 稿{f['vol']}/ADS{av}")
    if f["page"] and str(ap) and not str(ap).startswith(f["page"].split("\u2013")[0].strip("-")):
        notes.append(f"页 稿{f['page']}/ADS{ap}")
    if not (cand.get("doctypes") and "refereed" in (cand.get("doctypes") or [])):
        notes.append("非期刊正文(refereed) 类型=" + ",".join(cand.get("doctypes") or [])[:26])
    tag = "OK " if not notes else "CHK"
    if notes:
        bad.append((key, tag.strip(), "; ".join(notes)))
    else:
        ok += 1
    print(f"[{key:16s}] {tag} {cand.get('bibcode'):24s} {str(cand.get('pub'))[:26]:26s} "
          f"v{av:5s} p{str(ap):9s} {'; '.join(notes)}")
    time.sleep(0.35)

print(f"\n完全一致 {ok} / 未命中 {nf} / 共 {len(items)}")
print("\n== 需人工定夺 ==")
for k, kind, det in bad:
    print(f"  {k:16s} {kind:11s} {str(det)[:130]}")
