"""Reorder bibitems by first citation position."""
import re
import shutil
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
if "--scpma-legacy" not in sys.argv:
    sys.exit("本脚本按 SCPMA 编号制重排，会破坏 RAA 要求的字母序。"
             "确认要回到编号制请加 --scpma-legacy。")

TEX = ROOT / "docs/submission/raa_manuscript.tex"
BAK = TEX.with_name("raa_manuscript.tex.bak_before_reforder")

txt = TEX.read_text(encoding="utf-8")
m = re.search(r"\\begin\{thebibliography\}\{[^}]*\}", txt)
n = re.search(r"\\end\{thebibliography\}", txt)
if not m or not n:
    raise SystemExit("未找到 thebibliography 环境 —— 中止")
head, body, tail = txt[:m.end()], txt[m.end():n.start()], txt[n.start():]

pieces = re.split(r"(?=\s*\\bibitem\{)", body)
blocks, order_in_file = [], []
for p in pieces:
    if not p.strip():
        continue
    k = re.match(r"\s*\\bibitem\{([^}]+)\}", p)
    if not k:
        raise SystemExit(f"块首不是 bibitem，疑似 preamble 残留：{p[:60]!r}")
    blocks.append(p)
    order_in_file.append(k.group(1))

cited = []
for mm in re.finditer(r"\\cite[tp]?\{([^}]*)\}", head):
    for k in mm.group(1).split(","):
        k = k.strip()
        if k and k not in cited:
            cited.append(k)

missing = [k for k in cited if k not in order_in_file]
extra = [k for k in order_in_file if k not in cited]
if missing or extra:
    raise SystemExit(f"键集合不一致，中止。缺={missing} 多={extra}")

bykey = dict(zip(order_in_file, [b.strip() for b in blocks]))
new_body = "\n\n".join(bykey[k] for k in cited)

if sorted(re.findall(r"\\bibitem\{([^}]+)\}", new_body)) != sorted(order_in_file):
    raise SystemExit("重排后条目数/键发生变化 —— 中止")
norm = lambda s: re.sub(r"\s+", "", s)
if sorted(map(norm, blocks)) != sorted(map(norm, [bykey[k] for k in cited])):
    raise SystemExit("条目内容多重集合发生变化 —— 中止")

out = head + "\n\n" + new_body + "\n" + tail
seg = out[out.index("\\bibitem{"):out.index("\\end{thebibliography}")]
after = [b for b in re.split(r"(?=\\bibitem\{)", seg) if b.strip()]
if sorted(map(norm, after)) != sorted(map(norm, [b.strip() for b in blocks])):
    raise SystemExit("写盘前终检：条目内容与原始不一致 —— 中止")

if not BAK.exists():
    shutil.copy(TEX, BAK)
TEX.write_text(out, encoding="utf-8")
print(f"已重排 {len(cited)} 条；备份 -> {BAK.name}")
print("新顺序前 12 条:", cited[:12])
